function result = run_sign_pso(context, varargin)
%RUN_SIGN_PSO Search descriptor subsets using bounded, sign-encoded PSO.
%   CONTEXT fixes the compounds, medoid and objective J=-R*F*P. Named options
%   specify particle count, generations, seed and the search/mutation schedules.
%   Positions and velocities are nParticle-by-p; z=(position>0) defines a subset.
%   Position magnitudes are retained during the continuous update x=x+v.
%
%   RESULT contains the final subset, J/R/F/P and best-so-far scalar histories.
%   Selection compares K=(nR-1)*J with strict improvement; the first minimum
%   is retained on exact ties. Reported J and P use the mean-percentile scale.
%   All requested generations run. StallWindow assesses final-window change
%   after the search; it does not stop the optimizer. Caller RNG is restored.

%% Search options and context checks

parser = inputParser;
parser.FunctionName = mfilename;
addParameter(parser, 'PopulationSize', 50, @localPositiveInteger);
addParameter(parser, 'GenerationCount', 3000, @localPositiveInteger);
addParameter(parser, 'ScheduleHorizon', [], ...
    @(value) isempty(value) || localPositiveInteger(value));
addParameter(parser, 'RandomSeed', 0, @localNonnegativeInteger);
addParameter(parser, 'InertiaWeightStart', 0.80, @localNonnegativeScalar);
addParameter(parser, 'InertiaWeightEnd', 0.35, @localNonnegativeScalar);
addParameter(parser, 'CognitiveCoefficient', 1.49618, @localNonnegativeScalar);
addParameter(parser, 'SocialCoefficient', 1.49618, @localNonnegativeScalar);
addParameter(parser, 'ParticleMutationRateStart', 0.10, @localProbability);
addParameter(parser, 'ParticleMutationRateEnd', 0.03, @localProbability);
addParameter(parser, 'ParticleMutationRateFloor', 0.003, @localProbability);
addParameter(parser, 'FeatureMutationRateStart', 0.008, @localProbability);
addParameter(parser, 'FeatureMutationRateEnd', 0.001, @localProbability);
addParameter(parser, 'FeatureMutationRateFloor', 8e-5, @localProbability);
addParameter(parser, 'MutationSchedulePower', 2.0, @localPositiveFiniteScalar);
addParameter(parser, 'MutationQuenchStartFraction', 0.85, @localOpenUnitInterval);
addParameter(parser, 'MutationNoiseScale', 0.10, @localNonnegativeScalar);
addParameter(parser, 'VelocityLimit', 4.0, @localPositiveFiniteScalar);
addParameter(parser, 'PositionLimit', 4.0, @localPositiveFiniteScalar);
addParameter(parser, 'StallWindow', 200, @localPositiveInteger);
addParameter(parser, 'StallRelativeTolerance', 1e-4, @localPositiveFiniteScalar);
addParameter(parser, 'Verbose', true, @localLogicalScalar);
addParameter(parser, 'VerboseInterval', 25, @localPositiveInteger);
parse(parser, varargin{:});
options = parser.Results;

requiredFields = {'FeatureCount', 'ReferenceIndex', ...
    'TargetCountIncludingReference', 'TargetCountEvaluated', ...
    'BackgroundCount', 'DistanceEpsilon', 'ObjectiveDefinition', ...
    'PAggregation', 'PDefinition'};
for idxField = 1:numel(requiredFields)
    if ~isfield(context, requiredFields{idxField})
        error('run_sign_pso:InvalidContext', ...
            'Context is missing field %s.', requiredFields{idxField});
    end
end

expectedObjective = '-((D2-D1)/(D1+epsilon))*F*P';
if ~strcmp(context.ObjectiveDefinition, expectedObjective)
    error('run_sign_pso:ObjectiveMismatch', ...
        'Context does not specify the PSO-MDL relative-gap objective.');
end
if ~strcmp(context.PAggregation, 'mean_non_medoid')
    error('run_sign_pso:PAggregationMismatch', ...
        'The current fitness requires mean non-medoid proximity percentiles.');
end

particleCount = options.PopulationSize;
generationCount = options.GenerationCount;
if isempty(options.ScheduleHorizon)
    scheduleHorizon = generationCount;
else
    scheduleHorizon = options.ScheduleHorizon;
end
if scheduleHorizon > generationCount
    error('run_sign_pso:InvalidScheduleHorizon', ...
        'ScheduleHorizon must not exceed GenerationCount.');
end
featureCount = context.FeatureCount;

if options.InertiaWeightEnd > options.InertiaWeightStart
    error('run_sign_pso:InvalidInertiaSchedule', ...
        'InertiaWeightEnd must not exceed InertiaWeightStart.');
end
if options.ParticleMutationRateEnd > options.ParticleMutationRateStart
    error('run_sign_pso:InvalidParticleMutationSchedule', ...
        'ParticleMutationRateEnd must not exceed ParticleMutationRateStart.');
end
if options.FeatureMutationRateEnd > options.FeatureMutationRateStart
    error('run_sign_pso:InvalidFeatureMutationSchedule', ...
        'FeatureMutationRateEnd must not exceed FeatureMutationRateStart.');
end
if options.ParticleMutationRateFloor > options.ParticleMutationRateEnd
    error('run_sign_pso:InvalidParticleMutationFloor', ...
        'ParticleMutationRateFloor must not exceed ParticleMutationRateEnd.');
end
if options.FeatureMutationRateFloor > options.FeatureMutationRateEnd
    error('run_sign_pso:InvalidFeatureMutationFloor', ...
        'FeatureMutationRateFloor must not exceed FeatureMutationRateEnd.');
end
if options.StallWindow >= generationCount
    error('run_sign_pso:InvalidStallWindow', ...
        'StallWindow must be smaller than GenerationCount.');
end

%% Seeded initialization
% Generation 1 records the initial population; updates begin at generation 2.
savedRngState = rng;
cleanupRandomState = onCleanup(@() rng(savedRngState));
rng(options.RandomSeed, 'twister');

runTimer = tic;

population = 2 .* rand(particleCount, featureCount) - 1;
population = localAvoidExactZero(population);
velocity = zeros(particleCount, featureCount);

personalBestPosition = population;
[populationFitness, populationDiagnostics] = ...
    evaluate_pso_mdl_fitness(population, context);
personalBestFitness = populationFitness;
personalBestComparisonKey = populationDiagnostics.ComparisonKey;
personalBestD1 = populationDiagnostics.D1;
personalBestD2 = populationDiagnostics.D2;
personalBestSeparation = populationDiagnostics.Separation;
personalBestRelativeGap = populationDiagnostics.RelativeGap;
personalBestD2OverD1 = populationDiagnostics.D2OverD1;
personalBestP = populationDiagnostics.P;
personalBestF = populationDiagnostics.F;

[globalBestComparisonKey, globalBestIndex] = min(personalBestComparisonKey);
globalBestFitness = personalBestFitness(globalBestIndex);
if ~(isfinite(globalBestFitness) && globalBestFitness < 0)
    error('run_sign_pso:NoFeasibleInitialParticle', ...
        ['No initial particle has negative fitness with D2>D1 and D1>epsilon;' ...
         ' check data, epsilon, and fitness-validity rules.']);
end
globalBestPosition = personalBestPosition(globalBestIndex, :);

bestFitnessByGeneration = zeros(generationCount, 1);
bestComparisonKeyByGeneration = zeros(generationCount, 1);
selectedFeatureCountByGeneration = zeros(generationCount, 1);
bestD1ByGeneration = zeros(generationCount, 1);
bestD2ByGeneration = zeros(generationCount, 1);
bestRelativeGapByGeneration = zeros(generationCount, 1);
bestPByGeneration = zeros(generationCount, 1);

[currentInertiaWeight, currentParticleMutationRate, ...
    currentFeatureMutationRate, ~] = ...
    localScheduleValues(1);
recordGeneration(1);

if options.Verbose
    fprintf(['PSO-MDL relative-gap seed %d: generation 1/%d, J=%.9g, ' ...
        'F=%d, relative gap=%.6g, P=%.6g\n'], ...
        options.RandomSeed, generationCount, globalBestFitness, ...
        personalBestF(globalBestIndex), ...
        personalBestRelativeGap(globalBestIndex), ...
        personalBestP(globalBestIndex));
end

%% Particle dynamics, sign-flip mutation and best-so-far updates
for idxGeneration = 2:generationCount
    [currentInertiaWeight, currentParticleMutationRate, ...
        currentFeatureMutationRate, ~] = ...
        localScheduleValues(idxGeneration);
    % Draw independent learning coefficients for each particle-feature pair.
    r1 = rand(particleCount, featureCount);
    r2 = rand(particleCount, featureCount);
    globalBestMatrix = repmat(globalBestPosition, particleCount, 1);

    velocity = currentInertiaWeight .* velocity ...
        + options.CognitiveCoefficient .* r1 .* ...
            (personalBestPosition - population) ...
        + options.SocialCoefficient .* r2 .* ...
            (globalBestMatrix - population);

    velocity = min(max(velocity, -options.VelocityLimit), ...
        options.VelocityLimit);

    % Preserve continuous magnitudes; only their signs encode feature selection.
    population = population + velocity;
    population = min(max(population, -options.PositionLimit), ...
        options.PositionLimit);

    % Two-level mutation chooses particles first, then sites within those particles.
    mutationFlipCount = 0;
    mutatedParticleMask = rand(particleCount, 1) ...
        < currentParticleMutationRate;
    mutatedParticleIndex = find(mutatedParticleMask);
    if ~isempty(mutatedParticleIndex)
        mutatedSiteMask = rand(numel(mutatedParticleIndex), featureCount) ...
            < currentFeatureMutationRate;
        mutationFlipCount = nnz(mutatedSiteMask);

        selectedPosition = population(mutatedParticleIndex, :);
        if mutationFlipCount > 0
            currentValues = selectedPosition(mutatedSiteMask);
            currentValues = localAvoidExactZero(currentValues);

            % Positive jitter cannot undo the intended sign flip.
            jitterFactor = 1 + options.MutationNoiseScale .* ...
                randn(size(currentValues));
            jitterFactor = max(jitterFactor, 0.05);
            selectedPosition(mutatedSiteMask) = ...
                -currentValues .* jitterFactor;
            population(mutatedParticleIndex, :) = selectedPosition;

            selectedVelocity = velocity(mutatedParticleIndex, :);
            selectedVelocity(mutatedSiteMask) = 0;
            velocity(mutatedParticleIndex, :) = selectedVelocity;
        end
    end

    population = min(max(population, -options.PositionLimit), ...
        options.PositionLimit);
    population = localAvoidExactZero(population);

    [populationFitness, populationDiagnostics] = ...
        evaluate_pso_mdl_fitness(population, context);

    % Compare K directly, not a rescaled rounded J; keep strict tie behavior.
    improvedParticleIndex = find( ...
        populationDiagnostics.ComparisonKey < personalBestComparisonKey);
    personalBestPosition(improvedParticleIndex, :) = ...
        population(improvedParticleIndex, :);
    personalBestFitness(improvedParticleIndex) = ...
        populationFitness(improvedParticleIndex);
    personalBestComparisonKey(improvedParticleIndex) = ...
        populationDiagnostics.ComparisonKey(improvedParticleIndex);
    personalBestD1(improvedParticleIndex) = ...
        populationDiagnostics.D1(improvedParticleIndex);
    personalBestD2(improvedParticleIndex) = ...
        populationDiagnostics.D2(improvedParticleIndex);
    personalBestSeparation(improvedParticleIndex) = ...
        populationDiagnostics.Separation(improvedParticleIndex);
    personalBestRelativeGap(improvedParticleIndex) = ...
        populationDiagnostics.RelativeGap(improvedParticleIndex);
    personalBestD2OverD1(improvedParticleIndex) = ...
        populationDiagnostics.D2OverD1(improvedParticleIndex);
    personalBestP(improvedParticleIndex) = ...
        populationDiagnostics.P(improvedParticleIndex);
    personalBestF(improvedParticleIndex) = ...
        populationDiagnostics.F(improvedParticleIndex);

    [globalBestComparisonKey, globalBestIndex] = min(personalBestComparisonKey);
    globalBestFitness = personalBestFitness(globalBestIndex);
    globalBestPosition = personalBestPosition(globalBestIndex, :);
    recordGeneration(idxGeneration);

    if options.Verbose && (idxGeneration == 2 ...
            || mod(idxGeneration, options.VerboseInterval) == 0 ...
            || idxGeneration == generationCount)
        fprintf(['PSO-MDL seed %d: generation %d/%d, J=%.9g, ' ...
            'F=%d, R=%.6g, P=%.6g\n'], ...
            options.RandomSeed, idxGeneration, generationCount, ...
            globalBestFitness, personalBestF(globalBestIndex), ...
            personalBestRelativeGap(globalBestIndex), personalBestP(globalBestIndex));
    end
end

%% Post-run convergence assessment
% Relative change is (J_(G-window)-J_G)/max(abs(J_(G-window)),eps).
% This retrospective diagnostic never stops or steers the search.
objectiveTolerance = 1e-12 * max(1, max(abs(bestFitnessByGeneration)));
lastImprovementIndex = find( ...
    diff(bestFitnessByGeneration) < -objectiveTolerance, 1, 'last');
if isempty(lastImprovementIndex)
    lastImprovementGeneration = 1;
else
    lastImprovementGeneration = lastImprovementIndex + 1;
end

stallReferenceGeneration = generationCount - options.StallWindow;
stallReferenceFitness = bestFitnessByGeneration(stallReferenceGeneration);
finalWindowRelativeImprovement = ...
    (stallReferenceFitness - bestFitnessByGeneration(end)) ./ ...
    max(abs(stallReferenceFitness), eps);
convergedByStallCriterion = ...
    finalWindowRelativeImprovement < options.StallRelativeTolerance;

monotonicViolationCount = sum( ...
    diff(bestFitnessByGeneration) > objectiveTolerance);

if monotonicViolationCount ~= 0
    error('run_sign_pso:NonMonotonicBestHistory', ...
        'Global best-so-far fitness has %d monotonicity violations.', ...
        monotonicViolationCount);
end
if abs(bestFitnessByGeneration(end) - globalBestFitness) > objectiveTolerance
    error('run_sign_pso:FinalBestMismatch', ...
        'Final history value disagrees with globalBestFitness.');
end

%% Export the final model and scalar optimization histories
result = struct();
result.Version = 6;
result.MethodLabel = ...
    ['PSO-MDL sign-thresholded PSO with fixed-horizon time-varying inertia, ' ...
     'decaying sign-flip mutation with nonzero floor and bounded dynamics'];
result.ObjectiveDefinition = expectedObjective;
result.PAggregation = 'mean_non_medoid';
result.PDefinition = context.PDefinition;
result.PRange = [0, 100];
result.ComparisonKeyDefinition = 'K=-R*F*sum(p_non_medoid); J=K/(nRef-1)';
result.ComparisonKeyScale = context.TargetCountEvaluated;
result.ComparisonRule = 'strict improvement; first minimum on ties';
result.SelectionDefinition = 'position > 0';
result.SearchDefinition = ...
    ['continuous bounded latent positions; descriptor selected when ' ...
     'position>0; elementwise r1/r2; x=x+v; fixed-horizon linearly decreasing inertia; ' ...
     'quadratically decaying sign-flip mutation with a nonzero late-stage floor; ' ...
     'no per-generation magnitude reset'];
result.DistanceEpsilon = context.DistanceEpsilon;
result.D1ValidityRule = 'D1 > epsilon';
result.RandomSeed = options.RandomSeed;
result.PopulationSize = particleCount;
result.GenerationCount = generationCount;
result.ScheduleHorizon = scheduleHorizon;
result.InertiaWeightStart = options.InertiaWeightStart;
result.InertiaWeightEnd = options.InertiaWeightEnd;
result.CognitiveCoefficient = options.CognitiveCoefficient;
result.SocialCoefficient = options.SocialCoefficient;
result.ParticleMutationRateStart = options.ParticleMutationRateStart;
result.ParticleMutationRateEnd = options.ParticleMutationRateEnd;
result.ParticleMutationRateFloor = options.ParticleMutationRateFloor;
result.FeatureMutationRateStart = options.FeatureMutationRateStart;
result.FeatureMutationRateEnd = options.FeatureMutationRateEnd;
result.FeatureMutationRateFloor = options.FeatureMutationRateFloor;
result.MutationSchedulePower = options.MutationSchedulePower;
result.MutationQuenchStartFraction = ...
    options.MutationQuenchStartFraction;
result.MutationNoiseScale = options.MutationNoiseScale;
result.StallWindow = options.StallWindow;
result.StallRelativeTolerance = options.StallRelativeTolerance;
result.VelocityLimit = options.VelocityLimit;
result.PositionLimit = options.PositionLimit;
result.ReferenceIndex = context.ReferenceIndex;
result.TargetCount = context.TargetCountIncludingReference;
result.TargetCountEvaluated = context.TargetCountEvaluated;
result.BackgroundCount = context.BackgroundCount;
result.FeatureCount = featureCount;
result.BestPosition = globalBestPosition;
result.BestMask = globalBestPosition > 0;
result.BestFitness = globalBestFitness;
result.BestComparisonKey = globalBestComparisonKey;
result.BestD1 = personalBestD1(globalBestIndex);
result.BestD2 = personalBestD2(globalBestIndex);
result.BestSeparation = personalBestSeparation(globalBestIndex);
result.BestRelativeGap = personalBestRelativeGap(globalBestIndex);
result.BestD2OverD1 = personalBestD2OverD1(globalBestIndex);
result.BestP = personalBestP(globalBestIndex);
result.BestSelectedFeatureCount = personalBestF(globalBestIndex);
result.BestFitnessByGeneration = bestFitnessByGeneration;
result.BestComparisonKeyByGeneration = bestComparisonKeyByGeneration;
result.SelectedFeatureCountByGeneration = ...
    selectedFeatureCountByGeneration;
result.BestD1ByGeneration = bestD1ByGeneration;
result.BestD2ByGeneration = bestD2ByGeneration;
result.BestRelativeGapByGeneration = bestRelativeGapByGeneration;
result.BestPByGeneration = bestPByGeneration;
result.FinalWindowRelativeImprovement = finalWindowRelativeImprovement;
result.ConvergedByStallCriterion = convergedByStallCriterion;
result.LastImprovementGeneration = lastImprovementGeneration;
result.MonotonicViolationCount = monotonicViolationCount;
result.RuntimeSeconds = toc(runTimer);
result.CreatedAt = char(datetime('now', ...
    'Format', "yyyyMMdd'T'HHmmss"));

    function [wValue, particleRate, featureRate, progress] = ...
            localScheduleValues(generationIndex)
        % Progress is tied to the fixed schedule horizon, not fitness history.
        if scheduleHorizon <= 1
            progress = 1;
        else
            progress = min((generationIndex - 1) ./ ...
                (scheduleHorizon - 1), 1);
        end

        wValue = options.InertiaWeightStart + ...
            (options.InertiaWeightEnd - options.InertiaWeightStart) .* ...
            progress;

        decay = (1 - progress) .^ options.MutationSchedulePower;
        particleRate = options.ParticleMutationRateEnd + ...
            (options.ParticleMutationRateStart - ...
             options.ParticleMutationRateEnd) .* decay;
        featureRate = options.FeatureMutationRateEnd + ...
            (options.FeatureMutationRateStart - ...
             options.FeatureMutationRateEnd) .* decay;

        if progress >= options.MutationQuenchStartFraction
            quenchProgress = ...
                (progress - options.MutationQuenchStartFraction) ./ ...
                (1 - options.MutationQuenchStartFraction);
            quenchFactor = max(0, 1 - quenchProgress) .^ 2;
            particleRate = options.ParticleMutationRateFloor + ...
                (particleRate - options.ParticleMutationRateFloor) .* ...
                quenchFactor;
            featureRate = options.FeatureMutationRateFloor + ...
                (featureRate - options.FeatureMutationRateFloor) .* ...
                quenchFactor;
        end
    end

    function recordGeneration(generationIndex)
        bestFitnessByGeneration(generationIndex) = globalBestFitness;
        bestComparisonKeyByGeneration(generationIndex) = globalBestComparisonKey;
        selectedFeatureCountByGeneration(generationIndex) = ...
            personalBestF(globalBestIndex);
        bestD1ByGeneration(generationIndex) = ...
            personalBestD1(globalBestIndex);
        bestD2ByGeneration(generationIndex) = ...
            personalBestD2(globalBestIndex);
        bestRelativeGapByGeneration(generationIndex) = ...
            personalBestRelativeGap(globalBestIndex);
        bestPByGeneration(generationIndex) = ...
            personalBestP(globalBestIndex);
    end
end

function position = localAvoidExactZero(position)
zeroMask = position == 0;
if any(zeroMask, 'all')
    zeroCount = nnz(zeroMask);
    randomSigns = 2 .* (rand(zeroCount, 1) > 0.5) - 1;
    position(zeroMask) = randomSigns .* eps;
end
end

function tf = localPositiveInteger(value)
tf = isnumeric(value) && isscalar(value) && isfinite(value) ...
    && value == fix(value) && value >= 1;
end

function tf = localNonnegativeInteger(value)
tf = isnumeric(value) && isscalar(value) && isfinite(value) ...
    && value == fix(value) && value >= 0;
end

function tf = localFiniteScalar(value)
tf = isnumeric(value) && isscalar(value) && isfinite(value);
end

function tf = localNonnegativeScalar(value)
tf = localFiniteScalar(value) && value >= 0;
end

function tf = localPositiveFiniteScalar(value)
tf = localFiniteScalar(value) && value > 0;
end

function tf = localProbability(value)
tf = localFiniteScalar(value) && value >= 0 && value <= 1;
end

function tf = localOpenUnitInterval(value)
tf = localFiniteScalar(value) && value > 0 && value < 1;
end

function tf = localLogicalScalar(value)
tf = (islogical(value) && isscalar(value)) ...
    || (isnumeric(value) && isscalar(value) ...
        && isfinite(value) && (value == 0 || value == 1));
end

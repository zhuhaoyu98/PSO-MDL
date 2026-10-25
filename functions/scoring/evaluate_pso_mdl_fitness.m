function [fitnessValue, diagnostics] = ...
    evaluate_pso_mdl_fitness(populationPosition, context)
%EVALUATE_PSO_MDL_FITNESS Evaluate J = -R*F*P for a population of descriptor masks.
%   POPULATIONPOSITION is nParticle-by-p; z_j = (position_j > 0) selects columns.
%   CONTEXT is the fixed-medoid distance cache from build_fitness_context.
%
%   For each particle, d(x,z) = sqrt(sum_j z_j*(x_j-x_medoid,j)^2).
%   D1 = mean_t(d_t), t excluding the medoid; D2 = mean_b(d_b).
%   R = (D2-D1)/(D1+epsilon); F = sum_j z_j.
%   p_t = 100*count_b(d_b > d_t)/nB; P = sum_t(p_t)/(nR-1).
%   Equal distances do not contribute to p_t. P remains on the 0-100 scale:
%   it measures background-relative proximity, not absolute reference compactness.
%
%   Outputs J and DIAGNOSTICS fields are nParticle-by-1 vectors. Empty subsets
%   and D1<=epsilon receive J=0. All other objective terms retain their signs.
%   DIAGNOSTICS.ComparisonKey is K=-R*F*sum_t(p_t), evaluated before averaging.
%   J=K/(nR-1); PSO compares K directly to preserve floating-point tie decisions.

if ~(isnumeric(populationPosition) || islogical(populationPosition)) ...
        || ~isreal(populationPosition) ...
        || ~ismatrix(populationPosition)
    error('evaluate_pso_mdl_fitness:InvalidPopulation', ...
        'Particle positions must be a real two-dimensional numeric matrix or logical masks.');
end
if any(~isfinite(populationPosition), 'all')
    error('evaluate_pso_mdl_fitness:NonFinitePopulation', ...
        'Particle positions must not contain NaN or Inf.');
end

requiredFields = {'TargetSquaredDelta', 'BackgroundSquaredDelta', ...
    'TargetCountIncludingReference', 'TargetCountEvaluated', ...
    'BackgroundCount', 'FeatureCount', ...
    'DistanceEpsilon'};
for idxField = 1:numel(requiredFields)
    if ~isfield(context, requiredFields{idxField})
        error('evaluate_pso_mdl_fitness:InvalidContext', ...
            'Context is missing field %s.', requiredFields{idxField});
    end
end

countFields = {'TargetCountIncludingReference', 'TargetCountEvaluated', ...
    'BackgroundCount', 'FeatureCount'};
for fieldIndex = 1:numel(countFields)
    value = context.(countFields{fieldIndex});
    if ~isnumeric(value) || ~isreal(value) || ~isscalar(value) ...
            || ~isfinite(value) || value < 1 || value ~= fix(value)
        error('evaluate_pso_mdl_fitness:InvalidCounts', ...
            'Context counts must be finite positive integer scalars.');
    end
end
if context.TargetCountIncludingReference < 2 ...
        || context.TargetCountEvaluated ~= context.TargetCountIncludingReference - 1
    error('evaluate_pso_mdl_fitness:InvalidCounts', ...
        'Expected nRef>=2, nNonMedoid=nRef-1, nBackground>=1 and nFeature>=1.');
end
if ~isequal(size(context.TargetSquaredDelta), ...
        [context.TargetCountEvaluated, context.FeatureCount]) ...
        || ~isequal(size(context.BackgroundSquaredDelta), ...
        [context.BackgroundCount, context.FeatureCount])
    error('evaluate_pso_mdl_fitness:ContextSizeMismatch', ...
        'Squared-distance matrices must match the reference/background counts.');
end
if ~(isnumeric(context.DistanceEpsilon) ...
        && isscalar(context.DistanceEpsilon) ...
        && isfinite(context.DistanceEpsilon) ...
        && context.DistanceEpsilon > 0)
    error('evaluate_pso_mdl_fitness:InvalidEpsilon', ...
        'context.DistanceEpsilon must be a finite, strictly positive scalar.');
end

[particleCount, featureCount] = size(populationPosition);
if featureCount ~= context.FeatureCount
    error('evaluate_pso_mdl_fitness:FeatureCountMismatch', ...
        'Population has %d features; context requires %d.', ...
        featureCount, context.FeatureCount);
end

%% Descriptor selection and output allocation
% F counts selected columns, not the sum or magnitude of the latent positions.
selectedMask = populationPosition > 0;
selectedFeatureCount = sum(selectedMask, 2);
fitnessValue = zeros(particleCount, 1);
comparisonKey = zeros(particleCount, 1);
D1 = zeros(particleCount, 1);
D2 = zeros(particleCount, 1);
relativeGap = zeros(particleCount, 1);
D2OverD1 = zeros(particleCount, 1);
P = zeros(particleCount, 1);
separation = zeros(particleCount, 1);
validRelativeGapMask = false(particleCount, 1);

activeParticleIndex = find(selectedFeatureCount > 0);
if ~isempty(activeParticleIndex)
    % Columns now represent active particles: distances are (nR-1)-by-nActive
    % and nB-by-nActive. Reductions below average over compounds, not particles.
    activeMaskMatrix = double(selectedMask(activeParticleIndex, :)).';
    targetDistance = sqrt( ...
        context.TargetSquaredDelta * activeMaskMatrix);
    backgroundDistance = sqrt( ...
        context.BackgroundSquaredDelta * activeMaskMatrix);

    if ~isreal(targetDistance) || ~isreal(backgroundDistance) ...
            || any(~isfinite(targetDistance), 'all') ...
            || any(~isfinite(backgroundDistance), 'all')
        error('evaluate_pso_mdl_fitness:InvalidDistances', ...
            'Reference and background distances must be real and finite.');
    end

    activeD1 = mean(targetDistance, 1).';
    activeD2 = mean(backgroundDistance, 1).';
    activeP = zeros(numel(activeParticleIndex), 1);
    % Reference-specific proximity percentiles: use the strict upper tail of
    % the background distances. Accumulate before dividing by nR-1, not nR.
    for idxTarget = 1:context.TargetCountEvaluated
        activeP = activeP + sum( ...
            backgroundDistance > targetDistance(idxTarget, :), 1).' ...
            ./ context.BackgroundCount .* 100;
    end

    % The context already excludes the medoid; this is the arithmetic mean P.
    activePSum = activeP;
    activeP = activeP ./ context.TargetCountEvaluated;

    activeSeparation = activeD2 - activeD1;
    % Besides adding epsilon to the denominator, reject D1<=epsilon by the
    % existing zero-fitness rule. Do not replace this with max(D1,epsilon).
    activeValid = isfinite(activeD1) ...
        & activeD1 > context.DistanceEpsilon ...
        & isfinite(activeD2) & isfinite(activeP) ...
        & isfinite(activeSeparation);
    activeRelativeGap = zeros(numel(activeParticleIndex), 1);
    activeD2OverD1 = zeros(numel(activeParticleIndex), 1);
    activeFitness = zeros(numel(activeParticleIndex), 1);
    activeComparisonKey = zeros(numel(activeParticleIndex), 1);
    activeRelativeGap(activeValid) = ...
        activeSeparation(activeValid) ./ ...
        (activeD1(activeValid) + context.DistanceEpsilon);

    activeD2OverD1(activeValid) = ...
        activeD2(activeValid) ./ activeD1(activeValid);
    % K=(nR-1)*J has the same minimizer for a fixed reference set in exact
    % arithmetic. Compute K first: rescaling an already rounded J may merge ties.
    activeComparisonKey(activeValid) = -activeRelativeGap(activeValid) ...
        .* selectedFeatureCount(activeParticleIndex(activeValid)) ...
        .* activePSum(activeValid);
    activeFitness(activeValid) = activeComparisonKey(activeValid) ...
        ./ context.TargetCountEvaluated;

    D1(activeParticleIndex) = activeD1;
    D2(activeParticleIndex) = activeD2;
    relativeGap(activeParticleIndex) = activeRelativeGap;
    D2OverD1(activeParticleIndex) = activeD2OverD1;
    P(activeParticleIndex) = activeP;
    separation(activeParticleIndex) = activeSeparation;
    validRelativeGapMask(activeParticleIndex) = activeValid;
    fitnessValue(activeParticleIndex) = activeFitness;
    comparisonKey(activeParticleIndex) = activeComparisonKey;
end

if any(~isfinite(fitnessValue)) || any(~isfinite(comparisonKey))
    error('evaluate_pso_mdl_fitness:NonFiniteFitness', ...
        'Fitness evaluation produced NaN or Inf.');
end

% Legacy field RelativeGap denotes R; Separation denotes the raw D2-D1.
if nargout > 1
    diagnostics = struct();
    diagnostics.F = selectedFeatureCount;
    diagnostics.D1 = D1;
    diagnostics.D2 = D2;
    diagnostics.Separation = separation;
    diagnostics.RelativeGap = relativeGap;
    diagnostics.D2OverD1 = D2OverD1;
    diagnostics.P = P;
    diagnostics.PAggregation = 'mean_non_medoid';
    diagnostics.PRange = [0, 100];
    diagnostics.DistanceEpsilon = context.DistanceEpsilon;
    diagnostics.ValidRelativeGapMask = validRelativeGapMask;
    diagnostics.Fitness = fitnessValue;
    diagnostics.ComparisonKey = comparisonKey;
    diagnostics.ComparisonKeyScale = context.TargetCountEvaluated;
end
end

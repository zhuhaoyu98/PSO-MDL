function [fitnessValue, diagnostics] = ...
    evaluate_pso_mdl_fitness(populationPosition, context)
%EVALUATE_PSO_MDL_FITNESS Evaluate the batched mean-P objective J=-R*F*P.
%   Positions are particle-by-feature; position>0 selects a descriptor.
%   F counts selected descriptors. D1/D2 average Euclidean distances from the
%   fixed reference medoid to other references/background molecules.
%   R=(D2-D1)/(D1+epsilon). P averages 100*mean(d_background>d_t)
%   over nRef-1 non-medoid references, retaining the 0-100 percentage scale.
%   Equal distances are excluded. P is relative proximity, not absolute compactness.
%   Empty selections and D1<=epsilon retain the validated zero-fitness rule.
%   Outputs are particle-by-1 fitness and diagnostic vectors.
%   Diagnostics.ComparisonKey is K=-R*F*sum(p_t), computed before averaging.
%   K=(nRef-1)*J in exact arithmetic; PSO compares K to retain the historical
%   floating-point tie behavior. Returned fitness is J=K/(nRef-1), not K.

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
    % Sum unchanged strict-tail percentiles, then average over references.
    for idxTarget = 1:context.TargetCountEvaluated
        activeP = activeP + sum( ...
            backgroundDistance > targetDistance(idxTarget, :), 1).' ...
            ./ context.BackgroundCount .* 100;
    end

    % The medoid was excluded when the fixed context was built.
    activePSum = activeP;
    activeP = activeP ./ context.TargetCountEvaluated;

    activeSeparation = activeD2 - activeD1;
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
    % Preserve the historical multiplication order before normalizing J.
    % Multiplying an already rounded J cannot recover lost tie distinctions.
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

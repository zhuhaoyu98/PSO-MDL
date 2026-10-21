function context = build_fitness_context( ...
    targetData, backgroundData, distanceEpsilon)
%BUILD_FITNESS_CONTEXT Build the fixed-medoid context for J=-R*F*P.
%   Reference and background inputs are finite row-by-feature double matrices.
%   The medoid minimizes total Euclidean distance in the full retained space.
%   R=(D2-D1)/(D1+epsilon); P is the mean non-medoid background-relative
%   proximity percentile (0-100, strict background distance > reference distance).
%   Squared per-feature differences are cached without changing distance arithmetic.

arguments
    targetData double
    backgroundData double
    distanceEpsilon (1, 1) double = 1e-12
end

if ~ismatrix(targetData) || ~ismatrix(backgroundData) ...
        || ~isreal(targetData) || ~isreal(backgroundData)
    error('build_fitness_context:InvalidRank', ...
        'Reference and background inputs must be two-dimensional matrices.');
end
if size(targetData, 1) < 2 || isempty(backgroundData)
    error('build_fitness_context:InsufficientRows', ...
        'At least two references and one background molecule are required.');
end
if size(targetData, 2) ~= size(backgroundData, 2) ...
        || size(targetData, 2) < 1
    error('build_fitness_context:FeatureCountMismatch', ...
        'Reference and background must have equal, nonzero feature counts.');
end
if any(~isfinite(targetData), 'all') ...
        || any(~isfinite(backgroundData), 'all')
    error('build_fitness_context:NonFiniteInput', ...
        'Inputs must be finite, imputed, scaled numeric matrices.');
end
if ~(isfinite(distanceEpsilon) && distanceEpsilon > 0)
    error('build_fitness_context:InvalidEpsilon', ...
        'distanceEpsilon must be a finite, strictly positive scalar.');
end

% Determine the medoid once in the full retained descriptor space.
rowSquaredNorm = sum(targetData.^2, 2);
squaredDistance = rowSquaredNorm + rowSquaredNorm.' ...
    - 2 .* (targetData * targetData.');
squaredDistance = max(squaredDistance, 0);
[~, referenceIndex] = min(sum(sqrt(squaredDistance), 2));

referenceDrug = targetData(referenceIndex, :);
targetMask = true(size(targetData, 1), 1);
targetMask(referenceIndex) = false;
targetWithoutReference = targetData(targetMask, :);

context = struct();
context.Version = 2;
context.MethodLabel = 'PSO-MDL normalized relative-gap objective';
context.ObjectiveDefinition = ...
    '-((D2-D1)/(D1+epsilon))*F*P';
context.SearchDefinition = ...
    'continuous signed positions; descriptor selected when position>0';
context.ReferenceDefinition = ...
    'target medoid determined in the full retained descriptor space';
context.PAggregation = 'mean_non_medoid';
context.PDefinition = ...
    'mean_t(100*count_b(d_background_b>d_reference_t)/nBackground), t!=medoid';
context.PRange = [0, 100];
context.DistanceEpsilon = distanceEpsilon;
context.D1ValidityRule = 'D1 > epsilon';
context.ReferenceIndex = referenceIndex;
context.ReferenceDrug = referenceDrug;
context.TargetWithoutReference = targetWithoutReference;
context.TargetSquaredDelta = (targetWithoutReference - referenceDrug).^2;
context.BackgroundSquaredDelta = (backgroundData - referenceDrug).^2;
if any(~isfinite(context.TargetSquaredDelta), 'all') ...
        || any(~isfinite(context.BackgroundSquaredDelta), 'all')
    error('build_fitness_context:NonFiniteSquaredDelta', ...
        'Precomputed squared differences overflowed; check the input scale.');
end
context.TargetCountIncludingReference = size(targetData, 1);
context.TargetCountEvaluated = size(targetWithoutReference, 1);
context.BackgroundCount = size(backgroundData, 1);
context.FeatureCount = size(targetData, 2);
context.CreatedAt = char(datetime('now', ...
    'Format', "yyyyMMdd'T'HHmmss"));
end

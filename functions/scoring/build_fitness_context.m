function context = build_fitness_context( ...
    targetData, backgroundData, distanceEpsilon)
%BUILD_FITNESS_CONTEXT Define the reference origin and cache distance components.
%   TARGETDATA = X_R (nR-by-p) and BACKGROUNDDATA = X_B (nB-by-p) are finite,
%   already imputed/scaled matrices. Rows are compounds; columns are descriptors.
%   Target denotes the reference set, not a measured response or protein target.
%
%   The fixed medoid minimizes the sum of Euclidean distances to the nR
%   references in the full retained space. It is not refitted in each subset.
%   CONTEXT stores its row index, scaled coordinates and squared differences
%   for the nR-1 other references and all nB background compounds.
%   DISTANCEEPSILON is the positive constant in R=(D2-D1)/(D1+epsilon).

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

%% Fixed medoid in the full retained descriptor space
% Equivalent to pairwise Euclidean distances; clamp round-off below zero before
% the square root. MATLAB min selects the first row when total distances tie.
rowSquaredNorm = sum(targetData.^2, 2);
squaredDistance = rowSquaredNorm + rowSquaredNorm.' ...
    - 2 .* (targetData * targetData.');
squaredDistance = max(squaredDistance, 0);
[~, referenceIndex] = min(sum(sqrt(squaredDistance), 2));

% Exclude the zero self-distance from D1 and the reference-percentile average P.
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
%% Cached squared coordinate differences
% A subset mask selects columns in these matrices before taking square roots;
% distances are Euclidean, not mean squared or root-mean-square distances.
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

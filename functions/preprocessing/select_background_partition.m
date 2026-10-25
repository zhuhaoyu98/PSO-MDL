function [partitionTable, backgroundTable, holdoutTable, ...
    excludedTable, partitionMetadata, propertySummary] = ...
    select_background_partition( ...
    deduplicationAudit, candidateProperties, ...
    targetProperties, targetNames, propertyNames, options)
%SELECT_BACKGROUND_PARTITION Select background compounds and freeze membership.
%   CANDIDATEPROPERTIES is nCandidate-by-q and TARGETPROPERTIES is nReference-by-q,
%   with q matching properties in the order of PROPERTYNAMES. Rows must align
%   with the identity audit and reference names. Target denotes reference drugs.
%
%   Property-matched mode orders candidates by scaled RMS property distance to
%   each reference and allocates disjoint representatives by round-robin quotas.
%   Random mode samples without replacement using the configured seed.
%   Output partitions are mutually exclusive; background draw order is retained.

arguments
    deduplicationAudit table
    candidateProperties double
    targetProperties double
    targetNames (:, 1) string
    propertyNames (:, 1) string
    options.SelectionMode (1, 1) string = "property_matched"
    options.SizeMode (1, 1) string = "capped_fraction"
    options.FixedCount (1, 1) double ...
        {mustBeInteger, mustBeNonnegative} = 20000
    options.Fraction (1, 1) double ...
        {mustBeGreaterThan(options.Fraction, 0), ...
         mustBeLessThanOrEqual(options.Fraction, 1)} = 0.10
    options.MaximumCount (1, 1) double ...
        {mustBeInteger, mustBePositive} = 20000
    options.MinimumCount (1, 1) double ...
        {mustBeInteger, mustBeNonnegative} = 100
    options.MinimumHoldoutCount (1, 1) double ...
        {mustBeInteger, mustBePositive} = 100
    options.RandomSeed (1, 1) double ...
        {mustBeInteger, mustBeNonnegative} = 20260724
end

requiredVariables = [ ...
    "SourceRow", "DragonName", "StableIdentifier", ...
    "IsEligibleUnique", "DeduplicationStatus"];
assert(all(ismember(requiredVariables, ...
    string(deduplicationAudit.Properties.VariableNames))), ...
    'Deduplication audit is missing fields required for partitioning.');

selectionMode = lower(strtrim(options.SelectionMode));
sizeMode = lower(strtrim(options.SizeMode));
assert(ismember(selectionMode, ["property_matched", "random"]), ...
    'SelectionMode must be "property_matched" or "random".');
assert(ismember(sizeMode, ["fixed", "fraction", "capped_fraction"]), ...
    ['SizeMode must be "fixed", "fraction", or' ...
     '"capped_fraction"。']);

candidateCount = height(deduplicationAudit);
isEligible = logical(deduplicationAudit.IsEligibleUnique);
eligibleRows = find(isEligible);
eligibleCount = numel(eligibleRows);
maximumAllowedBackground = ...
    eligibleCount - options.MinimumHoldoutCount;
assert(maximumAllowedBackground >= options.MinimumCount, ...
    ['Only %d eligible unique candidates remain; cannot provide at least %d background and' ...
     ' %d holdout molecules. Check inputs or explicitly lower these limits.'], ...
    eligibleCount, options.MinimumCount, options.MinimumHoldoutCount);

% capped_fraction uses min(MaximumCount, round(Fraction*nEligible)), followed
% by the configured minimum background and minimum screening-count constraints.
switch sizeMode
    case "fixed"
        backgroundCount = options.FixedCount;
    case "fraction"
        backgroundCount = round(options.Fraction * eligibleCount);
    case "capped_fraction"
        backgroundCount = min(options.MaximumCount, ...
            round(options.Fraction * eligibleCount));
end
backgroundCount = max(backgroundCount, options.MinimumCount);
if backgroundCount > maximumAllowedBackground
    if sizeMode == "fixed"
        error('mdl:FixedBackgroundTooLarge', ...
            ['FixedCount=%d; reserving at least %d holdout molecules allows at most' ...
             ' %d background molecules.'], ...
            backgroundCount, options.MinimumHoldoutCount, ...
            maximumAllowedBackground);
    end
    backgroundCount = maximumAllowedBackground;
end
assert(backgroundCount > 0 && backgroundCount < eligibleCount, ...
    'Background count must be positive and smaller than eligible candidate count.');

nearestTargetIndex = nan(candidateCount, 1);
nearestTargetName = strings(candidateCount, 1);
nearestTargetPropertyDistance = nan(candidateCount, 1);
selectedForTargetIndex = zeros(backgroundCount, 1);
selectedDistance = nan(backgroundCount, 1);
propertyCenter = zeros(1, 0);
propertyScale = zeros(1, 0);
targetQuota = zeros(numel(targetNames), 1);
selectedCountByTarget = zeros(numel(targetNames), 1);

if selectionMode == "property_matched"
    assert(~isempty(propertyNames), ...
        'property_matched requires at least one matching property.');
    assert(size(candidateProperties, 1) == candidateCount, ...
        'Candidate-property row count disagrees with identity manifest.');
    assert(size(candidateProperties, 2) == numel(propertyNames), ...
        'Candidate-property column count disagrees with propertyNames.');
    assert(size(targetProperties, 1) == numel(targetNames) && ...
        size(targetProperties, 2) == numel(propertyNames), ...
        'Reference-property dimensions disagree with reference or property names.');
    assert(~isempty(targetNames), 'property_matched requires reference molecules.');

    eligibleProperties = candidateProperties(eligibleRows, :);
    propertyCenter = median(eligibleProperties, 1, 'omitmissing');
    assert(all(isfinite(propertyCenter)), ...
        'A matching property has no finite values in eligible candidates.');

    for propertyIndex = 1:numel(propertyNames)
        missingCandidate = ...
            ~isfinite(eligibleProperties(:, propertyIndex));
        missingTarget = ...
            ~isfinite(targetProperties(:, propertyIndex));
        eligibleProperties(missingCandidate, propertyIndex) = ...
            propertyCenter(propertyIndex);
        targetProperties(missingTarget, propertyIndex) = ...
            propertyCenter(propertyIndex);
    end

    % Matching uses median absolute deviation about the median, with SD fallback.
    % This is distinct from Step 01's mean-AD descriptor scaling for the objective.
    propertyScale = median(abs( ...
        eligibleProperties - propertyCenter), 1);
    for propertyIndex = 1:numel(propertyNames)
        if ~isfinite(propertyScale(propertyIndex)) || ...
                propertyScale(propertyIndex) <= sqrt(eps)
            fallbackScale = std( ...
                eligibleProperties(:, propertyIndex), 0, 1);
            assert(isfinite(fallbackScale) && ...
                fallbackScale > sqrt(eps), ...
                'Matching property %s has no usable variation in eligible candidates.', ...
                propertyNames(propertyIndex));
            propertyScale(propertyIndex) = fallbackScale;
        end
    end

    candidateStandardized = ...
        (eligibleProperties - propertyCenter) ./ propertyScale;
    targetStandardized = ...
        (targetProperties - propertyCenter) ./ propertyScale;
    targetCount = size(targetStandardized, 1);
    % RMS across q scaled properties: this selects background records only.
    % It is not D1/D2, which use Euclidean distances in the optimized subspace.
    distanceToTargets = zeros( ...
        eligibleCount, targetCount, 'single');
    for targetIndex = 1:targetCount
        distanceToTargets(:, targetIndex) = single(sqrt(mean( ...
            (candidateStandardized - ...
            targetStandardized(targetIndex, :)) .^ 2, 2)));
    end

    [nearestDistanceEligible, nearestIndexEligible] = ...
        min(distanceToTargets, [], 2);
    nearestTargetIndex(eligibleRows) = nearestIndexEligible;
    nearestTargetName(eligibleRows) = ...
        targetNames(nearestIndexEligible);
    nearestTargetPropertyDistance(eligibleRows) = ...
        double(nearestDistanceEligible);

    sortedPositionByTarget = zeros( ...
        eligibleCount, targetCount, 'uint32');
    for targetIndex = 1:targetCount

        [~, currentOrder] = sort( ...
            distanceToTargets(:, targetIndex), 'ascend');
        sortedPositionByTarget(:, targetIndex) = ...
            uint32(currentOrder);
    end

    % Balance reference-specific quotas and never assign one compound twice.
    baseQuota = floor(backgroundCount / targetCount);
    targetQuota = repmat(baseQuota, targetCount, 1);
    targetQuota(1:mod(backgroundCount, targetCount)) = ...
        targetQuota(1:mod(backgroundCount, targetCount)) + 1;

    selectedRowsInDrawOrder = zeros(backgroundCount, 1);
    selectedEligiblePosition = false(eligibleCount, 1);
    nextPositionByTarget = ones(targetCount, 1);
    selectedCount = 0;
    while selectedCount < backgroundCount
        madeSelection = false;
        for targetIndex = 1:targetCount
            if selectedCountByTarget(targetIndex) >= ...
                    targetQuota(targetIndex)
                continue;
            end
            while nextPositionByTarget(targetIndex) <= eligibleCount
                eligiblePosition = double(sortedPositionByTarget( ...
                    nextPositionByTarget(targetIndex), targetIndex));
                nextPositionByTarget(targetIndex) = ...
                    nextPositionByTarget(targetIndex) + 1;
                if ~selectedEligiblePosition(eligiblePosition)
                    selectedCount = selectedCount + 1;
                    sourceIndex = eligibleRows(eligiblePosition);
                    selectedRowsInDrawOrder(selectedCount) = sourceIndex;
                    selectedForTargetIndex(selectedCount) = targetIndex;
                    selectedDistance(selectedCount) = ...
                        distanceToTargets(eligiblePosition, targetIndex);
                    selectedEligiblePosition(eligiblePosition) = true;
                    selectedCountByTarget(targetIndex) = ...
                        selectedCountByTarget(targetIndex) + 1;
                    madeSelection = true;
                    break;
                end
            end
        end
        assert(madeSelection, ...
            'Cannot select %d property-matched background molecules under the uniqueness constraint.', ...
            backgroundCount);
    end
else

    previousRandomState = rng;
    cleanupRandomState = onCleanup(@() rng(previousRandomState));
    rng(options.RandomSeed, 'twister');
    randomPosition = randperm(eligibleCount, backgroundCount);
    selectedRowsInDrawOrder = eligibleRows(randomPosition);
    clear cleanupRandomState
end

% Only eligible, non-background records enter the subsequent candidate ranking.
isBackground = false(candidateCount, 1);
isBackground(selectedRowsInDrawOrder) = true;
isHoldout = isEligible & ~isBackground;
isExcluded = ~isEligible;
assert(~any(isBackground & isHoldout) && ...
    ~any(isBackground & isExcluded) && ...
    ~any(isHoldout & isExcluded) && ...
    all(isBackground | isHoldout | isExcluded), ...
    'Background, holdout, and exclusions are not mutually exclusive and exhaustive.');

partitionGroup = "Excluded_" + ...
    string(deduplicationAudit.DeduplicationStatus);
partitionGroup(isHoldout) = "ScreeningHoldout";
if selectionMode == "property_matched"
    partitionGroup(isBackground) = ...
        "BackgroundTrainingPropertyMatched";
else
    partitionGroup(isBackground) = "BackgroundTrainingRandom";
end

partitionTable = deduplicationAudit;
partitionTable.NearestTargetIndex = nearestTargetIndex;
partitionTable.NearestTargetName = nearestTargetName;
partitionTable.NearestTargetPropertyDistance = ...
    nearestTargetPropertyDistance;
partitionTable.PartitionGroup = partitionGroup;

backgroundTable = partitionTable(selectedRowsInDrawOrder, :);
backgroundTable.DrawOrder = (1:backgroundCount).';
backgroundTable.MatchedTargetIndex = selectedForTargetIndex;
backgroundTable.MatchedTargetName = strings(backgroundCount, 1);
hasMatchedTarget = selectedForTargetIndex > 0;
backgroundTable.MatchedTargetName(hasMatchedTarget) = ...
    targetNames(selectedForTargetIndex(hasMatchedTarget));
backgroundTable.MatchedTargetPropertyDistance = selectedDistance;
backgroundTable = movevars(backgroundTable, ...
    {'DrawOrder', 'MatchedTargetIndex', ...
     'MatchedTargetName', 'MatchedTargetPropertyDistance'}, ...
    'Before', 1);
holdoutTable = partitionTable(isHoldout, :);
excludedTable = partitionTable(isExcluded, :);

propertySummary = localBuildPropertySummary( ...
    candidateProperties, targetProperties, ...
    selectedRowsInDrawOrder, isEligible, propertyNames);

partitionMetadata = struct();
partitionMetadata.Version = 3;
partitionMetadata.Method = char(selectionMode);
partitionMetadata.BackgroundSizeMode = char(sizeMode);
partitionMetadata.BackgroundFraction = options.Fraction;
partitionMetadata.BackgroundFixedCountParameter = options.FixedCount;
partitionMetadata.BackgroundMaximumCount = options.MaximumCount;
partitionMetadata.BackgroundMinimumCount = options.MinimumCount;
partitionMetadata.MinimumHoldoutCount = options.MinimumHoldoutCount;
partitionMetadata.RandomSeed = options.RandomSeed;
partitionMetadata.CandidateRows = candidateCount;
partitionMetadata.EligibleUniqueCount = eligibleCount;
partitionMetadata.BackgroundCount = nnz(isBackground);
partitionMetadata.ScreeningHoldoutCount = nnz(isHoldout);
partitionMetadata.ExcludedCount = nnz(isExcluded);
partitionMetadata.MissingIdentifierExcludedCount = nnz( ...
    deduplicationAudit.IsMissingIdentifier);
partitionMetadata.ReferenceOverlapExcludedCount = nnz( ...
    deduplicationAudit.IsReferenceOverlap);
partitionMetadata.DuplicateNonRepresentativeExcludedCount = nnz( ...
    deduplicationAudit.DeduplicationStatus == ...
    "DuplicateNonRepresentativeExcluded");
partitionMetadata.BackgroundRowsInDrawOrder = ...
    selectedRowsInDrawOrder;
partitionMetadata.IsBackground = isBackground;
partitionMetadata.IsScreeningHoldout = isHoldout;
partitionMetadata.IsExcluded = isExcluded;
partitionMetadata.PropertyNames = propertyNames;
partitionMetadata.PropertyCenter = propertyCenter;
partitionMetadata.PropertyScale = propertyScale;
partitionMetadata.TargetNames = targetNames;
partitionMetadata.TargetQuota = targetQuota;
partitionMetadata.SelectedCountByTarget = selectedCountByTarget;
partitionMetadata.CreatedAt = char(datetime( ...
    'now', 'TimeZone', 'local', ...
    'Format', 'yyyy-MM-dd''T''HH:mm:ssXXX'));
end

function summaryTable = localBuildPropertySummary( ...
    candidateProperties, targetProperties, ...
    selectedRows, isEligible, propertyNames)

groupName = strings(0, 1);
propertyNameColumn = strings(0, 1);
sampleCount = zeros(0, 1);
meanValue = zeros(0, 1);
medianValue = zeros(0, 1);
standardDeviation = zeros(0, 1);
minimumValue = zeros(0, 1);
maximumValue = zeros(0, 1);

if isempty(candidateProperties) || isempty(propertyNames)
    summaryTable = table(groupName, propertyNameColumn, ...
        sampleCount, meanValue, medianValue, ...
        standardDeviation, minimumValue, maximumValue, ...
        'VariableNames', {'Group', 'Property', 'N', ...
        'Mean', 'Median', 'SD', 'Minimum', 'Maximum'});
    return;
end

groupLabels = ["Reference", "EligibleUniqueCandidatePool", ...
    "SelectedBackground"];
groupMatrices = {targetProperties, ...
    candidateProperties(isEligible, :), ...
    candidateProperties(selectedRows, :)};
for groupIndex = 1:numel(groupLabels)
    currentMatrix = groupMatrices{groupIndex};
    for propertyIndex = 1:numel(propertyNames)
        finiteValue = currentMatrix( ...
            isfinite(currentMatrix(:, propertyIndex)), propertyIndex);
        groupName(end + 1, 1) = groupLabels(groupIndex);
        propertyNameColumn(end + 1, 1) = ...
            propertyNames(propertyIndex);
        sampleCount(end + 1, 1) = numel(finiteValue);
        meanValue(end + 1, 1) = mean(finiteValue);
        medianValue(end + 1, 1) = median(finiteValue);
        standardDeviation(end + 1, 1) = std(finiteValue);
        minimumValue(end + 1, 1) = min(finiteValue);
        maximumValue(end + 1, 1) = max(finiteValue);
    end
end

summaryTable = table(groupName, propertyNameColumn, ...
    sampleCount, meanValue, medianValue, ...
    standardDeviation, minimumValue, maximumValue, ...
    'VariableNames', {'Group', 'Property', 'N', ...
    'Mean', 'Median', 'SD', 'Minimum', 'Maximum'});
end

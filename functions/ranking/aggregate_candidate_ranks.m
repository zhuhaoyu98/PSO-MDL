function result = aggregate_candidate_ranks(rankByCandidate, sourceRows)
%AGGREGATE_CANDIDATE_RANKS Summarize equally weighted single-model ranks.
%   RANKBYCANDIDATE is N-by-M with each column a permutation of 1:N.
%   SOURCEROWS is an N-by-1 unique source-row vector used only to resolve
%   equal mean ranks deterministically. No fitness weighting is applied.
%   Result fields contain ranks in input row order, the final ordering, and
%   cross-model uncertainty. P10/P90 use the original nearest order statistic,
%   not an interpolated percentile. No candidate distances are rescaled.
validateattributes(rankByCandidate, {'numeric'}, {'2d', 'real', 'finite', ...
    'integer', 'positive', 'nonempty'});
expectedHoldoutRows = size(rankByCandidate, 1);
modelCount = size(rankByCandidate, 2);
validateattributes(sourceRows, {'numeric'}, {'vector', 'real', 'finite', ...
    'integer', 'positive', 'numel', expectedHoldoutRows});
sourceRows = sourceRows(:);
assert(numel(unique(sourceRows)) == expectedHoldoutRows, ...
    'Source rows must be unique.');
for model = 1:modelCount
    assert(isequal(sort(double(rankByCandidate(:, model))), ...
        (1:expectedHoldoutRows).'), 'Each model must contain a complete rank permutation.');
end
rankDouble = double(rankByCandidate);
meanRank = mean(rankDouble, 2);
medianRank = median(rankDouble, 2);
rankSD = std(rankDouble, 0, 2);
bestRank = min(rankDouble, [], 2);
worstRank = max(rankDouble, [], 2);
sortedRanks = sort(rankDouble, 2, 'ascend');
rankP10 = sortedRanks(:, max(1, ceil(0.10 * modelCount)));
rankP90 = sortedRanks(:, max(1, ceil(0.90 * modelCount)));
clear sortedRanks
[~, meanRankOrder] = sortrows( ...
    [meanRank, double(sourceRows)], [1, 2]);
finalRank = zeros(expectedHoldoutRows, 1);
finalRank(meanRankOrder) = (1:expectedHoldoutRows).';
result = struct('DoubleRanks', rankDouble, 'MeanRank', meanRank, ...
    'MedianRank', medianRank, 'RankSD', rankSD, 'BestRank', bestRank, ...
    'WorstRank', worstRank, 'RankP10', rankP10, 'RankP90', rankP90, ...
    'Order', meanRankOrder, 'FinalRank', finalRank);
end

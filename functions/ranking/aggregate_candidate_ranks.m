function result = aggregate_candidate_ranks(rankByCandidate, sourceRows)
%AGGREGATE_CANDIDATE_RANKS Combine M models by mean ordinal rank.
%   RANKBYCANDIDATE is N-by-M: each column is a permutation of 1:N (1 is best).
%   SOURCEROWS is an N-by-1 unique record locator used only for exact rank ties.
%   meanRank_i = sum_m(rank_im)/M; neither J nor descriptor count weights models.
%
%   RESULT stores summaries in input row order, the final sorting order and
%   final ranks. P10/P90 are nearest-rank order statistics at ceil(0.10*M) and
%   ceil(0.90*M), not interpolated percentiles or confidence-interval limits.
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
% Aggregate across model columns. Distances from different subspaces are not
% averaged because their descriptor counts and distance scales may differ.
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
% A source-row tie-break makes the final rank permutation reproducible.
[~, meanRankOrder] = sortrows( ...
    [meanRank, double(sourceRows)], [1, 2]);
finalRank = zeros(expectedHoldoutRows, 1);
finalRank(meanRankOrder) = (1:expectedHoldoutRows).';
result = struct('DoubleRanks', rankDouble, 'MeanRank', meanRank, ...
    'MedianRank', medianRank, 'RankSD', rankSD, 'BestRank', bestRank, ...
    'WorstRank', worstRank, 'RankP10', rankP10, 'RankP90', rankP90, ...
    'Order', meanRankOrder, 'FinalRank', finalRank);
end

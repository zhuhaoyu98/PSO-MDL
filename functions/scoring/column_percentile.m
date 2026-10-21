function percentileValue = column_percentile(values, percentile)
%COLUMN_PERCENTILE Original linearly interpolated percentile along rows.
%   VALUES is observations-by-variables. PERCENTILE is a scalar in [0,100].
%   Used for optimization/overlap summaries, not candidate P10/P90 ranks.
%   Position 1+(n-1)*percentile/100 preserves the research reduction order.
values = sort(double(values), 1, 'ascend');
rowCount = size(values, 1);
if rowCount == 1
    percentileValue = values;
    return
end
position = 1 + (rowCount - 1) .* percentile ./ 100;
lowerIndex = floor(position);
upperIndex = ceil(position);
fraction = position - lowerIndex;
percentileValue = values(lowerIndex, :) .* (1 - fraction) ...
    + values(upperIndex, :) .* fraction;
end

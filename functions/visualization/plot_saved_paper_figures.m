function outputs = plot_saved_paper_figures(dataDirectory, outputDirectory, options)
%PLOT_SAVED_PAPER_FIGURES Step 04 helper for the frozen paper-result bundle.
%   OUTPUTS = PLOT_SAVED_PAPER_FIGURES() uses data/paper_results and writes
%   results/paper_figures. No PSO, candidate scoring, or ensemble resampling
%   is performed. Data are the real mean-P study summaries, not toy results.
%   Optional name-value arguments: Show=false, Overwrite=false, Dpi=600.
%   The drawing operations/styles follow the formal Step 04; only input
%   loading and output orchestration are independent of its large MAT files.
arguments
    dataDirectory (1,1) string = ""
    outputDirectory (1,1) string = ""
    options.Show (1,1) logical = false
    options.Overwrite (1,1) logical = false
    options.Dpi (1,1) double {mustBePositive,mustBeFinite} = 600
end
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
if strlength(dataDirectory)==0, dataDirectory=fullfile(root,'data','paper_results'); end
if strlength(outputDirectory)==0, outputDirectory=fullfile(root,'results','paper_figures'); end
dataPath=string(java.io.File(char(dataDirectory)).getCanonicalPath());
outPath=string(java.io.File(char(outputDirectory)).getCanonicalPath());
assert(~strcmpi(dataPath,outPath) && ~startsWith(lower(outPath),lower(dataPath)+filesep), ...
    'pso_mdl:UnsafeOutput','Figures must be written outside the input data directory.');
data=load_paper_results(dataDirectory);
meta=data.Metadata;
if ~isfolder(outputDirectory), mkdir(outputDirectory); end
figureDir=char(outputDirectory);
overwriteExistingFigures=options.Overwrite;
showFiguresOnScreen=options.Show;
figureResolutionDpi=options.Dpi;
figureFontName=meta.FigureFontName;
baseFontSizePt=meta.FigureFontSize;
fitnessDisplayScale=meta.FitnessDisplayScale;
regularFigureSizeCm=meta.RegularSizeCm;
tallFigureSizeCm=meta.TallSizeCm;
modelCount=meta.ModelCount;
generationCount=meta.GenerationCount;
step02FigureTable=data.History;
step03JaccardTable=data.TopKJaccard;
step03TopCandidateTable=data.Top20;
candidateFigureCount=height(data.Top20);
featureSelectionFrequency=data.Models.featureSelectionFrequency;
modelSummaryTable=data.Models.modelSummaryTable;
finalF=double(modelSummaryTable.SelectedFeatureCount(:));
finalRelativeGap=double(modelSummaryTable.RelativeGap(:));
rhoFvsR=data.Correlations.SpearmanRho(1);
ensemble=data.EnsembleSummary;
ensembleSubsetSizes=ensemble.EnsembleSize.';
ensembleSpearmanQ1=ensemble.SpearmanQ1;
ensembleSpearmanMedian=ensemble.MedianSpearman;
ensembleSpearmanQ3=ensemble.SpearmanQ3;
ensembleTopKQ1=ensemble.TopKJaccardQ1;
ensembleTopKMedian=ensemble.MedianTopKJaccard;
ensembleTopKQ3=ensemble.TopKJaccardQ3;
actualEnsembleTopK=meta.EnsembleTopK;

fig3aBase = fullfile(figureDir, 'FIG3A');
fig3bBase = fullfile(figureDir, 'FIG3B');
fig3cBase = fullfile(figureDir, 'FIG3C');
fig4aBase = fullfile(figureDir, 'FIG4A');
fig4bBase = fullfile(figureDir, 'FIG4B');
fig4dBase = fullfile(figureDir, 'FIG4D');
fig4eBase = fullfile(figureDir, 'FIG4E');
fig4cBase = fullfile(figureDir, 'FIG4C');

allFigureBases = {fig3aBase, fig3bBase, fig3cBase, fig4aBase, ...
    fig4bBase, fig4dBase, fig4eBase, fig4cBase};
for baseIndex = 1:numel(allFigureBases)
    localPrepareFigureOutput(allFigureBases{baseIndex}, overwriteExistingFigures);
end

visibilityOptions = {'off', 'on'};
figureVisibility = visibilityOptions{1 + double(showFiguresOnScreen)};

blue = [0, 114, 178] ./ 255;
lightBlue = [86, 180, 233] ./ 255;
vermillion = [213, 94, 0] ./ 255;
lightOrange = [0.96, 0.82, 0.64];
lightGray = [0.72, 0.72, 0.72];
darkGray = [0.30, 0.30, 0.30];

generation = double(step02FigureTable.PSOIteration);
fitnessMedian = double(step02FigureTable.MedianBestFitness);
fitnessQ1 = double(step02FigureTable.BestFitnessQ1);
fitnessQ3 = double(step02FigureTable.BestFitnessQ3);
descriptorMedian = double(step02FigureTable.MedianSelectedDescriptorCount);
descriptorQ1 = double(step02FigureTable.SelectedDescriptorCountQ1);
descriptorQ3 = double(step02FigureTable.SelectedDescriptorCountQ3);

topK = double(step03JaccardTable.TopK);
jaccardMedian = double(step03JaccardTable.MedianJaccard);
jaccardQ1 = double(step03JaccardTable.JaccardQ1);
jaccardQ3 = double(step03JaccardTable.JaccardQ3);

fig3a = figure('Visible', figureVisibility, 'Color', 'w', ...
    'Units', 'centimeters', 'Position', [2, 2, regularFigureSizeCm(1), regularFigureSizeCm(2)]);
ax = axes(fig3a);
hold(ax, 'on');
fill(ax, [generation; flipud(generation)], ...
    [fitnessQ1; flipud(fitnessQ3)] ./ fitnessDisplayScale, ...
    lightBlue, 'EdgeColor', 'none', 'FaceAlpha', 0.32, 'HandleVisibility', 'off');
plot(ax, generation, fitnessMedian ./ fitnessDisplayScale, ...
    'Color', blue, 'LineWidth', 1.55, 'HandleVisibility', 'off');
xlabel(ax, 'PSO generation');
ylabel(ax, 'Best-so-far objective (\times10^6)');
localStyleAxes(ax, figureFontName, baseFontSizePt);
ax.YAxis.Exponent = 0;
xlim(ax, [1, generationCount]);
hold(ax, 'off');
drawnow;
localExportFigure(fig3a, fig3aBase, figureResolutionDpi);

fig3b = figure('Visible', figureVisibility, 'Color', 'w', ...
    'Units', 'centimeters', 'Position', [2, 2, regularFigureSizeCm(1), regularFigureSizeCm(2)]);
ax = axes(fig3b);
hold(ax, 'on');
fill(ax, [generation; flipud(generation)], ...
    [descriptorQ1; flipud(descriptorQ3)], ...
    lightOrange, 'EdgeColor', 'none', 'FaceAlpha', 0.50, 'HandleVisibility', 'off');
plot(ax, generation, descriptorMedian, ...
    'Color', vermillion, 'LineWidth', 1.55, 'HandleVisibility', 'off');
xlabel(ax, 'PSO generation');
ylabel(ax, 'Number of selected descriptors');
localStyleAxes(ax, figureFontName, baseFontSizePt);
xlim(ax, [1, generationCount]);
ylim(ax, localRoundedDataRange([descriptorQ1; descriptorQ3], 100));
hold(ax, 'off');
drawnow;
localExportFigure(fig3b, fig3bBase, figureResolutionDpi);

fig3c = figure('Visible', figureVisibility, 'Color', 'w', ...
    'Units', 'centimeters', 'Position', [2, 2, regularFigureSizeCm(1), regularFigureSizeCm(2)]);
ax = axes(fig3c);
hold(ax, 'on');
scatter(ax, finalF, finalRelativeGap, 26, blue, 'filled', ...
    'MarkerFaceAlpha', 0.78, 'MarkerEdgeColor', 'w', 'LineWidth', 0.45);
xlabel(ax, 'Number of selected descriptors (F)');
ylabel(ax, 'Relative separation (R)');
localStyleAxes(ax, figureFontName, baseFontSizePt);
set(ax, 'XGrid', 'off', 'YGrid', 'off');
text(ax, 0.97, 0.96, sprintf('\\rho = %s%.3f\nn = %d', char(8722), abs(rhoFvsR), modelCount), ...
    'Units', 'normalized', 'HorizontalAlignment', 'right', ...
    'VerticalAlignment', 'top', 'FontName', figureFontName, ...
    'FontSize', baseFontSizePt, 'Color', darkGray, 'Interpreter', 'tex');
hold(ax, 'off');
drawnow;
localExportFigure(fig3c, fig3cBase, figureResolutionDpi);

fig4a = figure('Visible', figureVisibility, 'Color', 'w', ...
    'Units', 'centimeters', 'Position', [2, 2, regularFigureSizeCm(1), regularFigureSizeCm(2)]);
ax = axes(fig4a);
hold(ax, 'on');
sortedFeatureFrequency = sort(featureSelectionFrequency, 'descend');
frequencyOrder = (1:numel(sortedFeatureFrequency)).';
plot(ax, frequencyOrder, sortedFeatureFrequency, ...
    'Color', darkGray, 'LineWidth', 1.25);
xlabel(ax, 'Descriptor frequency rank');
ylabel(ax, 'Selection frequency');
localStyleAxes(ax, figureFontName, baseFontSizePt);
set(ax, 'YLim', [0, 1], 'YGrid', 'off');
xlim(ax, [1, numel(sortedFeatureFrequency)]);
hold(ax, 'off');
drawnow;
localExportFigure(fig4a, fig4aBase, figureResolutionDpi);

fig4b = figure('Visible', figureVisibility, 'Color', 'w', ...
    'Units', 'centimeters', 'Position', [2, 2, regularFigureSizeCm(1), regularFigureSizeCm(2)]);
ax = axes(fig4b);
hold(ax, 'on');
fill(ax, [topK; flipud(topK)], [jaccardQ1; flipud(jaccardQ3)], ...
    lightBlue, 'EdgeColor', 'none', 'FaceAlpha', 0.28, 'HandleVisibility', 'off');
plot(ax, topK, jaccardMedian, '-o', ...
    'Color', blue, 'MarkerFaceColor', blue, 'MarkerEdgeColor', 'w', ...
    'MarkerSize', 4.2, 'LineWidth', 1.35);
set(ax, 'XScale', 'log', 'YLim', [0, 1]);
if isscalar(topK)
    xlim(ax, [0.8 .* topK, 1.25 .* topK]);
else
    xlim(ax, [min(topK) ./ 1.15, max(topK) .* 1.15]);
end
xticks(ax, topK);
xticklabels(ax, localCompactIntegerLabels(topK));
xlabel(ax, 'Top-K threshold');
ylabel(ax, 'Pairwise Jaccard index');
localStyleAxes(ax, figureFontName, baseFontSizePt);
hold(ax, 'off');
drawnow;
localExportFigure(fig4b, fig4bBase, figureResolutionDpi);

fig4d = figure('Visible', figureVisibility, 'Color', 'w', ...
    'Units', 'centimeters', 'Position', [2, 2, regularFigureSizeCm(1), regularFigureSizeCm(2)]);
ax = axes(fig4d);
hold(ax, 'on');
fill(ax, [ensembleSubsetSizes(:); flipud(ensembleSubsetSizes(:))], ...
    [ensembleSpearmanQ1; flipud(ensembleSpearmanQ3)], ...
    lightBlue, 'EdgeColor', 'none', 'FaceAlpha', 0.28, 'HandleVisibility', 'off');
plot(ax, ensembleSubsetSizes, ensembleSpearmanMedian, '-o', ...
    'Color', blue, 'MarkerFaceColor', blue, 'MarkerEdgeColor', 'w', ...
    'MarkerSize', 4.2, 'LineWidth', 1.35);
set(ax, 'XScale', 'log');
xticks(ax, ensembleSubsetSizes);
xticklabels(ax, compose('%d', ensembleSubsetSizes));
xlim(ax, [min(ensembleSubsetSizes) ./ 1.15, max(ensembleSubsetSizes) .* 1.15]);
ylim(ax, [0.994, 1.000]);
yticks(ax, [0.994, 0.996, 0.998, 1.000]);
xlabel(ax, 'Ensemble size, M');
ylabel(ax, 'Spearman \rho between ensembles');
localStyleAxes(ax, figureFontName, baseFontSizePt);
hold(ax, 'off');
drawnow;
localExportFigure(fig4d, fig4dBase, figureResolutionDpi);

fig4e = figure('Visible', figureVisibility, 'Color', 'w', ...
    'Units', 'centimeters', 'Position', [2, 2, regularFigureSizeCm(1), regularFigureSizeCm(2)]);
ax = axes(fig4e);
hold(ax, 'on');
fill(ax, [ensembleSubsetSizes(:); flipud(ensembleSubsetSizes(:))], ...
    [ensembleTopKQ1; flipud(ensembleTopKQ3)], ...
    lightOrange, 'EdgeColor', 'none', 'FaceAlpha', 0.45, 'HandleVisibility', 'off');
plot(ax, ensembleSubsetSizes, ensembleTopKMedian, '-o', ...
    'Color', vermillion, 'MarkerFaceColor', vermillion, 'MarkerEdgeColor', 'w', ...
    'MarkerSize', 4.2, 'LineWidth', 1.35);
set(ax, 'XScale', 'log');
xticks(ax, ensembleSubsetSizes);
xticklabels(ax, compose('%d', ensembleSubsetSizes));
xlim(ax, [min(ensembleSubsetSizes) ./ 1.15, max(ensembleSubsetSizes) .* 1.15]);
ylim(ax, localUpperBoundedRange(min(ensembleTopKQ1), 1, 0.035));
xlabel(ax, 'Ensemble size, M');
ylabel(ax, sprintf('Top-%d Jaccard between ensembles', actualEnsembleTopK));
localStyleAxes(ax, figureFontName, baseFontSizePt);
hold(ax, 'off');
drawnow;
localExportFigure(fig4e, fig4eBase, figureResolutionDpi);

fig4c = figure('Visible', figureVisibility, 'Color', 'w', ...
    'Units', 'centimeters', 'Position', [2, 2, tallFigureSizeCm(1), tallFigureSizeCm(2)]);
ax = axes(fig4c);
hold(ax, 'on');
candidateY = (1:candidateFigureCount).';
rangeLegend = plot(ax, nan, nan, '-', 'Color', lightGray, 'LineWidth', 0.9);
percentileLegend = plot(ax, nan, nan, '-', 'Color', blue, 'LineWidth', 2.0);
meanLegend = scatter(ax, nan, nan, 24, vermillion, 'filled', ...
    'MarkerEdgeColor', 'w', 'LineWidth', 0.45);
for candidateIndex = 1:candidateFigureCount
    plot(ax, ...
        [step03TopCandidateTable.BestRankAcrossSeeds(candidateIndex), ...
         step03TopCandidateTable.WorstRankAcrossSeeds(candidateIndex)], ...
        [candidateY(candidateIndex), candidateY(candidateIndex)], ...
        '-', 'Color', lightGray, 'LineWidth', 0.75, 'HandleVisibility', 'off');
    plot(ax, ...
        [step03TopCandidateTable.RankP10AcrossSeeds(candidateIndex), ...
         step03TopCandidateTable.RankP90AcrossSeeds(candidateIndex)], ...
        [candidateY(candidateIndex), candidateY(candidateIndex)], ...
        '-', 'Color', blue, 'LineWidth', 2.0, 'HandleVisibility', 'off');
end
scatter(ax, step03TopCandidateTable.MeanRankAcrossSeeds, candidateY, 24, ...
    vermillion, 'filled', 'MarkerEdgeColor', 'w', 'LineWidth', 0.45);
set(ax, 'XScale', 'log', 'YDir', 'reverse', ...
    'YLim', [0.3, candidateFigureCount + 0.7]);
maximumRank = max(step03TopCandidateTable.WorstRankAcrossSeeds);
rankUpperLimit = max(10, 10 .^ ceil(log10(maximumRank)));
xlim(ax, [0.8, rankUpperLimit]);
yticks(ax, candidateY);
yticklabels(ax, compose('%d', step03TopCandidateTable.Rank));
xlabel(ax, 'Rank in individual PSO model');
ylabel(ax, 'Final ensemble rank');
localStyleAxes(ax, figureFontName, baseFontSizePt);
legend(ax, [meanLegend, percentileLegend, rangeLegend], ...
    {'Mean rank', '10th-90th percentile', 'Min-max'}, ...
    'Location', 'southoutside', 'Box', 'off', 'FontSize', 7);
hold(ax, 'off');
drawnow;
localExportFigure(fig4c, fig4cBase, figureResolutionDpi);


ids=["FIG3A";"FIG3B";"FIG3C";"FIG4A";"FIG4B";"FIG4C";"FIG4D";"FIG4E"];
outputs.Directory=outputDirectory;
outputs.Files=table(ids,fullfile(outputDirectory,ids+".png"), ...
    fullfile(outputDirectory,ids+".pdf"),'VariableNames',{'Figure','PNG','PDF'});
outputs.Source="saved mean-P model and aggregate statistics";
outputs.Handles=[fig3a fig3b fig3c fig4a fig4b fig4c fig4d fig4e];
if ~options.Show, close(outputs.Handles); outputs.Handles=[]; end
fprintf('Replotted FIG3A-C and FIG4A-E from saved summaries; no study computations rerun.\n');
end

function labels = localCompactIntegerLabels(values)
values = double(values(:));
labels = strings(size(values));
for valueIndex = 1:numel(values)
    if values(valueIndex) >= 1000 ...
            && mod(values(valueIndex), 1000) == 0
        labels(valueIndex) = sprintf('%gk', values(valueIndex) ./ 1000);
    else
        labels(valueIndex) = sprintf('%g', values(valueIndex));
    end
end
end

function yRange = localUpperBoundedRange(observedMinimum, upperBound, padding)
observedMinimum = double(observedMinimum);
lowerBound = max(-1, observedMinimum - padding);
if upperBound == 1
    lowerBound = max(0, lowerBound);
end
if lowerBound >= upperBound
    lowerBound = upperBound - max(padding, 0.01);
end
yRange = [lowerBound, upperBound];
end

function yRange = localRoundedDataRange(values, stepSize)
values = double(values(:));
values = values(isfinite(values));
assert(~isempty(values), 'Cannot set an axis range from empty data.');
lowerBound = floor(min(values) ./ stepSize) .* stepSize;
upperBound = ceil(max(values) ./ stepSize) .* stepSize;
if lowerBound == upperBound
    lowerBound = lowerBound - stepSize;
    upperBound = upperBound + stepSize;
end
yRange = [lowerBound, upperBound];
end

function localStyleAxes(axisHandle, fontName, fontSize)
set(axisHandle, ...
    'FontName', fontName, ...
    'FontSize', fontSize, ...
    'LineWidth', 0.8, ...
    'TickDir', 'out', ...
    'TickLength', [0.015, 0.015], ...
    'Box', 'off', ...
    'Layer', 'top', ...
    'XGrid', 'off', ...
    'YGrid', 'off', ...
    'XMinorGrid', 'off', ...
    'YMinorGrid', 'off');
end

function localPrepareFigureOutput(outputBase, overwriteExisting)
pngFile = [outputBase '.png'];
pdfFile = [outputBase '.pdf'];
if ~overwriteExisting
    assert(~isfile(pngFile) && ~isfile(pdfFile), ...
        ['Output already exists and overwrite is disabled: %s' newline ...
         'Use a new output directory, or explicitly enable Overwrite.'], outputBase);
end
end

function localExportFigure(figureHandle, outputBase, resolutionDpi)
exportgraphics(figureHandle, [outputBase '.png'], ...
    'Resolution', resolutionDpi);
exportgraphics(figureHandle, [outputBase '.pdf'], ...
    'ContentType', 'vector');
end



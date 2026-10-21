function data = load_paper_results(dataDirectory)
%LOAD_PAPER_RESULTS Validate and read the curated, mean-P study output bundle.
%   DATA = LOAD_PAPER_RESULTS() uses data/paper_results relative to this repo.
%   Only saved model outputs and aggregate statistics are loaded. No raw
%   descriptor table, chemical identity, or full candidate-rank matrix is used.
%   File digests detect accidental changes, not authenticity or legal rights.
arguments
    dataDirectory (1,1) string = ""
end
if strlength(dataDirectory)==0
    root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
    dataDirectory=fullfile(root,'data','paper_results');
end
assert(isfolder(dataDirectory),'pso_mdl:MissingPaperData', ...
    'Curated paper-result directory was not found.');
manifest=readtable(fullfile(dataDirectory,'provenance.tsv'), ...
    'FileType','text','Delimiter','\t','TextType','string');
for k=1:height(manifest)
    name=manifest.File(k);
    assert(~contains(name,["/","\",".."]),'pso_mdl:InvalidManifest','Expected a file basename.');
    file=fullfile(dataDirectory,name);
    assert(isfile(file),'pso_mdl:MissingPaperData','Missing paper artifact: %s',name);
    assert(strcmpi(compute_file_sha256(file),manifest.SHA256(k)), ...
        'pso_mdl:PaperDataHashMismatch','Paper artifact digest mismatch: %s',name);
end
data.Metadata=jsondecode(fileread(fullfile(dataDirectory,'metadata.json')));
assert(strcmp(data.Metadata.PAggregation,'mean_non_medoid') && ...
    isequal(data.Metadata.PRange(:),[0;100]), ...
    'pso_mdl:WrongPDefinition','Expected mean non-medoid P on the 0-100 scale.');
data.Models=load(fullfile(dataDirectory,'final_models.mat'));
data.Trajectories=load(fullfile(dataDirectory,'optimization_trajectories.mat'));
tableFiles={'History','optimization_history.tsv'; ...
    'ModelSummary','model_summary.tsv';'Correlations','model_correlations.tsv'; ...
    'FeatureFrequency','feature_frequency.tsv';'MaskJaccard','mask_jaccard.tsv'; ...
    'TopKJaccard','jaccard_at_k_summary.tsv'; ...
    'EnsembleSummary','independent_ensemble_summary.tsv'; ...
    'EnsembleRepeats','independent_ensemble_repeats.tsv'; ...
    'Top20','top20_rank_uncertainty.tsv'};
for k=1:size(tableFiles,1)
    data.(tableFiles{k,1})=readtable(fullfile(dataDirectory,tableFiles{k,2}), ...
        'FileType','text','Delimiter','\t','VariableNamingRule','preserve');
end
m=data.Models; h=data.Trajectories; meta=data.Metadata;
assert(isequal(size(m.bestMasks),[meta.ModelCount meta.FeatureCount]));
assert(isequal(size(h.J),[meta.ModelCount meta.GenerationCount]));
assert(isequal(m.seeds(:),meta.Seeds(:)) && isequal(h.seeds(:),m.seeds(:)));
assert(isequal(sum(m.bestMasks,1).',m.featureSelectionCount));
assert(isequal(m.featureSelectionCount/meta.ModelCount,m.featureSelectionFrequency));
assert(isequal(sum(m.bestMasks,2),m.modelSummaryTable.SelectedFeatureCount));
assert(all(isfinite(h.P),'all') && all(h.P>=0 & h.P<=100,'all'));
assert(max(abs(h.J+h.R.*h.F.*h.P),[],'all')<1e-8);
assert(isequal(h.F(:,end),m.modelSummaryTable.SelectedFeatureCount));
assert(isequal(h.J(:,end),m.modelSummaryTable.BestFitness));
assert(isequal(data.History.PSOIteration,(1:meta.GenerationCount).'));
assert(max(abs(data.History.MedianBestFitness-median(h.J,1).'))<1e-8);
assert(isequal(data.History.MedianSelectedDescriptorCount,median(h.F,1).'));
assert(height(data.Top20)==20 && isequal(data.Top20.Rank,(1:20).'));
assert(all(data.Top20.Top20CountAcrossSeeds/meta.ModelCount == ...
    data.Top20.Top20FrequencyAcrossSeeds));
assert(all(data.TopKJaccard.ModelPairCount==nchoosek(meta.ModelCount,2)));
assert(isequal(data.EnsembleSummary.EnsembleSize,meta.EnsembleSizes(:)));
end

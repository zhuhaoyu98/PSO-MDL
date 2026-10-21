function paths = pso_mdl_paths(cfg)
%PSO_MDL_PATHS Create the run-specific output layout below Output.Directory.
paths.Root = cfg.Output.Directory;
paths.Results = paths.Root;
paths.Cache = fullfile(paths.Root, 'cache');
paths.Audit = fullfile(paths.Root, 'audit');
paths.Figures = fullfile(paths.Root, 'figures');
paths.Status = fullfile(paths.Root, 'status');
names = fieldnames(paths);
for index = 1:numel(names)
    directory = paths.(names{index});
    if ~isfolder(directory), mkdir(directory); end
end
end

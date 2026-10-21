function cfg = resolve_pso_mdl_config(configInput)
%RESOLVE_PSO_MDL_CONFIG Resolve profile names or a complete configuration struct.
%   Relative input/output paths are resolved against the repository, never pwd.
%   Each Step adds functions/ before calling this resolver; source subfolders
%   and runtime requirements are initialized here without a separate setup.
root = fileparts(fileparts(mfilename('fullpath')));
folders = {'io', 'preprocessing', 'scoring', 'optimization', 'ranking', 'visualization'};
for index = 1:numel(folders)
    addpath(fullfile(root, 'functions', folders{index}));
end
if verLessThan('matlab', '24.1')
    error('pso_mdl:MATLABVersion', 'MATLAB R2024a or later is required.');
end
if isempty(which('mad')) || isempty(which('corr')) || isempty(ver('stats'))
    error('pso_mdl:MissingToolbox', ...
        'Statistics and Machine Learning Toolbox is required for mad and corr.');
end
if nargin < 1, configInput = "formal"; end
if isstruct(configInput)
    cfg = configInput;
else
    cfg = pso_mdl_config(configInput);
end
groups = {'Profile', 'Input', 'Output', 'Identity', 'Background', ...
    'Preprocessing', 'PSO', 'Stability', 'Ranking', 'Figures', 'IO'};
assert(isscalar(cfg) && all(isfield(cfg, groups)), ...
    'Use pso_mdl_config to construct a complete configuration.');
inputs = {'CandidateDescriptors', 'ReferenceDescriptors', ...
    'CandidateManifest', 'ReferenceManifest', 'MissingnessCache'};
for index = 1:numel(inputs)
    name = inputs{index};
    cfg.Input.(name) = localAbsolutePath(cfg.Input.(name), root);
end
cfg.Output.Directory = localAbsolutePath(cfg.Output.Directory, root);
assert(~isempty(cfg.Output.Directory), 'Output.Directory cannot be empty.');
assert(numel(cfg.PSO.Seeds) >= 2 && numel(unique(cfg.PSO.Seeds)) == numel(cfg.PSO.Seeds), ...
    'At least two unique PSO seeds are required.');
end

function resolved = localAbsolutePath(value, root)
value = string(value);
assert(isscalar(value) && ~ismissing(value), 'Each file path must be a scalar string.');
if strlength(value) == 0, resolved = ''; return; end
if ~isempty(regexp(char(value), '^(?:[A-Za-z]:[\\/]|/|\\\\)', 'once'))
    resolved = char(value);
else
    resolved = fullfile(root, char(value));
end
end

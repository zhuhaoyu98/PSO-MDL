function environment = configure_numerical_execution(cfg)
%CONFIGURE_NUMERICAL_EXECUTION Apply and record the training thread limit.
%   CFG.Execution specifies ThreadLimit and ExpectedBLASBranch. Older complete
%   configurations without this group use the current project's defaults.
%   The returned structure records the effective runtime, not just the request.
%   This function changes only the current MATLAB process's thread limit.
%   The caller is responsible for restoring its previous limit if required.
%   MKL branch selection is not forced after startup. An AVX2 request may be
%   unsupported on non-Intel CPUs; a mismatch is reported without claiming
%   bitwise equivalence to the historical NIRO computation environment.
defaults = pso_mdl_config("formal");
settings = defaults.Execution;
if isfield(cfg, 'Execution')
    names = fieldnames(settings);
    for index = 1:numel(names)
        if isfield(cfg.Execution, names{index})
            settings.(names{index}) = cfg.Execution.(names{index});
        end
    end
end
validateattributes(settings.ThreadLimit, {'numeric'}, ...
    {'real','finite','scalar','integer','positive'}, mfilename, 'ThreadLimit');
expected = string(settings.ExpectedBLASBranch);
assert(isscalar(expected) && ~ismissing(expected), ...
    'ExpectedBLASBranch must be a scalar, nonmissing string.');
previous = maxNumCompThreads;
maxNumCompThreads(settings.ThreadLimit);
actual = maxNumCompThreads;
assert(actual == settings.ThreadLimit, 'The requested thread limit was not applied.');
blas = version('-blas');
token = regexp(blas, 'CNR branch ([^)]+)', 'tokens', 'once');
branch = "unknown";
if ~isempty(token), branch = string(token{1}); end
matched = strlength(expected) == 0 || strcmpi(branch, expected);
environment = struct('MATLABVersion',version,'Release',version('-release'), ...
    'ComputerName',getenv('COMPUTERNAME'),'Processor',getenv('PROCESSOR_IDENTIFIER'), ...
    'BLAS',blas,'BLASBranch',branch,'ExpectedBLASBranch',expected, ...
    'BLASBranchMatchesExpected',matched,'RequestedThreadLimit',settings.ThreadLimit, ...
    'Threads',actual,'PreviousThreadLimit',previous, ...
    'MKL_CBWR',getenv('MKL_CBWR'),'MKL_NUM_THREADS',getenv('MKL_NUM_THREADS'), ...
    'OMP_NUM_THREADS',getenv('OMP_NUM_THREADS'));
fprintf('Numerical execution: %d computational threads; MKL CNR branch %s.\n', actual, branch);
if ~matched
    warning('pso_mdl:BLASBranchMismatch', ...
        ['MATLAB reports CNR branch "%s", not the historical NIRO branch "%s". ' ...
         'The thread limit was applied, but the MKL branch was not forced. ' ...
         'Exact cross-machine reproduction is not guaranteed.'], branch, expected);
end
end

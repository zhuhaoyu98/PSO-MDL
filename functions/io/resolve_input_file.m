function filePath = resolve_input_file(matches, description, isRequired)
%RESOLVE_INPUT_FILE Require one non-directory input, or allow an absent optional input.
%   MATCHES is a file path/pattern or the result of dir.
%   DESCRIPTION identifies the input in errors; paths are retained in diagnostics.
%   ISREQUIRED defaults to true. Multiple matches are always rejected.
if nargin < 3, isRequired = true; end
if ischar(matches) || isstring(matches)
    pattern = string(matches);
    assert(isscalar(pattern) && ~ismissing(pattern), 'Expected one file path or pattern.');
    description = sprintf('%s [%s]', description, pattern);
    matches = match_input_files(pattern);
end
if ~isempty(matches), matches = matches(~[matches.isdir]); end
if isempty(matches) && ~isRequired
    filePath = "";
    return;
end
if numel(matches) ~= 1
    error('mdl:InputMatchCount', ...
        ['%s requires one file (required=%d); found %d. ' ...
         'Check functions/pso_mdl_config.m and the required upstream Step outputs.'], ...
        description, isRequired, numel(matches));
end
filePath = fullfile(matches(1).folder, matches(1).name);
end

function matches = match_input_files(pattern)
%MATCH_INPUT_FILES Match a configured path; an empty optional path matches none.
if strlength(string(pattern)) == 0
    matches = struct('name', {}, 'folder', {}, 'date', {}, 'bytes', {}, ...
        'isdir', {}, 'datenum', {});
else
    matches = dir(pattern);
end
end

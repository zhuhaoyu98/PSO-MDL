function write_json(filePath, value)
%WRITE_JSON Export UTF-8 JSON with indentation and a final newline.
fileId = fopen(filePath, 'w', 'n', 'UTF-8');
assert(fileId >= 0, 'Cannot write JSON: %s', filePath);
cleanup = onCleanup(@() fclose(fileId));
fprintf(fileId, '%s\n', jsonencode(value, 'PrettyPrint', true));
end

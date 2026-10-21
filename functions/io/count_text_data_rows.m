function dataRowCount = count_text_data_rows(filePath)
%COUNT_TEXT_DATA_ROWS Count nonempty data records after a single header row.
%   FILEPATH names a UTF-8 text table. Blank records do not count as molecules.

arguments
    filePath (1, 1) string
end

assert(isfile(filePath), 'Text file not found: %s', filePath);
fileIdentifier = fopen(filePath, 'r', 'n', 'UTF-8');
assert(fileIdentifier >= 0, 'Cannot open text file: %s', filePath);
cleanupObject = onCleanup(@() fclose(fileIdentifier));

headerLine = fgetl(fileIdentifier);
assert(ischar(headerLine), 'Text file is empty: %s', filePath);

dataRowCount = 0;
while true
    currentLine = fgetl(fileIdentifier);
    if ~ischar(currentLine)
        break;
    end
    if ~isempty(currentLine)
        dataRowCount = dataRowCount + 1;
    end
end
end

function [selectedData, selectedNames, descriptorNames, audit] = ...
    read_descriptor_rows( ...
    filePath, selectedRows, expectedTotalRows, ...
    expectedDescriptorCount, chunkRows, progressInterval)
%READ_DESCRIPTOR_ROWS Stream selected molecule rows in their requested order.
%   SELECTEDROWS are unique 1-based physical data-row positions. Expected row
%   and feature counts are checked, including consecutive sequence numbers.
%   Outputs contain selected raw descriptors, names, ordered headers and audit.
%   SELECTEDDATA is numel(SELECTEDROWS)-by-p in requested row order, not sorted
%   source order. Its units are the original Dragon descriptor units.

arguments
    filePath (1, 1) string
    selectedRows (:, 1) double
    expectedTotalRows (1, 1) double {mustBeInteger, mustBePositive}
    expectedDescriptorCount (1, 1) double ...
        {mustBeInteger, mustBePositive}
    chunkRows (1, 1) double {mustBeInteger, mustBePositive} = 500
    progressInterval (1, 1) double ...
        {mustBeInteger, mustBePositive} = 25000
end

assert(~isempty(selectedRows), 'selectedRows must not be empty.');
assert(all(isfinite(selectedRows) & selectedRows == fix(selectedRows) & ...
    selectedRows >= 1 & selectedRows <= expectedTotalRows), ...
    'selectedRows must be finite integers in 1:%d.', expectedTotalRows);
assert(numel(unique(selectedRows)) == numel(selectedRows), ...
    'selectedRows must not contain duplicates.');

[sortedRows, originalPositions] = sort(selectedRows);
selectedCount = numel(selectedRows);

stream = descriptor_table_stream(filePath);
cleanupFile = onCleanup(@() stream.Close());
headerLine = stream.Header;
headerColumns = strsplit( ...
    headerLine, sprintf('\t'), 'CollapseDelimiters', false);

expectedTotalColumns = expectedDescriptorCount + 2;
assert(numel(headerColumns) == expectedTotalColumns, ...
    'Expected %d header columns; found %d: %s', ...
    expectedTotalColumns, numel(headerColumns), filePath);

descriptorNames = string(headerColumns(3:end));
assert(numel(unique(descriptorNames)) == expectedDescriptorCount, ...
    'Duplicate descriptor names in table: %s', filePath);

formatSpec = ['%f%q', repmat('%f', 1, expectedDescriptorCount)];

selectedData = nan(selectedCount, expectedDescriptorCount);
selectedNames = strings(selectedCount, 1);
selectedSequenceNumbers = nan(selectedCount, 1);

totalRowsRead = 0;
nextSelectedPosition = 1;
nextProgress = progressInterval;
firstName = "";
lastName = "";
timerId = tic;

while ~stream.AtEnd()
    data = stream.Read(formatSpec, chunkRows, ...
        'Delimiter', '\t', ...
        'MultipleDelimsAsOne', false, ...
        'CollectOutput', true, ...
        'EmptyValue', NaN, ...
        'ReturnOnError', false);

    if isempty(data) || isempty(data{1})
        if stream.AtEnd()
            break;
        end
        error('Parsing stopped near data row %d: %s', ...
            totalRowsRead + 1, filePath);
    end

    sequenceNumbers = data{1};
    moleculeNames = string(data{2});
    descriptorValues = data{3};
    rowsInChunk = size(descriptorValues, 1);

    assert(numel(sequenceNumbers) == rowsInChunk && ...
        numel(moleculeNames) == rowsInChunk, ...
        'Sequence-number, name, and descriptor row counts disagree: %s', filePath);
    assert(size(descriptorValues, 2) == expectedDescriptorCount, ...
        'Near row %d, parsed only %d descriptors; expected %d: %s', ...
        totalRowsRead + 1, size(descriptorValues, 2), ...
        expectedDescriptorCount, filePath);

    expectedSequence = ...
        (totalRowsRead + 1:totalRowsRead + rowsInChunk).';
    assert(isequal(sequenceNumbers, expectedSequence), ...
        'Nonconsecutive No. values near row %d: %s', ...
        totalRowsRead + 1, filePath);

    if totalRowsRead == 0
        firstName = moleculeNames(1);
    end
    lastName = moleculeNames(end);

    chunkLastGlobalRow = totalRowsRead + rowsInChunk;
    selectionStart = nextSelectedPosition;
    while nextSelectedPosition <= selectedCount && ...
            sortedRows(nextSelectedPosition) <= chunkLastGlobalRow
        nextSelectedPosition = nextSelectedPosition + 1;
    end
    selectionEnd = nextSelectedPosition - 1;

    if selectionEnd >= selectionStart
        sortedPositions = selectionStart:selectionEnd;
        globalRows = sortedRows(sortedPositions);
        localRows = globalRows - totalRowsRead;
        outputPositions = originalPositions(sortedPositions);

        selectedData(outputPositions, :) = ...
            descriptorValues(localRows, :);
        selectedNames(outputPositions) = moleculeNames(localRows);
        selectedSequenceNumbers(outputPositions) = ...
            sequenceNumbers(localRows);
    end

    totalRowsRead = chunkLastGlobalRow;
    if totalRowsRead >= nextProgress
        fprintf('  Scanned %d / %d rows (%.1f%%); extracted %d / %d requested rows\n', ...
            totalRowsRead, expectedTotalRows, ...
            100 * totalRowsRead / expectedTotalRows, ...
            nextSelectedPosition - 1, selectedCount);
        nextProgress = nextProgress + progressInterval;
    end
end

clear cleanupFile

assert(totalRowsRead == expectedTotalRows, ...
    'Expected to scan %d rows; scanned %d: %s', ...
    expectedTotalRows, totalRowsRead, filePath);
assert(nextSelectedPosition == selectedCount + 1, ...
    'Extracted only %d / %d requested rows: %s', ...
    nextSelectedPosition - 1, selectedCount, filePath);
assert(isequal(selectedSequenceNumbers, selectedRows), ...
    'Output row order or sequence numbers disagree with selectedRows: %s', filePath);
assert(~any(ismissing(selectedNames)), ...
    'Empty NAME in extracted rows: %s', filePath);

audit = struct();
audit.FilePath = char(filePath);
audit.FileBytes = stream.FileBytes;
audit.SourceFiles = stream.Files;
audit.PartsManifestFile = stream.ManifestFile;
audit.TotalRows = totalRowsRead;
audit.TotalColumns = expectedTotalColumns;
audit.DescriptorColumns = expectedDescriptorCount;
audit.SelectedRows = selectedCount;
audit.FirstName = char(firstName);
audit.LastName = char(lastName);
audit.ElapsedSeconds = toc(timerId);
audit.ScannedAt = char(datetime( ...
    'now', 'TimeZone', 'local', ...
    'Format', 'yyyy-MM-dd''T''HH:mm:ssXXX'));
end

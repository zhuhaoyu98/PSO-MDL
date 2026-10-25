function [selectedData, moleculeNames, selectedDescriptorNames, audit] = ...
    read_descriptor_columns( ...
    filePath, selectedDescriptorIndices, expectedTotalRows, ...
    expectedDescriptorCount, chunkRows, progressInterval)
%READ_DESCRIPTOR_COLUMNS Stream selected descriptor columns while preserving requested order.
%   Inputs specify 1-based unique descriptor indices, expected dimensions and
%   chunk/progress sizes. The two leading identifier columns are not descriptors.
%   Outputs retain all molecules and return requested columns in input order.
%   SELECTEDDATA is nCompound-by-numel(SELECTEDDESCRIPTORINDICES), in raw units.
%   Use descriptor indices 1:p, not text-column indices that include No. and NAME.

arguments
    filePath (1, 1) string
    selectedDescriptorIndices (:, 1) double
    expectedTotalRows (1, 1) double {mustBeInteger, mustBePositive}
    expectedDescriptorCount (1, 1) double ...
        {mustBeInteger, mustBePositive}
    chunkRows (1, 1) double {mustBeInteger, mustBePositive} = 2000
    progressInterval (1, 1) double ...
        {mustBeInteger, mustBePositive} = 25000
end

assert(~isempty(selectedDescriptorIndices), ...
    'selectedDescriptorIndices must not be empty.');
assert(all(isfinite(selectedDescriptorIndices) & ...
    selectedDescriptorIndices == fix(selectedDescriptorIndices) & ...
    selectedDescriptorIndices >= 1 & ...
    selectedDescriptorIndices <= expectedDescriptorCount), ...
    'Descriptor indices must be finite integers in 1:%d.', expectedDescriptorCount);
assert(numel(unique(selectedDescriptorIndices)) == ...
    numel(selectedDescriptorIndices), ...
    'selectedDescriptorIndices must not contain duplicates.');

selectedDescriptorIndices = selectedDescriptorIndices(:);
selectedCount = numel(selectedDescriptorIndices);
[selectedIndicesInFileOrder, ~] = sort(selectedDescriptorIndices);
[foundRequestedIndex, requestedOrderFromFileOrder] = ismember( ...
    selectedDescriptorIndices, selectedIndicesInFileOrder);
assert(all(foundRequestedIndex), ...
    'Cannot construct the requested descriptor-column ordering.');

stream = descriptor_table_stream(filePath);
cleanupFile = onCleanup(@() stream.Close());
headerLine = stream.Header;
headerColumns = strsplit( ...
    headerLine, sprintf('\t'), 'CollapseDelimiters', false);
expectedTotalColumns = expectedDescriptorCount + 2;
assert(numel(headerColumns) == expectedTotalColumns, ...
    'Expected %d header columns; found %d: %s', ...
    expectedTotalColumns, numel(headerColumns), filePath);

allDescriptorNames = string(headerColumns(3:end));
assert(numel(unique(allDescriptorNames)) == expectedDescriptorCount, ...
    'Duplicate descriptor names in table: %s', filePath);
selectedDescriptorNames = ...
    allDescriptorNames(selectedDescriptorIndices);
selectedDescriptorNames = selectedDescriptorNames(:);

descriptorFormats = repmat("%*f", expectedDescriptorCount, 1);
descriptorFormats(selectedDescriptorIndices) = "%f";
formatSpec = ['%f%q', char(join(descriptorFormats, ""))];

selectedData = nan(expectedTotalRows, selectedCount);
moleculeNames = strings(expectedTotalRows, 1);

totalRowsRead = 0;
nextProgress = progressInterval;
firstName = "";
lastName = "";
timerId = tic;

while ~stream.AtEnd()
    data = stream.Read(formatSpec, chunkRows, ...
        'Delimiter', '\t', ...
        'MultipleDelimsAsOne', false, ...
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
    currentNames = string(data{2});

    descriptorValuesInFileOrder = horzcat(data{3:end});
    descriptorValues = descriptorValuesInFileOrder( ...
        :, requestedOrderFromFileOrder);
    rowsInChunk = numel(sequenceNumbers);

    assert(numel(currentNames) == rowsInChunk, ...
        'Sequence-number and name row counts disagree: %s', filePath);
    assert(size(descriptorValues, 1) == rowsInChunk && ...
        size(descriptorValues, 2) == selectedCount, ...
        'Near row %d, failed to parse %d selected descriptors: %s', ...
        totalRowsRead + 1, selectedCount, filePath);

    expectedSequence = ...
        (totalRowsRead + 1:totalRowsRead + rowsInChunk).';
    assert(isequal(sequenceNumbers, expectedSequence), ...
        'Nonconsecutive No. values near row %d: %s', ...
        totalRowsRead + 1, filePath);

    if totalRowsRead == 0
        firstName = currentNames(1);
    end
    lastName = currentNames(end);

    outputRows = totalRowsRead + (1:rowsInChunk);
    selectedData(outputRows, :) = descriptorValues;
    moleculeNames(outputRows) = currentNames;
    totalRowsRead = totalRowsRead + rowsInChunk;

    if totalRowsRead >= nextProgress
        fprintf('  Scanned %d / %d rows (%.1f%%); retained %d descriptor columns\n', ...
            totalRowsRead, expectedTotalRows, ...
            100 * totalRowsRead / expectedTotalRows, selectedCount);
        nextProgress = nextProgress + progressInterval;
    end
end

clear cleanupFile

assert(totalRowsRead == expectedTotalRows, ...
    'Expected to scan %d rows; scanned %d: %s', ...
    expectedTotalRows, totalRowsRead, filePath);
assert(~any(ismissing(moleculeNames)), ...
    'Empty NAME in output: %s', filePath);

audit = struct();
audit.FilePath = char(filePath);
audit.FileBytes = stream.FileBytes;
audit.SourceFiles = stream.Files;
audit.PartsManifestFile = stream.ManifestFile;
audit.TotalRows = totalRowsRead;
audit.TotalColumns = expectedTotalColumns;
audit.InputDescriptorColumns = expectedDescriptorCount;
audit.SelectedDescriptorIndices = selectedDescriptorIndices;
audit.SelectedDescriptorNames = selectedDescriptorNames;
audit.SelectedDescriptorCount = selectedCount;
audit.FirstName = char(firstName);
audit.LastName = char(lastName);
audit.ElapsedSeconds = toc(timerId);
audit.ScannedAt = char(datetime( ...
    'now', 'TimeZone', 'local', ...
    'Format', 'yyyy-MM-dd''T''HH:mm:ssXXX'));

end

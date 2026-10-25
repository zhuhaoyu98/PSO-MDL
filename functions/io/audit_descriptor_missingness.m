function [missingCount, infCount, descriptorNames, audit] = ...
    audit_descriptor_missingness( ...
    filePath, expectedTotalRows, expectedDescriptorCount, ...
    chunkRows, progressInterval)
%AUDIT_DESCRIPTOR_MISSINGNESS Count NaN and Inf entries in a descriptor file by streaming rows.
%   Inputs specify the expected row/descriptor counts and chunk size. The two
%   leading columns must be a consecutive integer and a nonempty molecule name.
%   Outputs are 1-by-feature counts, ordered descriptor names and a read audit.
%   NaN and Inf are counted separately. This function neither imputes values
%   nor decides which descriptors are retained for modeling.

arguments
    filePath (1, 1) string
    expectedTotalRows (1, 1) double {mustBeInteger, mustBePositive}
    expectedDescriptorCount (1, 1) double ...
        {mustBeInteger, mustBePositive}
    chunkRows (1, 1) double {mustBeInteger, mustBePositive} = 500
    progressInterval (1, 1) double ...
        {mustBeInteger, mustBePositive} = 25000
end

[~, descriptorNames, headerAudit] = read_descriptor_header(filePath);
assert(numel(descriptorNames) == expectedDescriptorCount, ...
    'Descriptor header has %d features; expected %d: %s', ...
    numel(descriptorNames), expectedDescriptorCount, filePath);

stream = descriptor_table_stream(filePath);
cleanupObject = onCleanup(@() stream.Close());

formatSpec = ['%f%q', repmat('%f', 1, expectedDescriptorCount)];
missingCount = zeros(1, expectedDescriptorCount);
infCount = zeros(1, expectedDescriptorCount);
totalRowsRead = 0;
nextProgress = progressInterval;
firstName = "";
lastName = "";
timerIdentifier = tic;

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
        error('mdl:DragonParseStopped', ...
            'Parsing stopped near data row %d: %s', ...
            totalRowsRead + 1, filePath);
    end

    sequenceNumber = data{1};
    moleculeName = string(data{2});
    descriptorBlock = horzcat(data{3:end});
    rowsInBlock = numel(sequenceNumber);

    assert(size(descriptorBlock, 1) == rowsInBlock && ...
        size(descriptorBlock, 2) == expectedDescriptorCount, ...
        'Near row %d, failed to parse all %d descriptors: %s', ...
        totalRowsRead + 1, expectedDescriptorCount, filePath);
    expectedSequence = ...
        (totalRowsRead + 1:totalRowsRead + rowsInBlock).';
    assert(isequal(sequenceNumber, expectedSequence), ...
        'Nonconsecutive No. values near row %d: %s', ...
        totalRowsRead + 1, filePath);
    assert(numel(moleculeName) == rowsInBlock && ...
        ~any(ismissing(moleculeName) | strlength(moleculeName) == 0), ...
        'Empty NAME value near row %d: %s', ...
        totalRowsRead + 1, filePath);

    if totalRowsRead == 0
        firstName = moleculeName(1);
    end
    lastName = moleculeName(end);
    missingCount = missingCount + sum(isnan(descriptorBlock), 1);
    infCount = infCount + sum(isinf(descriptorBlock), 1);
    totalRowsRead = totalRowsRead + rowsInBlock;

    if totalRowsRead >= nextProgress
        fprintf('  Missingness audit: %d / %d rows (%.1f%%)\n', ...
            totalRowsRead, expectedTotalRows, ...
            100 * totalRowsRead / expectedTotalRows);
        nextProgress = nextProgress + progressInterval;
    end
end

clear cleanupObject
assert(totalRowsRead == expectedTotalRows, ...
    'Expected %d rows; read %d: %s', ...
    expectedTotalRows, totalRowsRead, filePath);

audit = headerAudit;
audit.TotalRows = totalRowsRead;
audit.FirstName = char(firstName);
audit.LastName = char(lastName);
audit.ElapsedSeconds = toc(timerIdentifier);
audit.AuditedAt = char(datetime( ...
    'now', 'TimeZone', 'local', ...
    'Format', 'yyyy-MM-dd''T''HH:mm:ssXXX'));
end

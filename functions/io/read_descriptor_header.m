function [identifierNames, descriptorNames, audit] = ...
    read_descriptor_header(filePath)
%READ_DESCRIPTOR_HEADER Read and validate a Dragon-compatible descriptor table header.
%   The first two tab-separated columns identify records. Remaining, nonempty,
%   unique names define the descriptor space and its order.

arguments
    filePath (1, 1) string
end

stream = descriptor_table_stream(filePath);
cleanupObject = onCleanup(@() stream.Close());
headerLine = stream.Header;
headerColumns = string(strsplit( ...
    headerLine, sprintf('\t'), 'CollapseDelimiters', false));
assert(numel(headerColumns) >= 3, ...
    'Descriptor table requires two identifier columns and at least one descriptor: %s', filePath);
assert(all(strlength(strtrim(headerColumns)) > 0), ...
    'Empty column name in descriptor table: %s', filePath);

identifierNames = headerColumns(1:2);
descriptorNames = headerColumns(3:end);
assert(numel(unique(descriptorNames)) == numel(descriptorNames), ...
    'Duplicate descriptor names in table: %s', filePath);

audit = struct();
audit.FilePath = char(filePath);
audit.FileBytes = stream.FileBytes;
audit.SourceFiles = stream.Files;
audit.PartsManifestFile = stream.ManifestFile;
audit.ExpectedPartSHA256 = stream.ExpectedHashes;
audit.IdentifierColumnNames = identifierNames;
audit.DescriptorCount = numel(descriptorNames);
audit.FirstDescriptor = char(descriptorNames(1));
audit.LastDescriptor = char(descriptorNames(end));
end

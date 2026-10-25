function stream = descriptor_table_stream(sourcePath)
%DESCRIPTOR_TABLE_STREAM Read one TXT or an ordered, manifested TXT bundle.
%   SOURCEPATH is a UTF-8 Dragon table or a directory containing
%   descriptor_parts.json. Each part repeats the same header and retains the
%   original global No. values. READ uses textscan options and returns full
%   logical chunks across part boundaries; no descriptor values are refitted.
%   Call CLOSE, normally through onCleanup, to release the current file.
%   Structural checks here cover part order, sizes, headers and row continuity;
%   full part SHA-256 verification is performed by Step 00 when enabled.
arguments
    sourcePath (1, 1) string
end

manifestPath = "";
expectedRows = NaN;
expectedHashes = "";
if isfolder(sourcePath)
    manifestPath = fullfile(sourcePath, 'descriptor_parts.json');
    assert(isfile(manifestPath), 'mdl:PartsManifestMissing', ...
        'Descriptor part directory requires descriptor_parts.json: %s', sourcePath);
    manifest = jsondecode(fileread(manifestPath));
    required = {'Version', 'PartCount', 'TotalRows', 'DescriptorCount', ...
        'SourceBytes', 'SourceSHA256', 'HeaderSHA256', 'Parts'};
    assert(isstruct(manifest) && all(isfield(manifest, required)), ...
        'mdl:PartsManifestInvalid', 'Incomplete descriptor parts manifest: %s', manifestPath);
    assert(isequal(manifest.Version, 1) && isstruct(manifest.Parts), ...
        'mdl:PartsManifestInvalid', 'Unsupported descriptor parts manifest.');
    parts = manifest.Parts(:);
    assert(all(isfield(parts, {'File', 'FirstRow', 'LastRow', 'Rows', 'Bytes', 'SHA256'})), ...
        'mdl:PartsManifestInvalid', 'Incomplete descriptor part records.');
    counts = [manifest.PartCount, manifest.TotalRows, manifest.DescriptorCount, ...
        manifest.SourceBytes, [parts.FirstRow], [parts.LastRow], [parts.Rows], [parts.Bytes]];
    assert(all(isfinite(counts) & counts >= 1 & counts == fix(counts)) && ...
        manifest.PartCount == numel(parts), ...
        'mdl:PartsManifestInvalid', 'Invalid descriptor part dimensions.');
    names = string({parts.File}).';
    expectedNames = compose('candidates_part_%03d.txt', (1:numel(parts)).');
    assert(isequal(names, expectedNames), 'mdl:PartsOrder', ...
        'Parts must be listed once in their original numbered order.');
    listed = dir(fullfile(sourcePath, '*.txt'));
    assert(isequal(sort(string({listed.name}).'), sort(names)), ...
        'mdl:PartsInventory', 'Missing or unlisted descriptor TXT parts: %s', sourcePath);
    expectedRows = double([parts.Rows]).';
    firstRows = double([parts.FirstRow]).';
    lastRows = double([parts.LastRow]).';
    assert(isequal(firstRows, [1; cumsum(expectedRows(1:end-1)) + 1]) && ...
        isequal(lastRows, cumsum(expectedRows)) && sum(expectedRows) == manifest.TotalRows, ...
        'mdl:PartsRows', 'Part row ranges are not consecutive.');
    expectedHashes = lower(string({parts.SHA256}).');
    allHashes = [expectedHashes; string(manifest.SourceSHA256); string(manifest.HeaderSHA256)];
    assert(all(cellfun(@(s) ~isempty(regexp(s, '^[0-9a-fA-F]{64}$', 'once')), ...
        cellstr(allHashes))), 'mdl:PartsManifestInvalid', 'Invalid SHA-256 records.');
    filePaths = fullfile(sourcePath, names);
else
    assert(isfile(sourcePath), 'mdl:DescriptorMissing', ...
        'Descriptor TXT or part directory not found: %s', sourcePath);
    filePaths = sourcePath;
end

partCount = numel(filePaths);
fileBytes = zeros(partCount, 1);
headerBytes = zeros(partCount, 1);
header = '';
for k = 1:partCount
    entry = dir(filePaths(k));
    assert(isscalar(entry) && ~entry.isdir, 'mdl:DescriptorMissing', ...
        'Descriptor file not found: %s', filePaths(k));
    fileBytes(k) = entry.bytes;
    [currentHeader, headerBytes(k), headerSHA] = localHeader(filePaths(k));
    if k == 1
        header = currentHeader;
    end
    assert(strcmp(currentHeader, header), 'mdl:PartsHeader', ...
        'Descriptor names or order differ in part: %s', filePaths(k));
    if strlength(manifestPath) > 0
        assert(entry.bytes == parts(k).Bytes && strcmpi(headerSHA, manifest.HeaderSHA256), ...
            'mdl:PartsIntegrity', 'Part size or header differs from its manifest: %s', filePaths(k));
    end
end
headerColumns = string(strsplit(header, sprintf('\t'), 'CollapseDelimiters', false));
assert(numel(headerColumns) >= 3 && all(strlength(strtrim(headerColumns)) > 0), ...
    'mdl:DescriptorHeader', 'Expected two identifiers and nonempty descriptor names: %s', sourcePath);
assert(numel(unique(headerColumns(3:end))) == numel(headerColumns) - 2, ...
    'mdl:DescriptorHeader', 'Duplicate descriptor names: %s', sourcePath);
if strlength(manifestPath) > 0
    assert(numel(headerColumns) - 2 == manifest.DescriptorCount && ...
        sum(fileBytes) - sum(headerBytes(2:end)) == manifest.SourceBytes, ...
        'mdl:PartsIntegrity', 'Part schema or reconstructed byte count disagrees with manifest.');
end

currentPart = 1;
fileId = -1;
rowsInPart = 0;
totalRows = 0;
openPart();
stream = struct('Header', header, 'Files', filePaths, ...
    'FileBytes', sum(fileBytes), 'PartBytes', fileBytes, ...
    'ManifestFile', manifestPath, 'ExpectedHashes', expectedHashes, ...
    'Read', @readChunk, 'AtEnd', @atEnd, 'Close', @closeStream);

    function openPart()
        if currentPart > partCount, return; end
        fileId = fopen(filePaths(currentPart), 'r', 'n', 'UTF-8');
        assert(fileId >= 0, 'mdl:DescriptorOpen', 'Cannot open: %s', filePaths(currentPart));
        currentHeader = fgetl(fileId);
        if ~ischar(currentHeader) || ~strcmp(currentHeader, header)
            closeStream();
            error('mdl:PartsHeader', 'Descriptor header changed while opening a part.');
        end
        rowsInPart = 0;
    end

    function data = readChunk(formatSpec, chunkRows, varargin)
        data = {};
        remaining = chunkRows;
        while remaining > 0 && currentPart <= partCount
            block = textscan(fileId, formatSpec, remaining, varargin{:});
            if isempty(block) || isempty(block{1})
                assert(feof(fileId), 'mdl:DragonParseStopped', ...
                    'Parsing stopped near global row %d: %s', totalRows + 1, filePaths(currentPart));
                nRows = 0;
            else
                nRows = size(block{1}, 1);
                assert(isequal(block{1}, (totalRows + (1:nRows)).'), ...
                    'mdl:DescriptorSequence', 'Nonconsecutive global No. values near row %d.', totalRows + 1);
                assert(all(cellfun(@(v) size(v, 1) == nRows, block)), ...
                    'mdl:DescriptorShape', 'Unequal parsed column lengths near row %d.', totalRows + 1);
                if isempty(data)
                    data = block;
                else
                    for column = 1:numel(data)
                        data{column} = [data{column}; block{column}]; %#ok<AGROW>
                    end
                end
                rowsInPart = rowsInPart + nRows;
                totalRows = totalRows + nRows;
                remaining = remaining - nRows;
                if strlength(manifestPath) > 0
                    assert(rowsInPart <= expectedRows(currentPart), 'mdl:PartsRows', ...
                        'Too many records in part: %s', filePaths(currentPart));
                end
            end
            if feof(fileId)
                if strlength(manifestPath) > 0
                    assert(rowsInPart == expectedRows(currentPart), 'mdl:PartsRows', ...
                        'Wrong record count in part: %s', filePaths(currentPart));
                end
                fclose(fileId);
                fileId = -1;
                currentPart = currentPart + 1;
                openPart();
            end
        end
    end

    function value = atEnd()
        value = currentPart > partCount;
    end

    function closeStream()
        if fileId >= 0
            fclose(fileId);
            fileId = -1;
        end
    end
end

function [header, byteCount, hash] = localHeader(filePath)
fileId = fopen(filePath, 'r', 'n', 'UTF-8');
assert(fileId >= 0, 'mdl:DescriptorOpen', 'Cannot open: %s', filePath);
cleanup = onCleanup(@() fclose(fileId));
header = fgetl(fileId);
assert(ischar(header), 'mdl:DescriptorHeader', 'Empty descriptor table: %s', filePath);
byteCount = ftell(fileId);
fseek(fileId, 0, 'bof');
rawHeader = fread(fileId, byteCount, '*uint8');
digest = java.security.MessageDigest.getInstance('SHA-256');
digest.update(typecast(rawHeader, 'int8'));
hash = lower(reshape(dec2hex(typecast(digest.digest(), 'uint8'), 2).', 1, []));
end

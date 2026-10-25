function [distanceByMask, moleculeNames, selectedDescriptorNames, audit] = ...
    score_candidates_streaming( ...
    filePath, selectedDescriptorIndices, imputationValues, ...
    normalizationScale, referenceVector, masks, ...
    expectedTotalRows, expectedDescriptorCount, ...
    chunkRows, progressInterval)
%SCORE_CANDIDATES_STREAMING Compute distances in fitted descriptor subspaces.
%   FILEPATH is a Dragon TXT or ordered part directory. Descriptor indices refer
%   to the raw columns after the two identifiers. MASKS is M-by-pUnion;
%   REFERENCEVECTOR is the scaled medoid (1-by-pUnion). Imputation values and
%   scaling divisors are fixed 1-by-pUnion rows from background preprocessing.
%
%   DISTANCEBYMASK is N-by-M: d_im=sqrt(sum_j z_mj*(x_ij-x_medoid,j)^2).
%   All N source rows are scored; Step 03 subsequently retains screening rows.
%   No per-compound columns are dropped, and no preprocessing is refitted.

arguments
    filePath (1, 1) string
    selectedDescriptorIndices (:, 1) double
    imputationValues (1, :) double
    normalizationScale (1, :) double
    referenceVector (1, :) double
    masks {mustBeNumericOrLogical}
    expectedTotalRows (1, 1) double {mustBeInteger, mustBePositive}
    expectedDescriptorCount (1, 1) double ...
        {mustBeInteger, mustBePositive}
    chunkRows (1, 1) double {mustBeInteger, mustBePositive} = 500
    progressInterval (1, 1) double ...
        {mustBeInteger, mustBePositive} = 25000
end

assert(~isempty(selectedDescriptorIndices), ...
    'selectedDescriptorIndices must not be empty.');
assert(all(isfinite(selectedDescriptorIndices) ...
    & selectedDescriptorIndices == fix(selectedDescriptorIndices) ...
    & selectedDescriptorIndices >= 1 ...
    & selectedDescriptorIndices <= expectedDescriptorCount), ...
    'Descriptor indices must be finite integers in 1:%d.', expectedDescriptorCount);
assert(numel(unique(selectedDescriptorIndices)) ...
    == numel(selectedDescriptorIndices), ...
    'selectedDescriptorIndices must not contain duplicates.');

selectedDescriptorIndices = selectedDescriptorIndices(:);
selectedCount = numel(selectedDescriptorIndices);
imputationValues = imputationValues(:).';
normalizationScale = normalizationScale(:).';
referenceVector = referenceVector(:).';
masks = logical(masks);

assert(numel(imputationValues) == selectedCount, ...
    'Frozen imputation-value count disagrees with selected features.');
assert(numel(normalizationScale) == selectedCount, ...
    'Frozen scale count disagrees with selected features.');
assert(numel(referenceVector) == selectedCount, ...
    'Reference-vector length disagrees with selected features.');
assert(size(masks, 2) == selectedCount && size(masks, 1) >= 1, ...
    'Model masks must be nModel-by-selectedCount.');
assert(all(sum(masks, 2) > 0), 'Every model mask must select at least one feature.');
assert(all(isfinite(imputationValues)), 'Frozen imputation values contain NaN or Inf.');
assert(all(isfinite(normalizationScale) & normalizationScale > 0), ...
    'Frozen normalization scales must be finite and positive.');
assert(all(isfinite(referenceVector)), 'Reference vector contains NaN or Inf.');

[selectedIndicesInFileOrder, ~] = sort(selectedDescriptorIndices);
[foundRequestedIndex, requestedOrderFromFileOrder] = ismember( ...
    selectedDescriptorIndices, selectedIndicesInFileOrder);
assert(all(foundRequestedIndex), 'Cannot construct the requested column ordering.');

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

modelCount = size(masks, 1);
distanceByMask = nan(expectedTotalRows, modelCount);
moleculeNames = strings(expectedTotalRows, 1);
imputedCountBySelectedFeature = zeros(selectedCount, 1);

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
            break
        end
        error('Parsing stopped near data row %d: %s', ...
            totalRowsRead + 1, filePath);
    end

    sequenceNumbers = data{1};
    currentNames = string(data{2});
    descriptorValuesInFileOrder = horzcat(data{3:end});
    descriptorValues = descriptorValuesInFileOrder( ...
        :, requestedOrderFromFileOrder);
    clear descriptorValuesInFileOrder data

    rowsInChunk = numel(sequenceNumbers);
    assert(numel(currentNames) == rowsInChunk, ...
        'Sequence-number and name row counts disagree: %s', filePath);
    assert(size(descriptorValues, 1) == rowsInChunk ...
        && size(descriptorValues, 2) == selectedCount, ...
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

    % Apply one common imputation/scaling transform to every compound.
    % Keeping the full selected schema makes distances comparable within a model.
    nonfiniteMask = ~isfinite(descriptorValues);
    if any(nonfiniteMask, 'all')
        imputedCountBySelectedFeature = ...
            imputedCountBySelectedFeature + sum(nonfiniteMask, 1).';
        [~, nonfiniteColumn] = find(nonfiniteMask);
        imputationColumn = imputationValues(:);
        descriptorValues(nonfiniteMask) = ...
            imputationColumn(nonfiniteColumn);
    end
    descriptorValues = descriptorValues ./ normalizationScale;
    assert(all(isfinite(descriptorValues), 'all'), ...
        'Frozen preprocessing left NaN or Inf near row %d.', ...
        totalRowsRead + 1);

    % Matrix multiplication selects each model's columns from the union.
    % The square root is taken after summing squared coordinate differences.
    squaredDelta = (descriptorValues - referenceVector).^2;
    currentDistance = sqrt(squaredDelta * double(masks).');

    outputRows = totalRowsRead + (1:rowsInChunk);
    distanceByMask(outputRows, :) = currentDistance;
    moleculeNames(outputRows) = currentNames;
    totalRowsRead = totalRowsRead + rowsInChunk;

    if totalRowsRead >= nextProgress
        fprintf(['  Scanned and scored %d / %d rows (%.1f%%);' ...
            ' %d union descriptors, %d models\n'], ...
            totalRowsRead, expectedTotalRows, ...
            100 * totalRowsRead / expectedTotalRows, ...
            selectedCount, modelCount);
        nextProgress = nextProgress + progressInterval;
    end
end

clear cleanupFile

assert(totalRowsRead == expectedTotalRows, ...
    'Expected to scan %d rows; scanned %d: %s', ...
    expectedTotalRows, totalRowsRead, filePath);
assert(~any(ismissing(moleculeNames)), ...
    'Empty NAME in output: %s', filePath);
assert(all(isfinite(distanceByMask), 'all'), ...
    'Output distances contain NaN or Inf.');

audit = struct();
audit.Version = 1;
audit.Method = 'streamed frozen preprocessing and medoid Euclidean distance';
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
audit.ModelCount = modelCount;
audit.ImputedCountBySelectedFeature = ...
    imputedCountBySelectedFeature;
audit.TotalImputedValues = sum(imputedCountBySelectedFeature);
audit.FirstName = char(firstName);
audit.LastName = char(lastName);
audit.ElapsedSeconds = toc(timerId);
audit.ScannedAt = char(datetime( ...
    'now', 'TimeZone', 'local', ...
    'Format', 'yyyy-MM-dd''T''HH:mm:ssXXX'));
end

function mustBeNumericOrLogical(value)
if ~(isnumeric(value) || islogical(value)) || ~ismatrix(value)
    error('masks must be a two-dimensional numeric or logical matrix.');
end
if any(value(:) ~= 0 & value(:) ~= 1)
    error('masks must contain only zeros and ones.');
end
end

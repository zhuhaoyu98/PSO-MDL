function [auditTable, summary] = ...
    deduplicate_candidate_manifest(manifestFile, options)
%DEDUPLICATE_CANDIDATE_MANIFEST Deduplicate a row-aligned identity manifest.
%   Exact stable identifiers are trimmed and uppercased. Representatives are
%   chosen by configured stable fields, then SourceRow for exact ties. Missing
%   identities and overlaps with the reference set follow the explicit options.
%   This function does not infer identities or standardize chemical structures.
%   Outputs are a row-by-row audit table and summary; SourceRow is an IO locator.
%   Exact full-identifier matching is not scaffold deduplication or salt/tautomer
%   standardization; the molecular identity policy is supplied by the manifest.

arguments
    manifestFile (1, 1) string
    options.IdentifierColumn (1, 1) string = "INCHIKEY"
    options.DragonNameColumn (1, 1) string = "DragonName"
    options.OriginalTitleColumn (1, 1) string = "OriginalTitle"
    options.ExternalIdColumn (1, 1) string = "IDNUMBER"
    options.SmilesColumn (1, 1) string = "SMILES"
    options.TieBreakColumns (1, :) string = ["IDNUMBER", "DragonName"]
    options.MissingIdentifierPolicy (1, 1) string = "exclude"
    options.ReferenceManifestFile (1, 1) string = ""
    options.ReferenceIdentifierColumn (1, 1) string = "INCHIKEY"
    options.ExcludeReferenceOverlap (1, 1) logical = true
end

assert(isfile(manifestFile), 'Candidate identity manifest not found: %s', manifestFile);
missingPolicy = lower(strtrim(options.MissingIdentifierPolicy));
assert(ismember(missingPolicy, ["exclude", "error"]), ...
    'MissingIdentifierPolicy must be "exclude" or "error".');

importOptions = detectImportOptions(manifestFile, ...
    'FileType', 'text', 'Delimiter', '\t', ...
    'TextType', 'string', 'VariableNamingRule', 'preserve');
manifest = readtable(manifestFile, importOptions);
variableNames = string(manifest.Properties.VariableNames);
assert(ismember(options.IdentifierColumn, variableNames), ...
    'Candidate manifest is missing identifier column %s.', options.IdentifierColumn);
assert(ismember(options.DragonNameColumn, variableNames), ...
    'Candidate manifest is missing descriptor-name column %s.', options.DragonNameColumn);

candidateCount = height(manifest);
assert(candidateCount > 0, 'Candidate manifest contains no data rows.');
sourceRow = (1:candidateCount).';
dragonName = localGetStringColumn( ...
    manifest, options.DragonNameColumn, candidateCount);
assert(all(~ismissing(dragonName) & strlength(strtrim(dragonName)) > 0), ...
    'Candidate manifest column %s contains empty names.', options.DragonNameColumn);

rawIdentifier = localGetStringColumn( ...
    manifest, options.IdentifierColumn, candidateCount);
stableIdentifier = upper(strtrim(rawIdentifier));
stableIdentifier(ismissing(stableIdentifier)) = "";
isMissingIdentifier = strlength(stableIdentifier) == 0;
if missingPolicy == "error" && any(isMissingIdentifier)
    error('mdl:MissingStableIdentifier', ...
        ['Candidate manifest has %d rows without %s. Supply stable identities,' ...
         ' or explicitly set MissingIdentifierPolicy to "exclude".'], ...
        nnz(isMissingIdentifier), options.IdentifierColumn);
end

referenceIdentifier = strings(0, 1);
if options.ExcludeReferenceOverlap
    assert(strlength(options.ReferenceManifestFile) > 0 && ...
        isfile(options.ReferenceManifestFile), ...
        ['ExcludeReferenceOverlap=true requires a valid' ...
         'ReferenceManifestFile。']);
    referenceOptions = detectImportOptions( ...
        options.ReferenceManifestFile, ...
        'FileType', 'text', 'Delimiter', '\t', ...
        'TextType', 'string', 'VariableNamingRule', 'preserve');
    referenceManifest = readtable( ...
        options.ReferenceManifestFile, referenceOptions);
    referenceVariables = string( ...
        referenceManifest.Properties.VariableNames);
    assert(ismember(options.ReferenceIdentifierColumn, ...
        referenceVariables), ...
        'Reference manifest is missing identifier column %s.', ...
        options.ReferenceIdentifierColumn);
    referenceIdentifier = upper(strtrim(string( ...
        referenceManifest.(char( ...
        options.ReferenceIdentifierColumn)))));
    referenceIdentifier(ismissing(referenceIdentifier)) = "";
    referenceIdentifier = unique( ...
        referenceIdentifier(strlength(referenceIdentifier) > 0));
    assert(~isempty(referenceIdentifier), ...
        'Reference manifest has no nonempty stable identifiers.');
end

isReferenceOverlap = false(candidateCount, 1);
if ~isempty(referenceIdentifier)
    isReferenceOverlap = ~isMissingIdentifier & ...
        ismember(stableIdentifier, referenceIdentifier);
end

% Select one representative per identifier before background assignment.
% Keep all rows in the audit so every exclusion maps to its original descriptor row.
nonmissingRows = find(~isMissingIdentifier);
identifierGroupSize = zeros(candidateCount, 1);
representativeSourceRow = nan(candidateCount, 1);
isRepresentative = false(candidateCount, 1);
duplicateGroupCount = 0;
uniqueCandidateIdentifierCount = 0;

if ~isempty(nonmissingRows)
    [uniqueIdentifier, ~, groupIndex] = unique( ...
        stableIdentifier(nonmissingRows), 'sorted');
    groupSize = accumarray(groupIndex, 1);
    identifierGroupSize(nonmissingRows) = groupSize(groupIndex);
    duplicateGroupCount = nnz(groupSize > 1);
    uniqueCandidateIdentifierCount = numel(uniqueIdentifier);

    sortTable = table(stableIdentifier(nonmissingRows), ...
        'VariableNames', {'StableIdentifier'});
    usedTieBreakColumns = strings(0, 1);
    for columnIndex = 1:numel(options.TieBreakColumns)
        currentColumn = options.TieBreakColumns(columnIndex);
        if ismember(currentColumn, variableNames)
            tieValue = upper(strtrim(localGetStringColumn( ...
                manifest, currentColumn, candidateCount)));
            tieValue(ismissing(tieValue)) = "";
            sortTable.(sprintf('TieBreak%02d', columnIndex)) = ...
                tieValue(nonmissingRows);
            usedTieBreakColumns(end + 1, 1) = currentColumn;
        end
    end
    sortTable.SourceRowFinalTieBreak = sourceRow(nonmissingRows);
    [~, sortOrder] = sortrows(sortTable);
    sortedSourceRows = nonmissingRows(sortOrder);
    sortedIdentifier = stableIdentifier(sortedSourceRows);
    isFirstInGroup = [true; ...
        sortedIdentifier(2:end) ~= sortedIdentifier(1:end-1)];
    representativeRows = sortedSourceRows(isFirstInGroup);
    isRepresentative(representativeRows) = true;

    [foundRepresentative, representativeLocation] = ismember( ...
        stableIdentifier(nonmissingRows), ...
        stableIdentifier(representativeRows));
    assert(all(foundRepresentative), ...
        'Cannot map stable identifiers to deduplication representatives.');
    representativeSourceRow(nonmissingRows) = ...
        representativeRows(representativeLocation);
else
    usedTieBreakColumns = strings(0, 1);
end

isEligibleUnique = isRepresentative & ...
    ~isMissingIdentifier & ~isReferenceOverlap;
exclusionReason = repmat( ...
    "EligibleUniqueRepresentative", candidateCount, 1);
exclusionReason(isMissingIdentifier) = ...
    "MissingStableIdentifierExcluded";
exclusionReason(isReferenceOverlap) = ...
    "ReferenceOverlapExcluded";
exclusionReason(~isMissingIdentifier & ...
    ~isReferenceOverlap & ~isRepresentative) = ...
    "DuplicateNonRepresentativeExcluded";

originalTitle = localGetStringColumn( ...
    manifest, options.OriginalTitleColumn, candidateCount);
externalId = localGetStringColumn( ...
    manifest, options.ExternalIdColumn, candidateCount);
canonicalSmiles = localGetStringColumn( ...
    manifest, options.SmilesColumn, candidateCount);
inputRowLabel = strings(candidateCount, 1);
if ismember("GlobalRow", variableNames)
    inputRowLabel = string(manifest.("GlobalRow"));
end

auditTable = table( ...
    sourceRow, dragonName, originalTitle, externalId, ...
    stableIdentifier, canonicalSmiles, inputRowLabel, ...
    identifierGroupSize, representativeSourceRow, ...
    isRepresentative, isMissingIdentifier, ...
    isReferenceOverlap, isEligibleUnique, exclusionReason, ...
    'VariableNames', { ...
    'SourceRow', 'DragonName', 'OriginalTitle', 'ExternalID', ...
    'StableIdentifier', 'CanonicalSMILES', 'InputRowLabel', ...
    'IdentifierGroupSize', 'RepresentativeSourceRow', ...
    'IsRepresentative', 'IsMissingIdentifier', ...
    'IsReferenceOverlap', 'IsEligibleUnique', ...
    'DeduplicationStatus'});

summary = struct();
summary.Version = 1;
summary.Method = ...
    'stable-identifier global deduplication before partitioning';
summary.ManifestFile = char(manifestFile);
summary.IdentifierColumn = char(options.IdentifierColumn);
summary.IdentifierNormalization = 'trim + uppercase';
summary.MissingIdentifierPolicy = char(missingPolicy);
summary.TieBreakColumnsUsed = cellstr(usedTieBreakColumns);
summary.SourceRowRole = [ ...
    'technical descriptor-file locator only; not the molecular ' ...
    'identity or duplicate criterion'];
summary.CandidateRows = candidateCount;
summary.UniqueCandidateIdentifierCount = ...
    uniqueCandidateIdentifierCount;
summary.DuplicateIdentifierGroupCount = duplicateGroupCount;
summary.DuplicateNonRepresentativeCount = nnz( ...
    exclusionReason == "DuplicateNonRepresentativeExcluded");
summary.MissingIdentifierExcludedCount = nnz(isMissingIdentifier);
summary.ReferenceIdentifierCount = numel(referenceIdentifier);
summary.ReferenceOverlapExcludedCount = nnz(isReferenceOverlap);
summary.EligibleUniqueCount = nnz(isEligibleUnique);
summary.CreatedAt = char(datetime( ...
    'now', 'TimeZone', 'local', ...
    'Format', 'yyyy-MM-dd''T''HH:mm:ssXXX'));
end

function value = localGetStringColumn(inputTable, columnName, rowCount)

if strlength(columnName) > 0 && ...
        ismember(columnName, string(inputTable.Properties.VariableNames))
    value = string(inputTable.(char(columnName)));
else
    value = strings(rowCount, 1);
end
value = value(:);
value(ismissing(value)) = "";
end

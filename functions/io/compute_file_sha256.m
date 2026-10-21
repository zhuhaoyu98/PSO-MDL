function hexadecimalDigest = compute_file_sha256(filePath)
%COMPUTE_FILE_SHA256 Hash a file in 8 MiB blocks using SHA-256.
%   FILEPATH is a scalar string. The output is a lowercase hexadecimal string.
%   File bytes are read by MATLAB and reinterpreted as signed Java bytes.

arguments
    filePath (1, 1) string
end

assert(isfile(filePath), 'File for SHA-256 calculation not found: %s', filePath);
messageDigest = java.security.MessageDigest.getInstance('SHA-256');
fileIdentifier = fopen(filePath, 'rb');
assert(fileIdentifier >= 0, 'Cannot read file for hashing: %s', filePath);
cleanupObject = onCleanup(@() fclose(fileIdentifier));
blockBytes = 8 * 1024 * 1024;

while ~feof(fileIdentifier)
    unsignedBlock = fread(fileIdentifier, blockBytes, '*uint8');
    if isempty(unsignedBlock)
        break;
    end

    signedBlock = typecast(unsignedBlock, 'int8');
    messageDigest.update(signedBlock, 0, numel(signedBlock));
end

rawDigest = int8(messageDigest.digest());
unsignedDigest = typecast(rawDigest, 'uint8');
hexadecimalDigest = lower(string(reshape(dec2hex( ...
    unsignedDigest, 2).', 1, [])));
clear cleanupObject
end

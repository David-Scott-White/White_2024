function createDirectoryIfNotExists(directoryPath)

% thanks chatGPT!

    % Check if the directory exists
    if exist(directoryPath, 'dir')
        disp(['Directory "', directoryPath, '" already exists.'])
    else
        % Create the directory
        mkdir(directoryPath);
        disp(['Directory "', directoryPath, '" created successfully.'])
    end
end
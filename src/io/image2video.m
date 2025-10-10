function image2video(imageFolder)
    % Collect image files (supports PNG and JPG)
    imageFiles = dir(fullfile(imageFolder, '*.png'));
    if isempty(imageFiles)
        imageFiles = dir(fullfile(imageFolder, '*.jpg'));
    end
    if isempty(imageFiles)
        error('No PNG or JPG images found in the specified folder.');
    end
    imageNames = {imageFiles.name};
    numFrames = numel(imageNames);

    % Read first image to get size
    sampleImg = imread(fullfile(imageFolder, imageNames{1}));
    [imgHeight, imgWidth, channels] = size(sampleImg);

    % User selects frame rate
    while true
        frameRate = input('Enter desired frame rate (frames/sec, e.g., 24, 30, 60): ');
        if isnumeric(frameRate) && frameRate > 0 && frameRate <= 120
            break;
        else
            disp('Please enter a valid frame rate between 1 and 120.');
        end
    end

    quality = 75; % Default quality
    while true
        % Estimate video size (approximate, for user decision)
        durationSec = numFrames / frameRate;
        % Simple estimate (uncompressed): size = frames * width * height * channels (bytes)
        bitsPerPixel = 24; % RGB, 8 bits per channel
        estimatedUncompressedSize = imgHeight * imgWidth * channels * numFrames / (1024^2); % MB
        % Crude compression factor based on quality (higher quality = smaller factor)
        compressionFactor = (100 - quality) / 100 * 15 + 1.5; % simple heuristic
        estimatedSizeMB = estimatedUncompressedSize / compressionFactor;

        fprintf('Proposed video:\n');
        fprintf(' - Resolution: %dx%d\n', imgWidth, imgHeight);
        fprintf(' - Frames: %d\n', numFrames);
        fprintf(' - Frame rate: %.1f fps\n', frameRate);
        fprintf(' - Quality: %d\n', quality);
        fprintf(' - Estimated size: %.2f MB\n', estimatedSizeMB);

        reply = input('Accept settings? [y/n, or enter new quality (0-100)]: ', 's');
        if strcmpi(reply, 'y')
            break;
        elseif strcmpi(reply, 'n')
            quality = input('Enter new quality (0-100): ');
        else
            qval = str2double(reply);
            if ~isnan(qval) && qval >= 0 && qval <= 100
                quality = qval;
            end
        end
    end

    % Create video file
    outVid = VideoWriter(fullfile(imageFolder, 'output_video.mp4'), 'MPEG-4');
    outVid.Quality = quality;
    outVid.FrameRate = frameRate;
    open(outVid);

    for i = 1:numFrames
        img = imread(fullfile(imageFolder, imageNames{i}));
        % Convert uint16 or uint32 to uint8 by scaling to 0-255
        if isa(img, 'uint16')
            img = uint8(double(img) / double(intmax('uint16')) * 255);
        elseif isa(img, 'uint32')
            img = uint8(double(img) / double(intmax('uint32')) * 255);
        end
        % Also handle grayscale images by replicating channels if needed
        if size(img,3) == 1
            img = repmat(img, [1 1 3]);
        end
        writeVideo(outVid, img);
    end

    close(outVid);
    disp('Video created successfully.');
end

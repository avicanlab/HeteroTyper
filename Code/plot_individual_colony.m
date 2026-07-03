%% 11.11.2025 - HeteroTyper Pipeline for Bright Plates

function plot_individual_colony(data)

    % Define Room Temperature Incubation Time
    incTime = 20;

    %% Interactive plotting of a single colony with rainbow segmentation and green border
    
    % Ask user for plate index
    plate_id = inputdlg('Plate index (1-N)?');
    plate_id = str2double(plate_id{1});
    
    if plate_id < 1 || plate_id > length(data.processed)
        error('Invalid plate index.');
    end
    
    % Ask user for colony number
    col_id = inputdlg('Colony number?');
    col_id = str2double(col_id{1});
    
    % Load plate data
    raw_image    = data.processed{plate_id}.img_final;
    mask_clean   = data.processed{plate_id}.colonies.debug.binary;
    region_props = data.processed{plate_id}.colonies.region_props;
    
    % Determine number of colonies
    % Determine number of colonies
    if istable(region_props)
        n_colonies = height(region_props);   % <-- use height for table
    elseif isstruct(region_props)
        n_colonies = numel(region_props);
    else
        error('Unexpected region_props format.');
    end
    
    % Check if colony index is valid
    if col_id < 1 || col_id > n_colonies
        error('Invalid colony number for this plate.');
    end
    
    
    % Rainbow segmentation mask
    if islogical(mask_clean)
        labeled_mask = bwlabel(mask_clean);
    else
        labeled_mask = mask_clean;
    end
    mask_rgb = label2rgb(labeled_mask,'jet','k','shuffle');
    
    % Get colony centroid and metrics
    if istable(region_props)
        centroid     = round(region_props.Centroid(col_id,:));
        area_val     = region_props.Area(col_id);
        ecc_val      = region_props.Eccentricity(col_id);
        circ_val     = region_props.Circularity(col_id);
        sol_val      = region_props.Solidity(col_id);
    elseif isstruct(region_props)
        centroid     = round(region_props(col_id).Centroid);
        area_val     = region_props(col_id).Area;
        ecc_val      = region_props(col_id).Eccentricity;
        circ_val     = region_props(col_id).Circularity;
        sol_val      = region_props(col_id).Solidity;
    end

    lagT_val  = colonies.lag_time(col_id);
    int_val   = colonies.timecourse_intensity_smoothed(end, col_id);
    size_val  = colonies.timecourse_size_smoothed(end, col_id);
    intS_val  = int_val ./ size_val;
    
    x = centroid(1); y = centroid(2);
    
    % Crop around colony
    [img_h, img_w, ~] = size(raw_image);
    crop_size = 150;
    x_min = max(1, x - crop_size); x_max = min(img_w, x + crop_size);
    y_min = max(1, y - crop_size); y_max = min(img_h, y + crop_size);
    
    cropped_raw  = raw_image(y_min:y_max, x_min:x_max, :);
    cropped_mask = mask_rgb(y_min:y_max, x_min:x_max, :);
    
    % Prepare grayscale base
    if size(cropped_raw,3) == 3
        gray_img = rgb2gray(cropped_raw);
    else
        gray_img = cropped_raw;
    end
    gray_img = imadjust(gray_img,[0.025 0.15],[]);
    gray_img_rgb = repmat(gray_img,[1 1 3]);
    
    % Green outline for the selected colony
    colony_mask = (labeled_mask == col_id);
    colony_crop_mask = colony_mask(y_min:y_max, x_min:x_max);
    perim = bwperim(colony_crop_mask);
    thick_border = imdilate(perim, strel('disk',1));
    
    r = gray_img_rgb(:,:,1);
    g = gray_img_rgb(:,:,2);
    b = gray_img_rgb(:,:,3);
    r(thick_border) = 0;
    g(thick_border) = 255;
    b(thick_border) = 0;
    gray_img_rgb(:,:,1) = r;
    gray_img_rgb(:,:,2) = g;
    gray_img_rgb(:,:,3) = b;
    
    % Blend with rainbow mask
    blended_img = imfuse(gray_img_rgb, cropped_mask, 'blend');
    
    % Add colony metrics
    text_lines = {
        sprintf('Lag=%.1f', lagT_val)
        sprintf('Size=%.2f', size_val)
        sprintf('Area=%.1f', area_val)
        sprintf('Int=%.1f', int_val)
        sprintf('IntSize=%.1f', intS_val)
        sprintf('Ecc=%.2f', ecc_val)
        sprintf('Cir=%.2f', circ_val)
        sprintf('Sol=%.2f', sol_val)
    };
    for t = 1:numel(text_lines)
        blended_img = insertText(blended_img, ...
            [size(blended_img,2)-5, 5 + 12*(t-1)], text_lines{t}, ...
            'FontSize', 10, 'TextColor', 'white', ...
            'BoxColor', 'black', 'BoxOpacity', 0, ...
            'AnchorPoint', 'RightTop');
    end
    
    % Show the result
    figure('Name',sprintf('Plate %d - Colony %d', plate_id,col_id));
    imshow(blended_img);

end

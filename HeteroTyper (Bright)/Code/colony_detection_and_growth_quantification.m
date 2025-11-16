%% 11.11.2025 - HeteroTyper Pipeline for Bright Plates

function colonies = colony_detection_and_growth_quantification(inp, position_name, fn, params, plot_flag)

    %% Define Border Radius
    b_radius = 0.85;

    % --- Ensure fn is a cell array of filenames ---
    if ischar(fn) || isstring(fn), fn = {fn};  end

    % === Extract time info directly from filenames ===
    try
        time_info = extract_time_info(fn);
        time_h = time_info.elapsed_time_h(:);
    catch
        warning('extract_time_info failed — using frame index as time.');
        time_h = (0:(numel(inp.img)-1))';
        time_info.elapsed_time_h = time_h;
    end

    n_frames = numel(inp.img);
    n_time   = numel(time_h);

    if n_time ~= n_frames
        warning('Time vector length (%d) != number of frames (%d). Trimming to min length.', n_time, n_frames);
        L = min(n_time, n_frames);
        time_h = time_h(1:L);
    else
        L = n_frames;
    end

    % ---------- Parameters ----------
    LoG_threshold  = params.LoG_thresh;
    size_threshold = params.size_thresh;
    ecc_threshold  = params.eccentricity_thresh;

    fprintf('colony_detection... %s (frames=%d, timepoints=%d used)\n', position_name, n_frames, L);

    % ---------- Use last frame to segment colonies ----------
    last_img = inp.img{end};
    if size(last_img,3) == 3
        last_img_gray = rgb2gray(last_img);
    else
        last_img_gray = last_img;
    end

    %% --- Detect colonies and get segmentation mask ---
    [mask, debug] = detect_colonies_in_img(last_img_gray, LoG_threshold, size_threshold, plot_flag);

    % if plot_flag == 1
    %     mask_plot = label2rgb(mask,'jet','k','shuffle');
    %     figure('Name', strcat(position_name,'-colony detection'));
    %     imshowpair(last_img_gray, mask_plot, 'blend');
    % end

    %% --- Ensure mask is labeled ---
    if islogical(mask), mask = bwlabel(mask); end
    if ~isequal(size(mask), size(last_img_gray))
        error('Mask and intensity image have different sizes!');
    end

    %% --- Compute region properties ---
    geom_props = regionprops('table', mask, last_img_gray, ...
        'PixelIdxList', 'Area', 'Eccentricity', 'Centroid', 'Circularity', ...
        'ConvexArea', 'ConvexHull', 'EulerNumber', 'MinorAxisLength', 'MajorAxisLength', ...
        'Perimeter', 'Solidity', 'MeanIntensity', 'MinIntensity', 'MaxIntensity');

    
    %% --- Direct pixel count (explicit area confirmation) ---
    labeledIDs = unique(mask);
    labeledIDs(labeledIDs == 0) = [];
    pixel_counts = zeros(numel(labeledIDs),1);
    for k = 1:numel(labeledIDs)
        pixel_counts(k) = nnz(mask == labeledIDs(k));
    end
    geom_props.Size_inPixels = pixel_counts;

    %% --- Border exclusion logic (using 85% radius) ---
    [H,W] = size(mask);
    [X,Y] = meshgrid(1:W,1:H);
    cx = W/2; cy = H/2;
    radius = min(cx,cy) * b_radius;
    dist_from_center = sqrt((X-cx).^2 + (Y-cy).^2);
    mask_valid = dist_from_center <= radius;

    border_touch = false(height(geom_props),1);
    for r = 1:height(geom_props)
        px = geom_props.PixelIdxList{r};
        if any(~mask_valid(px))
            border_touch(r) = true;
        end
    end

    %% --- Apply filters (size, eccentricity, border) ---
    pass_size   = geom_props.Area >= size_threshold;
    pass_shape  = geom_props.Eccentricity <= ecc_threshold;
    pass_border = ~border_touch;
    colony_ok   = pass_size & pass_shape & pass_border;

    % Keep label mapping consistent
    props_labels = labeledIDs(:);          % regionprops rows correspond to label IDs ascending
    ok_labels    = props_labels(colony_ok);

    mask_filtered = ismember(mask, ok_labels);

    % ---------- Build growth matrices over the FULL stack ----------
    nr_colonies = height(geom_props);
    T = L;                                  % effective number of usable frames
    sum_intensity = zeros(T, nr_colonies);
    colony_size   = zeros(T, nr_colonies);

    for t = 1:T
        img_t = inp.img{t};
        if size(img_t,3) == 3, img_t = rgb2gray(img_t); end

        % Foreground for size (same threshold logic as detection)
        bw_t = imbinarize(img_t, LoG_threshold);

        for j = 1:nr_colonies
            pix = geom_props.PixelIdxList{j};
            if isempty(pix), continue; end
            sum_intensity(t,j) = sum(img_t(pix));
            colony_size(t,j)   = sum(bw_t(pix));    % # of foreground pixels
        end
    end

    % Smooth along time
    sum_intensity_sm = medfilt1(sum_intensity(1:T,:), 3, [], 1);
    colony_size_sm   = medfilt1(colony_size(1:T,:),   3, [], 1);

    % ---------- Lag time, AUC, final size (use smoothed size) ----------
    lag_thr = params.lag_time_thresh;
    lag_ix  = zeros(1, nr_colonies);
    AUC     = zeros(1, nr_colonies);
    finalSize = zeros(1, nr_colonies);

    for j = 1:nr_colonies
        hit = find(colony_size_sm(:,j) > lag_thr, 1, 'first');
        if ~isempty(hit)
            lag_ix(j) = hit;
        else
            lag_ix(j) = T; % never crossed: assign last time
        end
        AUC(j)     = trapz(time_h, colony_size_sm(:,j));
        finalSize(j) = colony_size_sm(T,j);
    end

    % ---------- Save outputs ----------
    colonies.mask                  = mask;
    colonies.mask_filtered         = mask_filtered;
    colonies.debug                 = debug;

    geom_props.flag_colony_ok      = colony_ok(:);
    colonies.region_props          = geom_props;
    colonies.flag_colony_ok        = colony_ok(:);

    colonies.new = struct( ...
        'time_info',                      struct('elapsed_time_h', time_h(:)), ...
        'timecourse_intensity',           sum_intensity, ...
        'timecourse_intensity_smoothed',  sum_intensity_sm, ...
        'timecourse_size',                colony_size, ...
        'timecourse_size_smoothed',       colony_size_sm, ...
        'lag_time',                       time_h(lag_ix), ...
        'AUC',                            AUC, ...
        'final_col_size',                 finalSize);

    % --- Legacy aliases (keep old callers happy) ---
    colonies.time_info                       = colonies.new.time_info;
    colonies.timecourse_intensity            = colonies.new.timecourse_intensity;
    colonies.timecourse_intensity_smoothed   = colonies.new.timecourse_intensity_smoothed;
    colonies.timecourse_size                 = colonies.new.timecourse_size;
    colonies.timecourse_size_smoothed        = colonies.new.timecourse_size_smoothed;
    colonies.lag_time                        = colonies.new.lag_time;
    colonies.AUC                             = colonies.new.AUC;
    colonies.final_col_size                  = colonies.new.final_col_size;

    % Book-keeping / logging
    colonies.control_check.detected_colony_count  = max(mask(:));
    colonies.control_check.accepted_colony_count  = sum(colony_ok);
    colonies.control_check.regionprops_count      = height(geom_props);
    colonies.control_check.timestamp              = datestr(now);

    if plot_flag == 1
        mask_rgb = label2rgb(mask,'jet','k','shuffle');
        figure('Name', strcat(position_name,' - colony detection')); imshowpair(last_img_gray, mask_rgb, 'blend');
    end

    fprintf('[%s] colonies: detected=%d, accepted=%d\n', position_name, colonies.control_check.detected_colony_count, sum(colony_ok));

end



%% === Helper: detect colonies in single image ===
function [out, debug] = detect_colonies_in_img(img_in, LoG_threshold, size_threshold, plot_flag)

    fig_w = 1800;
    fig_h = 400;

    if size(img_in,3) == 3
        img_gray = rgb2gray(img_in);
    else
        img_gray = img_in;
    end

    %% Background flattening
    background = imgaussfilt(img_gray, 30); 
    img_flat = imsubtract(img_gray, background);
    img_norm = mat2gray(img_flat);

    img_contrast = imadjust(img_norm, [], [], 1); 
    % imadjust, adjusts the intensity values of the image to improve contrast
    % imadjust(img_norm, [], [], gamma); 
    % gamma = 1, (normalization placeholder) no gamma correction
    % gamma < 1, brightening mid-tones
    % gamma > 1, darkening mid-tones
    
    img_smooth = imgaussfilt(img_contrast, 1);
    % imgaussfilt(X, Y), applies a Gaussian blur to smooth the image (to reduce noise and small intensity variations)
    % Y, is the standard deviation (σ) of the Gaussian kernel
    % Smaller σ (e.g., 1) → light smoothing (fine details mostly preserved)
    % Larger σ (e.g., 5–10) → stronger blur (more detail lost, but noise reduced)

   
    % Thresholding + fill
    bw = imbinarize(img_smooth, LoG_threshold);
    bw = imfill(bw, 'holes');
    bw = bwareaopen(bw, size_threshold);

    %% Watershed segmentation
    D = imcomplement(bwdist(~bw));
    D2 = imhmin(D,1);
    L = watershed(D2);
    bw(~L) = 0;

    % Exclude border colonies (85% radius)
    [H,W] = size(bw);
    [X,Y] = meshgrid(1:W,1:H);
    cx = W/2; cy = H/2;
    radius = min(cx,cy) * b_radius;
    dist_from_center = sqrt((X-cx).^2 + (Y-cy).^2);
    mask_valid = dist_from_center <= radius;

    labeled = bwlabel(bw);
    stats = regionprops(labeled, 'PixelIdxList');
    remove_idx = false(1,numel(stats));
    for i = 1:numel(stats)
        if any(~mask_valid(stats(i).PixelIdxList))
            remove_idx(i) = true;
        end
    end
    keep_idx = find(~remove_idx);
    bw = ismember(labeled, keep_idx);

    % Final labeled mask
    cc = bwconncomp(bw);
    out = labelmatrix(cc);

    % Debug images and save into the data under "debug" parameter
    debug.gray = img_gray;
    debug.background = background;
    debug.flat = img_flat;
    debug.norm = img_norm;
    debug.contrast = img_contrast;
    debug.smooth = img_smooth;
    debug.binary = bw;
    debug.segmented = label2rgb(out,'jet','k','shuffle');
    debug.mask_clean = out;

    if plot_flag == 1
        images = {debug.gray, debug.background, debug.flat, debug.norm, ...
                  debug.contrast, debug.smooth, debug.binary, debug.segmented};
        titles = {'1. Grayscale','2. Background','3. Flattened','4. Normalized', ...
                  '5. Contrasted','6. Smoothed','7. Binary','8. Final'};
        figure('Name','Colony Detection Debug','Position',[100 100 fig_w fig_h]);
        for k = 1:8
            subplot(2,4,k);
            imshow(images{k},[]);
            hold on;
            if k >= 7, viscircles([W/2, H/2], radius, 'Color','w','LineWidth',0.3); end
            title(titles{k});
        end
    end
end
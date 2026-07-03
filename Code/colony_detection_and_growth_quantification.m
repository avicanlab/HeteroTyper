%% 11.11.2025 - HeteroTyper Pipeline for Bright Plates
%%
%% SINGLE UNIFIED PIPELINE — no redundant filtering:
%%
%%  Stage 1 [detect_colonies_in_img]:
%%    Background subtraction → binarize → split merged → raw labeled mask
%%    NO border removal here. Just detect everything.
%%
%%  Stage 2 [this function — filters]:
%%    From raw mask, compute regionprops, then apply filters:
%%      F1:  size   >= size_thresh
%%      F2:  eccentricity <= ecc_thresh  (removes extreme elongated noise)
%%      F2b: near-border annular zone (border_zone_inner to b_radius):
%%           stricter eccentricity = params.border_ecc_thresh (default 0.55)
%%           zone inner edge = params.border_zone_inner (default 0.75)
%%      F3:  ALL pixels inside the b_radius plate circle  (border rule)
%%    → mask_filtered: result after F1+F2+F2b+F3
%%
%%  Stage 3: mask_filtered used directly for growth quantification.
%%
%%  Figures (NEW ORDER):
%%    _1_debug    : 8-panel preprocessing steps (raw mask in panel 8)
%%    _2_detection: raw mask blended with gray image (all detected, pre-filter)
%%    _3_cleanup  : white=kept  green=removed  (swapped from old step 4)
%%    _4_segmented: mask_filtered overlaid on gray image  (swapped from old step 3)
%%
%%  Merged-colony splitting (tunable):
%%    params.split_h_thresh  — minimum valley depth in the distance transform
%%                             between two adjacent colony peaks to trigger a
%%                             split. Seeds placed at each colony's thickest
%%                             interior point (medial axis), not at centroids —
%%                             so heterogeneous colony sizes are handled fairly.
%%                             Increase to split less (default 2).

function colonies = colony_detection_and_growth_quantification(inp, position_name, fn, params, plot_flag)

    %% ---- Parameter extraction -----------------------------------------------
    b_radius = 0.85;
    if isfield(params, 'border_range') && ~isempty(params.border_range)
        b_radius = params.border_range;
    end

    border_zone_inner = 0.75;
    if isfield(params, 'border_zone_inner') && ~isempty(params.border_zone_inner)
        border_zone_inner = params.border_zone_inner;
    end

    border_ecc_thresh = 0.55;
    if isfield(params, 'border_ecc_thresh') && ~isempty(params.border_ecc_thresh)
        border_ecc_thresh = params.border_ecc_thresh;
    end

    % split_h_thresh: height of the valley in the distance transform that
    % must exist between two peaks for them to be treated as separate colonies.
    % Purely topology-driven — works correctly when colony sizes are heterogeneous
    % because each connected object seeds from its own distance-transform ridge,
    % not from a centroid neighbourhood. Increase to split less. Default 2.
    split_h_thresh = 2;
    if isfield(params, 'split_h_thresh') && ~isempty(params.split_h_thresh)
        split_h_thresh = params.split_h_thresh;
    end

    if ischar(fn) || isstring(fn), fn = {fn}; end

    %% ---- Time info ----------------------------------------------------------
    try
        time_info = extract_time_info(fn);
        time_h    = time_info.elapsed_time_h(:);
    catch
        warning('extract_time_info failed — using frame index as time.');
        time_h = (0:(numel(inp.img)-1))';
    end

    n_frames = numel(inp.img);
    n_time   = numel(time_h);
    if n_time ~= n_frames
        warning('Time vector (%d) != frames (%d). Trimming.', n_time, n_frames);
        time_h = time_h(1:min(n_time,n_frames));
    end
    L = numel(time_h);

    LoG_threshold  = params.LoG_thresh;
    size_threshold = params.size_thresh;
    ecc_threshold  = params.eccentricity_thresh;

    fprintf('colony_detection... %s (frames=%d, timepoints=%d)\n', position_name, n_frames, L);

    %% ---- Last frame to grayscale --------------------------------------------
    last_img = inp.img{end};
    if size(last_img,3) == 3
        last_img_gray = rgb2gray(last_img);
    else
        last_img_gray = last_img;
    end

    %% -----------------------------------------------------------------------
    %% Stage 1: Raw detection + merged-colony splitting
    %% -----------------------------------------------------------------------
    [mask_raw, debug] = detect_colonies_in_img(last_img_gray, LoG_threshold, ...
                                               size_threshold, split_h_thresh, params);
    if islogical(mask_raw), mask_raw = bwlabel(mask_raw); end

    %% -----------------------------------------------------------------------
    %% Stage 2: regionprops on raw mask
    %% -----------------------------------------------------------------------
    geom_props = regionprops('table', mask_raw, last_img_gray, ...
        'PixelIdxList', 'Area', 'Eccentricity', 'Centroid', 'Circularity', ...
        'ConvexArea', 'ConvexHull', 'EulerNumber', 'MinorAxisLength', ...
        'MajorAxisLength', 'Perimeter', 'Solidity', ...
        'MeanIntensity', 'MinIntensity', 'MaxIntensity');

    labeledIDs = unique(mask_raw);
    labeledIDs(labeledIDs == 0) = [];
    pixel_counts = arrayfun(@(id) nnz(mask_raw == id), labeledIDs);
    geom_props.Size_inPixels = pixel_counts(:);

    %% -----------------------------------------------------------------------
    %% Stage 2: Filtering (F1, F2, F2b near-border zone, F3)
    %% -----------------------------------------------------------------------
    [H, W] = size(mask_raw);
    [X, Y] = meshgrid(1:W, 1:H);

    if isfield(params,'plate_center_current') && ~any(isnan(params.plate_center_current))
        cx          = params.plate_center_current(1);
        cy          = params.plate_center_current(2);
        plate_r_raw = params.plate_r_current;
    else
        cx          = W/2;
        cy          = H/2;
        plate_r_raw = min(cx,cy);
    end
    radius = plate_r_raw * b_radius;

    dist_from_center = sqrt((X - cx).^2 + (Y - cy).^2);
    circle_pxl       = dist_from_center <= radius;   % accepted plate disk

    n_reg = height(geom_props);
    pass_size            = geom_props.Area >= size_threshold;
    pass_ecc             = geom_props.Eccentricity <= ecc_threshold;
    pass_border          = true(n_reg, 1);
    pass_near_border_ecc = true(n_reg, 1);

    for k = 1:n_reg
        px = geom_props.PixelIdxList{k};

        % F3: every pixel must be inside the plate circle
        if any(~circle_pxl(px))
            pass_border(k) = false;
            continue;
        end

        % F2b: centroid in near-border annular zone → stricter eccentricity
        cx_col = geom_props.Centroid(k,1);
        cy_col = geom_props.Centroid(k,2);
        d_col  = sqrt((cx_col - cx)^2 + (cy_col - cy)^2);
        if d_col > (plate_r_raw * border_zone_inner)
            if geom_props.Eccentricity(k) > border_ecc_thresh
                pass_near_border_ecc(k) = false;
            end
        end
    end

    colony_ok = pass_size & pass_ecc & pass_near_border_ecc & pass_border;

    % Diagnostic: which filter removes how many colonies (exclusive cascade)
    fail_size_only        = ~pass_size;
    fail_ecc_only         = pass_size  & ~pass_ecc;
    fail_near_border_only = pass_size  &  pass_ecc & ~pass_near_border_ecc;
    fail_border_only      = pass_size  &  pass_ecc &  pass_near_border_ecc & ~pass_border;
    fprintf('  [filter] detected=%d\n', n_reg);
    fprintf('    removed by F1 size        : %d  (area < %d px)\n',               sum(fail_size_only),        size_threshold);
    fprintf('    removed by F2 eccentricity: %d  (ecc > %.2f)\n',                 sum(fail_ecc_only),         ecc_threshold);
    fprintf('    removed by F2b near-border: %d  (ecc > %.2f in annular zone)\n', sum(fail_near_border_only), border_ecc_thresh);
    fprintf('    removed by F3 border      : %d  (centroid/pixels outside plate)\n', sum(fail_border_only));
    fprintf('    final kept                : %d\n', sum(colony_ok));

    props_labels  = labeledIDs(:);
    ok_labels     = props_labels(colony_ok);
    mask_filtered = double(mask_raw) .* double(ismember(mask_raw, ok_labels));

    clean_ids = find(colony_ok);
    nr_clean  = numel(clean_ids);

    %% -----------------------------------------------------------------------
    %% Stage 3: Build growth matrices — filtered colonies only
    %% -----------------------------------------------------------------------
    T             = L;
    sum_intensity = zeros(T, nr_clean);
    colony_size   = zeros(T, nr_clean);

    for t = 1:T
        img_t = inp.img{t};
        if size(img_t,3) == 3, img_t = rgb2gray(img_t); end
        bw_t = imbinarize(img_t, LoG_threshold);
        for j = 1:nr_clean
            pix = geom_props.PixelIdxList{clean_ids(j)};
            if isempty(pix), continue; end
            sum_intensity(t,j) = sum(img_t(pix));
            colony_size(t,j)   = sum(bw_t(pix));
        end
    end

    sum_intensity_sm = medfilt1(sum_intensity, 3, [], 1);
    colony_size_sm   = medfilt1(colony_size,   3, [], 1);

    lag_thr = params.lag_time_thresh;
    lag_ix  = zeros(1, nr_clean);
    AUC     = zeros(1, nr_clean);
    finalSz = zeros(1, nr_clean);
    for j = 1:nr_clean
        hit = find(colony_size_sm(:,j) > lag_thr, 1, 'first');
        if ~isempty(hit), lag_ix(j) = hit; else, lag_ix(j) = T; end
        AUC(j)     = trapz(time_h, colony_size_sm(:,j));
        finalSz(j) = colony_size_sm(T,j);
    end

    %% -----------------------------------------------------------------------
    %% Save outputs
    %% -----------------------------------------------------------------------
    colonies.mask               = mask_raw;
    colonies.mask_filtered      = mask_filtered;
    colonies.mask_clean         = mask_filtered;
    colonies.mask_clean_labeled = mask_filtered;
    colonies.debug              = debug;

    geom_props.flag_colony_ok   = colony_ok(:);
    colonies.region_props       = geom_props;
    colonies.region_props_clean = geom_props(clean_ids, :);
    colonies.flag_colony_ok     = colony_ok(:);

    colonies.new = struct( ...
        'time_info',                     struct('elapsed_time_h', time_h(:)), ...
        'timecourse_intensity',          sum_intensity, ...
        'timecourse_intensity_smoothed', sum_intensity_sm, ...
        'timecourse_size',               colony_size, ...
        'timecourse_size_smoothed',      colony_size_sm, ...
        'lag_time',                      time_h(lag_ix), ...
        'AUC',                           AUC, ...
        'final_col_size',                finalSz);

    colonies.time_info                     = colonies.new.time_info;
    colonies.timecourse_intensity          = colonies.new.timecourse_intensity;
    colonies.timecourse_intensity_smoothed = colonies.new.timecourse_intensity_smoothed;
    colonies.timecourse_size               = colonies.new.timecourse_size;
    colonies.timecourse_size_smoothed      = colonies.new.timecourse_size_smoothed;
    colonies.lag_time                      = colonies.new.lag_time;
    colonies.AUC                           = colonies.new.AUC;
    colonies.final_col_size                = colonies.new.final_col_size;

    colonies.control_check.detected_colony_count = max(mask_raw(:));
    colonies.control_check.accepted_colony_count = nr_clean;
    colonies.control_check.clean_colony_count    = nr_clean;
    colonies.control_check.regionprops_count     = n_reg;
    colonies.control_check.timestamp             = datestr(now);

    %% -----------------------------------------------------------------------
    %% Figures
    %%   _2_detection : raw detected regions blended with grayscale (pre-filter)
    %%   _3_cleanup   : white=kept  green=removed  (NOW step 3, was step 4)
    %%   _4_segmented : mask_filtered overlaid on gray  (NOW step 4, was step 3)
    %% -----------------------------------------------------------------------
    if plot_flag == 1

        % _2_detection --------------------------------------------------------
        figure('Name', strcat(position_name, ' - colony detection'));
        imshowpair(last_img_gray, label2rgb(mask_raw, 'jet', 'k', 'shuffle'), 'blend');
        title(sprintf('%s — detected: %d', position_name, max(mask_raw(:))), ...
            'Interpreter', 'none', 'FontSize', 18, 'FontWeight', 'bold');

        % _3_cleanup: colour-coded by which filter removed each colony
        %   white   = kept
        %   red     = F1: failed size  (area < size_thresh)
        %   magenta = F2: failed global eccentricity
        %   yellow  = F2b: failed near-border eccentricity
        %   green   = F3: failed border (centroid/pixels outside plate)
        [Hv, Wv] = size(mask_raw);
        vis = zeros(Hv, Wv, 3, 'uint8');

        % White: kept
        kept_mask = mask_filtered > 0;
        vis(:,:,1) = vis(:,:,1) + uint8(kept_mask)*255;
        vis(:,:,2) = vis(:,:,2) + uint8(kept_mask)*255;
        vis(:,:,3) = vis(:,:,3) + uint8(kept_mask)*255;

        % Colour-code each removed colony by the filter that caught it
        for k = 1:n_reg
            if colony_ok(k), continue; end
            pxm = uint8(mask_raw == props_labels(k));
            if     fail_size_only(k)                         % red
                vis(:,:,1) = vis(:,:,1) + pxm*220;
            elseif fail_ecc_only(k)                          % magenta
                vis(:,:,1) = vis(:,:,1) + pxm*220;
                vis(:,:,3) = vis(:,:,3) + pxm*220;
            elseif fail_near_border_only(k)                  % yellow
                vis(:,:,1) = vis(:,:,1) + pxm*220;
                vis(:,:,2) = vis(:,:,2) + pxm*220;
            elseif fail_border_only(k)                       % green
                vis(:,:,2) = vis(:,:,2) + pxm*200;
            end
        end

        figure('Name', strcat(position_name, ' - colony cleanup'));
        imshow(vis);
        title(sprintf('%s  white=kept | red=size  magenta=ecc  yellow=near-border  green=border  (%d→%d)', ...
            position_name, max(mask_raw(:)), nr_clean), 'Interpreter','none', 'FontSize', 18, 'FontWeight', 'bold');

        % _4_segmented (step 4 — was step 3) — gray plate + coloured overlay --
        % Convert grayscale to uint8 RGB background
        if isa(last_img_gray,'uint8')
            gray3 = repmat(last_img_gray, [1 1 3]);
        else
            gray3 = repmat(im2uint8(mat2gray(last_img_gray)), [1 1 3]);
        end
        col_overlay = label2rgb(uint16(mask_filtered), 'jet', [0 0 0], 'shuffle');
        col_mask    = repmat(mask_filtered > 0, [1 1 3]);
        blend_img   = gray3;
        blend_img(col_mask) = uint8( ...
            0.40 * double(gray3(col_mask)) + ...
            0.60 * double(col_overlay(col_mask)));
        figure('Name', strcat(position_name, ' - segmented'));
        imshow(blend_img);
        title(sprintf('%s — final colonies (%d)', position_name, nr_clean), ...
            'Interpreter', 'none', 'FontSize', 18, 'FontWeight', 'bold');
    end

    fprintf('[%s] detected=%d  final=%d\n', position_name, max(mask_raw(:)), nr_clean);
end


%% =========================================================================
%% detect_colonies_in_img
%% Stage 1: background subtraction → binarize → edge-topology watershed
%%          to split merged/touching colonies → label.
%% NO border removal — happens in the main function.
%%
%% Splitting strategy — topology-driven via imextendedmax (NOT erosion):
%%
%%   The distance transform D = bwdist(~bw) assigns each foreground pixel
%%   its distance to the nearest background pixel.  Each colony's interior
%%   forms a ridge in D whose peak height equals approximately the colony
%%   radius.  Where two colonies touch, D forms a saddle/valley near zero.
%%
%%   imextendedmax(D, split_h_thresh) finds every local maximum of D that
%%   rises at least split_h_thresh above its surrounding "collar" — a purely
%%   topological test that is INDEPENDENT of the spatial distance between
%%   peaks.  A small colony next to a large one retains its own seed as long
%%   as the valley between them is deep enough.  This is exactly what fails
%%   with erosion-based seeding, where the large colony's wider peak
%%   dominates and the small colony's seed is swallowed.
%%
%%   A per-component safety guarantee ensures every connected blob gets at
%%   least one seed even if split_h_thresh is set higher than any valley.
%%
%% params.split_h_thresh (default 2):
%%   Minimum "depth" of the valley in D between two peaks for a split to
%%   occur.  Geometrically: the contact neck between two colonies must have
%%   a distance-transform value at least split_h_thresh px lower than each
%%   colony's interior peak.
%%   • Increase (e.g. 4–8) → fewer, more conservative splits.
%%   • Decrease (e.g. 1)   → splits even lightly touching colonies.
%% =========================================================================
function [out, debug] = detect_colonies_in_img(img_in, LoG_threshold, size_threshold, split_h_thresh, params)

    if size(img_in,3) == 3
        img_gray = rgb2gray(img_in);
    else
        img_gray = img_in;
    end

    %% Background estimation
    if isfield(params,'plate_r_current') && ~isnan(params.plate_r_current) && params.plate_r_current > 0
        bg_sigma = round(params.plate_r_current / 3);
    else
        bg_sigma = 40;
    end
    bg_sigma = max(bg_sigma, 30);
    bg_sigma = min(bg_sigma, 120);
    fprintf('  [detect] bg_sigma=%d  split_h_thresh=%.1f\n', bg_sigma, split_h_thresh);

    background   = imgaussfilt(img_gray, bg_sigma);
    img_flat     = imsubtract(img_gray, background);
    img_norm     = mat2gray(img_flat);
    img_contrast = imadjust(img_norm, [], [], 1);
    img_smooth   = imgaussfilt(img_contrast, 1);

    %% Threshold → fill → area open
    bw = imbinarize(img_smooth, LoG_threshold);
    bw = imfill(bw, 'holes');
    bw = bwareaopen(bw, size_threshold);

    %% -----------------------------------------------------------------
    %% Edge-topology watershed — split merged/touching colonies
    %% -----------------------------------------------------------------
    % Distance transform: value = distance to nearest background pixel.
    % Peaks correspond to the thickest interior point of each colony,
    % i.e., the point farthest from the colony edge — NOT the centroid.
    D = bwdist(~bw);

    % Find extended maxima: connected regions of D that are local maxima
    % AND stand at least split_h_thresh above their surrounding terrain.
    % Unlike erosion, this respects each colony's own topology regardless
    % of how close or how differently sized neighbouring colonies are.
    h       = max(0.5, split_h_thresh);
    markers = imextendedmax(D, h) & bw;

    % Safety: guarantee at least one seed per connected component of bw.
    % If split_h_thresh is too large, some objects may receive no marker.
    cc_bw = bwconncomp(bw);
    for q = 1:cc_bw.NumObjects
        pix_q = cc_bw.PixelIdxList{q};
        if ~any(markers(pix_q))
            [~, peak_local] = max(D(pix_q));
            markers(pix_q(peak_local)) = true;
        end
    end

    % Marker-controlled watershed on the inverted distance transform.
    % Watershed lines (L_ws == 0) become the boundaries between split colonies.
    D_neg  = imcomplement(D);
    D_neg2 = imimposemin(D_neg, markers);
    L_ws   = watershed(D_neg2);

    bw_split = bw;
    bw_split(L_ws == 0) = 0;                         % cut at watershed ridges
    bw_split = bwareaopen(bw_split, size_threshold);  % discard tiny split fragments

    %% Label — no border removal here
    cc  = bwconncomp(bw_split);
    out = labelmatrix(cc);

    %% Debug struct
    debug.gray       = img_gray;
    debug.background = background;
    debug.flat       = img_flat;
    debug.norm       = img_norm;
    debug.contrast   = img_contrast;
    debug.smooth     = img_smooth;
    debug.binary     = bw_split;
    debug.segmented  = label2rgb(out, 'jet', 'k', 'shuffle');
    debug.mask_clean = out;

    %% Debug figure
    if isfield(params, 'plate_center_current') && ~any(isnan(params.plate_center_current))
        b_radius = 0.85;
        if isfield(params,'border_range'), b_radius = params.border_range; end
        cx_p   = params.plate_center_current(1);
        cy_p   = params.plate_center_current(2);
        radius = params.plate_r_current * b_radius;
    else
        [Hg, Wg] = size(img_gray);
        cx_p = Wg/2; cy_p = Hg/2;
        radius = min(cx_p,cy_p) * 0.85;
    end

    images = {debug.gray, debug.background, debug.flat, debug.norm, ...
              debug.contrast, debug.smooth, debug.binary, debug.segmented};
    titles = {'1. Grayscale','2. Background','3. Flattened','4. Normalized', ...
              '5. Contrasted','6. Smoothed','7. Binary (split)','8. Raw detected'};

    % ---- Debug figure: 3x3 grid, tile 9 empty ------------------------------
    % All spacing values are PERCENTAGES of the figure size (0-100).
    n_cols      = 3;   
    n_rows      = 3;
    pad_l_pct   = 1;   % left   figure margin  %
    pad_r_pct   = 1;   % right  figure margin  %
    pad_t_pct   = 1.5;   % top    figure margin  %
    pad_b_pct   = 1;   % bottom figure margin  %
    col_gap_pct = 1.5; % horizontal gap BETWEEN columns  %
    row_gap_pct = 1;   % vertical   gap BETWEEN rows     %

    % Convert to normalised [0-1] units for axes placement
    pad_l   = pad_l_pct   / 100;
    pad_r   = pad_r_pct   / 100;
    pad_t   = pad_t_pct   / 100;
    pad_b   = pad_b_pct   / 100;
    col_gap = col_gap_pct / 100;
    row_gap = row_gap_pct / 100;

    ax_w = (1 - pad_l - pad_r - (n_cols-1)*col_gap) / n_cols;
    ax_h = (1 - pad_t - pad_b - (n_rows-1)*row_gap) / n_rows;

    fig_dbg = figure('Name','Colony Detection Debug','Position',[50 50 1260 1260]);
    axes_dbg = gobjects(1, 8);
    for k = 1:8
        col_k = mod(k-1, n_cols);          % 0-based column index
        row_k = floor((k-1) / n_cols);     % 0-based row index (0 = top)
        x0 = pad_l + col_k * (ax_w + col_gap);
        y0 = 1 - pad_t - (row_k+1)*ax_h - row_k*row_gap;  % bottom-left y
        axes_dbg(k) = axes('Position', [x0, y0, ax_w, ax_h]); 
        imshow(images{k}, []);
        hold on;
        if k >= 7
            viscircles([cx_p, cy_p], radius, 'Color','w','LineWidth',0.5);
        end
        title(titles{k}, 'FontSize', 13, 'FontWeight', 'bold');
        axes_dbg(k).TitleFontSizeMultiplier = 1;
    end
    % tile 9 intentionally left blank
end
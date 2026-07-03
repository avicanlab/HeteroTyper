%% 11.11.2025 - HeteroTyper Pipeline for Bright Plates
%% Per-plate center detection: each plate gets its own median center computed
%% across a sample of its own image stack (not a global cross-plate median).
%%
%% Logic:
%%   1) For each plate, sample up to max_qc_frames images from its stack.
%%   2) Detect circle center + radius in every sampled frame via Hough.
%%   3) Compute per-plate median center and median radius from those frames.
%%   4) Flag each sampled frame: GREEN if within 10% radius deviation from
%%      the plate's own median, RED if outside.
%%   5) The plate's median center/radius is what gets used downstream for
%%      border masking in colony detection — NOT a cross-plate global.
%%   6) Summary figure shows one subplot per plate with:
%%        Blue  = per-plate median circle (what will be used)
%%        Green = individual frame detections within 10% of median
%%        Red   = individual frame detections deviating > 10%

function data = find_generic_plate_center(data)

fn        = data.metadata.fn;
nr_plates = data.params.nr_plates;

% How many frames to sample per plate for center estimation.
% Default 5 for speed; set to Inf to use all frames.
max_qc_frames = get_param(data.params, 'plate_center_max_qc_frames', 5);

% Deviation threshold: fraction of per-plate median radius
dev_thresh = get_param(data.params, 'plate_center_dev_thresh', 0.10);

img_display  = cell(nr_plates, 1);   % first frame image for display only
center_med   = nan(nr_plates, 2);    % per-plate median center  [cx cy]
radius_med   = nan(nr_plates, 1);    % per-plate median radius
frame_centers = cell(nr_plates, 1); % all sampled frame centers
frame_radii  = cell(nr_plates, 1);  % all sampled frame radii
frame_good   = cell(nr_plates, 1);  % logical: within dev_thresh of median

%% Per-plate loop
for i = 1:nr_plates

    t1     = strcat(fn(i).folder, '\', fn(i).name);
    fn_tmp = dir(t1);
    fn_tmp = fn_tmp(~[fn_tmp.isdir]);
    fn_tmp = fn_tmp(~startsWith({fn_tmp.name}, '.'));

    if isempty(fn_tmp)
        warning('[plate %d] No image files found. Skipping.', i);
        continue;
    end

    % Choose which frames to sample
    n_files   = numel(fn_tmp);
    frame_idx = choose_frames(n_files, max_qc_frames);
    n_sampled = numel(frame_idx);

    fprintf('\n[plate %d] %s — sampling %d/%d frames\n', ...
        i, fn(i).name, n_sampled, n_files);

    fc = nan(n_sampled, 2);
    fr = nan(n_sampled, 1);

    for q = 1:n_sampled
        fpath = fullfile(fn_tmp(frame_idx(q)).folder, fn_tmp(frame_idx(q)).name);
        im    = imread(fpath);

        % Store first frame for display
        if q == 1, img_display{i} = im; end

        [fc(q,:), fr(q)] = find_center_individual(im);
        fprintf('  frame %d: center=(%.1f,%.1f) r=%.1f\n', ...
            frame_idx(q), fc(q,1), fc(q,2), fr(q));
    end

    % Per-plate median center and radius
    cm = median(fc, 1, 'omitnan');
    rm = median(fr,    'omitnan');

    % Flag each sampled frame
    deviations = sqrt(sum((fc - cm).^2, 2)) ./ rm;
    fg = deviations <= dev_thresh;

    center_med(i,:)   = cm;
    radius_med(i)     = rm;
    frame_centers{i}  = fc;
    frame_radii{i}    = fr;
    frame_good{i}     = fg;

    n_good = sum(fg);
    fprintf('  [plate %d] median center=(%.1f,%.1f) r=%.1f | %d/%d frames within %.0f%%\n', ...
        i, cm(1), cm(2), rm, n_good, n_sampled, dev_thresh*100);
end

%% Store per-plate results in params
% center_median and radius are now [nr_plates x 2] and [nr_plates x 1]
data.params.center_median  = center_med;
data.params.radius         = radius_med;
data.params.plate_r        = radius_med;
data.params.frame_centers  = frame_centers;
data.params.frame_radii    = frame_radii;
data.params.frame_good     = frame_good;

%% Summary plot
% -------------------------------------------------------------------------
% Layout: tiledlayout with tight padding so images fill the figure.
% Title:  PosXXXX extracted from folder name (same regex as main script).
% Draw order: blue (median/consensus) first, then green or red ON TOP so
%             the quality colour is never hidden behind the blue circle.
% Colour rule: ALL frames ok  -> blue + green overlay
%              ANY frame bad  -> blue + red overlay
% -------------------------------------------------------------------------
[best_rows, best_cols] = best_subplot_grid(nr_plates);

% Manual axes positioning — independent control over column gap vs row gap.
% All spacing values are PERCENTAGES of the figure size (0-100).
pad_l_pct   = 1;   % left   figure margin  %
pad_r_pct   = 1;   % right  figure margin  %
pad_t_pct   = 1.5;   % top    figure margin  %  (extra room for bold titles)
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

ax_w = (1 - pad_l - pad_r - (best_cols-1)*col_gap) / best_cols;
ax_h = (1 - pad_t - pad_b - (best_rows-1)*row_gap) / best_rows;

tile_px  = 420;
fig_w_px = tile_px * best_cols;
fig_h_px = tile_px * best_rows;
fig = figure('Name', 'find plate centers', ...
             'Position', [50, 50, fig_w_px, fig_h_px]);

for i = 1:nr_plates
    col_i = mod(i-1, best_cols);
    row_i = floor((i-1) / best_cols);
    x0 = pad_l + col_i * (ax_w + col_gap);
    y0 = 1 - pad_t - (row_i+1)*ax_h - row_i*row_gap;
    axes('Position', [x0, y0, ax_w, ax_h]); 

    if isempty(img_display{i}), continue; end
    imshow(img_display{i}); hold on;

    fc = frame_centers{i};
    fr = frame_radii{i};
    fg = frame_good{i};
    cm = center_med(i,:);
    rm = radius_med(i);

    % ── Step 1: draw blue median circle first (bottom layer) ────────────
    if ~isnan(rm)
        viscircles(cm, rm * 0.9, 'Color', 'b', 'LineWidth', 2.5);
    end

    % ── Step 2: draw green or red ON TOP of blue (top layer) ────────────
    % All frames within tolerance → green; any outlier → red.
    % Use the median circle radius so the overlay sits exactly on the blue.
    if ~isnan(rm)
        all_ok   = all(fg);
        ovl_col  = ternary(all_ok, [0 0.85 0], [0.95 0.1 0.1]);  % green / red
        viscircles(cm, rm * 0.9, 'Color', ovl_col, 'LineWidth', 1.2);
    end

    % ── Title: use Position column from metadata, fall back to PosXXXX / folder name ─
    pos_label = '';
    if isfield(data.params, 'position_labels') && ~isempty(data.params.position_labels)
        labels = data.params.position_labels;
        if i <= numel(labels)
            % handles both cell arrays of strings and string arrays
            pos_label = char(labels{i});
        end
    end
    if isempty(pos_label)
        pos_label = regexp(fn(i).name, 'Pos\d{4}', 'match', 'once');
    end
    if isempty(pos_label)
        pos_label = fn(i).name;
    end
    pct_good = round(100 * sum(fg) / numel(fg));
    title(sprintf('%s  (%d%% ok)', pos_label, pct_good), ...
          'FontSize', 13, 'FontWeight', 'bold', 'Interpreter', 'none');
end

%% Save figure
if isfield(data.params, 'output_path') && ~isempty(data.params.output_path)
    if ~exist(data.params.output_path, 'dir'), mkdir(data.params.output_path); end
    out_png = fullfile(data.params.output_path, '00_plate_center_summary.png');
    saveas(fig, out_png);
    fprintf('\n[save] Plate center summary -> %s\n', out_png);
end

end


%% =========================================================================
%% find_center_individual — Hough primary, blob fallback
%% =========================================================================
function [center_pos, radius] = find_center_individual(img)

img_gray = rgb2gray(img);
[H, W]   = size(img_gray);

scale     = min(512 / min(H,W), 1.0);
img_small = imresize(img_gray, scale);
[Hs, Ws]  = size(img_small);
shortS    = min(Hs, Ws);

r_lo = round(shortS * 0.25);
r_hi = round(shortS * 0.85);

try
    [cb, rb] = imfindcircles(img_small, [r_lo r_hi], ...
        'ObjectPolarity', 'bright', 'Sensitivity', 0.93, 'EdgeThreshold', 0.04);
    [cd, rd] = imfindcircles(img_small, [r_lo r_hi], ...
        'ObjectPolarity', 'dark',   'Sensitivity', 0.93, 'EdgeThreshold', 0.04);
    all_c = [cb; cd];
    all_r = [rb; rd];
catch
    all_c = []; all_r = [];
end

if ~isempty(all_c)
    img_mid = [Ws/2, Hs/2];
    n_cands = size(all_c, 1);
    scores  = zeros(n_cands, 1);
    r_range = max(all_r) - min(all_r);

    for k = 1:n_cands
        c = all_c(k,:);
        r = all_r(k);
        centre_pen   = norm(c - img_mid) / shortS;
        r_bonus      = ternary(r_range > 0, (r - min(all_r)) / r_range, 0);
        margin       = min([c(1), c(2), Ws-c(1), Hs-c(2)]);
        boundary_pen = max(0, 1 - margin / (shortS * 0.15));
        scores(k)    = -1.5*centre_pen + 1.2*r_bonus - 2.0*boundary_pen;
    end

    [~, ix]    = max(scores);
    center_pos = all_c(ix,:) ./ scale;
    radius     = all_r(ix)   ./ scale;
    return;
end

%% Blob fallback
img_bin   = ~imbinarize(img_small);
im_filled = imfill(img_bin, 'holes');
im_filled = bwareaopen(im_filled, round(Hs*Ws*0.05));
rp        = regionprops(im_filled, 'Area', 'Centroid', 'EquivDiameter');

if isempty(rp)
    center_pos = [W/2, H/2];
    radius     = min(H,W) / 2;
    return;
end

[~, ix]    = max([rp.Area]);
center_pos = rp(ix).Centroid ./ scale;
radius     = (rp(ix).EquivDiameter / 2) ./ scale;

end


%% =========================================================================
%% choose_frames — evenly spaced sample indices
%% =========================================================================
function frame_idx = choose_frames(n_files, max_frames)
if isinf(max_frames) || max_frames >= n_files
    frame_idx = 1:n_files;
    return;
end
max_frames = max(1, round(max_frames));
frame_idx  = unique(round(linspace(1, n_files, max_frames)));
end


%% =========================================================================
%% best_subplot_grid — auto grid closest to 4:3
%% =========================================================================
function [best_rows, best_cols] = best_subplot_grid(n)
target_ratio = 4/3;
max_waste    = max(1, ceil(n * 0.10));
best_rows = 1; best_cols = n;
best_err  = abs(n - target_ratio);
for nr = 1:n
    nc = ceil(n/nr);
    if nc < nr, break; end
    if nc*nr - n > max_waste, continue; end
    err = abs(nc/nr - target_ratio);
    if err < best_err
        best_err = err; best_rows = nr; best_cols = nc;
    end
end
end


%% =========================================================================
%% ternary helper
%% =========================================================================
function out = ternary(cond, a, b)
if cond, out = a; else, out = b; end
end


%% =========================================================================
%% get_param helper
%% =========================================================================
function val = get_param(params, name, default_val)
if isfield(params, name) && ~isempty(params.(name))
    val = params.(name);
else
    val = default_val;
end
end
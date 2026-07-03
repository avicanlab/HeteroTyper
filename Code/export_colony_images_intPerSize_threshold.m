%% HeteroTyper Pipeline for Bright Plates
% Exports colony images split into "valid" and "filtered" categories
% based on intensity-per-size (IntPerSize) thresholds chosen interactively.
% Uses the same colony population as plot_combined_samples — filtering
% is read from ht.groups; no re-filtering from raw data.
%
% Requires preprocess_pipeline_data(data) to have been run first.
%
%% USAGE:
%   export_colony_img_intPerSize_threshold(data, ht);

function export_colony_img_intPerSize_threshold(data, ht)
    HT_FLOG = -1;

    if ~isstruct(ht) || ~isfield(ht, 'groups')
        error(['Second input must be the ht struct produced by preprocess_pipeline_data.\n' ...
               'Run:  preprocess_pipeline_data(data)  first.']);
    end

    p      = ht.params;
    labels = ht.labels;

    % Images go into a dedicated subfolder
    img_root = fullfile(p.out_dir, 'Colony Images - IntensitySize threshold');

    % ------------------------------------------------------------------
    %  Open log file
    % ------------------------------------------------------------------
    timestamp_log = datestr(now, 'yyyymmdd_HHMMSS');
    log_path = fullfile(p.out_dir, ...
                        sprintf('export_colony_img_intPerSize_threshold_%s.txt', timestamp_log));
    HT_FLOG = fopen(log_path, 'w');
    if HT_FLOG == -1, warning('Could not open log file: %s', log_path); end
    ht_fprintf(HT_FLOG, 'Log file: %s\n', log_path);

    % ------------------------------------------------------------------
    %  Ask user how to define the IntPerSize threshold(s)
    % ------------------------------------------------------------------
    fprintf('\n--- IntPerSize threshold mode ---\n');
    fprintf('   How would you like to define the threshold(s)?\n');
    fprintf('   [1] Lower bound only   — export colonies with IntPerSize < L\n');
    fprintf('   [2] Upper bound only   — export colonies with IntPerSize > H\n');
    fprintf('   [3] Both bounds        — export colonies with IntPerSize < L  OR  > H\n\n');

    while true
        raw_mode = strtrim(input('   Enter mode (1, 2, or 3): ', 's'));
        mode = str2double(raw_mode);
        if ismember(mode, [1 2 3]), break; end
        fprintf('   WARNING  Please enter 1, 2, or 3.\n');
    end

    IntS_thr_L = NaN;
    IntS_thr_H = NaN;

    if mode == 1 || mode == 3
        while true
            raw = strtrim(input('   Enter lower bound L (colonies with IntPerSize < L are "filtered"): ', 's'));
            val = str2double(raw);
            if ~isnan(val) && val >= 0, IntS_thr_L = val; break; end
            fprintf('   WARNING  Please enter a non-negative number.\n');
        end
    end

    if mode == 2 || mode == 3
        while true
            raw = strtrim(input('   Enter upper bound H (colonies with IntPerSize > H are "filtered"): ', 's'));
            val = str2double(raw);
            if ~isnan(val) && val >= 0
                if mode == 3 && val <= IntS_thr_L
                    fprintf('   WARNING  H must be greater than L (%.4g).\n', IntS_thr_L);
                    continue;
                end
                IntS_thr_H = val;
                break;
            end
            fprintf('   WARNING  Please enter a non-negative number.\n');
        end
    end

    % Build threshold description for log and folder names
    switch mode
        case 1, thr_desc = sprintf('IntPerSize_lt_%.2f', IntS_thr_L);
        case 2, thr_desc = sprintf('IntPerSize_gt_%.2f', IntS_thr_H);
        case 3, thr_desc = sprintf('IntPerSize_lt_%.2f_or_gt_%.2f', IntS_thr_L, IntS_thr_H);
    end

    fprintf('\n   -> Mode: %s\n\n', thr_desc);
    ht_fprintf(HT_FLOG, 'Threshold mode: %s\n\n', thr_desc);

    crop_size        = 150;
    border_thickness = 1;

    for g = 1:length(ht.groups)
        grp          = ht.groups(g);
        sample_group = labels{g};

        % ------------------------------------------------------------------
        %  Iterate plates in the same order as preprocess_pipeline_data.
        % ------------------------------------------------------------------
        slot = 0;

        for j = 1:length(grp.plate_indices)
            plate_idx = grp.plate_indices(j);

            colonies = data.processed{plate_idx}.colonies;
            if ~isfield(colonies, 'new') || ~isfield(colonies.new, 'lag_time')
                continue;
            end
            n_here = length(colonies.new.lag_time);
            if n_here == 0, continue; end

            idx_range = slot + (1:n_here);
            slot = slot + n_here;

            % Plate identifier
            pos_raw = data.metadata.original.Position(plate_idx);
            if iscell(pos_raw), plate_no = pos_raw{1}; else, plate_no = pos_raw; end

            % Slice ht.groups arrays for this plate
            size_p = grp.size(idx_range);
            lag_p  = grp.lag_time(idx_range);
            int_p  = grp.intensity(idx_range);
            ips_p  = grp.int_per_size(idx_range);
            cen_p  = grp.centroid(idx_range, :);

            % Additional metrics for annotation
            mint_p = grp.mean_intensity(idx_range);
            peri_p = grp.perimeter(idx_range);
            circ_p = grp.circularity(idx_range);
            ecc_p  = grp.eccentricity(idx_range);
            sol_p  = grp.solidity(idx_range);

            % ------------------------------------------------------------------
            %  Base population: same as plot_combined_samples (non-NaN size)
            %  plus finite IntPerSize (required for threshold comparison).
            % ------------------------------------------------------------------
            base_ok = ~isnan(size_p) & isfinite(ips_p);

            % Apply user-chosen threshold(s) to identify "filtered" colonies
            switch mode
                case 1, filter_mask = base_ok & (ips_p <  IntS_thr_L);
                case 2, filter_mask = base_ok & (ips_p >  IntS_thr_H);
                case 3, filter_mask = base_ok & (ips_p <  IntS_thr_L | ips_p > IntS_thr_H);
            end
            % "Valid" = passed base filter AND NOT filtered by IntPerSize
            valid_mask = base_ok & ~filter_mask;

            n_valid    = sum(valid_mask);
            n_filtered = sum(filter_mask);

            if n_valid + n_filtered == 0
                ht_fprintf(HT_FLOG, 'Plate %s (%s): no colonies to export\n', ...
                           num2str(plate_no), sample_group);
                continue;
            end

            % Load image data only when there is something to export
            raw_image = data.processed{plate_idx}.img_final;
            [img_h, img_w, ~] = size(raw_image);

            % Build labeled mask for green border overlay
            mask_bin = data.processed{plate_idx}.colonies.debug.binary;
            if islogical(mask_bin)
                labeled_mask = bwlabel(mask_bin);
            else
                labeled_mask = mask_bin;
            end
            mask_rgb = label2rgb(labeled_mask, 'jet', 'k', 'shuffle');

            % Map ht.groups colony indices back to region_props ok_idx
            props  = colonies.region_props;
            ok_idx = find(props.flag_colony_ok == 1);

            category_sets = { ...
                'Filtered colonies', find(filter_mask); ...
                'Valid colonies',    find(valid_mask)   ...
            };

            for cat_i = 1:size(category_sets, 1)
                cat_name = category_sets{cat_i, 1};
                ci_list  = category_sets{cat_i, 2};
                if isempty(ci_list), continue; end

                for k = 1:length(ci_list)
                    ci = ci_list(k);

                    cx = round(cen_p(ci, 1));
                    cy = round(cen_p(ci, 2));
                    if isnan(cx) || isnan(cy), continue; end

                    x_min = max(1, cx-crop_size);  x_max = min(img_w, cx+crop_size);
                    y_min = max(1, cy-crop_size);  y_max = min(img_h, cy+crop_size);

                    cropped_raw  = raw_image(y_min:y_max, x_min:x_max, :);
                    cropped_mask = mask_rgb(y_min:y_max,  x_min:x_max, :);

                    gray_rgb = repmat(make_gray(cropped_raw), [1 1 3]);

                    % Green border around the specific colony
                    % Map the local colony index back to the original region_props index
                    if ci <= length(ok_idx)
                        actual_idx  = ok_idx(ci);
                        colony_mask = (labeled_mask == actual_idx);
                        colony_crop = colony_mask(y_min:y_max, x_min:x_max);
                        thick_border = imdilate(bwperim(colony_crop), strel('disk', border_thickness));
                        gray_rgb(:,:,1) = set_channel(gray_rgb(:,:,1), thick_border, 0);
                        gray_rgb(:,:,2) = set_channel(gray_rgb(:,:,2), thick_border, 255);
                        gray_rgb(:,:,3) = set_channel(gray_rgb(:,:,3), thick_border, 0);
                    end

                    blended = imfuse(gray_rgb, cropped_mask, 'blend');

                    text_lines = { ...
                        sprintf('Lag=%.1f',     lag_p(ci)), ...
                        sprintf('Size=%.2f',    size_p(ci)), ...
                        sprintf('Int=%.1f',     int_p(ci)), ...
                        sprintf('MeanInt=%.1f', mint_p(ci)), ...
                        sprintf('IntSize=%.1f', ips_p(ci)), ...
                        sprintf('Cir=%.2f',     circ_p(ci)), ...
                        sprintf('Ecc=%.2f',     ecc_p(ci)), ...
                        sprintf('Sol=%.2f',     sol_p(ci)), ...
                        sprintf('Peri=%.2f',    peri_p(ci)) ...
                    };
                    for t = 1:numel(text_lines)
                        blended = insertText(blended, ...
                            [size(blended,2)-5, 5+12*(t-1)], text_lines{t}, ...
                            'FontSize', 10, 'TextColor', 'white', ...
                            'BoxColor', 'black', 'BoxOpacity', 0, 'AnchorPoint', 'RightTop');
                    end

                    out_dir = fullfile(img_root, thr_desc, cat_name, sample_group);
                    if ~exist(out_dir, 'dir'), mkdir(out_dir); end
                    fname = sprintf('IntSize_%.2f_%s_Plate%s_Colony%d.png', ...
                                    ips_p(ci), sample_group, num2str(plate_no), ci);
                    imwrite(blended, fullfile(out_dir, fname));
                end
            end

            ht_fprintf(HT_FLOG, 'Plate %s (%s): valid=%d, filtered=%d\n', ...
                       num2str(plate_no), sample_group, n_valid, n_filtered);
        end
    end

    ht_fprintf(HT_FLOG, '\nExport complete. Images saved to: %s\n', img_root);
    if HT_FLOG ~= -1, fclose(HT_FLOG); fprintf('Log saved: %s\n', log_path); end
end


function gray = make_gray(img)
    if size(img,3) == 3, gray = rgb2gray(img); else, gray = img; end
    gray = imadjust(gray, [0.025 0.15], []);
end


function ch = set_channel(ch, mask, val)
    ch(mask) = val;
end


function ht_fprintf(HT_FLOG, varargin)
    fprintf(varargin{:});
    if ~isempty(HT_FLOG) && HT_FLOG ~= -1
        fprintf(HT_FLOG, varargin{:});
    end
end
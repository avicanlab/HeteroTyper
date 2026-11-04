%% HeteroTyper Pipeline
% Exports cropped colony images (blended raw + segmentation) for colonies
% with appearance time above a user-defined threshold.
% Uses the same colony population as plot_combined_samples — filtering
% is read from ht.groups; no re-filtering from raw data.
%
% Requires preprocess_pipeline_data(data) to have been run first.
%
%% USAGE:
%   export_colony_images(data, ht);

function export_colony_images(data, ht)
    HT_FLOG = -1;

    if ~isstruct(ht) || ~isfield(ht, 'groups')
        error(['Second input must be the ht struct produced by preprocess_pipeline_data.\n' ...
               'Run:  preprocess_pipeline_data(data)  first.']);
    end

    p      = ht.params;
    labels = ht.labels;

    % Images go into a dedicated subfolder
    img_root = fullfile(p.out_dir, 'Colony Images');

    % ------------------------------------------------------------------
    %  Open log file
    % ------------------------------------------------------------------
    timestamp_log = datestr(now, 'yyyymmdd_HHMMSS');
    log_path = fullfile(p.out_dir, sprintf('export_colony_images_%s.txt', timestamp_log));
    HT_FLOG = fopen(log_path, 'w');
    if HT_FLOG == -1, warning('Could not open log file: %s', log_path); end
    ht_fprintf(HT_FLOG, 'Log file: %s\n', log_path);

    % ------------------------------------------------------------------
    %  Ask user for appearance-time threshold
    % ------------------------------------------------------------------
    fprintf('\n--- Appearance-time threshold for colony image export ---\n');
    fprintf('   Only colonies with appearance_time > threshold will be exported.\n');
    while true
        raw = strtrim(input('   Enter appearance-time threshold in hours [e.g. 38]: ', 's'));
        appearance_threshold = str2double(raw);
        if ~isnan(appearance_threshold) && appearance_threshold >= 0, break; end
        fprintf('   WARNING  Please enter a non-negative number.\n');
    end
    fprintf('   -> Appearance-time threshold: %.4g h\n\n', appearance_threshold);
    ht_fprintf(HT_FLOG, 'Appearance-time threshold: %.4g h\n\n', appearance_threshold);

    crop_size = 150;

    total_plates = sum(cellfun(@(g) length(g.plate_indices), num2cell(ht.groups)));
    plate_count  = 0;

    for g = 1:length(ht.groups)
        grp          = ht.groups(g);
        sample_group = labels{g};

        % ------------------------------------------------------------------
        %  Iterate plates in the same order as preprocess_pipeline_data,
        %  slicing ht.groups arrays by the number of colonies per plate.
        %  This is the only correct way to get per-colony plate identity
        %  without re-filtering from data.processed.
        % ------------------------------------------------------------------
        slot = 0;

        for j = 1:length(grp.plate_indices)
            plate_idx = grp.plate_indices(j);
            plate_count = plate_count + 1;

            colonies = data.processed{plate_idx}.colonies;
            if ~isfield(colonies, 'new') || ~isfield(colonies.new, 'appearance_time')
                continue;
            end
            n_here = length(colonies.new.appearance_time);
            if n_here == 0, continue; end

            idx_range = slot + (1:n_here);
            slot = slot + n_here;

            % Plate identifier
            pos_raw = data.metadata.original.Position(plate_idx);
            if iscell(pos_raw), plate_no = pos_raw{1}; else, plate_no = pos_raw; end

            % Slice ht.groups arrays for this plate's colonies
            appearance_p  = grp.appearance_time(idx_range);
            size_p = grp.size(idx_range);
            int_p  = grp.intensity(idx_range);
            cen_p  = grp.centroid(idx_range, :);

            % Additional metrics for annotation
            area_p  = grp.area(idx_range);
            mint_p  = grp.mean_intensity(idx_range);
            ips_p   = grp.int_per_size(idx_range);
            peri_p  = grp.perimeter(idx_range);
            circ_p  = grp.circularity(idx_range);
            ecc_p   = grp.eccentricity(idx_range);
            sol_p   = grp.solidity(idx_range);

            % ------------------------------------------------------------------
            %  Select colonies: same population as plot_combined_samples
            %  (non-NaN size = passed size filter) PLUS appearance > threshold.
            %  No additional ecc_threshold filter — matches plot_combined_samples.
            % ------------------------------------------------------------------
            target = find(~isnan(size_p) & isfinite(appearance_p) & (appearance_p > appearance_threshold));

            if isempty(target)
                ht_fprintf(HT_FLOG, 'Plate %s (%s): 0 colonies above threshold\n', ...
                           num2str(plate_no), sample_group);
                continue;
            end

            % Load image data for this plate
            raw_image  = data.processed{plate_idx}.img_final;
            mask_clean = data.processed{plate_idx}.colonies.debug.segmented;
            [img_h, img_w, ~] = size(raw_image);

            n_exported = 0;
            for k = 1:length(target)
                ci = target(k);

                cx = round(cen_p(ci, 1));
                cy = round(cen_p(ci, 2));
                if isnan(cx) || isnan(cy), continue; end

                x_min = max(1, cx-crop_size);  x_max = min(img_w, cx+crop_size);
                y_min = max(1, cy-crop_size);  y_max = min(img_h, cy+crop_size);

                cropped_raw  = raw_image(y_min:y_max, x_min:x_max, :);
                cropped_mask = mask_clean(y_min:y_max, x_min:x_max, :);

                gray_img = make_gray(cropped_raw);
                blended  = imfuse(repmat(gray_img, [1 1 3]), cropped_mask, 'blend');

                text_lines = { ...
                    sprintf('Appearance=%.1f',     appearance_p(ci)), ...
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

                out_dir = fullfile(img_root, sample_group);
                if ~exist(out_dir, 'dir'), mkdir(out_dir); end
                fname = sprintf('%s_Plate%s_Colony%d.png', sample_group, num2str(plate_no), ci);
                imwrite(blended, fullfile(out_dir, fname));
                n_exported = n_exported + 1;
            end

            ht_fprintf(HT_FLOG, 'Plate %s (%s): exported %d / %d colonies above threshold\n', ...
                       num2str(plate_no), sample_group, n_exported, length(target));
        end
    end

    ht_fprintf(HT_FLOG, '\nExport complete. Images saved to: %s\n', img_root);
    if HT_FLOG ~= -1, fclose(HT_FLOG); fprintf('Log saved: %s\n', log_path); end
end


function gray = make_gray(img)
    if size(img,3) == 3, gray = rgb2gray(img); else, gray = img; end
    gray = imadjust(gray, [0.025 0.15], []);
end


function ht_fprintf(HT_FLOG, varargin)
    fprintf(varargin{:});
    if ~isempty(HT_FLOG) && HT_FLOG ~= -1
        fprintf(HT_FLOG, varargin{:});
    end
end
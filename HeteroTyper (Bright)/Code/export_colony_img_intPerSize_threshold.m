%% 11.11.2025 - HeteroTyper Pipeline for Bright Plates

function data = export_colony_img_intPerSize_threshold(data)

    nr_plates = length(data.processed);
    
    % Define Room Temperature Incubation Time
    incTime = 20;
    
    % Define min-max colony numbers
    min_col = 10;
    max_col = 650;

    %% --- Define Parameters/Thresholds ---
    Size_threshold = 100;
    Ecc_threshold  = 0.7;
    IntS_thr_L = 25;
    IntS_thr_H = 75;
    
    crop_size = 150;

    file_path = 'D:\Gizem\HeteroTyper\Test_Bright\Test_output\Colony Images - IntensitySize threshold\';            
    
    %% --- Extract plates with valid growth quantification ---
    for i = 1:nr_plates
        growth_available(i,1) = data.processed{i}.growth_quant;
    
        % --- Skip plates with too many or too few colonies ---
        % manually adapt if needed: 
        n_col = length(find(data.processed{i}.colonies.region_props.flag_colony_ok));
        if n_col > max_col || n_col < min_col
            growth_available(i,1) = 0;
        end
    end
    
    ix_growth = find(growth_available);
    
    %% split things by strain/site 
    ix_3H = find(strcmpi(data.metadata.original.Time,'3H')); 
    ix_7H = find(strcmpi(data.metadata.original.Time,'7H')); 
    ix_24H = find(strcmpi(data.metadata.original.Time,'24H'));  
    ix_48H = find(strcmpi(data.metadata.original.Time,'48H'));
    labels_type = {'3H','7H','24H','48H'}; 
    
    ix{1} = intersect(ix_3H,ix_growth);
    ix{2} = intersect(ix_7H,ix_growth);
    ix{3} = intersect(ix_24H,ix_growth);
    ix{4} = intersect(ix_48H,ix_growth);
    
    
    % ---------- PROGRESS ----------
    total_tasks = sum(cellfun(@length, ix));
    task_count = 0;

    %% === MAIN LOOP ===
    for i = 1:length(ix)
        ix_tmp = ix{i};
    
        for j = 1:length(ix_tmp)
            task_count = task_count + 1;
    
            sample_group   = labels_type{i};
            plate_no       = data.metadata.original.Position(ix_tmp(j));

            colonies = data.processed{ix_tmp(j)}.colonies;
            region_props = colonies.region_props;
    
            %% --- Extract image and colony data ---
            raw_image    = data.processed{ix_tmp(j)}.img_final;
            mask_clean   = colonies.debug.binary;

            col_flag        = props.flag_colony_ok;
            lag_time        = colonies.lag_time + incTime; 
            col_size        = colonies.timecourse_size_smoothed(end,:)'; 
            col_int         = colonies.timecourse_intensity_smoothed(end,:)';
            col_area        = props.Area(:);
            mean_int        = props.MeanIntensity(:); 
            col_peri        = props.Perimeter(:); 
            col_circ        = props.Circularity(:); 
            col_ecc         = props.Eccentricity(:); 
            col_sol         = props.Solidity(:); 
            centroids       = props.Centroid(:); 
            
            int_per_size    = col_int ./ col_size;

            
            %% --- Apply filtering: enforce size & eccentricity ---
            valid_size_ecc  = (col_flag == 1) & (col_size >= Size_threshold & col_ecc < Ecc_threshold) ;                % Apply "Size" & "Eccentricity" thresholds here!
            valid_intSize   = (isfinite(int_per_size) & (int_per_size >= IntS_thr_H | int_per_size <= IntS_thr_L));     % Apply "Intensity per size" thresholds here!
    
            valid_mask      = valid_size_ecc & valid_intSize;
            filtered_mask   = ~filter_mask;
            col_flag        = filter_mask;
    
            %% --- Prepare category indices ---
            valid_idx       = find(valid_mask);
            filtered_idx    = find(filtered_mask);
    
            category_sets = {
                'Filtered colonies', filtered_idx;
                'Valid colonies', valid_idx
            };
    
            %% --- RGB label mask ---
            if islogical(mask_clean)
                labeled_mask = bwlabel(mask_clean);
            else
                labeled_mask = mask_clean;
            end
            mask_rgb = label2rgb(labeled_mask, 'jet', 'k', 'shuffle');
    
            %% --- Get centroid data ---
            if istable(region_props)
                centroid_data = region_props.Centroid;
            elseif isstruct(region_props)
                centroid_data = {region_props.Centroid}';
            else
                warning('Unexpected region_props format (Plate %d)', plate_no);
                continue;
            end
    
            %% --- Image geometry ---
            [img_h, img_w, ~] = size(raw_image);
            plate_center = [img_w/2, img_h/2];
            plate_radius = min(img_w, img_h) * 0.85; % usable area
    
            %% --- Export loop per category ---
            for cat_i = 1:size(category_sets,1)
                category_name = category_sets{cat_i,1};
                idx_set = category_sets{cat_i,2};
    
                if isempty(idx_set)
                    continue;
                end
    
                border_thickness = 1; % adjustable outline thickness
    
                for colony_id = idx_set'
                    if colony_id > numel(col_flag) || ~col_flag(colony_id)
                        continue;
                    end
    
                    % --- Get centroid safely ---
                    if iscell(centroid_data)
                        centroid = round(centroid_data{colony_id});
                    else
                        centroid = round(centroid_data(colony_id, :));
                    end
                    if numel(centroid) ~= 2 || any(isnan(centroid))
                        continue;
                    end
                    x = centroid(1); y = centroid(2);
    
                    % --- Get colony mask & border check ---
                    colony_mask = (labeled_mask == colony_id);
                    [yy, xx] = find(colony_mask);
                    dist_px = sqrt((xx - plate_center(1)).^2 + (yy - plate_center(2)).^2);
                    if any(dist_px > plate_radius) || nnz(colony_mask) == 0
                        col_flag(colony_id) = 0;
                        continue;
                    end
    
                    % --- Crop colony ---
                    x_min = max(1, x - crop_size);
                    x_max = min(img_w, x + crop_size);
                    y_min = max(1, y - crop_size);
                    y_max = min(img_h, y + crop_size);
    
                    cropped_raw  = raw_image(y_min:y_max, x_min:x_max, :);
                    cropped_mask = mask_rgb(y_min:y_max, x_min:x_max, :);
    
                    % --- Prepare grayscale image ---
                    if size(cropped_raw,3) == 3
                        gray_img = rgb2gray(cropped_raw);
                    else
                        gray_img = cropped_raw;
                    end
                    gray_img = imadjust(gray_img,[0.025 0.15],[]);
                    gray_img_rgb = repmat(gray_img,[1 1 3]);
    
                    % --- Compute thicker green border ---
                    colony_crop_mask = colony_mask(y_min:y_max, x_min:x_max);
                    perim = bwperim(colony_crop_mask);
                    thick_border = imdilate(perim, strel('disk', border_thickness));
    
                    % --- Apply green border ---
                    r = gray_img_rgb(:,:,1);
                    g = gray_img_rgb(:,:,2);
                    b = gray_img_rgb(:,:,3);
                    r(thick_border) = 0;
                    g(thick_border) = 255;
                    b(thick_border) = 0;
                    gray_img_rgb(:,:,1) = r;
                    gray_img_rgb(:,:,2) = g;
                    gray_img_rgb(:,:,3) = b;
    
                    % --- Blend with rainbow mask ---
                    blended_img = imfuse(gray_img_rgb, cropped_mask, 'blend');
    
                    % --- Fetch metrics ---
                    lagT_val         = lag_time(colony_id);
                    size_val         = col_size(colony_id);
                    area_val         = col_area(colony_id);
                    int_val          = col_int(colony_id);
                    meanInt_val      = mean_int(colony_id);
                    intS_val         = int_val ./ size_val;
                    peri_val         = col_peri(colony_id);
                    circ_val         = col_circ(colony_id);
                    ecc_val          = col_ecc(colony_id);
                    sol_val          = col_sol(colony_id);
    
                    % --- Overlay text ---
                    text_lines = {
                        sprintf('Lag=%.1f', lagT_val)
                        sprintf('Size=%.2f', size_val)
                        %sprintf('Area=%.1f', area_val)
                        sprintf('Int=%.1f', int_val)
                        sprintf('MeanInt=%.1f', meanInt_val)
                        sprintf('IntSize=%.1f', intS_val)
                        sprintf('Cir=%.2f', circ_val)
                        sprintf('Ecc=%.2f', ecc_val)
                        sprintf('Sol=%.2f', sol_val)
                        sprintf('Peri=%.2f', peri_val)
                    };
                    
                    for t = 1:numel(text_lines)
                        blended_img = insertText(blended_img, ...
                            [size(blended_img,2)-5, 5 + 12*(t-1)], text_lines{t}, ...
                            'FontSize', 10, 'TextColor', 'white', ...
                            'BoxColor', 'black', 'BoxOpacity', 0, ...
                            'AnchorPoint', 'RightTop');
                    end
    
                    % --- Save image ---
                    output_dir = fullfile(file_path, category_name, sample_group);
                    if ~exist(output_dir, 'dir')
                        mkdir(output_dir);
                    end
                    filename = sprintf('IntSize_%.2f_%s_Plate%d_Colony%d.png', ...
                                       int_per_size_val, sample_group, plate_no, colony_id);
                    imwrite(blended_img, fullfile(output_dir, filename));
                end
            end
    
            fprintf('Plate %d - Final colonies exported: %d\n', ...
                    plate_no, sum(col_flag));
        end
    end

end

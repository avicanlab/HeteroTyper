%% 11.11.2025 - HeteroTyper Pipeline for Bright Plates

function data = export_colony_images(data) 

    nr_plates = length(data.processed); 

    % Define Room Temperature Incubation Time
    incTime = 20;
    
    % Define min-max colony numbers
    min_col = 10;
    max_col = 650;


    % Define filter parameters
    size_threshold = 100; 
    ecc_threshold  = 0.70;
    lag_threshold  = 38;

    % adjust window size as needed
    crop_size = 150; 

    file_path = 'D:\Gizem\HeteroTyper\Test_Bright\Test_output\Colony Images\';

    %% extract plates where we did quantify growth parameters 
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
    
    ix1 = ix_3H; 
    ix1 = intersect(ix1,ix_growth); 
    ix2 = ix_7H; 
    ix2 = intersect(ix2,ix_growth); 
    ix3 = ix_24H; 
    ix3 = intersect(ix3,ix_growth); 
    ix4 = ix_48H; 
    ix4 = intersect(ix4,ix_growth);

    ix{1} = ix1; 
    ix{2} = ix2; 
    ix{3} = ix3; 
    ix{4} = ix4; 
    
    % ---------- PROGRESS BAR SETUP ----------
    total_tasks = sum(cellfun(@length, ix));  % total number of plates
    task_count = 0;

    for i = 1:length(ix) 

        ix_tmp = ix{i}; 
        mismatch_report = [];
        group_imgs = {};
        group_labels = {};
    
        for j = 1:length(ix_tmp) 
            % increment global task count
            task_count = task_count + 1;

            %% ----- BEGIN COLONY IMAGE EXTRACTION -----

            sample_group = labels_type{i};
            plate_no = data.metadata.original.Position(ix_tmp(j)); % plate number

            raw_image = data.processed{ix_tmp(j)}.img_final; % assumes .img_final exists (original plate image)
            mask_clean = data.processed{ix_tmp(j)}.colonies.debug.segmented; % segmentation mask

            colonies = data.processed{ix_tmp(j)}.colonies.new;
            props = data.processed{ix_tmp(j)}.colonies.region_props;
            
            % Colony data
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

            
            % Find colonies with lag time > threshold and flagged as OK
            target_colonies = find(col_flag & (lag_time > lag_threshold) & (col_size > size_threshold) & (col_ecc < ecc_threshold));
            
            if ~isempty(target_colonies)
                % Get the indices of colonies that passed initial segmentation
                n_colonies = size(props.Centroid,1);
                valid_idx = target_colonies(target_colonies <= n_colonies);

                % Generate color-coded mask (like in plot_individual_plate)
                mask_rgb = mask_clean;
            
                for k = 1:length(valid_idx)
                    colony_id = valid_idx(k);
            
                    % Get centroid
                    centroid = round(props.Centroid(colony_id, :));
                    x = centroid(1);
                    y = centroid(2);
            
                    % Define bounds (make sure we don’t exceed image limits)
                    [img_h,img_w,~] = size(raw_image);
                    x_min = max(1, x - crop_size);
                    x_max = min(size(raw_image, 2), x + crop_size);
                    y_min = max(1, y - crop_size);
                    y_max = min(size(raw_image, 1), y + crop_size);
            
                    % Crop raw and mask images
                    cropped_raw = raw_image(y_min:y_max, x_min:x_max, :);
                    cropped_mask = mask_rgb(y_min:y_max, x_min:x_max, :);
            
                    % Convert raw image to brightened grayscale (like your original figure)
                    if size(cropped_raw,3) == 3
                        gray_img = rgb2gray(cropped_raw);
                    else
                        gray_img = cropped_raw;
                    end
                    gray_img = imadjust(gray_img,[0.025 0.15],[]);   % match original adjustment
                    gray_img_rgb = repmat(gray_img,[1 1 3]);          % convert back to RGB for blending
            
                    % Blend segmentation mask with brightened gray colony
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
                    
                    % save in memory
                    group_imgs{end+1} = blended_img;
                    group_labels{end+1} = sprintf('P%d-C%d', plate_no, colony_id);


                    % Save
                    output_dir = fullfile(file_path, sample_group);
                    if ~exist(output_dir, 'dir')
                        mkdir(output_dir);
                    end
                    filename = sprintf('%s_Plate%d_Colony%d.png', sample_group, plate_no, colony_id);
                    fullpath = fullfile(output_dir, filename);

                    imwrite(blended_img, fullpath);
                end
            end


            %% ----- END COLONY IMAGE EXTRACTION -----

            % ---------- UPDATE PROGRESS ----------
            fprintf('\rProcessing overall: %d / %d (%.1f%%)', ...
                task_count, total_tasks, (task_count/total_tasks)*100);
            if task_count == total_tasks
                fprintf('\n'); % newline when done
            end

            
        end 
         
    end 

end

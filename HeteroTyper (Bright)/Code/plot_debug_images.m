%% 11.11.2025 - HeteroTyper Pipeline for Bright Plates 

% Run:   plot_debug_images(data_name);

function data = plot_debug_images(data)

    nr_plates = length(data.processed); 
    
    % Define min-max colony numbers
    min_col = 10;
    max_col = 650;

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

    fprintf('Number of sample groups:  %d\n', length(ix));

    for i = 1:length(ix)
        ix_tmp = ix{i};
        for j = 1:length(ix_tmp)
            % Extract debug struct for this plate/position
            debug = data.processed{ix_tmp(j)}.colonies.debug;
        
            % Convert image to grayscale safely
            raw_image = debug.gray;
            img_segm  = debug.segmented;
    
            % Convert raw image to brightened grayscale (like your original figure)
            if size(raw_image,3) == 3
                gray = rgb2gray(raw_image);
                else
                    gray = raw_image;
            end
            gray = imadjust(gray,[0.025 0.15],[]);   % match original adjustment
            gray_rgb = repmat(gray,[1 1 3]);          % convert back to RGB for blending
            
            mask_rgb = img_segm;
    
            % Blend segmentation mask with brightened gray colony
            blended_img = imfuse(gray_rgb, mask_rgb, 'blend');
    
    
            
            images = {debug.gray, debug.background, debug.flat, debug.norm, ...
                      debug.contrast, debug.smooth, debug.binary, blended_img};
            titles = {'1. Grayscale','2. Background','3. Background Flattened','4. Normalized', ...
                      '5. Contrasted','6. Smoothed','7. Binary Mask (Filtered)','8. Final Segmentation'};
        
            % Create figure
            fig_w = 2000;  
            fig_h = 800; 
            figure('Name',['Colony Detection Debug - ' num2str(ix_tmp(j))],'Position',[100 100 fig_w fig_h]);
            rows = 2; 
            cols = 4;
            gap = 0.01;      % 1% gap between subplots
            margins = 0.02;  % 2% outer margin
        
            % Calculate tile size
            w = (1 - (cols+1)*gap - 2*margins) / cols;
            h = (1 - (rows+1)*gap - 2*margins) / rows;
        
            % Plot loop
            for k = 1:(rows*cols)
                [row, col] = ind2sub([rows, cols], k);
        
                xpos = margins + (col-1)*(w+gap);
                ypos = 1 - (row)*(h+gap) - margins + gap;
        
                ax = axes('Position',[xpos ypos w h]); 
                imshow(images{k},[]);
                title(titles{k}, 'Interpreter','none');
            end
        end
    
        
    end

end


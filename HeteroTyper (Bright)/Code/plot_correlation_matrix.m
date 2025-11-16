%% 11.11.2025 - HeteroTyper Pipeline for Bright Plates

function data = plot_correlation_matrix(data) 

    nr_plates = length(data.processed);

    % Define Room Temperature Incubation Time
    incTime = 20; 
    
    % Define min-max colony numbers
    min_col = 10;
    max_col = 650;


    % Define filter parameters
    size_threshold = 100; 
    ecc_threshold = 0.70;


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


    %% Parameters to plot
    group_labels = {'3H','7H','24H','48H'};
    corr_param_names = {'LagTime','FinalSize','Area','Intensity','MeanIntensity','IntensityPerSize','Perimeter','Circularity','Eccentricity','Solidity'};

    % === Define colors from HEX ===
    hex2rgb = @(hex) sscanf(hex(2:end),'%2x%2x%2x',[1 3])/255;
    
    % blue   = hex2rgb('#2D54CB'); % blue
    blue   = hex2rgb('#2A50A1'); % blue
    white  = hex2rgb('#ffffff'); % white
    red    = hex2rgb('#aa0024'); % red
    
    % === Build colormap ===
    n = 256; 
    n_half = n/2;
    
    % -1 → 0 = blue → white
    blue_to_white = [linspace(blue(1), white(1), n_half)', ...
                     linspace(blue(2), white(2), n_half)', ...
                     linspace(blue(3), white(3), n_half)'];
    
    % 0 → +1 = white → red
    white_to_red = [linspace(white(1), red(1), n_half)', ...
                    linspace(white(2), red(2), n_half)', ...
                    linspace(white(3), red(3), n_half)'];
    
    % Full colormap
    custom_cmap = [blue_to_white; white_to_red];


    figure('Name','Correlation Matrices by Timepoint');
    for i = 1:length(ix) 

        ix_tmp = ix{i}; 

        % Initialize combined arrays for this group
        combined_lag_time  = [];
        combined_size      = [];
        combined_area      = [];
        combined_int       = [];
        combined_mean_int  = [];
        combined_peri     = [];
        combined_solid     = [];
        combined_circ      = [];
        combined_ecc       = [];
        combined_intSize   = [];


        % Initialize tracker before the for j=1:length(ix_tmp) loop
        mismatch_report = [];  
    
        for j = 1:length(ix_tmp) 

            idx_plate = ix_tmp(j);
            colonies  = data.processed{idx_plate}.colonies;
            props     = colonies.region_props;

            % Retrieve parameters
            lag_time   = colonies.new.lag_time(:) + incTime;  % Incubation time at room temperature (20 hours) included here
            col_size     = colonies.new.timecourse_size_smoothed(end,:)';
            col_area     = props.Area(:);
            col_int      = colonies.new.timecourse_intensity_smoothed(end,:)';
            col_meanInt  = props.MeanIntensity(:); 
            col_peri      = props.Perimeter(:);
            col_sol      = props.Solidity(:);
            col_circ     = props.Circularity(:);
            col_ecc      = props.Eccentricity(:);
            col_flag     = props.flag_colony_ok(:);

            % Generate combined sample group tables 
            % Final colony size and eccentricity filtered, and "non-infinite" intensity per size data will be used!
            ints_per_size    = col_int ./ col_size;

            % Keep only valid entries (nonzero, finite, flagged colonies)
            valid_idx_intS = (col_flag == 1) & isfinite(ints_per_size) & (col_size > size_threshold) & (col_ecc < ecc_threshold);
            combined_intSize     = [combined_intSize; ints_per_size(valid_idx_intS)];
            combined_lag_time    = [combined_lag_time; lag_time(valid_idx_intS)];
            combined_size        = [combined_size; col_size(valid_idx_intS)];
            combined_area        = [combined_area; col_area(valid_idx_intS)];
            combined_int         = [combined_int; col_int(valid_idx_intS)];
            combined_mean_int    = [combined_mean_int, col_meanInt(valid_idx_intS)];
            combined_peri        = [combined_peri; col_peri(valid_idx_intS)];
            combined_circ        = [combined_circ; col_circ(valid_idx_intS)];
            combined_ecc         = [combined_ecc; col_ecc(valid_idx_intS)];
            combined_solid       = [combined_solid; col_sol(valid_idx_intS)];

            
        end 
        
        % Build data matrix for this group
        param_matrix = [combined_lag_time, ...
                        combined_size, ...
                        combined_area, ...
                        combined_int, ...
                        combined_mean_int, ...
                        combined_intSize, ...
                        combined_peri, ...
                        combined_circ, ...
                        combined_ecc, ...
                        combined_solid];
    
        % Remove invalid rows
        param_matrix = param_matrix(all(isfinite(param_matrix),2), :);
    
        % Compute correlation
        R = corr(param_matrix, 'Type', 'Pearson');
    
        % Plot in subplot
        subplot(2,2,i);
        h = heatmap(corr_param_names, corr_param_names, R, ...
            'Colormap', custom_cmap, ...
            'ColorLimits', [-1 1], ...
            'CellLabelFormat','%.2f');
    
        title([group_labels{i} ' Correlation']);
        set(gca,'FontSize',12);

    end 
end

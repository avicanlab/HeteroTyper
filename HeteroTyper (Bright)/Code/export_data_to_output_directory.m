%% 11.11.2025 - HeteroTyper Pipeline for Bright Plates 

function data = export_data_to_output_directory(data) 

    nr_plates = length(data.processed); 

    % Define Room Temperature Incubation Time
    incTime = 20;
    
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
    

    % ---------- PROGRESS BAR SETUP ----------
    total_tasks = sum(cellfun(@length, ix));  % total number of plates
    task_count = 0;


    % Define filter parameters
    size_threshold = 100; 
    ecc_threshold = 0.70;

    
    %% === EXPORT FOR R: Combined Morphology Table (All Parameters) ===
    % Define output folder
    output_dir = 'D:\Gizem\HeteroTyper\Test_Bright\Test_output\';
    if ~exist(output_dir, 'dir')
        mkdir(output_dir);
    end
    
    % Initialize one master table for all groups
    allGroupsTable_Bright     = table();
    allGroupsTable_Bright_2   = table();
    
    % Loop over all sample groups (3H, 7H, 24H, 48H)
    for g = 1:length(labels_type)
        group_label = labels_type{g};
        group_indices = ix{g};
    
        % Initialize combined arrays for this group
        combined_name               = [];
        combined_folder             = [];
        combined_plate_group        = [];
        combined_lag_time_group     = [];
        combined_size_group         = [];
        combined_area_group         = [];
        combined_int_group          = [];
        combined_perim_group        = [];
        combined_circ_group         = [];
        combined_ecc_group          = [];
        combined_solid_group        = [];
        combined_centroid_x         = [];
        combined_centroid_y         = [];


        combined_name_2               = [];
        combined_folder_2             = [];
        combined_plate_group_2        = [];
        combined_lag_time_group_2     = [];
        combined_size_group_2         = [];
        combined_area_group_2         = [];
        combined_int_group_2          = [];
        combined_int_per_size_2       = [];
        combined_perim_group_2        = [];
        combined_circ_group_2         = [];
        combined_ecc_group_2          = [];
        combined_solid_group_2        = [];
        combined_centroid_x_2         = [];
        combined_centroid_y_2         = [];
    
    
        % Loop through plates belonging to this group
        for j = 1:length(group_indices)
            task_count = task_count + 1;

            idx_plate = group_indices(j);
            colonies  = data.processed{idx_plate}.colonies;
            props     = colonies.region_props;
            plate_no  = data.metadata.original.Position(task_count);
            f_folder  = data.metadata.fn(task_count).folder;
            f_name    = data.metadata.fn(task_count).name;
    
            % Retrieve parameters
            lag_time   = colonies.new.lag_time(:) + incTime;  % Incubation time at room temperature (20 hours) included here
            col_size   = colonies.new.timecourse_size_smoothed(end,:)';
            col_area   = props.Area(:);
            col_int    = colonies.new.timecourse_intensity_smoothed(end,:)';
            col_per    = props.Perimeter(:);
            col_sol    = props.Solidity(:);
            col_circ   = props.Circularity(:);
            col_ecc    = props.Eccentricity(:);
            col_flag   = props.flag_colony_ok(:);
            col_cent_x = props.Centroid(:,1); 
            col_cent_y = props.Centroid(:,2); 
    
            % Filter valid colonies (use same logic as before)
            valid_idx = (col_flag == 1) & isfinite(col_size) & (col_size > size_threshold) & (col_ecc < ecc_threshold);

            % Generate combined sample group tables 
            % Final colony size and eccentricity filtered data will be used!
            combined_name           = [combined_name; repmat({f_name}, sum(valid_idx), 1)];
            combined_folder         = [combined_folder; repmat({f_folder}, sum(valid_idx), 1)];
            combined_plate_group    = [combined_plate_group; repmat(plate_no, sum(valid_idx), 1)];
            combined_lag_time_group = [combined_lag_time_group; lag_time(valid_idx)];
            combined_size_group     = [combined_size_group; col_size(valid_idx)];
            combined_area_group     = [combined_area_group; col_area(valid_idx)];
            combined_int_group      = [combined_int_group; col_int(valid_idx)];
            combined_perim_group    = [combined_perim_group; col_per(valid_idx)];
            combined_circ_group     = [combined_circ_group; col_circ(valid_idx)];
            combined_ecc_group      = [combined_ecc_group; col_ecc(valid_idx)];
            combined_solid_group    = [combined_solid_group; col_sol(valid_idx)];
            combined_centroid_x     = [combined_centroid_x; col_cent_x(valid_idx)];
            combined_centroid_y     = [combined_centroid_y; col_cent_y(valid_idx)];
    


            % Generate combined sample group tables 
            % Final colony size and eccentricity filtered, and "non-infinite" intensity per size data will be used!
            ints_per_size    = col_int ./ col_size;
            % Keep only valid entries (nonzero, finite, flagged colonies)
            valid_idx_intS = (col_flag == 1) & isfinite(ints_per_size) & (col_size > size_threshold) & (col_ecc < ecc_threshold);
            combined_name_2               = [combined_name_2; repmat({f_name}, sum(valid_idx_intS), 1)];
            combined_folder_2             = [combined_folder_2; repmat({f_folder}, sum(valid_idx_intS), 1)];
            combined_plate_group_2        = [combined_plate_group_2; repmat(plate_no, sum(valid_idx_intS), 1)];
            combined_lag_time_group_2     = [combined_lag_time_group_2; lag_time(valid_idx_intS)];
            combined_size_group_2         = [combined_size_group_2; col_size(valid_idx_intS)];
            combined_area_group_2         = [combined_area_group_2; col_area(valid_idx_intS)];
            combined_int_group_2          = [combined_int_group_2; col_int(valid_idx_intS)];
            combined_int_per_size_2       = [combined_int_per_size_2; ints_per_size(valid_idx_intS)];
            combined_perim_group_2        = [combined_perim_group_2; col_per(valid_idx_intS)];
            combined_circ_group_2         = [combined_circ_group_2; col_circ(valid_idx_intS)];
            combined_ecc_group_2          = [combined_ecc_group_2; col_ecc(valid_idx_intS)];
            combined_solid_group_2        = [combined_solid_group_2; col_sol(valid_idx_intS)];
            combined_centroid_x_2         = [combined_centroid_x_2; col_cent_x(valid_idx_intS)];
            combined_centroid_y_2         = [combined_centroid_y_2; col_cent_y(valid_idx_intS)];

        end
    
        % Equalize vector lengths (precaution)
        n_rows = min([numel(combined_lag_time_group), numel(combined_size_group), ...
                      numel(combined_int_group), numel(combined_perim_group), ...
                      numel(combined_solid_group), numel(combined_circ_group), numel(combined_ecc_group)]);
    
        n_rows_2 = min([numel(combined_lag_time_group_2), numel(combined_int_per_size_2)]);
    
        % Create table for current group
        groupTable = table( ...
            repmat({group_label}, n_rows, 1), ...
            combined_name(1:n_rows), ...
            combined_folder(1:n_rows), ...
            combined_plate_group(1:n_rows), ...
            combined_lag_time_group(1:n_rows), ...
            combined_size_group(1:n_rows), ...
            combined_area_group(1:n_rows), ...
            combined_int_group(1:n_rows), ...
            combined_perim_group(1:n_rows), ...
            combined_circ_group(1:n_rows), ...
            combined_ecc_group(1:n_rows), ...
            combined_solid_group(1:n_rows), ...
            combined_centroid_x(1:n_rows), ...
            combined_centroid_y(1:n_rows), ...
            'VariableNames', {'SampleGroup','Name','Folder','PlateNo','LagTime','FinalSize','Area','Intensity','Perimeter','Circularity','Eccentricity','Solidity','Centroid_X','Centroid_Y'});
    
        groupTable_2 = table( ...
            repmat({group_label}, n_rows_2, 1), ...
            combined_name_2(1:n_rows_2), ...
            combined_folder_2(1:n_rows_2), ...
            combined_plate_group_2(1:n_rows_2), ...
            combined_lag_time_group_2(1:n_rows_2), ...
            combined_size_group_2(1:n_rows_2), ...
            combined_area_group_2(1:n_rows_2), ...
            combined_int_group_2(1:n_rows_2), ...
            combined_int_per_size_2(1:n_rows_2), ...
            combined_perim_group_2(1:n_rows_2), ...
            combined_circ_group_2(1:n_rows_2), ...
            combined_ecc_group_2(1:n_rows_2), ...
            combined_solid_group_2(1:n_rows_2), ...
            combined_centroid_x_2(1:n_rows_2), ...
            combined_centroid_y_2(1:n_rows_2), ...
            'VariableNames', {'SampleGroup','Name','Folder','PlateNo','LagTime','FinalSize','Area','Intensity','IntPerSize','Perimeter','Circularity','Eccentricity','Solidity','Centroid_X','Centroid_Y'});
    
        % Combine into master table
        allGroupsTable_Bright   = [allGroupsTable_Bright; groupTable];
        allGroupsTable_Bright_2 = [allGroupsTable_Bright_2; groupTable_2];
    end
    
    %% -- Export as one file --
    filename = fullfile(output_dir, '\Morphology_AllParameters_Bright.xlsx');
    writetable(allGroupsTable_Bright, filename);
    fprintf('✅ Exported combined morphology data for all groups to:\n   %s\nTotal colonies: %d\n', filename, height(allGroupsTable_Bright));
    
    filename_2 = fullfile(output_dir, '\Morphology_AllParameters_Bright_2.xlsx');
    writetable(allGroupsTable_Bright_2, filename_2);
    fprintf('✅ Exported combined morphology data for all groups to:\n   %s\nTotal colonies: %d\n', filename_2, height(allGroupsTable_Bright_2));


    %% -- Save to MATLAB workspace --
    assignin('base', 'allGroupsTable_Bright', allGroupsTable_Bright);
    mat_file1 = fullfile(output_dir, 'allGroupsTable_Bright.mat');
    save(mat_file1, 'allGroupsTable_Bright')

    assignin('base', 'allGroupsTable_2_Bright', allGroupsTable_Bright_2);
    mat_file2 = fullfile(output_dir, 'allGroupsTable_Bright_2.mat');
    save(mat_file2, 'allGroupsTable_Bright_2')
end

%% HeteroTyper Pipeline for Bright Plates
% Exports combined morphology tables to Excel and .mat files.
% Uses the same colony population as plot_combined_samples — all data
% is read directly from ht.groups; no re-filtering from raw data.
%
% Requires preprocess_pipeline_data(data) to have been run first.
%
%% USAGE:
%   export_data_to_output_directory(data, ht);

function export_data_to_output_directory(data, ht)
    HT_FLOG = -1;

    if ~isstruct(ht) || ~isfield(ht, 'groups')
        error(['Second input must be the ht struct produced by preprocess_pipeline_data.\n' ...
               'Run:  preprocess_pipeline_data(data)  first.']);
    end

    p      = ht.params;
    labels = ht.labels;

    % Output goes into a dedicated subfolder so the main output directory
    % is not cluttered by data files.
    output_dir = fullfile(p.out_dir, 'Exported Data');
    if ~exist(output_dir, 'dir'), mkdir(output_dir); end

    % ------------------------------------------------------------------
    %  Open log file
    % ------------------------------------------------------------------
    timestamp_log = datestr(now, 'yyyymmdd_HHMMSS');
    log_path = fullfile(output_dir, ...
                        sprintf('export_data_to_output_directory_%s.txt', timestamp_log));
    HT_FLOG = fopen(log_path, 'w');
    if HT_FLOG == -1
        warning('Could not open log file: %s', log_path);
        HT_FLOG = -1;
    end
    ht_fprintf(HT_FLOG, 'Log file: %s\n', log_path);
    ht_fprintf(HT_FLOG, 'Output directory: %s\n\n', output_dir);

    allGroupsTable_Bright   = table();
    allGroupsTable_Bright_2 = table();

    for g = 1:length(ht.groups)
        grp          = ht.groups(g);
        group_label  = labels{g};
        group_indices = grp.plate_indices;

        % ------------------------------------------------------------------
        %  Build per-colony metadata arrays (plate number, filename, folder).
        %  We iterate plates in the same order preprocess_pipeline_data did,
        %  counting how many colonies came from each plate by checking which
        %  ht.groups entries are non-NaN for that plate's slot range.
        %  This is the only correct way — no re-filtering from data.processed.
        % ------------------------------------------------------------------
        n_total = grp.n_colonies;
        colony_plate_no  = cell(n_total, 1);
        colony_name      = cell(n_total, 1);
        colony_folder    = cell(n_total, 1);

        slot = 0;
        for j = 1:length(group_indices)
            plate_idx = group_indices(j);
            colonies  = data.processed{plate_idx}.colonies;

            % Count how many colonies from this plate entered the accumulator
            if isfield(colonies, 'new') && isfield(colonies.new, 'lag_time')
                n_here = length(colonies.new.lag_time);
            else
                continue;
            end
            if n_here == 0, continue; end

            % Plate identifier
            pos_raw = data.metadata.original.Position(plate_idx);
            if iscell(pos_raw), plate_no = pos_raw{1}; else, plate_no = pos_raw; end

            % Filename / folder
            if isfield(data.metadata, 'fn') && length(data.metadata.fn) >= plate_idx
                f_name   = data.metadata.fn(plate_idx).name;
                f_folder = data.metadata.fn(plate_idx).folder;
            else
                f_name   = sprintf('plate_%d', plate_idx);
                f_folder = '';
            end

            idx_range = slot + (1:n_here);
            for ci = idx_range
                colony_plate_no{ci} = plate_no;
                colony_name{ci}     = f_name;
                colony_folder{ci}   = f_folder;
            end
            slot = slot + n_here;
        end

        % ------------------------------------------------------------------
        %  Table 1: all colonies that passed the size filter
        %  (non-NaN size = passed; same population as plot_combined_samples)
        % ------------------------------------------------------------------
        valid = ~isnan(grp.size);
        n_v   = sum(valid);

        groupTable = table( ...
            repmat({group_label}, n_v, 1), ...
            colony_name(valid), ...
            colony_folder(valid), ...
            colony_plate_no(valid), ...
            grp.lag_time(valid), ...
            grp.size(valid), ...
            grp.area(valid), ...
            grp.intensity(valid), ...
            grp.perimeter(valid), ...
            grp.circularity(valid), ...
            grp.eccentricity(valid), ...
            grp.solidity(valid), ...
            grp.centroid(valid, 1), ...
            grp.centroid(valid, 2), ...
            'VariableNames', {'SampleGroup','Name','Folder','PlateNo', ...
                              'LagTime','FinalSize','Area','Intensity', ...
                              'Perimeter','Circularity','Eccentricity','Solidity', ...
                              'Centroid_X','Centroid_Y'});

        % ------------------------------------------------------------------
        %  Table 2: colonies that also have a finite int_per_size
        % ------------------------------------------------------------------
        valid2 = valid & isfinite(grp.int_per_size);
        n_v2   = sum(valid2);

        groupTable_2 = table( ...
            repmat({group_label}, n_v2, 1), ...
            colony_name(valid2), ...
            colony_folder(valid2), ...
            colony_plate_no(valid2), ...
            grp.lag_time(valid2), ...
            grp.size(valid2), ...
            grp.area(valid2), ...
            grp.intensity(valid2), ...
            grp.int_per_size(valid2), ...
            grp.perimeter(valid2), ...
            grp.circularity(valid2), ...
            grp.eccentricity(valid2), ...
            grp.solidity(valid2), ...
            grp.centroid(valid2, 1), ...
            grp.centroid(valid2, 2), ...
            'VariableNames', {'SampleGroup','Name','Folder','PlateNo', ...
                              'LagTime','FinalSize','Area','Intensity','IntPerSize', ...
                              'Perimeter','Circularity','Eccentricity','Solidity', ...
                              'Centroid_X','Centroid_Y'});

        ht_fprintf(HT_FLOG, 'Group %s: %d valid colonies, %d with IntPerSize\n', ...
                   group_label, n_v, n_v2);

        allGroupsTable_Bright   = [allGroupsTable_Bright;   groupTable];   %#ok<AGROW>
        allGroupsTable_Bright_2 = [allGroupsTable_Bright_2; groupTable_2]; %#ok<AGROW>
    end

    % --- Write Excel ---
    fn1 = fullfile(output_dir, sprintf('Morphology_AllParameters_%s.xlsx', timestamp_log));
    writetable(allGroupsTable_Bright, fn1);
    ht_fprintf(HT_FLOG, '\nExported Table 1 (size-filtered) to:\n  %s\n  Total colonies: %d\n', ...
               fn1, height(allGroupsTable_Bright));

    fn2 = fullfile(output_dir, sprintf('Morphology_AllParameters_IntPerSize_%s.xlsx', timestamp_log));
    writetable(allGroupsTable_Bright_2, fn2);
    ht_fprintf(HT_FLOG, 'Exported Table 2 (size-filtered + finite IntPerSize) to:\n  %s\n  Total colonies: %d\n', ...
               fn2, height(allGroupsTable_Bright_2));

    % --- Save .mat ---
    mat1 = fullfile(output_dir, sprintf('allGroupsTable_%s.mat', timestamp_log));
    save(mat1, 'allGroupsTable_Bright', 'allGroupsTable_Bright_2');
    ht_fprintf(HT_FLOG, 'Saved .mat file to:\n  %s\n', mat1);

    assignin('base', 'allGroupsTable_Bright',   allGroupsTable_Bright);
    assignin('base', 'allGroupsTable_Bright_2', allGroupsTable_Bright_2);

    ht_fprintf(HT_FLOG, '\nExport complete.\n');
    if HT_FLOG ~= -1, fclose(HT_FLOG); fprintf('Log saved: %s\n', log_path); end
end


function ht_fprintf(HT_FLOG, varargin)
    fprintf(varargin{:});
    if ~isempty(HT_FLOG) && HT_FLOG ~= -1
        fprintf(HT_FLOG, varargin{:});
    end
end
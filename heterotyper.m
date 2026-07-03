%% 18.05.2026 - HeteroTyper Computational Pipeline

function data = heterotyper()

addpath('.\Code');

%% =========================================================================
%% INTERACTIVE PARAMETER SETUP
%% =========================================================================

fprintf('\n');
fprintf('=============================================================\n');
fprintf('         HeteroTyper — Interactive Parameter Setup          \n');
fprintf('=============================================================\n\n');

%% --- 0. Target folder (inp) ---------------------------------------------
while true
    inp = strtrim(input('[ 0 ]  Enter path to the target image folder\n       (e.g. D:\\MyExperiment\\Images): ', 's'));
    if ~isempty(inp) && exist(inp, 'dir')
        break;
    end
    if isempty(inp)
        fprintf('  [!] Path cannot be empty.\n');
    else
        fprintf('  [!] Folder not found: %s\n      Please check the path and try again.\n', inp);
    end
end
assignin('base', 'inp', inp);   % save to base workspace so downstream scripts can access it

%% --- 1. Number of plates ------------------------------------------------
while true
    nr_plates_str = input('[ 1 ]  Enter number of plates in this experiment (e.g. 12, 20, 60, etc.): ', 's');
    nr_plates_val = str2double(strtrim(nr_plates_str));
    if ~isnan(nr_plates_val) && nr_plates_val > 0 && floor(nr_plates_val) == nr_plates_val
        break;
    end
    fprintf('  [!] Please enter a positive integer (e.g. 12).\n');
end

%% --- 2. Min colony number -----------------------------------------------
while true
    min_col_str = input('[ 2 ]  Enter minimum colony count for time-course analysis (default: 5): ', 's');
    min_col_str = strtrim(min_col_str);
    if isempty(min_col_str)
        min_col_val = 5;
        break;
    end
    min_col_val = str2double(min_col_str);
    if ~isnan(min_col_val) && min_col_val >= 0 && floor(min_col_val) == min_col_val
        break;
    end
    fprintf('  [!] Please enter a non-negative integer (e.g., 5, 10, etc.).\n');
end

%% --- 3. Max colony number -----------------------------------------------
while true
    max_col_str = input('[ 3 ]  Enter maximum colony count for time-course analysis (default: 500): ', 's');
    max_col_str = strtrim(max_col_str);
    if isempty(max_col_str)
        max_col_val = 500;
        break;
    end
    max_col_val = str2double(max_col_str);
    if ~isnan(max_col_val) && max_col_val > min_col_val && floor(max_col_val) == max_col_val
        break;
    end
    fprintf('  [!] Please enter an integer greater than the minimum (%d).\n', min_col_val);
end

%% --- 4. XLSX metadata file path -----------------------------------------
while true
    xlsx_path = strtrim(input('[ 4 ]  Enter full path to XLSX metadata file\n       (e.g. D:\\MyExp\\MetaData.xlsx): ', 's'));
    if exist(xlsx_path, 'file')
        break;
    end
    fprintf('  [!] File not found: %s\n      Please check the path and try again.\n', xlsx_path);
end

%% --- 5. Sheet name in the XLSX file -------------------------------------
sheet_name = strtrim(input('[ 5 ]  Enter the sheet name to load from that file (e.g. Sheet1, Experiment_1, etc.): ', 's'));
if isempty(sheet_name)
    sheet_name = 1;   % default to first sheet
    fprintf('  [i] No sheet name entered — loading first sheet.\n');
end

%% --- 6. Metadata column for subpanel titles -----------------------------
%%     First load the table so we can show the user the available columns.
fprintf('\n  Reading metadata file to retrieve available columns...\n');
try
    meta_tmp = readtable(xlsx_path, 'sheet', sheet_name, 'VariableNamingRule', 'preserve');
    col_names = meta_tmp.Properties.VariableNames;
    fprintf('\n  Columns found in metadata:\n');
    for c = 1:numel(col_names)
        fprintf('    [%2d]  %s\n', c, col_names{c});
    end
    fprintf('\n');
    while true
        pos_col_input = strtrim(input('[ 6 ]  Enter the column name to use for subpanel titles (type exactly as shown above): ', 's'));
        if ismember(pos_col_input, col_names)
            metadata_position_col = pos_col_input;
            break;
        end
        % Allow numeric index as shortcut
        idx = str2double(pos_col_input);
        if ~isnan(idx) && idx >= 1 && idx <= numel(col_names)
            metadata_position_col = col_names{idx};
            fprintf('  [i] Selected column: %s\n', metadata_position_col);
            break;
        end
        fprintf('  [!] Column "%s" not found. Type the name exactly as listed, or its index number.\n', pos_col_input);
    end
catch ME
    fprintf('  [!] Could not read metadata file to list columns: %s\n', ME.message);
    fprintf('      Falling back to manual entry.\n');
    meta_tmp = [];
    metadata_position_col = strtrim(input('[ 6 ]  Enter the metadata column name for subpanel titles (e.g. Position): ', 's'));
    if isempty(metadata_position_col)
        metadata_position_col = 'Position';
        fprintf('  [i] No input — defaulting to "Position".\n');
    end
end

%% --- 7. Output path for figures -----------------------------------------
fprintf('\n  Figure saving options:\n');
fprintf('    0  —  Disable figure saving\n');
fprintf('    1  —  Enable figure saving  (you will define the output folder path)\n');
fprintf('\n');
while true
    save_choice = strtrim(input('[ 7 ]  Select an option (0 or 1): ', 's'));
    if strcmp(save_choice, '0')
        output_path = '';
        fprintf('  [i] Figure saving disabled.\n');
        break;
    elseif strcmp(save_choice, '1')
        while true
            output_path = strtrim(input('\n       Enter output folder path for saved figures\n       (e.g. D:\\MyExperiment\\Results\\Exp_37C): ', 's'));
            if isempty(output_path)
                fprintf('  [!] Path cannot be empty when saving is enabled.\n');
            elseif any(ismember(output_path, {'<', '>', '"', '|', '?', '*'}))
                fprintf('  [!] Path contains invalid characters. Please try again.\n');
            else
                break;   % path accepted — folder will be created below if needed
            end
        end
        break;
    else
        fprintf('  [!] Please enter 0 or 1.\n');
    end
end

%% --- Summary printout ---------------------------------------------------
fprintf('\n-------------------------------------------------------------\n');
fprintf('  Parameters confirmed:\n');
fprintf('    Target folder       : %s\n',   inp);
fprintf('    Nr. plates          : %d\n',   nr_plates_val);
fprintf('    Min colony count    : %d\n',   min_col_val);
fprintf('    Max colony count    : %d\n',   max_col_val);
fprintf('    Metadata XLSX       : %s\n',   xlsx_path);
if isnumeric(sheet_name)
    fprintf('    Sheet               : (first sheet)\n');
else
    fprintf('    Sheet               : %s\n', sheet_name);
end
fprintf('    Subpanel title col  : %s\n',   metadata_position_col);
if isempty(output_path)
    fprintf('    Output path         : (disabled)\n');
else
    fprintf('    Output path         : %s\n', output_path);
end
fprintf('-------------------------------------------------------------\n\n');

%% =========================================================================
%% ASSIGN COLLECTED VALUES TO data.params
%% =========================================================================

%% Some parameters of interest
data.params.target_folder       = inp;
data.params.border_range        = 0.85;   % only consider inner 85% of the plate

data.params.LoG_thresh          = 0.09;  % intensity threshold to detect foreground pixels
data.params.size_thresh         = 100;   % minimal colony size in pixel
data.params.eccentricity_thresh = 0.70;  % maximal eccentricity to be considered a colony

%% Near-border zone filter (annular zone from border_zone_inner*R to border_range*R)
%% Colonies whose centroid falls in this zone must pass the stricter eccentricity.
data.params.border_zone_inner = 0.75;   % inner edge of the stricter zone (fraction of plate radius)
data.params.border_ecc_thresh = 0.55;   % stricter eccentricity limit for near-border colonies

%% Merged-colony splitting
data.params.split_h_thresh  = 2;

data.params.lag_time_thresh = 100;  % minimal pixel size for lag time threshold
data.params.early_DT_range  = 5;    % fold-change range for early doubling time

%% — User-defined parameters (collected above) —
data.params.min_colony_nr          = min_col_val;
data.params.max_colony_nr          = max_col_val;
data.params.nr_plates              = nr_plates_val;
data.params.metadata_position_col  = metadata_position_col;
data.params.output_path            = output_path;

%% =========================================================================
%% REST OF PIPELINE (unchanged)
%% =========================================================================

%% Create output folder if needed
if ~isempty(data.params.output_path) && ~exist(data.params.output_path, 'dir')
    mkdir(data.params.output_path);
    fprintf('Created output folder: %s\n', data.params.output_path);
end

%% Start log file
if ~isempty(data.params.output_path)
    log_file = fullfile(data.params.output_path, ...
        ['run_log_' datestr(now,'yyyymmdd_HHMMSS') '.txt']);
    diary(log_file);
    diary on;
    fprintf('=== HeteroTyper run started: %s ===\n', datestr(now));
end

fn = dir(data.params.target_folder);
fn = fn(3:end,:);
data.metadata.fn = fn;

plot_flag              = 1;
plot_flag_segmentation = 1;
plot_flag_growth       = 1;

%% Load metadata table (already read above; reuse if successful)
if ~isempty(meta_tmp)
    data.metadata.original = meta_tmp;
else
    data.metadata.original = readtable(xlsx_path, ...
        'sheet', sheet_name, 'VariableNamingRule', 'preserve');
end

%% Pass Position labels to params so find_generic_plate_center can use them
pos_col = data.params.metadata_position_col;
if ismember(pos_col, data.metadata.original.Properties.VariableNames)
    data.params.position_labels = data.metadata.original.(pos_col);
else
    data.params.position_labels = {};   % find_generic_plate_center will fall back to folder name
    warning('Metadata table has no ''%s'' column — plate titles will use folder names.', pos_col);
end

%% Find generic plate center
data = find_generic_plate_center(data);

%% Save plate center summary figure
if ~isempty(data.params.output_path)
    fig_centers = findobj('Type','figure','Name','find plate centers');
    if ~isempty(fig_centers)
        out_file = fullfile(data.params.output_path, 'Plate_center_summary.png');
        exportgraphics(fig_centers(1), out_file, 'Resolution', 150);
        fprintf('[save] Plate center summary -> %s\n', out_file);
        close(fig_centers);
    end
end

%% Loop through each folder (one per position / time-lapse)
for i = 1:length(fn)

    tmp = [];

    % ---- Extract short position ID: "Pos0003" from folder name ------------
    % Looks for PosXXXX (exactly 4 digits). Falls back to full folder name.
    pos_token = regexp(fn(i).name, 'Pos\d{4}', 'match', 'once');
    if isempty(pos_token)
        pos_token = fn(i).name;   % fallback: use full folder name
    end
    data.params.sample_name = pos_token;   % used for figure titles etc.

    t1 = strcat(fn(i).folder, '\', fn(i).name);

    % ---- Inject per-plate center FIRST, before any processing ----
    if ~isnan(data.params.radius(i)) && data.params.radius(i) > 0
        data.params.plate_r_current      = data.params.radius(i);
        data.params.plate_center_current = data.params.center_median(i,:);
    else
        valid_r = data.params.radius(~isnan(data.params.radius));
        valid_c = data.params.center_median(~any(isnan(data.params.center_median),2),:);
        data.params.plate_r_current      = median(valid_r);
        data.params.plate_center_current = median(valid_c, 1);
        fprintf('  [plate %d] WARNING: center detection failed, using cross-plate median\n', i);
    end
    fprintf('  [plate %d] center=(%.1f,%.1f) r=%.1f\n', i, ...
        data.params.plate_center_current(1), data.params.plate_center_current(2), ...
        data.params.plate_r_current);

    % Load first+last frames then full stack
    tmp.raw  = load_image_first_last(t1, fn(i).name);
    tmp.raw2 = load_image_stack(t1, fn(i).name);

    % Background subtraction and crop
    tmp.processed = post_processing_image_stack(tmp.raw,  data.params, plot_flag);

    % Detect colonies on full stack (last frame used for segmentation)
    tmp_full      = post_processing_image_stack(tmp.raw2, data.params, plot_flag);
    tmp.processed = tmp_full;

    tmp.processed.colonies = colony_detection_and_growth_quantification( ...
        tmp_full, fn(i).name, tmp.raw2.filename, data.params, plot_flag_segmentation);

    %% -----------------------------------------------------------------------
    %% Save figures with SHORT filename:  PosXXXX_Step1_debug.png  etc.
    %%
    %%   Step 1: Colony Detection Debug
    %%   Step 2: colony detection  (raw, pre-filter)
    %%   Step 3: colony cleanup    (white=kept  green=removed)
    %%   Step 4: segmented         (mask on gray overlay)
    %% -----------------------------------------------------------------------
    if ~isempty(data.params.output_path)

        % Step 1 — 8-panel debug
        fig_debug = findobj('Type','figure','Name','Colony Detection Debug');
        if ~isempty(fig_debug)
            out_f = fullfile(data.params.output_path, [pos_token '_Step1_debug.png']);
            exportgraphics(fig_debug(1), out_f, 'Resolution', 150);
            fprintf('  [save] Debug     -> %s\n', [pos_token '_Step1_debug.png']);
            close(fig_debug);
        end

        % Step 2 — raw detection overlay
        fig_overlay = findobj('Type','figure','Name', strcat(fn(i).name,' - colony detection'));
        if ~isempty(fig_overlay)
            out_f = fullfile(data.params.output_path, [pos_token '_Step2_detection.png']);
            exportgraphics(fig_overlay(1), out_f, 'Resolution', 150);
            fprintf('  [save] Detection -> %s\n', [pos_token '_Step2_detection.png']);
            close(fig_overlay);
        end

        % Step 3 — cleanup (white=kept, green=removed)
        fig_clean = findobj('Type','figure','Name', strcat(fn(i).name,' - colony cleanup'));
        if ~isempty(fig_clean)
            out_f = fullfile(data.params.output_path, [pos_token '_Step3_cleanup.png']);
            exportgraphics(fig_clean(1), out_f, 'Resolution', 150);
            fprintf('  [save] Cleanup   -> %s\n', [pos_token '_Step3_cleanup.png']);
            close(fig_clean(1));
        end

        % Step 4 — segmented on gray overlay
        fig_seg = findobj('Type','figure','Name', strcat(fn(i).name,' - segmented'));
        if ~isempty(fig_seg)
            out_f = fullfile(data.params.output_path, [pos_token '_Step4_segmented.png']);
            exportgraphics(fig_seg(1), out_f, 'Resolution', 150);
            fprintf('  [save] Segmented -> %s\n', [pos_token '_Step4_segmented.png']);
            close(fig_seg);
        end
    end

    %% flag_colony_ok is already set inside colony_detection_and_growth_quantification
    if size(tmp.processed.colonies.region_props, 1) == 0
        tmp.processed.colonies.flag_colony_ok = [];
        tmp.processed.colonies.mask_clean     = tmp.processed.colonies.mask;
    end

    %% Only quantify time courses if final colony count is within range
    nr_ok_colonies = sum(tmp.processed.colonies.flag_colony_ok == 1);

    if (nr_ok_colonies > data.params.min_colony_nr) && (nr_ok_colonies < data.params.max_colony_nr)
        tmp2.raw                         = load_image_stack(t1, fn(i).name);
        tmp2.processed                   = post_processing_image_stack(tmp2.raw, data.params, plot_flag);
        tmp2.processed.colonies          = tmp.processed.colonies;
        tmp2.processed.time_info         = extract_time_info(tmp2.raw.filename);
        tmp.processed.time_info          = tmp2.processed.time_info;
        data.processed{i,1}.time         = tmp.processed.time_info;
        data.processed{i,1}.growth_quant = 1;
    else
        data.processed{i,1}.growth_quant = 0;
    end

    data.processed{i,1}.colonies   = tmp.processed.colonies;
    data.processed{i,1}.center_pos = tmp.processed.center_pos;
    data.processed{i,1}.radius     = tmp.processed.radius;
    data.processed{i,1}.img_final  = tmp.processed.img{end};

end

%% Extract early doubling time
data = extract_early_growth(data);

%% Close log file
if ~isempty(data.params.output_path)
    fprintf('\n=== HeteroTyper run finished: %s ===\n', datestr(now));
    diary off;
end

end
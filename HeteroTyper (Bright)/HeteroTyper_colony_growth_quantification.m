%% 11.11.2025 - HeteroTyper Pipeline for Bright Plates

% test pipeline to analyze colony imaging time courses 
% 1. analysis performed separately for each colony
% 2. assumption: all images for given position are in the same folder and
% are sequential

function data = HeteroTyper_colony_growth_quantification(inp)

addpath('.\Code');

%% some parameters of interest
data.params.target_folder = inp; % where the pipeline expects the individual folders, one per position
data.params.border_range = 0.85; % only consider inner 85% of the plate

data.params.LoG_thresh = 0.09; % intensity threshold to detect foreground pixels
data.params.size_thresh = 100; % minimal colony size in pixel

data.params.eccentricity_thresh = 0.70; % maximal eccentricity to be considered a colony

data.params.lag_time_thresh = 100; % minimal pixel size for lag time threshold
data.params.early_DT_range = 5; % fold-change range for early double time: DT = time(DT_range*lag_time_thresh) - time(lag_time_thresh)

% only plates whose colony counts are between these two numbers will be
% analyzed over time
data.params.min_colony_nr = 10; % minimal number of clean colonies needed in plate
data.params.max_colony_nr = 650; % maximal number of clean colonies allowed in plate


fn = dir(data.params.target_folder);
fn = fn(3:end,:);
data.metadata.fn = fn;

plot_flag = 1;

plot_flag_preprocessing = 0;   % for pre-processing
plot_flag_segmentation = 1;    % for colony segmentation
plot_flag_growth = 1;

%% Load manual colony count (this part you will have to remove if you analyze other data sets
data.metadata.original = readtable('D:\Gizem\HeteroTyper\Metadata.xlsx','sheet','Meta_Bright');

%% find generic plate center
data = find_generic_plate_center(data);

%% Loop through each folder within target folder (assuming that time stacks are stored together in a folder for each position)

for i = 1:length(fn)
    
    tmp = [];
    
    % load image stack
    t1 = strcat(fn(i).folder,'\',fn(i).name); % get pointer
    tmp.raw = load_image_first_last(t1,fn(i).name);
    tmp.raw2 = load_image_stack(t1,fn(i).name);
    
    % postprocessing of images (background subtract and crop)
    tmp.processed = post_processing_image_stack(tmp.raw, data.params, plot_flag);

    % detect colonies
    tmp_full = post_processing_image_stack(tmp.raw2, data.params, plot_flag);
    tmp.processed = tmp_full; % keep downstream expectations unchanged
    tmp.processed.colonies = colony_detection_and_growth_quantification(tmp_full, fn(i).name, tmp.raw2.filename, data.params, plot_flag_segmentation);

    
    if(size(tmp.processed.colonies.region_props,1) > 0)
        % flag outlier colonies
        [tmp.processed.colonies.flag_colony_ok, tmp.processed.colonies.mask_clean] = colony_artifact_flag(tmp.processed.img{end},tmp.processed.colonies,tmp.processed.radius,tmp.processed.center_pos,data.params,plot_flag);
    else
        tmp.processed.colonies.flag_colony_ok = [];
        tmp.processed.colonies.mask_clean = tmp.processed.colonies.mask;
        tmp.processed.colonies.debug = tmp.processed.colonies.debug;
    end
    
    %% to avoid useless calculations: only quantify time courses if plate is suitable 
    nr_ok_colonies = length(tmp.processed.colonies.flag_colony_ok);
    
    if(nr_ok_colonies > data.params.min_colony_nr)&&(nr_ok_colonies < data.params.max_colony_nr)
        % load full image stack
        tmp2.raw = load_image_stack(t1,fn(i).name);
        % process images
        tmp2.processed = post_processing_image_stack(tmp2.raw, data.params, plot_flag);
        % get colony mask
        tmp2.processed.colonies = tmp.processed.colonies;
        tmp2.processed.time_info = extract_time_info(tmp2.raw.filename);
        tmp.processed.time_info = tmp2.processed.time_info;
       
        % extract colony parameter distributions
        tmp.processed.colonies.growth_params = growth_feat_distribution(tmp.processed.colonies ,tmp.processed.time_info.elapsed_time_h, data.params, fn(i).name, plot_flag_growth);
        data.processed{i,1}.time = tmp.processed.time_info;
        data.processed{i,1}.growth_quant = 1;
    else
        data.processed{i,1}.growth_quant = 0;
    end
    
    data.processed{i,1}.colonies = tmp.processed.colonies;
    
    data.processed{i,1}.center_pos = tmp.processed.center_pos;
    data.processed{i,1}.radius = tmp.processed.radius;
    
    data.processed{i,1}.img_final = tmp.processed.img{end};
    
end

%% extract early doubling time
data = extract_early_growth(data);

end
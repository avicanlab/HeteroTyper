%% 11.11.2025 - HeteroTyper Pipeline for Dark Plates 

function data = plot_lagTime_vs_morphology(data) 

    nr_plates = length(data.processed); 
    
    % Define maximum lag time
    max_lag = 55;

    % Define Room Temperature Incubation Time
    incTime = 24;
    
    % Define min-max colony numbers
    min_col = 10;
    max_col = 350;

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
    
    
    %% split things by organ
    ix_MLN = find(strcmpi(data.metadata.original.Organ,'MLN')); 
    ix_Spleen = find(strcmpi(data.metadata.original.Organ,'Spleen')); 
    ix_Liver = find(strcmpi(data.metadata.original.Organ,'Liver'));  
    labels_type = {'MLN','Spleen','Liver'}; 
    
    ix1 = ix_MLN; 
    ix1 = intersect(ix1,ix_growth); 
    ix2 = ix_Spleen; 
    ix2 = intersect(ix2,ix_growth); 
    ix3 = ix_Liver; 
    ix3 = intersect(ix3,ix_growth); 
    
    ix{1} = ix1; 
    ix{2} = ix2; 
    ix{3} = ix3; 
    
    colors = [0.37 0.21 0.65;0.12 0.69 0.70;0.87 0.71 0]; 

    %% Parameters to plot
    param_names = {'FinalSize','Area','Intensity','MeanIntensity','IntensityPerSize','Perimeter','Circularity','Eccentricity','Solidity'};


    % Define filter parameters
    size_threshold = 200; 
    ecc_threshold = 0.70;

    
    figure('Name','Combined plots'); 
    for i = 1:length(ix) 

        ix_tmp = ix{i}; 
        group_imgs = {};
        group_labels = {};

        % Generate empty combined arrays
        combined_lag_time      = []; 
        combined_size          = []; 
        combined_int           = []; 
        combined_int_per_size  = [];
        combined_mean_int      = []; 
        combined_area          = []; 
        combined_solidity      = []; 
        combined_circularity   = []; 
        combined_perimeter     = []; 
        combined_eccentricity  = []; 
        combined_centroid      = []; 

        time = [];
            
        % Initialize tracker before the for j=1:length(ix_tmp) loop
        mismatch_report = [];  
    
        for j = 1:length(ix_tmp) 
            colonies = data.processed{ix_tmp(j)}.colonies.new;
            props    = data.processed{ix_tmp(j)}.colonies.region_props;


            col_size        = colonies.timecourse_size_smoothed(end,:)'; 
            col_int         = colonies.timecourse_intensity_smoothed(end,:)';
            lag_time        = colonies.lag_time + incTime; 
            col_area        = props.Area(:);
            mean_int        = props.MeanIntensity(:); 
            perimeter       = props.Perimeter(:); 
            circularity     = props.Circularity(:); 
            eccentricity    = props.Eccentricity(:); 
            solidity        = props.Solidity(:); 
            centroids       = props.Centroid(:); 
        
            time            = colonies.time_info.elapsed_time_h + incTime; 
            col_flag        = props.flag_colony_ok(:);

            % Generate combined sample group tables 
            % Final colony size and eccentricity filtered, and "non-infinite" intensity per size data will be used!
            ints_per_size    = col_int ./ col_size;


            % Keep only valid entries (nonzero, finite, flagged colonies)
            valid_idx_intS = (col_flag == 1) & isfinite(ints_per_size) & (col_size > size_threshold) & (eccentricity < ecc_threshold);
            
            % Generate combined tables for growth parameters
            combined_lag_time      = [combined_lag_time; lag_time(valid_idx_intS)];
            combined_size          = [combined_size; col_size(valid_idx_intS)];
            combined_area          = [combined_area; col_area(valid_idx_intS)];
            combined_int           = [combined_int; col_int(valid_idx_intS)];
            combined_mean_int      = [combined_mean_int; mean_int(valid_idx_intS)];
            combined_int_per_size  = [combined_int_per_size; ints_per_size(valid_idx_intS)];
            combined_perimeter     = [combined_perimeter; perimeter(valid_idx_intS)];
            combined_circularity   = [combined_circularity; circularity(valid_idx_intS)];
            combined_eccentricity  = [combined_eccentricity; eccentricity(valid_idx_intS)];
            combined_solidity      = [combined_solidity; solidity(valid_idx_intS)];
            combined_centroid      = [combined_centroid; centroids(valid_idx_intS, :)];
           
            
          end 
            
            % Define y-axis limits for each parameter
            max_size        = max(combined_size);
            max_area        = max(combined_area);
            max_int         = max(combined_int);
            max_mean_int    = max(combined_mean_int);
            max_intSize     = max(combined_int_per_size);
            max_peri        = max(combined_perimeter);
            % max_circ        = max(combined_circularity);
            % max_ecc         = max(combined_eccentricity);
            % max_col         = max(combined_solidity);

            %% Define y-axis limits for each parameter
            y_limits.FinalSize         = [0 max_size];
            y_limits.Area              = [0 max_area];
            y_limits.Intensity         = [0 max_int];
            y_limits.MeanIntensity     = [0 max_mean_int];
            y_limits.IntensityPerSize  = [0 max_intSize];
            y_limits.Perimeter         = [0 max_peri];
            y_limits.Circularity       = [0 1];
            y_limits.Eccentricity      = [0 1];
            y_limits.Solidity          = [0 1];



            % Prepare parameter data
            param_data = {combined_size, combined_area, ...
                          combined_int, combined_mean_int, combined_int_per_size,...
                          combined_perimeter, combined_circularity, ...
                          combined_eccentricity, combined_solidity};

            
            % Initialize axis handle storage
            ax_handles = cell(1, length(param_names));
        
            % Plot scatter for each morphology vs LagTime
            for p = 1:length(param_names)
                subplot(length(ix), length(param_names), (i-1)*length(param_names) + p);
                scatter(combined_lag_time, param_data{p}, 10, colors(i,:), 'filled','MarkerFaceAlpha',0.4);
                
                if i == 1
                    title(param_names{p}, 'Interpreter','none');
                end
                if p == 1
                    ylabel(labels_type{i}, 'Interpreter','none');
                end
                xlabel('LagTime (h)');
                xlim([incTime max_lag]); 
                xticks(incTime:4:max_lag);
                set(gca,'FontSize',12);

                % Apply y-limits if defined
                if isfield(y_limits, param_names{p})
                    ylim(y_limits.(param_names{p}));
                end

                % store axis handle for this parameter
                ax_handles{p} = [ax_handles{p}, gca];
            end

    end 


end

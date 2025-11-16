%% 11.11.2025 - HeteroTyper Pipeline for Dark Plates

function data = plot_combined_samples_growth_curves(data) 

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
    

    % Following parameters (max_val, x_val, etc.) must be defined after identifying the maximum values for each parameter. 
    % For these parameters to be defined, code must be run once.
    max_val         = 10800; 
    xval            = [0:50:max_val]; % colony size bins 
    xval2           = [0:0.5:max_lag]; % lag time bins 
    
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
        

    % ---------- PROGRESS BAR SETUP ----------
    total_tasks = sum(cellfun(@length, ix));  % total number of plates
    task_count = 0;


    % Define filter parameters
    size_threshold = 200; 
    ecc_threshold = 0.70;
    
    
    figure('Name','Combined plots'); 
    for i = 1:length(ix) 

        ix_tmp = ix{i}; 
        group_imgs = {};
        group_labels = {};
      
        % Generate empty combined arrays
        combined_size          = []; 
        combined_size2         = []; % generate combined colony size matrix for growth curves 
        combined_eccentricity  = []; 

        % Initialize tracker before the for j=1:length(ix_tmp) loop
        mismatch_report = [];  
    
        for j = 1:length(ix_tmp)  
            % increment global task count
            task_count = task_count + 1;

            colonies = data.processed{ix_tmp(j)}.colonies.new;
            props    = data.processed{ix_tmp(j)}.colonies.region_props;


            col_size        = colonies.timecourse_size_smoothed(end,:)'; 
            col_size2       = colonies.timecourse_size_smoothed;            % Size values for all time points included (growth curves)
            col_ecc         = props.Eccentricity(:); 
        
            t               = colonies.time_info.elapsed_time_h + incTime; 
            col_flag        = props.flag_colony_ok(:);


            % Keep only valid entries (nonzero, finite, flagged colonies)
            valid_idx = (col_flag == 1) & (col_size > size_threshold) & (col_ecc < ecc_threshold);
            
            % Generate combined tables for growth parameters
            combined_size2         = [combined_size2; col_size2(:,valid_idx)];
        
        end 
            
            
            % Print mismatch summary after loop
            if ~isempty(mismatch_report)
                fprintf('\n⚠️ Colony count mismatches detected:\n');
                fprintf(' PlateIndex\tcol_size\tEccentricity\tflag_colony_ok\n');
                for k = 1:size(mismatch_report,1)
                    fprintf(' %d\t\t%d\t\t%d\t\t%d\n', ...
                        mismatch_report{k,1}, mismatch_report{k,2}, mismatch_report{k,3}, mismatch_report{k,4});
                end
            else
                fprintf('\n✅ No mismatches between col_size, eccentricity, and flag_colony_ok.\n');
            end

            
            fprintf('Size (Over Time) -- Max: %.2f\n', max(combined_size2(end,:)'));
            max_val_y = max(combined_size2(end,:)') + 100;
            t_max = max(t);

           
            % First column: Growth Curves 
            subplot(length(ix), 1, 1 * (i - 1) + 1),... 
            plot(t,combined_size2,'-'); 
            hold on; 
            axis([incTime t_max 0 max_val_y]);
            xlabel('Time (h)'); 
            ylabel('Colony Size (px)'); 
            set(gca, 'FontSize', 13); 

    end 

end

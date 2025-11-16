%% 11.11.2025 - HeteroTyper Pipeline for Dark Plates

function data = plot_combined_samples(data) 

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
    xval3           = [0:0.5:10]; % early DT bins 
    xval_area       = [0:50:10800];
    xval_int        = [0:0.1e5:9.6e5];
    xval_mean_int   = [0:0.05:100.05];
    xval_intSize    = [0:0.5:105];
    xval_peri       = [0:10:750];
    xval_solid      = [0:0.01:1.01];
    xval_circ       = [0:0.01:1.01];
    xval_eccen      = [0:0.01:1.01];
    
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

        combined_flag = []; 
        lag_time_hist_table = []; 
        size_hist_table = []; 
        int_hist_table = []; 
        combined_counts = []; 

        % Initialize tracker before the for j=1:length(ix_tmp) loop
        mismatch_report = [];  
    
        for j = 1:length(ix_tmp) 
            % increment global task count
            task_count = task_count + 1;

            colonies = data.processed{ix_tmp(j)}.colonies.new;
            props    = data.processed{ix_tmp(j)}.colonies.region_props;


            col_size        = colonies.timecourse_size_smoothed(end,:)'; 
            col_size2       = colonies.timecourse_size_smoothed;            % Size values for all time points included (growth curves)
            col_int         = colonies.timecourse_intensity_smoothed(end,:)';
            lag_time        = colonies.lag_time + incTime; 
            col_area        = props.Area(:);
            mean_int        = props.MeanIntensity(:); 
            col_peri        = props.Perimeter(:); 
            col_circ        = props.Circularity(:); 
            col_ecc         = props.Eccentricity(:); 
            col_sol         = props.Solidity(:); 
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
            combined_perimeter     = [combined_perimeter; col_peri(valid_idx_intS)];
            combined_circularity   = [combined_circularity; col_circ(valid_idx_intS)];
            combined_eccentricity  = [combined_eccentricity; col_ecc(valid_idx_intS)];
            combined_solidity      = [combined_solidity; col_sol(valid_idx_intS)];
            combined_centroid      = [combined_centroid; centroids(valid_idx_intS, :)]; 


            % === PROBABILITY DENSITY FUNCTION ===

            % Lag Time
            [counts, edges] = histcounts(combined_lag_time, xval2); % Compute histogram 
            if (j == 1) 
            combined_counts = counts; 
            else 
                combined_counts = combined_counts + counts; 
            end 
            bin_starts = edges(1:end-1); % Compute bin start points 
            total_cells = length(combined_lag_time); % Total number of cells 
            relative_freq = counts / total_cells; % Relative frequency 
            lag_time_hist_table = table(counts(:), bin_starts(:), repmat(total_cells, length(counts), 1), relative_freq(:), 'VariableNames', {'Count', 'BinStart', 'TotalCells', 'RelativeFrequency'}); % Create table 
            %disp(lag_time_hist_table); 
            

            % Size 
            [counts2, edges2] = histcounts(combined_size, xval); % Compute histogram 
            if (j == 1) 
                combined_counts2 = counts2; 
            else 
                combined_counts2 = combined_counts2 + counts2; 
            end 
            bin_starts2 = edges2(1:end-1); % Compute bin start points
            total_cells2 = length(combined_size); % Total number of cells
            relative_freq2 = counts2 / total_cells2; % Relative frequency
            size_hist_table = table(counts2(:), bin_starts2(:), repmat(total_cells2, length(counts2), 1), relative_freq2(:), 'VariableNames', {'Count', 'BinStart', 'TotalCells', 'RelativeFrequency'}); % Create table
            %disp(size_hist_table);
            

            % Intensity
            [counts3, edges3] = histcounts(combined_int, xval_int); % Compute histogram
            if (j == 1)
                combined_counts3 = counts3; 
            else 
                combined_counts3 = combined_counts3 + counts3; 
            end 
            bin_starts3 = edges3(1:end-1); % Compute bin start points
            total_cells3 = length(combined_int); % Total number of cells
            relative_freq3 = counts3 / total_cells3; % Relative frequency
            int_hist_table = table(counts3(:), bin_starts3(:), repmat(total_cells3, length(counts3), 1), relative_freq3(:), 'VariableNames', {'Count', 'BinStart', 'TotalCells', 'RelativeFrequency'}); % Create table
            % disp(int_hist_table); 
            

            % Intensity Per Size 
            [counts4, edges4] = histcounts(combined_int_per_size, xval_intSize); % Compute histogram 
            if (j == 1) 
                combined_counts4 = counts4; 
            else 
                combined_counts4 = combined_counts4 + counts4; 
            end 
            bin_starts4 = edges4(1:end-1); % Compute bin start points
            total_cells4 = length(combined_int_per_size); % Total number of cells
            relative_freq4 = counts4 / total_cells4; % Relative frequency
            intSize_hist_table = table(counts4(:), bin_starts4(:), repmat(total_cells4, length(counts4), 1), relative_freq4(:), 'VariableNames', {'Count', 'BinStart', 'TotalCells', 'RelativeFrequency'}); % Create table
            % disp(intSize_hist_table);
            

            % Perimeter
            [counts5, edges5] = histcounts(combined_perimeter, xval_peri); % Compute histogram
            if (j == 1)
                combined_counts5 = counts5; 
            else 
                combined_counts5 = combined_counts5 + counts5; 
            end 
            bin_starts5 = edges5(1:end-1); % Compute bin start points
            total_cells5 = length(combined_perimeter); % Total number of cells
            relative_freq5 = counts5 / total_cells5; % Relative frequency
            perim_hist_table = table(counts5(:), bin_starts5(:), repmat(total_cells5, length(counts5), 1), relative_freq5(:), 'VariableNames', {'Count', 'BinStart', 'TotalCells', 'RelativeFrequency'}); % Create table
            % disp(perim_hist_table);


            % Solidity
            [counts6, edges6] = histcounts(combined_solidity, xval_solid); % Compute histogram
            if (j == 1)
                combined_counts6 = counts6; 
            else 
                combined_counts6 = combined_counts6 + counts6; 
            end 
            bin_starts6 = edges6(1:end-1); % Compute bin start points
            total_cells6 = length(combined_solidity); % Total number of cells
            relative_freq6 = counts6 / total_cells6; % Relative frequency
            solid_hist_table = table(counts6(:), bin_starts6(:), repmat(total_cells6, length(counts6), 1), relative_freq6(:), 'VariableNames', {'Count', 'BinStart', 'TotalCells', 'RelativeFrequency'}); % Create table
            % disp(solid_hist_table); 


            % Circularity
            [counts7, edges7] = histcounts(combined_circularity, xval_circ); % Compute histogram 
            if (j == 1) 
                combined_counts7 = counts7; 
            else 
                combined_counts7 = combined_counts7 + counts7; 
            end 
            bin_starts7 = edges7(1:end-1); % Compute bin start points
            total_cells7 = length(combined_circularity); % Total number of cells
            relative_freq7 = counts7 / total_cells7; % Relative frequency
            circ_hist_table = table(counts7(:), bin_starts7(:), repmat(total_cells7, length(counts7), 1), relative_freq7(:), 'VariableNames', {'Count', 'BinStart', 'TotalCells', 'RelativeFrequency'}); % Create table
            % disp(circ_hist_table);


            % Eccentricity
            [counts8, edges8] = histcounts(combined_eccentricity, xval_eccen); % Compute histogram
            if (j == 1)
                combined_counts8 = counts8; 
            else 
                combined_counts8 = combined_counts8 + counts8; 
            end 
            bin_starts8 = edges8(1:end-1); % Compute bin start points
            total_cells8 = length(combined_eccentricity); % Total number of cells
            relative_freq8 = counts8 / total_cells8; % Relative frequency
            eccen_hist_table = table(counts8(:), bin_starts8(:), repmat(total_cells8, length(counts8), 1), relative_freq8(:), 'VariableNames', {'Count', 'BinStart', 'TotalCells', 'RelativeFrequency'}); % Create table
            % disp(eccen_hist_table);


            % Area
            [counts9, edges9] = histcounts(combined_area, xval_area); % Compute histogram
            if (j == 1)
                combined_counts9 = counts9; 
            else 
                combined_counts9 = combined_counts9 + counts9; 
            end 
            bin_starts9 = edges9(1:end-1); % Compute bin start points
            total_cells9 = length(combined_area); % Total number of cells
            relative_freq9 = counts9 / total_cells9; % Relative frequency
            area_hist_table = table(counts9(:), bin_starts9(:), repmat(total_cells9, length(counts9), 1), relative_freq9(:), 'VariableNames', {'Count', 'BinStart', 'TotalCells', 'RelativeFrequency'}); % Create table
            % disp(area_hist_table);


            % Mean Gray Intensity
            [counts10, edges10] = histcounts(combined_mean_int, xval_mean_int); % Compute histogram
            if (j == 1)
                combined_counts10 = counts10; 
            else 
                combined_counts10 = combined_counts10 + counts10; 
            end 
            bin_starts10 = edges10(1:end-1); % Compute bin start points
            total_cells10 = length(combined_mean_int); % Total number of cells
            relative_freq10 = counts10 / total_cells10; % Relative frequency
            meanInt_hist_table = table(counts10(:), bin_starts10(:), repmat(total_cells10, length(counts10), 1), relative_freq10(:), 'VariableNames', {'Count', 'BinStart', 'TotalCells', 'RelativeFrequency'}); % Create table
            % disp(meanInt_hist_table);
        end 
            

            %% === Export "LAG TIMES" to MATLAB base workspace ===
            varName = sprintf('combined_lag_time_%s', labels_type{i});
            assignin('base', varName, combined_lag_time);


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


            
            fprintf('Lag Time       -- Min: %.2f, Max: %.2f\n', min(combined_lag_time), max(combined_lag_time));
            fprintf('Size           -- Min: %.2f, Max: %.2f\n', min(combined_size), max(combined_size));
            fprintf('Area           -- Min: %.2f, Max: %.2f\n', min(combined_area), max(combined_area));
            fprintf('Intensity      -- Min: %.2f, Max: %.2f\n', min(combined_int), max(combined_int));
            fprintf('Mean Intensity -- Min: %.2f, Max: %.2f\n', min(combined_mean_int), max(combined_mean_int));
            fprintf('IntSize        -- Min: %.2f, Max: %.2f\n', min(combined_int_per_size), max(combined_int_per_size));
            fprintf('Perimeter      -- Min: %.2f, Max: %.2f\n', min(combined_perimeter), max(combined_perimeter));
            fprintf('Circularity    -- Min: %.2f, Max: %.2f\n', min(combined_circularity), max(combined_circularity));
            fprintf('Eccentricity   -- Min: %.2f, Max: %.2f\n', min(combined_eccentricity), max(combined_eccentricity));
            fprintf('Solidity       -- Min: %.2f, Max: %.2f\n', min(combined_solidity), max(combined_solidity));



            % Lag Time 
            subplot(length(ix), 10, 10 * (i - 1) + 1),... 
            xval_plot = 0.25:0.5:54.75; 
            bar(xval_plot,combined_counts./sum(combined_counts),'FaceColor',colors(i,:),'EdgeColor','none','FaceAlpha',1); 
            hold on; 
            plot(xval_plot, combined_counts./sum(combined_counts),'Color', 'k', 'LineWidth', 0.3);  % 'k' for black line 
            axis([incTime inf 0 max(combined_counts./sum(combined_counts))]); 
            xticks(incTime:4:56); 
            line([median(combined_lag_time) median(combined_lag_time)],[0 max(combined_counts./sum(combined_counts))],'Color','r','LineWidth',2);
            set(gca, 'FontSize', 11); 

            % Final Colony Size
            subplot(length(ix), 10, 10 * (i - 1) + 2),... 
            xval_plot2 = 0:50:10750; 
            bar(xval_plot2,combined_counts2./sum(combined_counts2),'FaceColor',colors(i,:),'EdgeColor','none','FaceAlpha',1); 
            hold on; 
            plot(xval_plot2, combined_counts2./sum(combined_counts2),'Color', 'k', 'LineWidth', 0.3); % 'k' for black line 
            axis([0 inf 0 max(combined_counts2./sum(combined_counts2))]); 
            line([median(combined_size) median(combined_size)],[0 max(combined_counts2./sum(combined_counts2))],'Color','r','LineWidth',3); 
            set(gca, 'FontSize', 11); 

            % Area
            subplot(length(ix), 10, 10 * (i - 1) + 3),... 
            xval_plot3 = 0:50:10750; 
            bar(xval_plot3,combined_counts9./sum(combined_counts9),'FaceColor',colors(i,:),'EdgeColor','none','FaceAlpha',1); 
            hold on; 
            plot(xval_plot3, combined_counts9./sum(combined_counts9),'Color', 'k', 'LineWidth', 0.3);  % 'k' for black line 
            axis([0 inf 0 max(combined_counts9./sum(combined_counts9))]); 
            line([median(combined_area) median(combined_area)],[0 max(combined_counts9./sum(combined_counts9))],'Color','r','LineWidth',1.5);
            set(gca, 'FontSize', 11); 

            % Intensity
            subplot(length(ix), 10, 10 * (i - 1) + 4),...
            xval_plot4 = 0:0.1e5:9.5e5; 
            bar(xval_plot4,combined_counts3./sum(combined_counts3),'FaceColor',colors(i,:),'EdgeColor','none','FaceAlpha',1); 
            hold on; 
            plot(xval_plot4, combined_counts3./sum(combined_counts3),'Color', 'k', 'LineWidth', 0.3); % 'k' for black line 
            axis([0 inf 0 max(combined_counts3./sum(combined_counts3))]); 
            line([median(combined_int) median(combined_int)],[0 max(combined_counts3./sum(combined_counts3))],'Color','r','LineWidth',1.5); 
            set(gca, 'FontSize', 11); 

            % Mean Gray Intensity
            subplot(length(ix), 10, 10 * (i - 1) + 5),... 
            xval_plot5 = 0:0.05:100;  
            bar(xval_plot5,combined_counts10./sum(combined_counts10),'FaceColor',colors(i,:),'EdgeColor','none','FaceAlpha',1); 
            hold on; 
            plot(xval_plot5, combined_counts10./sum(combined_counts10),'Color', 'k', 'LineWidth', 0.3); % 'k' for black line 
            axis([16 inf 0 max(combined_counts10./sum(combined_counts10))]); 
            line([median(combined_mean_int) median(combined_mean_int)],[0 max(combined_counts10./sum(combined_counts10))],'Color','r','LineWidth',1.5); 
            set(gca, 'FontSize', 11); 

            % Intensity Per Size
            subplot(length(ix), 10, 10 * (i - 1) + 6),... 
            xval_plot6 = 0.25:0.5:104.75; 
            bar(xval_plot6,combined_counts4./sum(combined_counts4),'FaceColor',colors(i,:),'EdgeColor','none','FaceAlpha',1); 
            hold on; 
            plot(xval_plot6, combined_counts4./sum(combined_counts4),'Color', 'k', 'LineWidth', 0.3);  % 'k' for black line 
            axis([0 inf 0 max(combined_counts4./sum(combined_counts4))]);  
            line([median(combined_int_per_size) median(combined_int_per_size)],[0 max(combined_counts4./sum(combined_counts4))],'Color','r','LineWidth',1.5);
            set(gca, 'FontSize', 11); 

            % Perimeter
            subplot(length(ix), 10, 10 * (i - 1) + 7),... 
            xval_plot7 = 0:10:740; 
            bar(xval_plot7,combined_counts5./sum(combined_counts5),'FaceColor',colors(i,:),'EdgeColor','none','FaceAlpha',1); 
            hold on; 
            plot(xval_plot7, combined_counts5./sum(combined_counts5),'Color', 'k', 'LineWidth', 0.3); % 'k' for black line _V2
            axis([0 inf 0 max(combined_counts5./sum(combined_counts5))]);  
            line([median(combined_perimeter) median(combined_perimeter)],[0 max(combined_counts5./sum(combined_counts5))],'Color','r','LineWidth',1.5); 
            set(gca, 'FontSize', 11); 

            % Circularity
            subplot(length(ix), 10, 10 * (i - 1) + 8),... 
            xval_plot8 = 0:0.01:1; 
            bar(xval_plot8,combined_counts7./sum(combined_counts7),'FaceColor',colors(i,:),'EdgeColor','none','FaceAlpha',1); 
            hold on; 
            plot(xval_plot8, combined_counts7./sum(combined_counts7),'Color', 'k', 'LineWidth', 0.3); % 'k' for black line 
            axis([0 inf 0 max(combined_counts7./sum(combined_counts7))]); 
            line([median(combined_circularity) median(combined_circularity)],[0 max(combined_counts6./sum(combined_counts6))],'Color','r','LineWidth',1.5); 
            set(gca, 'FontSize', 11); 

            % Eccentricity
            subplot(length(ix), 10, 10 * (i - 1) + 9),... 
            xval_plot9 = 0:0.01:1; 
            bar(xval_plot9,combined_counts8./sum(combined_counts9),'FaceColor',colors(i,:),'EdgeColor','none','FaceAlpha',1); 
            hold on; 
            plot(xval_plot9, combined_counts8./sum(combined_counts8),'Color', 'k', 'LineWidth', 0.3); % 'k' for black line 
            axis([0 inf 0 max(combined_counts8./sum(combined_counts8))]); 
            line([median(combined_eccentricity) median(combined_eccentricity)],[0 max(combined_counts8./sum(combined_counts8))],'Color','r','LineWidth',1.5); 
            set(gca, 'FontSize', 11); 

            % Solidity
            subplot(length(ix), 10, 10 * (i - 1) + 10),... 
            xval_plot10 = 0:0.01:1; 
            bar(xval_plot10,combined_counts6./sum(combined_counts6),'FaceColor',colors(i,:),'EdgeColor','none','FaceAlpha',1); 
            hold on; 
            plot(xval_plot10, combined_counts6./sum(combined_counts6),'Color', 'k', 'LineWidth', 0.3);  % 'k' for black line 
            axis([0 inf 0 max(combined_counts6./sum(combined_counts6))]); 
            line([median(combined_solidity) median(combined_solidity)],[0 max(combined_counts6./sum(combined_counts6))],'Color','r','LineWidth',1.5);
            set(gca, 'FontSize', 11);
    
    end 

end

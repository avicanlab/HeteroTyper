%% 11.11.2025 - HeteroTyper Pipeline for Dark Plates 

function data = plot_doublingTime(data)

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
    
    %% ---------------- MAIN LOOP ----------------
    figure('Name','Doubling Time Plots');
    size_threshold = 200;
    ecc_threshold = 0.70;
    combined_doublingT = [];
    
    for iGroup = 1:length(ix)
        ix_tmp = ix{iGroup};
        combined_doublingT_group = [];
    
        for j = 1:length(ix_tmp)
            plateIdx = ix_tmp(j);
            colonies = data.processed{plateIdx}.colonies.new;
            props    = data.processed{plateIdx}.colonies.region_props;
    
            % Use the correct field (not region_props)
            col_flag    = props.flag_colony_ok(:);
            col_ecc     = props.Eccentricity(:);
            col_size    = colonies.timecourse_size_smoothed(end,:)';
            timecourse  = colonies.timecourse_size_smoothed; % [time x colonies]
            t           = colonies.time_info.elapsed_time_h(:) + incTime;
    
            % Filter by flag and eccentricity
            valid_colonies = find(col_flag == 1 & col_size > size_threshold & col_ecc < ecc_threshold);
            doublingTime = nan(size(col_flag));
    
            for c = valid_colonies'
                size_curve = timecourse(:, c);
                if all(isnan(size_curve)) || all(size_curve == 0)
                    continue;
                end
    
                % --- Find first point where size ≥ 100 px
                idx_start = find(size_curve >= size_threshold, 1, 'first');
                if isempty(idx_start)
                    continue;
                end
    
                % Dynamic doubling target: 2× size at threshold-crossing
                start_size = size_curve(idx_start);
                target_size = 2 * start_size;
                t_start = t(idx_start);
    
                % --- Find first time ≥ 2× start size
                idx_end = find(size_curve >= target_size & (1:length(size_curve))' > idx_start, 1, 'first');
                if isempty(idx_end)
                    continue;
                end
                t_end = t(idx_end);
    
                % Doubling time
                doublingTime(c) = t_end - t_start;
                combined_doublingT_group = [combined_doublingT_group; doublingTime(c)];
            end
    
            % Store per-plate result
            data.processed{plateIdx}.colonies.new.doublingTime = doublingTime;
        end
    
        % === Histogram per group ===
        subplot(length(ix), 1, iGroup);
        validT = combined_doublingT_group(~isnan(combined_doublingT_group));
        if ~isempty(validT)
            histogram(validT, 'BinWidth', 0.5, 'FaceColor', colors(iGroup,:), ...
                'EdgeColor', 'none', 'FaceAlpha', 1);
            hold on;
            medVal = median(validT);
            line([medVal medVal], [0 max(histcounts(validT, 'BinWidth', 0.5))], ...
                'Color', 'r', 'LineWidth', 1.5);
            title(sprintf('%s (n=%d colonies) | Median = %.2f h', ...
                labels_type{iGroup}, length(validT), medVal));
            xlabel('Doubling Time (h)');
            ylabel('Frequency');
            set(gca, 'FontSize', 11);
        else
            title(sprintf('%s - No valid colonies', labels_type{iGroup}));
        end
    
        % combined_doublingT = [combined_doublingT; combined_doublingT_group];
    
        % %% ---------------- COMBINED PDF ----------------
        % subplot(length(ix), 1, length(ix));
        % valid_combined = combined_doublingT(~isnan(combined_doublingT));
        % 
        % if ~isempty(valid_combined)
        %     [counts, edges] = histcounts(valid_combined, 'BinWidth', 0.5);
        %     bin_centers = edges(1:end-1) + diff(edges)/2;
        %     relative_freq = counts / sum(counts); % Probability density
        %     bar(bin_centers, relative_freq, 'FaceColor', [0.5 0.5 0.5], ...
        %         'EdgeColor', 'none', 'FaceAlpha', 0.7);
        %     hold on;
        %     plot(bin_centers, smooth(relative_freq, 5), 'k-', 'LineWidth', 1.5);
        %     medVal = median(valid_combined);
        %     line([medVal medVal], [0 max(relative_freq)], 'Color', 'r', 'LineWidth', 1.5);
        %     title(sprintf('Combined PDF | n=%d | Median = %.2f h', ...
        %         length(valid_combined), medVal));
        %     xlabel('Doubling Time (h)');
        %     ylabel('Probability Density');
        %     set(gca, 'FontSize', 11);
        % end
    end
    
    
    
    %% ---------------- FINAL SUMMARY ----------------
    assignin('base','combined_doublingT',combined_doublingT);
    fprintf('\n✅ Doubling Time Computation Finished\n');
    fprintf('Total valid colonies: %d | Global Median DT = %.2f h\n', ...
        length(combined_doublingT), nanmedian(combined_doublingT));

end

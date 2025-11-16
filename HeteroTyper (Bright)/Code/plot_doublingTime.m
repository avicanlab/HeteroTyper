%% 11.11.2025 - HeteroTyper Pipeline for Bright Plates 

function data = plot_doublingTime(data)

    nr_plates = length(data.processed);
    
    % Define maximum lag time
    max_lag = 52;

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
    
    %% ---------------- DEFINE GROUPS ----------------
    ix_3H = find(strcmpi(data.metadata.original.Time,'3H'));
    ix_7H = find(strcmpi(data.metadata.original.Time,'7H'));
    ix_24H = find(strcmpi(data.metadata.original.Time,'24H'));
    ix_48H = find(strcmpi(data.metadata.original.Time,'48H'));
    
    labels_type = {'3H','7H','24H','48H'};
    
    ix{1} = intersect(ix_3H, ix_growth);
    ix{2} = intersect(ix_7H, ix_growth);
    ix{3} = intersect(ix_24H, ix_growth);
    ix{4} = intersect(ix_48H, ix_growth);
    
    colors = [0.37 0.21 0.65;0.12 0.69 0.70;0.87 0.71 0;0.90 0.52 0.1];

    
    %% ---------------- MAIN LOOP ----------------
    figure('Name','Doubling Time Plots');
    size_threshold = 100;
    ecc_threshold = 0.70;
    combined_doublingT = [];
    
    for iGroup = 1:length(ix)
        ix_tmp = ix{iGroup};
        combined_doublingT_group = [];
    
        for j = 1:length(ix_tmp)
            plateIdx = ix_tmp(j);
            colonies = data.processed{plateIdx}.colonies;
    
            col_flag = colonies.region_props.flag_colony_ok(:);
            col_ecc  = colonies.region_props.Eccentricity(:);
            col_size = colonies.new.timecourse_size_smoothed(end,:)';
            timecourse = colonies.new.timecourse_size_smoothed; % [time x colonies]
            t = data.processed{plateIdx}.colonies.new.time_info.elapsed_time_h(:) + incTime;
    
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

end

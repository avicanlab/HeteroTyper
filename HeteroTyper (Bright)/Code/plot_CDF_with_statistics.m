%% 11.11.2025 - HeteroTyper Pipeline for Bright Plates

%% --- RUNNING THE FUNCTION ---
% 1. (initially) Run "plot_combined_samples.m" to generate lag time tables into the MATLAB workspace
% 2. Run ->   lag_times = {combined_lag_time_3H, combined_lag_time_7H, combined_lag_time_24H, combined_lag_time_48H};
% 3. Run ->   group_labels = {'3H','7H','24H','48H'};
% 4. Run ->   colors = [0.161, 0.741, 0.584; 0.67, 0.74, 0.098; 0.941, 0.153, 0.42; 0.91, 0.627, 0.02];
% 5. Run ->   plot_CDF_with_statistics(lag_times, group_labels, colors);

function results = plot_CDF_with_statistics(lag_times, group_labels, colors)

    figure('Name','Lag Time Survival','Color','w'); hold on;
    results = struct();
    
    nGroups = length(lag_times);
    survData = cell(nGroups,1);

    for g = 1:nGroups
        lt = lag_times{g};
        if isempty(lt), continue; end

        % Truncate lag times
        lt = fix(lt * 10) / 10;

        % Sort and build empirical CDF
        lt = sort(lt(:));
        [xu, ~, ~] = unique(lt);
        n = length(lt);
        F = arrayfun(@(v) sum(lt <= v)/n, xu);
        S = 1 - F;  % survival function

        % Compute slope (finite difference)
        dS = diff(S) ./ diff(xu);   % delta survival / delta time
        maxSlope = max(abs(dS));
        meanSlope = mean(abs(dS));

        % Plot
        plot(xu, S, '-o', ...
            'LineWidth', 2, ...
            'MarkerSize', 2, ...
            'MarkerFaceColor', 'none', ...
            'Color', colors(g,:), ...
            'DisplayName', group_labels{g});

        % Save results
        results(g).group = group_labels{g};
        results(g).x = xu;
        results(g).survival = S;
        results(g).median = median(lt);
        results(g).q25 = prctile(lt,25);
        results(g).q75 = prctile(lt,75);
        results(g).n = n;
        results(g).slope.max = maxSlope;
        results(g).slope.mean = meanSlope;
        results(g).min_lag = min(lag_times{g});
        results(g).median_lag = median(lag_times{g});
        results(g).max_lag = max(lag_times{g});
        survData{g} = lt;
    end

    xlabel('Lag time (h)');
    ylabel('Survival (1 - CDF)');
    legend('show','Location','northeast');
    grid on; 
    box on;



    figure('Name','Lag Time Survival - Log Scale','Color','w'); hold on;
    results = struct();
    
    nGroups = length(lag_times);
    survData = cell(nGroups,1);

    for g = 1:nGroups
        lt = lag_times{g};
        if isempty(lt), continue; end

        % Truncate lag times
        lt = fix(lt * 10) / 10;

        % Sort and build empirical CDF
        lt = sort(lt(:));
        [xu, ~, ~] = unique(lt);
        n = length(lt);
        F = arrayfun(@(v) sum(lt <= v)/n, xu);
        S = 1 - F;  % survival function

        % Compute slope (finite difference)
        dS = diff(S) ./ diff(xu);   % delta survival / delta time
        maxSlope = max(abs(dS));
        meanSlope = mean(abs(dS));

        % Plot
        plot(xu, S, '-o', ...
            'LineWidth', 2, ...
            'MarkerSize', 2, ...
            'MarkerFaceColor', 'none', ...
            'Color', colors(g,:), ...
            'DisplayName', group_labels{g});

        % Save results
        results(g).group = group_labels{g};
        results(g).x = xu;
        results(g).survival = S;
        results(g).median = median(lt);
        results(g).q25 = prctile(lt,25);
        results(g).q75 = prctile(lt,75);
        results(g).n = n;
        results(g).slope.max = maxSlope;
        results(g).slope.mean = meanSlope;
        results(g).min_lag = min(lag_times{g});
        results(g).median_lag = median(lag_times{g});
        results(g).max_lag = max(lag_times{g});
        survData{g} = lt;
    end

    xlabel('Lag time (h)');
    ylabel('Survival (1 - CDF)');
    legend('show','Location','northeast');
    set(gca,'YScale','log');
    grid on; 
    box on;


    %% --- Pairwise statistical comparisons ---
    fprintf('\n=== Pairwise Statistics ===\n');
    for i = 1:nGroups
        for j = i+1:nGroups
            % Kolmogorov-Smirnov test
            [~, pKS] = kstest2(survData{i}, survData{j});
            
            % Log-rank test (if Bioinformatics toolbox available)
            censor1 = zeros(size(survData{i}));
            censor2 = zeros(size(survData{j}));
            try
                [~, pLR] = logrank(survData{i}, survData{j}, censor1, censor2);
            catch
                pLR = NaN;
            end

            fprintf('%s vs %s: KS p=%.4f, Log-rank p=%.4f\n', ...
                group_labels{i}, group_labels{j}, pKS, pLR);
        end
    end
    drawnow;

    %% --- Display slope summary ---
    fprintf('\n=== Slope Summary ===\n');
    for g = 1:nGroups
        fprintf('%s: min=%.2f h, median=%.2f h, max=%.2f h, max slope=%.4f, mean slope=%.4f\n', ...
            results(g).group, results(g).min_lag, results(g).median_lag, results(g).max_lag, ...
            results(g).slope.max, results(g).slope.mean);
    end

end

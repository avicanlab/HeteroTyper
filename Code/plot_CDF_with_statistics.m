%% HeteroTyper Pipeline for Bright Plates
% Plots empirical survival functions (1 - CDF) for lag time:
%   Figure 1 — linear y-axis
%   Figure 2 — log y-axis  (same S values, different display scale)
%   Figure 3 — standalone mean-slope bar chart
%
% Slope = dS/dt, computed once from the raw survival values.
% The log y-axis is a display-only transformation and does not
% affect slope computation.
%
% Requires preprocess_pipeline_data(data) to have been run first.
%
%% USAGE:
%   results = plot_CDF_with_statistics(ht);

function results = plot_CDF_with_statistics(ht)
    HT_FLOG = -1;  % log file handle (opened after params are loaded)


    % ------------------------------------------------------------------
    %  Open log file — mirrors everything printed to the terminal
    % ------------------------------------------------------------------
    timestamp_log = datestr(now, 'yyyymmdd_HHMMSS');
    log_path = fullfile(ht.params.out_dir, ...
                        sprintf('plot_CDF_with_statistics_%s.txt', timestamp_log));
    HT_FLOG = fopen(log_path, 'w');
    if HT_FLOG == -1
        warning('Could not open log file: %s', log_path);
        HT_FLOG = -1;
    end
    fprintf('Log file: %s\n', log_path);

    % ------------------------------------------------------------------
    %  Validate input
    % ------------------------------------------------------------------
    if ~isstruct(ht) || ~isfield(ht, 'groups')
        error(['Input must be the ht struct produced by preprocess_pipeline_data.\n' ...
               'Run:  preprocess_pipeline_data(data)  first.']);
    end

    % ------------------------------------------------------------------
    %  Ask user which sample groups to include
    % ------------------------------------------------------------------
    sel_idx = ask_group_selection(HT_FLOG, ht);

    nGroups      = length(sel_idx);
    lag_times    = cell(nGroups, 1);
    group_labels = cell(nGroups, 1);
    colors       = zeros(nGroups, 3);

    for k = 1:nGroups
        g               = sel_idx(k);
        lag_times{k}    = ht.groups(g).lag_time;
        group_labels{k} = ht.labels{g};
        colors(k, :)    = ht.colors(g, :);
    end

    ht_fprintf(HT_FLOG, '\n--- CDF plot parameters ---\n');
    ht_fprintf(HT_FLOG, '  Groups selected : %s\n', strjoin(group_labels, ', '));
    ht_fprintf(HT_FLOG, '  Output directory: %s\n', ht.params.out_dir);
    ht_fprintf(HT_FLOG, '---------------------------\n\n');

    % ------------------------------------------------------------------
    %  Compute survival curves and slopes once
    %  (log y-axis is display-only — S values and slopes are the same)
    % ------------------------------------------------------------------
    survData   = cell(nGroups, 1);
    xu_all     = cell(nGroups, 1);
    S_all      = cell(nGroups, 1);
    mean_slopes = nan(nGroups, 1);
    results    = struct();

    for g = 1:nGroups
        lt = lag_times{g};
        if isempty(lt), continue; end

        [xu, S, n, dS] = build_survival(lt);

        xu_all{g}   = xu;
        S_all{g}    = S;
        survData{g} = lt;

        results(g).group      = group_labels{g};
        results(g).x          = xu;
        results(g).survival   = S;
        % NaN-safe: lag_time may contain NaN entries for right-censored
        % colonies (set by preprocess). prctile already ignores NaN;
        % median does not, so use 'omitnan'.
        results(g).median     = median(lt, 'omitnan');
        results(g).q25        = prctile(lt, 25);
        results(g).q75        = prctile(lt, 75);
        results(g).n          = n;
        results(g).min_lag    = min(lt, [], 'omitnan');
        results(g).median_lag = median(lt, 'omitnan');
        results(g).max_lag    = max(lt, [], 'omitnan');
        results(g).slope.max  = max(abs(dS));
        results(g).slope.mean = mean(abs(dS));

        mean_slopes(g) = results(g).slope.mean;
    end

    % ------------------------------------------------------------------
    %  Figure 1 — Linear y-axis
    % ------------------------------------------------------------------
    hfig_lin = figure('Name','Lag Time Survival - Linear', ...
                      'Color','w', 'Units','inches', 'Position',[1 5 6 4.5]);
    ax_lin = axes('Position',[0.12 0.13 0.82 0.80]);
    hold(ax_lin, 'on');
    for g = 1:nGroups
        if isempty(xu_all{g}), continue; end
        plot(ax_lin, xu_all{g}, S_all{g}, '-o', ...
             'LineWidth',2, 'MarkerSize',2, ...
             'MarkerFaceColor','none', ...
             'Color',colors(g,:), ...
             'DisplayName',group_labels{g});
    end
    xlabel(ax_lin, 'Lag time (h)',       'FontWeight','bold');
    ylabel(ax_lin, 'Survival (1 - CDF)', 'FontWeight','bold');
    legend(ax_lin, 'show', 'Location','northeast');
    set(ax_lin, 'FontSize',11, 'Box','on');
    grid(ax_lin, 'on');

    % ------------------------------------------------------------------
    %  Figure 2 — Log y-axis  (same data, display scale only)
    % ------------------------------------------------------------------
    hfig_log = figure('Name','Lag Time Survival - Log Scale', ...
                      'Color','w', 'Units','inches', 'Position',[8 5 6 4.5]);
    ax_log = axes('Position',[0.12 0.13 0.82 0.80]);
    hold(ax_log, 'on');
    for g = 1:nGroups
        if isempty(xu_all{g}), continue; end
        plot(ax_log, xu_all{g}, S_all{g}, '-o', ...
             'LineWidth',2, 'MarkerSize',2, ...
             'MarkerFaceColor','none', ...
             'Color',colors(g,:), ...
             'DisplayName',group_labels{g});
    end
    xlabel(ax_log, 'Lag time (h)',       'FontWeight','bold');
    ylabel(ax_log, 'Survival (1 - CDF)', 'FontWeight','bold');
    legend(ax_log, 'show', 'Location','northeast');
    set(ax_log, 'YScale','log', 'FontSize',11, 'Box','on');
    grid(ax_log, 'on');

    % ------------------------------------------------------------------
    %  Figure 3 — Standalone mean-slope bar chart
    % ------------------------------------------------------------------
    hfig_bar = figure('Name','Mean Slope', ...
                      'Color','w', 'Units','inches', 'Position',[1 0.5 4 4]);
    ax_bar = axes('Position',[0.18 0.12 0.72 0.80]);
    draw_slope_bars(ax_bar, mean_slopes, group_labels, colors);

    drawnow;

    % ------------------------------------------------------------------
    %  Terminal — Pairwise statistics
    % ------------------------------------------------------------------
    ht_fprintf(HT_FLOG, '\n=== Pairwise Statistics ===\n');
    for i = 1:nGroups
        for j = i+1:nGroups
            if isempty(survData{i}) || isempty(survData{j}), continue; end
            [~, pKS] = kstest2(survData{i}, survData{j});
            try
                [~, pLR] = logrank(survData{i}, survData{j}, ...
                                   zeros(size(survData{i})), ...
                                   zeros(size(survData{j})));
            catch
                pLR = NaN;
            end
            ht_fprintf(HT_FLOG, '  %s vs %s:  KS p = %s,  Log-rank p = %s\n', ...
                group_labels{i}, group_labels{j}, fmt_p(pKS, '%.4f'), fmt_p(pLR, '%.4f'));
        end
    end

    % ------------------------------------------------------------------
    %  Terminal — Slope & range summary
    %  One summary: slope = dS/dt, identical for both display scales
    % ------------------------------------------------------------------
    ht_fprintf(HT_FLOG, '\n=== Slope & Range Summary  (dS/dt) ===\n');
    ht_fprintf(HT_FLOG, '  NOTE: the log y-axis is a display-only transformation.\n');
    ht_fprintf(HT_FLOG, '  Survival values S are unchanged, so slopes are the same\n');
    ht_fprintf(HT_FLOG, '  whether read from the linear or log figure.\n\n');
    ht_fprintf(HT_FLOG, '  %-6s  %8s  %10s  %8s  %9s  %10s\n', ...
            'Group', 'min (h)', 'median (h)', 'max (h)', 'max slope', 'mean slope');
    ht_fprintf(HT_FLOG, '  %s\n', repmat('-', 1, 60));
    for g = 1:nGroups
        if ~isfield(results(g),'group') || isempty(results(g).group), continue; end
        ht_fprintf(HT_FLOG, '  %-6s  %8.2f  %10.2f  %8.2f  %9.4f  %10.4f\n', ...
            results(g).group, ...
            results(g).min_lag, results(g).median_lag, results(g).max_lag, ...
            results(g).slope.max, results(g).slope.mean);
    end

    % ------------------------------------------------------------------
    %  Save all three figures (.fig, .pdf, .png)
    % ------------------------------------------------------------------
    out_dir   = ht.params.out_dir;
    timestamp = datestr(now, 'yyyymmdd_HHMMSS');
    grp_str   = strjoin(group_labels, '-');

    save_figure(HT_FLOG, hfig_lin, out_dir, sprintf('CDF_linear_%s_%s',    grp_str, timestamp));
    save_figure(HT_FLOG, hfig_log, out_dir, sprintf('CDF_log_%s_%s',       grp_str, timestamp));
    save_figure(HT_FLOG, hfig_bar, out_dir, sprintf('CDF_slope_%s_%s',     grp_str, timestamp));

    % Close log file
    if HT_FLOG ~= -1
        fclose(HT_FLOG);
        fprintf('Log saved: %s\n', log_path);
        HT_FLOG = -1;
    end

end  % end main function


%% ================================================================
%  LOCAL FUNCTION: build_survival
%  Builds empirical survival S = 1 - F from a lag time vector.
%  Returns unique x values, S, n, and finite-difference slope dS.
%  Slope = dS/dt — rate of change of survival fraction per hour.
% ================================================================
function [xu, S, n, dS] = build_survival(lt)
    lt         = lt(~isnan(lt));        % drop NaN (right-censored / size-rejected)
    lt         = fix(lt * 10) / 10;
    lt         = sort(lt(:));
    [xu, ~, ~] = unique(lt);
    n          = length(lt);
    F          = arrayfun(@(v) sum(lt <= v) / n, xu);
    S          = 1 - F;
    if length(xu) > 1
        dS = diff(S) ./ diff(xu);
    else
        dS = 0;
    end

end


%% ================================================================
%  LOCAL FUNCTION: draw_slope_bars
%  Standalone bar chart: group colours, no edge, clean axes, no box.
% ================================================================
function draw_slope_bars(ax, mean_slopes, group_labels, colors)

    nGroups = length(mean_slopes);
    valid   = ~isnan(mean_slopes);

    hold(ax, 'on');
    for g = 1:nGroups
        if ~valid(g), continue; end
        bar(ax, g, mean_slopes(g), 0.6, ...
            'FaceColor', colors(g,:), ...
            'EdgeColor', 'none', ...
            'FaceAlpha', 1);
    end

    y_max = max(mean_slopes(valid));
    if isempty(y_max) || y_max == 0, y_max = 0.1; end

    % Tight axis: 2% breathing gap, ~5 clean ticks.
    [y_ceil, tick_step] = nice_axis(y_max);

    ylim(ax, [0, y_ceil]);
    yticks(ax, 0 : tick_step : y_ceil);
    xlim(ax, [0.4, nGroups + 0.6]);
    set(ax, ...
        'XTick',              1:nGroups, ...
        'XTickLabel',         group_labels, ...
        'FontSize',           11, ...
        'Box',                'off', ...
        'TickDir',            'out', ...
        'XTickLabelRotation', 0);
    ylabel(ax, 'Mean Slope', 'FontWeight','bold');
    grid(ax, 'off');

end


%% ================================================================
%  LOCAL FUNCTION: nice_axis
%  Returns a tight (y_ceil, tick_step) pair for the slope bar axis.
%  Adds a 2% breathing gap above y_max, then finds the smallest
%  clean step that fits ~5 ticks within that ceiling.
% ================================================================
function [y_ceil, tick_step] = nice_axis(y_max)
    if y_max <= 0
        y_ceil    = 0.1;
        tick_step = 0.02;
        return;
    end
    raw_ceil   = y_max * 1.02;          % 2% gap above the tallest bar
    raw_step   = raw_ceil / 5;          % target ~5 ticks
    mag        = 10^floor(log10(raw_step));
    candidates = [1 2 2.5 5 10] * mag;
    for step = candidates
        y_ceil = ceil(raw_ceil / step) * step;
        if round(y_ceil / step) <= 6    % accept up to 6 ticks
            tick_step = step;
            return;
        end
    end
    % Fallback
    tick_step = candidates(end);
    y_ceil    = ceil(raw_ceil / tick_step) * tick_step;
end


%% ================================================================
%  LOCAL FUNCTION: ask_group_selection
%  Presents available group labels and returns selected indices.
% ================================================================
function sel_idx = ask_group_selection(HT_FLOG, ht)

    all_labels = ht.labels;
    n_avail    = length(all_labels);

    ht_fprintf(HT_FLOG, '\nAvailable sample groups:\n');
    ht_fprintf(HT_FLOG, '  +-----+--------+-------------+\n');
    ht_fprintf(HT_FLOG, '  | idx | label  | n colonies  |\n');
    ht_fprintf(HT_FLOG, '  +-----+--------+-------------+\n');
    for k = 1:n_avail
        ht_fprintf(HT_FLOG, '  | %3d | %-6s | %11d |\n', ...
                k, all_labels{k}, ht.groups(k).n_colonies);
    end
    ht_fprintf(HT_FLOG, '  +-----+--------+-------------+\n');
    ht_fprintf(HT_FLOG, '  Enter group labels separated by spaces or commas (e.g.  3h 7h)\n');
    ht_fprintf(HT_FLOG, '  or type  ALL  to include all groups.\n');

    while true
        raw    = strtrim(input('  Your selection: ', 's'));

        if strcmpi(raw, 'all')
            sel_idx = 1:n_avail;
            ht_fprintf(HT_FLOG, '  -> All %d groups selected.\n', n_avail);
            return;
        end

        tokens  = regexp(strtrim(raw), '[,\s]+', 'split');
        tokens  = tokens(~cellfun(@isempty, tokens));

        sel_idx = zeros(1, length(tokens));
        ok      = true;
        for t = 1:length(tokens)
            match = find(strcmpi(all_labels, tokens{t}));
            if isempty(match)
                ht_fprintf(HT_FLOG, '  WARNING  "%s" is not a recognised group label. Try again.\n', tokens{t});
                ok = false;
                break;
            end
            sel_idx(t) = match(1);
        end

        if ~ok, continue; end

        [~, ui] = unique(sel_idx, 'stable');
        sel_idx = sel_idx(ui);
        ht_fprintf(HT_FLOG, '  -> Selected groups: %s\n', strjoin(all_labels(sel_idx), ', '));
        return;
    end

end


%% ================================================================
%  LOCAL FUNCTION: save_figure
%  Saves a figure as .fig, .pdf (vector), and .png (300 dpi).
% ================================================================
function save_figure(HT_FLOG, hfig, out_dir, base_name)
    fig_path = fullfile(out_dir, base_name);
    savefig(hfig,  [fig_path '.fig']);
    exportgraphics(hfig, [fig_path '.pdf'], 'ContentType','vector');
    exportgraphics(hfig, [fig_path '.png'], 'Resolution',300);
    ht_fprintf(HT_FLOG, '  Saved: %s  (.fig / .pdf / .png)\n', base_name);
end


%% ================================================================
%  LOCAL FUNCTION: ht_fprintf
%  Writes to terminal and to the log file handle passed as HT_FLOG.
% ================================================================
function ht_fprintf(HT_FLOG, varargin)
    fprintf(varargin{:});
    if ~isempty(HT_FLOG) && HT_FLOG ~= -1
        fprintf(HT_FLOG, varargin{:});
    end
end


%% ================================================================
%  LOCAL FUNCTION: fmt_p
%  Text display for a p-value. kstest2 / logrank can legitimately
%  underflow to exact 0.0 in double precision for large n / a strong
%  distributional difference — the true value is just too small to
%  represent, not actually zero — shown as "<1e-300" instead.
% ================================================================
function s = fmt_p(p, fmt)
    if nargin < 2, fmt = '%.4g'; end
    if isnan(p)
        s = 'NaN';
    elseif p == 0
        s = '<1e-300';
    else
        s = sprintf(fmt, p);
    end
end
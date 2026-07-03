%% HeteroTyper Pipeline for Bright Plates
% Plots combined growth curves per sample group.
% Each line represents one individual colony, coloured distinctly.
% Requires preprocess_pipeline_data(data) to have been run first.
%
%% USAGE:
%   plot_combined_samples_growth_curves(ht);

function plot_combined_samples_growth_curves(ht)
    HT_FLOG = -1;  % log file handle (opened after params are loaded)

    check_ht(ht);

    p        = ht.params;
    n_groups = length(ht.groups);
    labels   = ht.labels;

    % ------------------------------------------------------------------
    %  Interactive group selection — mirrors plot_combined_samples
    % ------------------------------------------------------------------
    sel_groups = ask_group_selection(ht);
    n_sel      = length(sel_groups);

    % ------------------------------------------------------------------
    %  Figure layout — manual position grid so labels never clip
    % ------------------------------------------------------------------
    panel_w  = 5.0;   panel_h  = 2.0;
    gap_y    = 0.55;
    margin_l = 0.85;  margin_r = 0.20;
    margin_b = 0.70;  margin_t = 0.20;

    fig_w = margin_l + panel_w + margin_r;
    fig_h = margin_b + n_sel*panel_h + (n_sel-1)*gap_y + margin_t;


    % ------------------------------------------------------------------
    %  Open log file — mirrors everything printed to the terminal
    % ------------------------------------------------------------------
    timestamp_log = datestr(now, 'yyyymmdd_HHMMSS');
    log_path = fullfile(p.out_dir, ...
                        sprintf('plot_combined_samples_growth_curves_%s.txt', timestamp_log));
    HT_FLOG = fopen(log_path, 'w');
    if HT_FLOG == -1
        warning('Could not open log file: %s', log_path);
        HT_FLOG = -1;
    end
    fprintf('Log file: %s\n', log_path);

    hfig = figure('Name','Combined Growth Curves', ...
                  'Color','w', ...
                  'Units','inches', ...
                  'Position',[1 1 fig_w fig_h], ...
                  'PaperUnits','inches', ...
                  'PaperSize',[fig_w fig_h], ...
                  'PaperPosition',[0 0 fig_w fig_h]);

    global_y_max = 0;
    global_t_max = 0;
    for gi = 1:n_sel
        grp = ht.groups(sel_groups(gi));
        if isempty(grp.size_timecourse), continue; end
        global_y_max = max(global_y_max, max(grp.size_timecourse(end,:), [], 'omitnan'));
        global_t_max = max(global_t_max, max(grp.time_h));
    end
    y_ceil = global_y_max + 100;

    % ------------------------------------------------------------------
    %  Plot loop — one row per selected group
    % ------------------------------------------------------------------
    for gi = 1:n_sel
        g   = sel_groups(gi);
        grp = ht.groups(g);

        % Axes position (row 1 = top)
        left   = margin_l / fig_w;
        bottom = (margin_b + (n_sel-gi)*(panel_h+gap_y)) / fig_h;
        ax = axes('Position',[left, bottom, panel_w/fig_w, panel_h/fig_h]); %#ok<LAXES>

        if isempty(grp.size_timecourse)
            ht_fprintf(HT_FLOG, 'Group %s: no valid timecourse data, skipping.\n', labels{g});
            text(0.5, 0.5, 'No data', 'Units','normalized', ...
                 'HorizontalAlignment','center', 'FontSize',11);
            axis off;
            continue;
        end

        t_vec   = grp.time_h;
        size_tc = grp.size_timecourse;
        n_col   = size(size_tc, 2);

        ht_fprintf(HT_FLOG, 'Group %s — %d colonies, max final size = %.2f px\n', ...
                labels{g}, n_col, max(size_tc(end,:), [], 'omitnan'));

        % Per-colony colours: cycle through a shuffled colormap so each
        % line is visually distinct (mirrors the original jet-shuffle mask)
        col_colors = colony_colors(n_col);

        hold(ax, 'on');
        for c = 1:n_col
            plot(ax, t_vec, size_tc(:,c), '-', ...
                 'Color',     col_colors(c,:), ...
                 'LineWidth', 0.6);
        end

        axis(ax, [p.incTime, global_t_max, 0, y_ceil]);
        ylabel(ax, 'Colony size (px)', 'FontSize',11, 'FontWeight','bold');
        set(ax, 'FontSize',10, 'TickDir','out', 'Box','off');

        % x-label only on bottom row; group label always shown
        if gi == n_sel
            xlabel(ax, 'Time (h)', 'FontSize',11, 'FontWeight','bold');
        end

        % Group label inside panel — top-left corner
        xl = xlim(ax);
        yl = ylim(ax);
        text(ax, xl(1) + 0.02*(xl(2)-xl(1)), yl(2)*0.96, labels{g}, ...
             'FontSize',11, 'FontWeight','bold', 'Color','k', ...
             'VerticalAlignment','top');
    end

    % ------------------------------------------------------------------
    %  Save figure (.fig, .pdf, .png)
    % ------------------------------------------------------------------
    out_dir   = p.out_dir;

    

    timestamp = datestr(now, 'yyyymmdd_HHMMSS');
    base_name = sprintf('GrowthCurves_%s', timestamp);
    fig_path  = fullfile(out_dir, base_name);

    savefig(hfig,  [fig_path '.fig']);
    exportgraphics(hfig, [fig_path '.pdf'], 'ContentType','vector');
    exportgraphics(hfig, [fig_path '.png'], 'Resolution',300);

    ht_fprintf(HT_FLOG, '\nFigure saved:\n  %s.fig\n  %s.pdf\n  %s.png\n', ...
            fig_path, fig_path, fig_path);

    % Close log file
    if HT_FLOG ~= -1
        fclose(HT_FLOG);
        fprintf('Log saved: %s\n', log_path);
        HT_FLOG = -1;
    end

end


%% ================================================================
%  LOCAL FUNCTION: ask_group_selection
%  Interactive prompt — mirrors ask_selections in plot_combined_samples
%  but only asks for groups (no feature selection needed here).
% ================================================================
function sel_groups = ask_group_selection(ht)
    all_labels = ht.labels;
    n_avail    = length(all_labels);

    fprintf('\n============================================================\n');
    fprintf('  plot_combined_samples_growth_curves — Group Selection\n');
    fprintf('============================================================\n');
    fprintf('\n  Available groups:\n');
    fprintf('  +-----+--------+-------------+---------------------+\n');
    fprintf('  | idx | label  | n colonies  | timecourse cols     |\n');
    fprintf('  +-----+--------+-------------+---------------------+\n');
    for k = 1:n_avail
        grp  = ht.groups(k);
        n_tc = size(grp.size_timecourse, 2);
        fprintf('  | %3d | %-6s | %11d | %19d |\n', ...
                k, all_labels{k}, grp.n_colonies, n_tc);
    end
    fprintf('  +-----+--------+-------------+---------------------+\n');
    fprintf('\n  Enter group labels (e.g.  3h 24h)  or  ALL\n');
    fprintf('============================================================\n');

    while true
        raw_g = strtrim(input('  Groups: ', 's'));
        if strcmpi(raw_g, 'all')
            sel_groups = 1:n_avail;
            fprintf('  -> All %d groups selected.\n\n', n_avail);
            return;
        end
        tokens = regexp(raw_g, '[,\s]+', 'split');
        tokens = tokens(~cellfun(@isempty, tokens));
        sel = zeros(1, length(tokens));
        ok  = true;
        for t = 1:length(tokens)
            match = find(strcmpi(all_labels, tokens{t}));
            if isempty(match)
                fprintf('  WARNING  "%s" not recognised. Available: %s\n', ...
                        tokens{t}, strjoin(all_labels, ', '));
                ok = false; break;
            end
            sel(t) = match(1);
        end
        if ~ok, continue; end
        [~, ui]    = unique(sel, 'stable');
        sel_groups = sel(ui);
        fprintf('  -> Groups: %s\n\n', strjoin(all_labels(sel_groups), ', '));
        return;
    end
end


%% ================================================================
%  LOCAL FUNCTION: colony_colors
%  Returns an [n x 3] RGB matrix of visually distinct per-colony
%  colours by shuffling the 'lines' colormap (matches the
%  jet-shuffle colouring used in the original segmentation masks).
% ================================================================
function C = colony_colors(n)
    if n == 0
        C = zeros(0, 3);
        return;
    end
    % Base palette: tile the 'lines' colormap to cover n colonies
    base = lines(min(n, 256));
    if n > size(base,1)
        reps = ceil(n / size(base,1));
        base = repmat(base, reps, 1);
    end
    base = base(1:n, :);

    % Shuffle with a fixed seed so colours are reproducible across runs
    % but still look random (no two adjacent colonies share a colour)
    rng(42);
    idx = randperm(n);
    C   = base(idx, :);
end


%% ================================================================
%  LOCAL FUNCTION: check_ht
% ================================================================
function check_ht(ht)
    if ~isstruct(ht) || ~isfield(ht, 'groups')
        error(['Input must be the ht struct produced by preprocess_pipeline_data.\n' ...
               'Run:  preprocess_pipeline_data(data)  first.']);
    end
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
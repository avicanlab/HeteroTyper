%% HeteroTyper Pipeline
%  plot_lagTime_vs_growthRate  —  Lag time vs maximum Gompertz growth rate.
%
%  Scatter plot of lag time (h) vs mu_max (px/h) per selected group.
%  All groups in one figure; each panel coloured by log10(local_density)
%  using a viridis palette. A LOESS smooth line (red) is overlaid per panel
%  with its degree-2 polynomial equation. Significance is reported at
%  p < 0.05 (* / ** / *** / ns).
%
%  Gompertz fitting and all statistics are pre-computed in STEP 6d of
%  preprocess_pipeline_data — this script reads and visualises only.
%
%  USAGE:
%    plot_lagTime_vs_growthRate(ht)

function plot_lagTime_vs_growthRate(ht)
    HT_FLOG = -1;

    check_ht(ht);

    p      = ht.params;
    labels = ht.labels;
    colors = ht.colors;

    % ------------------------------------------------------------------
    %  Open log file
    % ------------------------------------------------------------------
    timestamp_log = datestr(now, 'yyyymmdd_HHMMSS');
    log_path = fullfile(p.out_dir, ...
        sprintf('plot_lagTime_vs_growthRate_%s.txt', timestamp_log));
    HT_FLOG = fopen(log_path, 'w');
    if HT_FLOG == -1
        warning('Could not open log file: %s', log_path);
        HT_FLOG = -1;
    end
    fprintf('Log file: %s\n', log_path);

    % ------------------------------------------------------------------
    %  Guard: growth field must exist
    % ------------------------------------------------------------------
    if ~isfield(ht.groups(1), 'growth')
        error(['plot_lagTime_vs_growthRate: ht.groups(g).growth not found.\n' ...
               'Re-run preprocess_pipeline_data (STEP 6d) first.']);
    end

    % ------------------------------------------------------------------
    %  Group selection
    % ------------------------------------------------------------------
    sel_idx  = ask_group_selection(HT_FLOG, ht);
    n_groups = length(sel_idx);

    ht_fprintf(HT_FLOG, '\n--- Lag Time vs Growth Rate ---\n');
    ht_fprintf(HT_FLOG, '  Groups : %s\n', strjoin(labels(sel_idx), ', '));
    ht_fprintf(HT_FLOG, '  Output : %s\n', p.out_dir);
    ht_fprintf(HT_FLOG, '-------------------------------\n\n');

    % ------------------------------------------------------------------
    %  Collect pooled data for shared axis limits
    % ------------------------------------------------------------------
    all_mu  = [];
    all_lag = [];
    for k = 1:n_groups
        g  = sel_idx(k);
        ok = ht.groups(g).growth.fit_ok & isfinite(ht.groups(g).growth.mu_max);
        all_mu  = [all_mu;  ht.groups(g).growth.mu_max(ok)];  %#ok<AGROW>
        all_lag = [all_lag; ht.groups(g).lag_time(ok)];        %#ok<AGROW>
    end

    if isempty(all_mu)
        ht_fprintf(HT_FLOG, 'No successful Gompertz fits found.\n');
        if HT_FLOG ~= -1, fclose(HT_FLOG); end
        return;
    end

    y_upper = ceil(max(all_mu) / 100) * 100;
    y_step  = 100;
    x_lo    = floor(min(all_lag) / 4) * 4;
    x_hi    = ceil(max(all_lag)  / 4) * 4;

    % ------------------------------------------------------------------
    %  Figure layout — grid of panels as square as possible.
    %  Colorbar is a dedicated axes so panels are never squeezed.
    % ------------------------------------------------------------------
    n_cols = ceil(sqrt(n_groups));
    n_rows = ceil(n_groups / n_cols);

    panel_w  = 3.5;   panel_h  = 3.2;
    gap_x    = 0.65;  gap_y    = 0.75;
    margin_l = 1.00;  margin_r = 0.20;
    margin_b = 0.85;  margin_t = 0.45;
    cb_w     = 0.20;  cb_gap   = 0.45;

    grid_w = n_cols*panel_w + (n_cols-1)*gap_x;
    grid_h = n_rows*panel_h + (n_rows-1)*gap_y;
    fig_w  = margin_l + grid_w + cb_gap + cb_w + margin_r;
    fig_h  = margin_b + grid_h + margin_t;

    hfig = figure('Name','Lag Time vs Growth Rate', ...
                  'Color','w', ...
                  'Units','inches', ...
                  'Position',[1 1 fig_w fig_h], ...
                  'PaperUnits','inches', ...
                  'PaperSize',[fig_w fig_h], ...
                  'PaperPosition',[0 0 fig_w fig_h], ...
                  'Visible','off');

    fs_tick  = 12;
    fs_label = 13;
    fs_annot = 13;
    fs_stat  = 11;

    cmap = viridis_cmap(256);

    % ------------------------------------------------------------------
    %  Log table header
    % ------------------------------------------------------------------
    ht_fprintf(HT_FLOG, '%-10s  %6s  %10s  %10s  %4s  %12s\n', ...
        'Group','n_fit','Spearman_r','Spearman_p','sig','wt_Spearman');
    ht_fprintf(HT_FLOG, '%s\n', repmat('-',1,58));

    % Pre-compute global log10-density range across ALL selected groups so
    % every panel uses the same colour scale and the colorbar is accurate.
    all_log_dens = [];
    for k = 1:n_groups
        g_pre = sel_idx(k);
        ok_pre = ht.groups(g_pre).growth.fit_ok & ...
                 isfinite(ht.groups(g_pre).growth.mu_max) & ...
                 isfinite(ht.groups(g_pre).lag_time);
        d_pre = ht.groups(g_pre).growth.local_density(ok_pre);
        if ~isempty(d_pre) && any(d_pre > 0)
            all_log_dens = [all_log_dens; log10(max(d_pre, 1))]; %#ok<AGROW>
        end
    end
    if isempty(all_log_dens)
        d_min_global = 0;  d_max_global = 1;
    else
        d_min_global = min(all_log_dens);
        d_max_global = max(all_log_dens);
        if d_max_global <= d_min_global
            d_max_global = d_min_global + 1;
        end
    end

    % ------------------------------------------------------------------
    %  Plot loop — one panel per group
    % ------------------------------------------------------------------
    for k = 1:n_groups
        g   = sel_idx(k);
        grp = ht.groups(g);

        ok     = grp.growth.fit_ok & isfinite(grp.growth.mu_max) & isfinite(grp.lag_time);
        lag_g  = grp.lag_time(ok);
        mu_g   = grp.growth.mu_max(ok);
        dens_g = grp.growth.local_density(ok);
        n_fit  = sum(ok);

        % Panel position
        col_k  = mod(k-1, n_cols);
        row_k  = floor((k-1) / n_cols);
        left   = (margin_l + col_k*(panel_w+gap_x)) / fig_w;
        bottom = (margin_b + (n_rows-1-row_k)*(panel_h+gap_y)) / fig_h;
        ax = axes('Position',[left, bottom, panel_w/fig_w, panel_h/fig_h]); %#ok<LAXES>
        hold(ax, 'on');

        % ── Scatter coloured by log10(local_density) ─────────────────
        % Colour each point by log10(local_density) mapped through the
        % global range so all panels share the same viridis scale.
        if ~isempty(dens_g) && any(dens_g > 0)
            log_dens  = log10(max(dens_g, 1));
            d_range   = max(d_max_global - d_min_global, eps);
            d_norm    = (log_dens - d_min_global) / d_range;
            d_norm    = max(0, min(1, d_norm));
            c_idx     = max(1, min(256, round(d_norm*255) + 1));
            pt_colors = cmap(c_idx, :);
        else
            pt_colors = repmat(colors(g,:), n_fit, 1);
        end

        scatter(ax, lag_g, mu_g, 27, pt_colors, 'filled', 'MarkerFaceAlpha', 0.80);

        % ── LOESS smooth (red) + polynomial equation of fitted curve ────
        loess_eq = '';
        if n_fit >= 10
            [lag_sort, si] = sort(lag_g);
            mu_sort        = mu_g(si);
            mu_loess       = loess_smooth(lag_sort, mu_sort, 0.75);
            mu_loess       = max(0, min(y_upper, mu_loess));
            plot(ax, lag_sort, mu_loess, '-', 'Color','r', 'LineWidth',1.5);

            % Fit degree-2 polynomial to the LOESS curve to get a concise equation
            % polyfit works on the smoothed values (not raw scatter) so it
            % captures the non-linear trend rather than noise.
            p_coef = polyfit(lag_sort, mu_loess, 2);
            a2 = p_coef(1);  a1 = p_coef(2);  a0 = p_coef(3);
            % Format: y = a2 x^2 + a1 x + a0  (sign-aware, compact)
            loess_eq = sprintf('y = %.3g x^2 %+.3g x %+.3g', a2, a1, a0);
        end

        % ── Load stats for annotation ────────────────────────────────
        stats_g = grp.growth;

        % ── Axis limits and ticks ─────────────────────────────────────
        xlim(ax, [x_lo, x_hi]);
        xticks(ax, x_lo:8:x_hi);
        ylim(ax, [0, y_upper]);
        yticks(ax, 0:y_step:y_upper);

        % ── Grid (major + one midpoint minor tick) ────────────────────
        set(ax, 'XGrid','on','YGrid','on', ...
                'XMinorGrid','on','YMinorGrid','on', ...
                'GridLineStyle','--', ...
                'GridColor',[0.72 0.72 0.72],'GridAlpha',1.0, ...
                'MinorGridLineStyle','--', ...
                'MinorGridColor',[0.87 0.87 0.87],'MinorGridAlpha',1.0);
        maj_x = xticks(ax);
        if numel(maj_x) >= 2
            ax.XAxis.MinorTickValues = (maj_x(1:end-1) + maj_x(2:end)) / 2;
        end
        maj_y = yticks(ax);
        if numel(maj_y) >= 2
            ax.YAxis.MinorTickValues = (maj_y(1:end-1) + maj_y(2:end)) / 2;
        end

        % ── Statistics annotation (top-right) ────────────────────────
        % Significance marker based on p < 0.05
        if ~isnan(stats_g.spearman_p) && stats_g.spearman_p < 0.001
            sig_str = '***';
        elseif ~isnan(stats_g.spearman_p) && stats_g.spearman_p < 0.01
            sig_str = '**';
        elseif ~isnan(stats_g.spearman_p) && stats_g.spearman_p < 0.05
            sig_str = '*';
        else
            sig_str = 'ns';
        end

        ann_lines = { ...
            sprintf('n = %d',                      n_fit), ...
            sprintf('\\rho = %.3f  p = %s (%s)', ...
                stats_g.spearman_rho, fmt_p(stats_g.spearman_p, '%.3g'), sig_str), ...
            sprintf('R^2 = %.3f',                  stats_g.spearman_rho^2), ...
            sprintf('\\rho_w = %.3f',              stats_g.wt_spearman), ...
            loess_eq ...   % LOESS polynomial equation (red) — empty if n<10
        };

        x_stat     = x_hi - 0.03*(x_hi - x_lo);
        y_stat     = y_upper * 0.98;
        y_step_ann = y_upper * 0.09;
        n_black    = 4;   % first 4 lines black (stats); last line red (LOESS eq)
        for ai = 1:length(ann_lines)
            if isempty(ann_lines{ai}), continue; end
            txt_color = [0 0 0];
            if ai > n_black
                txt_color = [0.75 0.10 0.10];   % red — matches LOESS line
            end
            text(ax, x_stat, y_stat - (ai-1)*y_step_ann, ann_lines{ai}, ...
                 'FontSize', fs_stat, 'Color', txt_color, ...
                 'HorizontalAlignment','right', 'VerticalAlignment','top', ...
                 'Interpreter','tex');
        end

        % ── Group label (top-left, group colour) ─────────────────────
        text(ax, x_lo + 0.03*(x_hi-x_lo), y_upper*0.98, labels{g}, ...
             'FontSize',fs_annot, 'FontWeight','bold', 'Color',colors(g,:), ...
             'HorizontalAlignment','left', 'VerticalAlignment','top', ...
             'Interpreter','none');

        % ── Axis labels (bottom row / first col only) ─────────────────
        if row_k == n_rows-1 || k == n_groups
            xlabel(ax, 'Lag time (h)', 'FontSize',fs_label, 'FontWeight','bold');
        end
        if col_k == 0
            ylabel(ax, 'Maximum growth rate (px/h)', ...
                   'FontSize',fs_label, 'FontWeight','bold');
        else
            set(ax, 'YTickLabel', {});
        end

        % ── Mirror ticks (top + right frame) ─────────────────────────
        set(ax, 'FontSize',fs_tick, 'TickDir','in', 'Box','on', ...
                'XMinorTick','off', 'YMinorTick','off');
        ax2 = axes('Position', ax.Position, ...
                   'XAxisLocation','top', 'YAxisLocation','right', ...
                   'Color','none', 'XColor',ax.XColor, 'YColor',ax.YColor, ...
                   'FontSize',fs_tick, 'TickDir','in', 'Box','off', ...
                   'XTick',ax.XTick, 'YTick',ax.YTick, ...
                   'XTickLabel',{}, 'YTickLabel',{}, ...
                   'XLim',ax.XLim, 'YLim',ax.YLim); %#ok<LAXES>
        uistack(ax2, 'bottom');
        axes(ax); %#ok<LAXES>

        % ── Log row ───────────────────────────────────────────────────
        ht_fprintf(HT_FLOG, '%-10s  %6d  %10.4f  %10s  %4s  %12.4f\n', ...
            labels{g}, n_fit, stats_g.spearman_rho, fmt_p(stats_g.spearman_p, '%.4g'), ...
            sig_str, stats_g.wt_spearman);
    end

    % ------------------------------------------------------------------
    %  Standalone colorbar — right of grid, vertically centred
    % ------------------------------------------------------------------
    cb_left   = (margin_l + grid_w + cb_gap) / fig_w;
    cb_bottom = margin_b / fig_h;
    cb_height = grid_h / fig_h;
    ax_cb = axes('Position',[cb_left, cb_bottom, cb_w/fig_w, cb_height], ...
                 'Visible','off'); %#ok<LAXES>
    colormap(ax_cb, cmap);
    clim(ax_cb, [d_min_global, d_max_global]);
    cb = colorbar(ax_cb, 'eastoutside');
    cb.Position    = [cb_left + 0.005, cb_bottom, 0.018, cb_height];
    cb.Label.String   = 'Colonies (log_{10})';
    cb.Label.FontSize = fs_label;
    cb.FontSize       = fs_tick;
    cb.TickDirection  = 'out';

    % ------------------------------------------------------------------
    %  Save (figure hidden until all files written)
    % ------------------------------------------------------------------
    grp_str   = strjoin(labels(sel_idx), '-');
    base_name = sprintf('LagTime_vs_GrowthRate_%s_%s', grp_str, timestamp_log);
    fig_path  = fullfile(p.out_dir, base_name);

    fprintf('Saving figure files...\n');
    savefig(hfig, [fig_path '.fig']);
    exportgraphics(hfig, [fig_path '.pdf'], 'ContentType','vector');
    exportgraphics(hfig, [fig_path '.png'], 'Resolution',300);

    ht_fprintf(HT_FLOG, '\nFigure saved:\n  %s.fig\n  %s.pdf\n  %s.png\n', ...
        fig_path, fig_path, fig_path);

    if HT_FLOG ~= -1
        fclose(HT_FLOG);
        fprintf('Log saved: %s\n', log_path);
    end

    set(hfig, 'Visible', 'on');
    drawnow;
end


%% ================================================================
%  LOCAL: loess_smooth
%  Fast LOESS via sub-sampled evaluation grid + linear interpolation.
%  Evaluates the weighted fit on n_eval evenly-spaced grid points
%  (default 200), then interpolates back to the original x locations.
%  This reduces complexity from O(N^2) to O(N * n_eval) ≈ O(N).
%  x must be sorted ascending.
%% ================================================================
function ys = loess_smooth(x, y, f, n_eval)
    if nargin < 4, n_eval = 200; end
    n      = length(x);
    h      = max(1, round(f * n));   % half-window in number of points
    x_grid = linspace(x(1), x(end), n_eval);
    y_grid = zeros(n_eval, 1);

    for gi = 1:n_eval
        x0 = x_grid(gi);
        % Distance from all data points to this grid node
        d  = abs(x - x0);
        % Find the h nearest neighbours
        [d_sort, ~] = sort(d);
        d_max = d_sort(min(h, n)) + eps;
        % Tricube weights
        u = min(d / d_max, 1);
        w = (1 - u.^3).^3;
        % Weighted linear fit
        sw  = sum(w);
        swx = sum(w .* x);
        swy = sum(w .* y);
        swxx= sum(w .* x.^2);
        swxy= sum(w .* x .* y);
        det = sw*swxx - swx^2;
        if abs(det) < 1e-12
            y_grid(gi) = swy / max(sw, eps);
        else
            b1 = (swxx*swy - swx*swxy) / det;
            b2 = (sw*swxy  - swx*swy)  / det;
            y_grid(gi) = b1 + b2*x0;
        end
    end

    % Interpolate grid values back to original x locations
    ys = interp1(x_grid, y_grid, x, 'linear', 'extrap');
    ys = ys(:);
end


%% ================================================================
%  LOCAL: viridis_cmap
%  11-point pchip interpolation of the viridis palette.
%  Matches R scale_color_viridis_c exactly.
%% ================================================================
function C = viridis_cmap(n)
    ctrl = [0.267, 0.005, 0.329;
            0.283, 0.141, 0.458;
            0.254, 0.265, 0.530;
            0.207, 0.372, 0.553;
            0.164, 0.471, 0.558;
            0.128, 0.567, 0.551;
            0.135, 0.659, 0.518;
            0.267, 0.749, 0.441;
            0.478, 0.821, 0.318;
            0.741, 0.873, 0.150;
            0.993, 0.906, 0.144];
    t_ctrl = linspace(0,1,size(ctrl,1));
    t_out  = linspace(0,1,n);
    C = [interp1(t_ctrl,ctrl(:,1),t_out,'pchip')', ...
         interp1(t_ctrl,ctrl(:,2),t_out,'pchip')', ...
         interp1(t_ctrl,ctrl(:,3),t_out,'pchip')'];
    C = max(0, min(1,C));
end


%% ================================================================
%  LOCAL: fmt_p  — format p-value for annotation/log text.
%  Spearman p (via corr()) can legitimately underflow to exact 0.0 in
%  double precision for large n / strong correlation — the true value
%  is just too small to represent, not actually zero — so that case is
%  shown as "<1e-300" rather than a misleading literal 0.
%% ================================================================
function s = fmt_p(p, fmt)
    if nargin < 2, fmt = '%.3g'; end
    if isnan(p)
        s = 'NaN';
    elseif p == 0
        s = '<1e-300';
    else
        s = sprintf(fmt, p);
    end
end


%% ================================================================
%  LOCAL: ask_group_selection
%% ================================================================
function sel_idx = ask_group_selection(HT_FLOG, ht)
    all_labels = ht.labels;
    n_avail    = length(all_labels);

    ht_fprintf(HT_FLOG, '\n============================================================\n');
    ht_fprintf(HT_FLOG, '  plot_lagTime_vs_growthRate -- Group Selection\n');
    ht_fprintf(HT_FLOG, '============================================================\n');
    ht_fprintf(HT_FLOG, '  +-----+--------+-------+----------+\n');
    ht_fprintf(HT_FLOG, '  | idx | label  | n_col | n_fitted |\n');
    ht_fprintf(HT_FLOG, '  +-----+--------+-------+----------+\n');
    for k = 1:n_avail
        n_fit = sum(ht.groups(k).growth.fit_ok);
        ht_fprintf(HT_FLOG, '  | %3d | %-6s | %5d | %8d |\n', ...
            k, all_labels{k}, ht.groups(k).n_colonies, n_fit);
    end
    ht_fprintf(HT_FLOG, '  +-----+--------+-------+----------+\n');
    ht_fprintf(HT_FLOG, '  Enter group labels (e.g. 3h 24h) or ALL:\n');

    while true
        raw = strtrim(input('  Your selection: ', 's'));
        if strcmpi(raw, 'all')
            sel_idx = 1:n_avail;
            ht_fprintf(HT_FLOG, '  -> All %d groups selected.\n', n_avail);
            return;
        end
        tokens  = regexp(raw, '[,\s]+', 'split');
        tokens  = tokens(~cellfun(@isempty, tokens));
        sel_idx = zeros(1, length(tokens));
        ok = true;
        for t = 1:length(tokens)
            match = find(strcmpi(all_labels, tokens{t}));
            if isempty(match)
                ht_fprintf(HT_FLOG, '  WARNING: "%s" not recognised. Try again.\n', tokens{t});
                ok = false; break;
            end
            sel_idx(t) = match(1);
        end
        if ~ok, continue; end
        [~, ui]  = unique(sel_idx,'stable');
        sel_idx  = sel_idx(ui);
        ht_fprintf(HT_FLOG, '  -> Selected: %s\n', strjoin(all_labels(sel_idx), ', '));
        return;
    end
end


%% ================================================================
%  LOCAL: check_ht
%% ================================================================
function check_ht(ht)
    if ~isstruct(ht) || ~isfield(ht,'groups')
        error(['Input must be the ht struct from preprocess_pipeline_data.\n' ...
               'Run:  preprocess_pipeline_data(data)  first.']);
    end
end


%% ================================================================
%  LOCAL: ht_fprintf
%% ================================================================
function ht_fprintf(HT_FLOG, varargin)
    fprintf(varargin{:});
    if ~isempty(HT_FLOG) && HT_FLOG ~= -1
        fprintf(HT_FLOG, varargin{:});
    end
end
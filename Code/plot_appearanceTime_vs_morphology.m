%% HeteroTyper Pipeline for Bright Plates
% Scatter plots of lag time vs selected morphological features,
% for selected sample groups.
% Requires preprocess_pipeline_data(data) to have been run first.
%
% USAGE:
%   plot_lagTime_vs_morphology(ht)

function plot_lagTime_vs_morphology(ht)
    HT_FLOG = -1;  % log file handle (opened after params are loaded)

    check_ht(ht);

    p      = ht.params;
    labels = ht.labels;
    colors = ht.colors;

    % ------------------------------------------------------------------
    %  Open log file — mirrors everything printed to the terminal
    % ------------------------------------------------------------------
    timestamp_log = datestr(now, 'yyyymmdd_HHMMSS');
    log_path = fullfile(p.out_dir, ...
                        sprintf('plot_lagTime_vs_morphology_%s.txt', timestamp_log));
    HT_FLOG = fopen(log_path, 'w');
    if HT_FLOG == -1
        warning('Could not open log file: %s', log_path);
        HT_FLOG = -1;
    end
    fprintf('Log file: %s\n', log_path);


    % ------------------------------------------------------------------
    %  All available morphological features
    % ------------------------------------------------------------------
    all_param_names  = {'Final Size','Area','Intensity','Mean Intensity', ...
                        'Intensity / Size','Perimeter','Circularity', ...
                        'Eccentricity','Solidity'};
    all_param_fields = {'size','area','intensity','mean_intensity', ...
                        'int_per_size','perimeter','circularity', ...
                        'eccentricity','solidity'};

    % Global y-limits per feature (consistent across groups)
    all_ylims = {[0, ht.global.size];   [0, ht.global.area]; ...
                 [0, ht.global.int];    [0, ht.global.mean_int]; ...
                 [0, ht.global.int_size]; [0, ht.global.peri]; ...
                 [0, 1]; [0, 1]; [0, 1]};

    % ------------------------------------------------------------------
    %  User prompts
    % ------------------------------------------------------------------
    sel_groups = ask_group_selection(HT_FLOG, ht);
    sel_params = ask_param_selection(HT_FLOG, all_param_names);

    n_groups = length(sel_groups);
    n_params = length(sel_params);

    param_names  = all_param_names(sel_params);
    param_fields = all_param_fields(sel_params);
    param_ylims  = all_ylims(sel_params);

    ht_fprintf(HT_FLOG, '\n--- Lag Time vs Morphology parameters ---\n');
    ht_fprintf(HT_FLOG, '  Groups   : %s\n', strjoin(labels(sel_groups), ', '));
    ht_fprintf(HT_FLOG, '  Features : %s\n', strjoin(param_names, ', '));
    ht_fprintf(HT_FLOG, '  Output   : %s\n', p.out_dir);
    ht_fprintf(HT_FLOG, '-----------------------------------------\n\n');

    % ------------------------------------------------------------------
    %  Figure layout — one panel per feature, all groups overlaid
    % ------------------------------------------------------------------
    panel_w  = 3.2;   panel_h  = 2.8;
    gap_x    = 0.85;   % wider gap: every panel has its own y-axis tick labels
    margin_l = 0.90;  margin_r = 1.60;   % right margin holds legend
    margin_b = 0.80;  margin_t = 0.35;

    fig_w = margin_l + n_params*panel_w + (n_params-1)*gap_x + margin_r;
    fig_h = margin_b + panel_h + margin_t;

    hfig = figure('Name','Lag Time vs Morphology', ...
                  'Color','w', ...
                  'Units','inches', ...
                  'Position',[1 1 fig_w fig_h], ...
                  'PaperUnits','inches', ...
                  'PaperSize',[fig_w fig_h], ...
                  'PaperPosition',[0 0 fig_w fig_h]);

    fs_tick  = 11;
    fs_label = 12;

    % Pre-create one axes per feature so we can overlay all groups
    ax_all = gobjects(1, n_params);
    ytick_offset = 0.55;   % inches reserved per panel for y-axis tick labels
    for pi = 1:n_params
        % Each panel (including those beyond the first) needs its own
        % left offset to accommodate y-axis tick labels.
        panel_left_offset = (pi-1) * (panel_w + gap_x);
        left = (margin_l + panel_left_offset) / fig_w;
        ax_all(pi) = axes('Position', ...   %#ok<LAXES>
            [left, margin_b/fig_h, panel_w/fig_w, panel_h/fig_h]);
        hold(ax_all(pi), 'on');
    end

    % ------------------------------------------------------------------
    %  Check whether morph_corr was computed (requires patched preprocess)
    % ------------------------------------------------------------------
    has_corr = isfield(ht.groups(1), 'morph_corr');
    if ~has_corr
        warning(['plot_lagTime_vs_morphology: ht.groups(g).morph_corr not found.\n' ...
                 'Re-run preprocess_pipeline_data to compute lag-morphology correlations.\n' ...
                 'Figures will be produced without stats annotations.']);
    end

    % ------------------------------------------------------------------
    %  Plot loop — iterate features in outer loop, groups in inner loop
    %  so all groups land on the same axes per feature
    % ------------------------------------------------------------------
    fs_stat = 9;   % annotation font size

    for pi = 1:n_params
        ax    = ax_all(pi);
        field = param_fields{pi};

        % ---- Log header for this feature ----
        ht_fprintf(HT_FLOG, '\n  Feature: %s\n', param_names{pi});
        ht_fprintf(HT_FLOG, '  %-10s  %8s  %8s  %10s  %8s  %10s  %s\n', ...
            'Group','r','R^2','p_pearson','rho','p_spearman','Equation');
        ht_fprintf(HT_FLOG, '  %s\n', repmat('-',1,76));

        for gi = 1:n_groups
            g   = sel_groups(gi);
            grp = ht.groups(g);

            lag  = grp.lag_time;
            feat = grp.(field);
            ok   = isfinite(lag) & isfinite(feat);

            scatter(ax, lag(ok), feat(ok), 12, colors(g,:), ...
                    'filled', 'MarkerFaceAlpha', 0.30, ...
                    'DisplayName', labels{g});

            % ---- Retrieve stats from preprocess ----
            if has_corr && isfield(grp.morph_corr, field)
                mc = grp.morph_corr.(field);
            else
                mc = struct('r',NaN,'r_sq',NaN,'p_pearson',NaN, ...
                            'rho',NaN,'p_spearman',NaN, ...
                            'slope',NaN,'intercept',NaN,'eq_str','');
            end

            % ---- 10% darker colour for regression line ----
            dark_col = max(0, colors(g,:) * 0.90);

            % ---- Linear regression line (spans actual data range only) ----
            if isfinite(mc.slope) && isfinite(mc.intercept)
                x_data_lo = min(lag(ok));
                x_data_hi = max(lag(ok));
                x_reg = linspace(x_data_lo, x_data_hi, 100);
                y_reg = mc.intercept + mc.slope * x_reg;
                % Clamp to y-axis limits so line doesn't escape the panel
                yl_now = param_ylims{pi};
                y_reg  = max(yl_now(1), min(yl_now(2), y_reg));
                plot(ax, x_reg, y_reg, '-', ...
                     'Color',     dark_col, ...
                     'LineWidth', 1.6, ...
                     'HandleVisibility', 'off');   % keep out of legend
            end

            % ---- Significance stars (Spearman, p < 0.05) ----
            % Spearman chosen: robust to non-normality and outliers,
            % appropriate for nested colony data from agar plates.
            sig_rho = star_str(mc.p_spearman);

            % ---- Rebuild equation string using slope/intercept ----
            if isfinite(mc.slope) && isfinite(mc.intercept)
                if mc.intercept >= 0
                    eq_fig = sprintf('y = %.3g x + %.3g', mc.slope, mc.intercept);
                else
                    eq_fig = sprintf('y = %.3g x - %.3g', mc.slope, abs(mc.intercept));
                end
            else
                eq_fig = '';
            end

            % ---- Figure annotation: rho, R^2, p-value + equation ----
            if isnan(mc.rho)
                ann = sprintf('%s  \rho=NaN  R^2=NaN  p=NaN', labels{g});
            elseif isempty(eq_fig)
                ann = sprintf('%s  \rho=%.3f  R^2=%.3f  p=%s %s', ...
                    labels{g}, mc.rho, mc.r_sq, fmt_p(mc.p_spearman), sig_rho);
            else
                ann = sprintf('%s  \rho=%.3f  R^2=%.3f  p=%s %s\n%s', ...
                    labels{g}, mc.rho, mc.r_sq, fmt_p(mc.p_spearman), sig_rho, eq_fig);
            end

            % Place annotation stacked per group — step down by group index
            yl      = param_ylims{pi};
            y_range = yl(2) - yl(1);
            y_ann   = yl(2) - (gi - 0.5) * y_range * 0.16;
            x_ann   = p.max_lag - 0.02*(p.max_lag - p.incTime);

            text(ax, x_ann, y_ann, ann, ...
                 'Color',               dark_col, ...
                 'FontSize',            fs_stat, ...
                 'FontWeight',          'bold', ...
                 'HorizontalAlignment', 'right', ...
                 'VerticalAlignment',   'middle', ...
                 'Interpreter',         'tex', ...
                 'Clipping',            'on');

            % ---- Log row ----
            ht_fprintf(HT_FLOG, '  %-10s  %8.4f  %8.4f  %10s  %8.4f  %10s  %s\n', ...
                labels{g}, mc.r, mc.r_sq, fmt_p(mc.p_pearson, '%.4g'), ...
                mc.rho, fmt_p(mc.p_spearman, '%.4g'), mc.eq_str);
        end

        xlim(ax, [p.incTime, p.max_lag]);
        xticks(ax, p.incTime:4:p.max_lag);
        ylim(ax, param_ylims{pi});

        ylabel(ax, param_names{pi}, ...
               'FontSize',fs_label, 'FontWeight','bold', ...
               'Interpreter','none');
        xlabel(ax, 'Lag Time (h)', 'FontSize',fs_label, 'FontWeight','bold');

        title(ax, param_names{pi}, ...
              'FontSize',fs_label, 'FontWeight','bold', ...
              'Interpreter','none');

        set(ax, 'FontSize',fs_tick, 'TickDir','out', 'Box','off');
        hold(ax, 'off');
    end

    % ------------------------------------------------------------------
    %  Single shared legend — placed in the right margin
    % ------------------------------------------------------------------
    lg = legend(ax_all(end), labels(sel_groups), ...
                'Location','eastoutside', ...
                'FontSize', fs_label, ...
                'Box','off', ...
                'Interpreter','none');
    % Nudge legend into the right margin
    lg.Position(1) = (margin_l + n_params*(panel_w+gap_x) - gap_x + 0.15) / fig_w;
    lg.Position(2) = (margin_b + panel_h/2) / fig_h - lg.Position(4)/2;

    % ------------------------------------------------------------------
    %  Save (.fig, .pdf, .png)
    % ------------------------------------------------------------------
    out_dir   = p.out_dir;
    timestamp = datestr(now, 'yyyymmdd_HHMMSS');
    grp_str   = strjoin(labels(sel_groups), '-');
    base_name = sprintf('LagTime_vs_Morphology_%s_%s', grp_str, timestamp);
    fig_path  = fullfile(out_dir, base_name);

    savefig(hfig,  [fig_path '.fig']);
    exportgraphics(hfig, [fig_path '.pdf'], 'ContentType','vector');
    exportgraphics(hfig, [fig_path '.png'], 'Resolution',300);

    ht_fprintf(HT_FLOG, 'Figure saved:\n  %s.fig\n  %s.pdf\n  %s.png\n', ...
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
% ================================================================
function sel_idx = ask_group_selection(HT_FLOG, ht)

    all_labels = ht.labels;
    n_avail    = length(all_labels);

    ht_fprintf(HT_FLOG, '\nAvailable sample groups:\n');
    ht_fprintf(HT_FLOG, '  +-----+--------+-------------+\n');
    ht_fprintf(HT_FLOG, '  | idx | label  | n colonies  |\n');
    ht_fprintf(HT_FLOG, '  +-----+--------+-------------+\n');
    for k = 1:n_avail
        ht_fprintf(HT_FLOG, '  | %3d | %-6s | %11d |\n', k, all_labels{k}, ht.groups(k).n_colonies);
    end
    ht_fprintf(HT_FLOG, '  +-----+--------+-------------+\n');
    ht_fprintf(HT_FLOG, '  Enter group labels separated by spaces or commas  (e.g.  3h 24h)\n');
    ht_fprintf(HT_FLOG, '  or type  ALL  to include all groups.\n');

    while true
        raw = strtrim(input('  Your selection: ', 's'));
        if strcmpi(raw, 'all')
            sel_idx = 1:n_avail;
            ht_fprintf(HT_FLOG, '  -> All %d groups selected.\n', n_avail);
            return;
        end
        tokens = regexp(strtrim(raw), '[,\s]+', 'split');
        tokens = tokens(~cellfun(@isempty, tokens));
        sel_idx = zeros(1, length(tokens));
        ok = true;
        for t = 1:length(tokens)
            match = find(strcmpi(all_labels, tokens{t}));
            if isempty(match)
                ht_fprintf(HT_FLOG, '  WARNING  "%s" not recognised. Try again.\n', tokens{t});
                ok = false; break;
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
%  LOCAL FUNCTION: ask_param_selection
% ================================================================
function sel_idx = ask_param_selection(HT_FLOG, param_names)

    n = length(param_names);
    ht_fprintf(HT_FLOG, '\nAvailable morphological features to plot against lag time:\n');
    ht_fprintf(HT_FLOG, '  +-----+---------------------+\n');
    for k = 1:n
        ht_fprintf(HT_FLOG, '  | %3d | %-19s |\n', k, param_names{k});
    end
    ht_fprintf(HT_FLOG, '  +-----+---------------------+\n');
    ht_fprintf(HT_FLOG, '  Enter feature numbers separated by commas  (e.g.  1,3,5)\n');
    ht_fprintf(HT_FLOG, '  or type  ALL  to include all %d features.\n', n);

    while true
        raw = strtrim(input('  Your selection: ', 's'));
        if strcmpi(raw, 'all')
            sel_idx = 1:n;
            ht_fprintf(HT_FLOG, '  -> All %d features selected.\n', n);
            return;
        end
        tokens = regexp(strtrim(raw), '[,\s]+', 'split');
        nums   = str2double(tokens);
        if any(isnan(nums)) || any(nums < 1) || any(nums > n) || any(nums ~= floor(nums))
            ht_fprintf(HT_FLOG, '  WARNING  Invalid input. Enter integers 1-%d or ALL.\n', n);
            continue;
        end
        [~, ui] = unique(nums, 'stable');
        sel_idx = nums(ui);
        ht_fprintf(HT_FLOG, '  -> Selected features: %s\n', strjoin(param_names(sel_idx), ', '));
        return;
    end
end


%% ================================================================
%  LOCAL FUNCTION: star_str
%  Returns significance stars for a p-value (threshold p < 0.05).
% ================================================================
function s = star_str(p)
    if isnan(p),       s = '';
    elseif p < 0.001,  s = '***';
    elseif p < 0.01,   s = '**';
    elseif p < 0.05,   s = '*';
    else,              s = 'ns';
    end
end


%% ================================================================
%  LOCAL FUNCTION: fmt_p
%  Text display for a p-value. Pearson/Spearman p (via corr()) can
%  legitimately underflow to exact 0.0 in double precision for large
%  n / strong correlation — the true value is just too small to
%  represent, not actually zero — shown as "<1e-300" instead.
% ================================================================
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
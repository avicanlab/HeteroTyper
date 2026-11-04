%% HeteroTyper Pipeline
%  plot_appearanceTime_vs_growthRate  —  Appearance time vs maximum Gompertz growth rate
%  for any set of sample groups (replaces plot_appearanceTime_vs_growthRate,
%  plot_appearanceTime_vs_growthRate_6s_stress and
%  plot_appearanceTime_vs_growthRate_6srna_growth).
%
%  Scatter plot of appearance time (h) vs mu_max (px/h) per selected group.
%
%  INTERACTIVE CHOICES (in this order)
%    1. Sample groups           labels / index numbers / ALL
%    2. Display order           panel order, as numbers of the selected
%                               groups (Enter = keep selection order)
%    3. Statistical comparisons all possible pairs of the selected groups
%                               are listed and numbered; enter the pair
%                               numbers to include (e.g. 1 4 6, ranges 1-3,
%                               ALL; Enter = none)
%
%  FIGURE 1 (local colony density):
%    One panel per group; every colony coloured by log10(local_density)
%    on a shared viridis scale. A LOESS smooth line (red) is overlaid per
%    panel.
%
%  FIGURE 2 (individual colonies by sample group):
%    Identical panel grid, axes and LOESS line as Figure 1, but every
%    colony is an individual point in its sample-group colour (ht.colors)
%    with alpha 0.4. No colourbar.
%
%  FIGURE 3 (pairwise overlay — one panel per chosen comparison):
%    The two groups of the pair overlaid; colonies in their sample-group
%    colour with alpha 0.25, one LOESS line per group in a darker shade of
%    the group colour, legend. Grid of up to 3 columns. Shared x-axis
%    (same as Figures 1-2); each panel has its own y-axis limit (0 to a
%    rounded value above the highest mu_max of the pair, ~15 % headroom for
%    the legend), so panels are NOT directly comparable by height.
%    Subtitle = significance stars of the two comparisons below (p_BH).
%    Skipped when no comparison is chosen.
%
%  PAIRWISE COMPARISONS (chosen pairs only)
%    mu_max : Mann–Whitney U test (two-sided, normal approximation with tie
%             and continuity correction) of mu_max between the two groups;
%             effect size = rank-biserial r (>0: first group higher) and the
%             median ratio.
%    rho    : comparison of the two Spearman rho (appearance time vs mu_max, from
%             STEP 6d) by Fisher z with the Bonett–Wright variance
%             (1 + rho^2/2) / (n - 3); is the appearance–growth relationship
%             equally strong in both groups?
%    Each test is Benjamini–Hochberg adjusted across the chosen pairs.
%    Colonies of all positions are pooled and treated as independent;
%    large n makes small differences significant.
%
%  Per-group statistics (n, Spearman rho / p / stars, rho^2, the density-
%  weighted rho_w and the LOESS polynomial) and the comparisons are written
%  to the log and to AppearanceTime_vs_GrowthRate_stats_<groups>_<ts>.xlsx.
%
%  Gompertz fitting and the per-group statistics are pre-computed in STEP
%  6d of preprocess_pipeline_data. No toolbox function is needed here.
%
%  USAGE:
%    plot_appearanceTime_vs_growthRate(ht)

function plot_appearanceTime_vs_growthRate(ht)
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
        sprintf('plot_appearanceTime_vs_growthRate_%s.txt', timestamp_log));
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
        error(['plot_appearanceTime_vs_growthRate: ht.groups(g).growth not found.\n' ...
               'Re-run preprocess_pipeline_data (STEP 6d) first.']);
    end

    % ------------------------------------------------------------------
    %  Group selection, display order, comparisons
    % ------------------------------------------------------------------
    sel_idx  = ask_group_selection(HT_FLOG, ht);
    sel_idx  = ask_order(HT_FLOG, sel_idx, labels);
    n_groups = length(sel_idx);
    pairs    = ask_pairs(HT_FLOG, labels(sel_idx));   % [P x 2] indices k into G (display order)

    ht_fprintf(HT_FLOG, '\n--- Appearance Time vs Growth Rate ---\n');
    ht_fprintf(HT_FLOG, '  Groups      : %s\n', strjoin(labels(sel_idx), ', '));
    if isempty(pairs)
        ht_fprintf(HT_FLOG, '  Comparisons : none\n');
    else
        ht_fprintf(HT_FLOG, '  Comparisons : %s\n', strjoin(pair_labels(labels(sel_idx), pairs), ', '));
    end
    ht_fprintf(HT_FLOG, '  Output      : %s\n', p.out_dir);
    ht_fprintf(HT_FLOG, '-------------------------------\n\n');

    % ------------------------------------------------------------------
    %  Collect pooled data for shared axis limits
    % ------------------------------------------------------------------
    all_mu  = [];
    all_appearance = [];
    for k = 1:n_groups
        g  = sel_idx(k);
        ok = ht.groups(g).growth.fit_ok & isfinite(ht.groups(g).growth.mu_max);
        all_mu  = [all_mu;  ht.groups(g).growth.mu_max(ok)];  %#ok<AGROW>
        all_appearance = [all_appearance; ht.groups(g).appearance_time(ok)];        %#ok<AGROW>
    end

    if isempty(all_mu)
        ht_fprintf(HT_FLOG, 'No successful Gompertz fits found.\n');
        if HT_FLOG ~= -1, fclose(HT_FLOG); end
        return;
    end

    y_upper = ceil(max(all_mu) / 100) * 100;
    y_step  = 100;
    x_lo    = floor(min(all_appearance) / 4) * 4;
    x_hi    = ceil(max(all_appearance)  / 4) * 4;

    % ------------------------------------------------------------------
    %  Per-group data + statistics (computed once; written to the log and
    %  the xlsx table)
    % ------------------------------------------------------------------
    G = struct('g', {}, 'appearance', {}, 'mu', {}, 'dens', {}, 'n_fit', {}, ...
               'appearance_sort', {}, 'mu_loess', {}, 'p_coef', {}, 'loess_eq', {}, 'sig', {});
    for k = 1:n_groups
        g   = sel_idx(k);
        grp = ht.groups(g);
        ok  = grp.growth.fit_ok & isfinite(grp.growth.mu_max) & isfinite(grp.appearance_time);
        G(k).g     = g;
        G(k).appearance   = grp.appearance_time(ok);
        G(k).mu    = grp.growth.mu_max(ok);
        G(k).dens  = grp.growth.local_density(ok);
        G(k).n_fit = sum(ok);
        G(k).p_coef   = [NaN NaN NaN];
        G(k).loess_eq = '';
        if G(k).n_fit >= 10
            [G(k).appearance_sort, si] = sort(G(k).appearance);
            mu_sort       = G(k).mu(si);
            G(k).mu_loess = max(0, min(y_upper, loess_smooth(G(k).appearance_sort, mu_sort, 0.75)));
            % degree-2 polynomial fitted to the LOESS curve (concise equation)
            G(k).p_coef   = polyfit(G(k).appearance_sort, G(k).mu_loess, 2);
            G(k).loess_eq = sprintf('y = %.3g x^2 %+.3g x %+.3g', G(k).p_coef);
        end
        G(k).sig = sig_stars(grp.growth.spearman_p);
    end

    % ------------------------------------------------------------------
    %  Pairwise comparisons (chosen pairs only), BH across pairs per test
    % ------------------------------------------------------------------
    C = struct('i', {}, 'j', {}, 'label', {}, 'n1', {}, 'n2', {}, 'med1', {}, 'med2', {}, ...
               'U', {}, 'z_mw', {}, 'p_mw', {}, 'p_mw_bh', {}, 'r_rb', {}, ...
               'rho1', {}, 'rho2', {}, 'z_rho', {}, 'p_rho', {}, 'p_rho_bh', {});
    for r = 1:size(pairs, 1)
        i = pairs(r, 1);  j = pairs(r, 2);
        st_i = ht.groups(G(i).g).growth;  st_j = ht.groups(G(j).g).growth;
        [U, z_mw, p_mw, r_rb] = mann_whitney(G(i).mu, G(j).mu);
        [z_rho, p_rho] = compare_spearman(st_i.spearman_rho, G(i).n_fit, ...
                                          st_j.spearman_rho, G(j).n_fit);
        C(end+1) = struct('i', i, 'j', j, ...
            'label', sprintf('%s vs %s', labels{G(i).g}, labels{G(j).g}), ...
            'n1', G(i).n_fit, 'n2', G(j).n_fit, ...
            'med1', median_or_nan(G(i).mu), 'med2', median_or_nan(G(j).mu), ...
            'U', U, 'z_mw', z_mw, 'p_mw', p_mw, 'p_mw_bh', NaN, 'r_rb', r_rb, ...
            'rho1', st_i.spearman_rho, 'rho2', st_j.spearman_rho, ...
            'z_rho', z_rho, 'p_rho', p_rho, 'p_rho_bh', NaN); %#ok<AGROW>
    end
    if ~isempty(C)
        q1 = bh_adjust([C.p_mw]);   q2 = bh_adjust([C.p_rho]);
        for r = 1:numel(C), C(r).p_mw_bh = q1(r);  C(r).p_rho_bh = q2(r); end
    end

    % ------------------------------------------------------------------
    %  Figure layout — grid of panels as square as possible.
    %  Colorbar is a dedicated axes so panels are never squeezed.
    % ------------------------------------------------------------------
    n_cols = ceil(sqrt(n_groups));
    n_rows = ceil(n_groups / n_cols);

    panel_w  = 3.5;   panel_h  = 3.2;
    gap_x    = 0.65;  gap_y    = 0.75;
    margin_l = 1.00;  margin_r = 0.30;
    margin_b = 0.85;  margin_t = 0.45;
    cb_w     = 0.18;  cb_gap   = 0.35;   % colourbar bar width / gap to grid (in)
    cb_lbl_w = 1.00;                     % room for colourbar tick labels + title (in)

    grid_w = n_cols*panel_w + (n_cols-1)*gap_x;
    grid_h = n_rows*panel_h + (n_rows-1)*gap_y;
    fig_h  = margin_b + grid_h + margin_t;
    % Figure 1 adds the colourbar (bar + its labels) to the right of the
    % grid; Figure 2 uses the identical panel grid without a colourbar.
    fig_w_all = [margin_l + grid_w + cb_gap + cb_w + cb_lbl_w, ...
                 margin_l + grid_w + margin_r];

    fs_tick  = 12;
    fs_label = 13;
    fs_annot = 13;

    cmap = viridis_cmap(256);

    % ------------------------------------------------------------------
    %  Log tables
    % ------------------------------------------------------------------
    ht_fprintf(HT_FLOG, '%-14s  %6s  %10s  %10s  %4s  %12s  %s\n', ...
        'Group','n_fit','Spearman_r','Spearman_p','sig','wt_Spearman','LOESS poly (deg 2)');
    ht_fprintf(HT_FLOG, '%s\n', repmat('-',1,100));
    for k = 1:n_groups
        st = ht.groups(G(k).g).growth;
        ht_fprintf(HT_FLOG, '%-14s  %6d  %10.4f  %10s  %4s  %12.4f  %s\n', ...
            labels{G(k).g}, G(k).n_fit, st.spearman_rho, fmt_p(st.spearman_p, '%.4g'), ...
            G(k).sig, st.wt_spearman, G(k).loess_eq);
    end
    ht_fprintf(HT_FLOG, '\n');

    if ~isempty(C)
        ht_fprintf(HT_FLOG, 'PAIRWISE COMPARISONS (p_BH: Benjamini–Hochberg across the %d chosen pairs, per test)\n', numel(C));
        ht_fprintf(HT_FLOG, '%-32s  %9s  %9s  %10s  %4s  %7s  %7s  %7s  %10s  %4s\n', ...
            'Pair', 'med mu_1', 'med mu_2', 'p_BH (MW)', 'sig', 'r_rb', 'rho_1', 'rho_2', 'p_BH (rho)', 'sig');
        ht_fprintf(HT_FLOG, '%s\n', repmat('-',1,112));
        for r = 1:numel(C)
            c = C(r);
            ht_fprintf(HT_FLOG, '%-32s  %9.3g  %9.3g  %10s  %4s  %7.3f  %7.3f  %7.3f  %10s  %4s\n', ...
                c.label, c.med1, c.med2, fmt_p(c.p_mw_bh, '%.4g'), stars_na(c.p_mw_bh), c.r_rb, ...
                c.rho1, c.rho2, fmt_p(c.p_rho_bh, '%.4g'), stars_na(c.p_rho_bh));
        end
        ht_fprintf(HT_FLOG, ['  MW: Mann–Whitney U on mu_max (r_rb = rank-biserial r, >0: first group higher).\n' ...
                             '  rho: Fisher z comparison of the two Spearman rho (Bonett–Wright variance).\n' ...
                             '  Colonies pooled over positions and treated as independent.\n\n']);
    end

    % Pre-compute global log10-density range across ALL selected groups so
    % every panel uses the same colour scale and the colorbar is accurate.
    all_log_dens = [];
    for k = 1:n_groups
        d_pre = G(k).dens;
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
    %  Figures 1 and 2 with the identical panel grid:
    %    1 = colonies coloured by local colony density (default)
    %    2 = colonies coloured by sample group, alpha 0.4
    % ------------------------------------------------------------------
    fig_names = {'Appearance Time vs Growth Rate - local colony density', ...
                 'Appearance Time vs Growth Rate - colonies by sample group', ...
                 'Appearance Time vs Growth Rate - pairwise overlay'};
    fig_stems = {'AppearanceTime_vs_GrowthRate', 'AppearanceTime_vs_GrowthRate_byGroup', ...
                 'AppearanceTime_vs_GrowthRate_pairs_overlay'};
    grp_str   = regexprep(strjoin(labels(sel_idx), '-'), '\s+', '');

    for fig_mode = 1:2
        fig_w = fig_w_all(fig_mode);
        hfig = figure('Name',fig_names{fig_mode}, 'Color','w', 'Units','inches', ...
                      'Position',[1 1 fig_w fig_h], 'PaperUnits','inches', ...
                      'PaperSize',[fig_w fig_h], 'PaperPosition',[0 0 fig_w fig_h], ...
                      'Visible','off');

        for k = 1:n_groups
            g = G(k).g;

            % Panel position
            col_k  = mod(k-1, n_cols);
            row_k  = floor((k-1) / n_cols);
            left   = (margin_l + col_k*(panel_w+gap_x)) / fig_w;
            bottom = (margin_b + (n_rows-1-row_k)*(panel_h+gap_y)) / fig_h;
            ax = axes('Position',[left, bottom, panel_w/fig_w, panel_h/fig_h]);
            hold(ax, 'on');

            if fig_mode == 1
                % Figure 1: colour = log10(local_density) on the shared viridis scale
                dens_g = G(k).dens;
                if ~isempty(dens_g) && any(dens_g > 0)
                    d_norm    = (log10(max(dens_g, 1)) - d_min_global) / max(d_max_global - d_min_global, eps);
                    c_idx     = max(1, min(256, round(max(0, min(1, d_norm))*255) + 1));
                    pt_colors = cmap(c_idx, :);
                else
                    pt_colors = repmat(colors(g,:), G(k).n_fit, 1);
                end
                scatter(ax, G(k).appearance, G(k).mu, 27, pt_colors, 'filled', 'MarkerFaceAlpha', 0.80);
            else
                % Figure 2: individual colonies in the sample-group colour, alpha 0.4
                scatter(ax, G(k).appearance, G(k).mu, 27, colors(g,:), 'filled', ...
                        'MarkerFaceAlpha', 0.40, 'MarkerEdgeColor', 'none');
            end

            % LOESS smooth (red); its equation is in the log / xlsx table
            if G(k).n_fit >= 10
                plot(ax, G(k).appearance_sort, G(k).mu_loess, '-', 'Color','r', 'LineWidth',1.5);
            end

            style_panel(ax, x_lo, x_hi, y_upper, y_step, fs_tick);

            % Group label (top-left, group colour)
            text(ax, x_lo + 0.03*(x_hi-x_lo), y_upper*0.98, labels{g}, ...
                 'FontSize',fs_annot, 'FontWeight','bold', 'Color',colors(g,:), ...
                 'HorizontalAlignment','left', 'VerticalAlignment','top', ...
                 'Interpreter','none');

            if row_k == n_rows-1 || k == n_groups
                xlabel(ax, 'Appearance time (h)', 'FontSize',fs_label, 'FontWeight','bold');
            end
            if col_k == 0
                ylabel(ax, 'Maximum growth rate (px/h)', 'FontSize',fs_label, 'FontWeight','bold');
            else
                set(ax, 'YTickLabel', {});
            end
            add_mirror_axes(ax, fs_tick);
        end

        % Figure 1 only: standalone colorbar right of the grid
        if fig_mode == 1
            cb_left   = (margin_l + grid_w + cb_gap) / fig_w;
            cb_bottom = margin_b / fig_h;
            cb_height = grid_h / fig_h;
            ax_cb = axes('Position',[cb_left, cb_bottom, cb_w/fig_w, cb_height], 'Visible','off');
            colormap(ax_cb, cmap);
            clim(ax_cb, [d_min_global, d_max_global]);
            cb = colorbar(ax_cb, 'eastoutside');
            cb.Units          = 'normalized';
            cb.Position       = [cb_left, cb_bottom, cb_w/fig_w, cb_height];
            cb.Label.String   = 'Colonies (log_{10})';
            cb.Label.FontSize = fs_label;
            cb.FontSize       = fs_tick;
            cb.TickDirection  = 'out';
        end

        save_figure(hfig, p.out_dir, fig_stems{fig_mode}, grp_str, timestamp_log, ...
                    fig_names{fig_mode}, HT_FLOG);
    end   % fig_mode

    % ------------------------------------------------------------------
    %  FIGURE 3 — one panel per chosen comparison, the two groups overlaid,
    %  grid of up to 3 columns. Every colony in its sample-group colour,
    %  alpha 0.25; one LOESS line per group (SHOW_LOESS_FIG3 = false hides
    %  them). Subtitle: stars of the two comparisons (p_BH).
    % ------------------------------------------------------------------
    SHOW_LOESS_FIG3 = true;

    if isempty(C)
        ht_fprintf(HT_FLOG, '\nFigure 3 skipped: no pairwise comparison chosen.\n');
    else
        n_p3     = numel(C);
        n_c3     = min(n_p3, 3);              % grid: up to 3 columns
        n_r3     = ceil(n_p3 / n_c3);
        gap_x3   = gap_x + 0.25;              % room for each panel's own y tick labels
        gap_y3   = gap_y + 0.35;              % room for title + subtitle
        note_h   = 0.40;                      % footnote below the grid (in)
        margin_t3 = margin_t + 0.35;
        fig_w  = margin_l + n_c3*panel_w + (n_c3-1)*gap_x3 + margin_r;
        fig_h3 = note_h + margin_b + n_r3*panel_h + (n_r3-1)*gap_y3 + margin_t3;
        hfig = figure('Name',fig_names{3}, 'Color','w', 'Units','inches', ...
                      'Position',[1 1 fig_w fig_h3], 'PaperUnits','inches', ...
                      'PaperSize',[fig_w fig_h3], 'PaperPosition',[0 0 fig_w fig_h3], ...
                      'Visible','off');
        for pp = 1:n_p3
            c  = C(pp);
            kk = [c.i, c.j];
            c3 = mod(pp-1, n_c3);  r3 = floor((pp-1) / n_c3);
            left   = (margin_l + c3*(panel_w+gap_x3)) / fig_w;
            bottom = (note_h + margin_b + (n_r3-1-r3)*(panel_h+gap_y3)) / fig_h3;
            ax = axes('Position',[left, bottom, panel_w/fig_w, panel_h/fig_h3]);
            hold(ax, 'on');
            h_leg = gobjects(0);
            for k = kk
                g = G(k).g;
                scatter(ax, G(k).appearance, G(k).mu, 27, colors(g,:), 'filled', ...
                        'MarkerFaceAlpha', 0.25, 'MarkerEdgeColor', 'none');
                % legend proxy at full opacity (alpha-0.25 dots are too pale in a legend)
                h_leg(end+1) = plot(ax, NaN, NaN, 'o', 'LineStyle','none', 'MarkerSize', 7, ...
                                    'MarkerFaceColor', colors(g,:), 'MarkerEdgeColor', colors(g,:), ...
                                    'DisplayName', labels{g}); %#ok<AGROW>
            end
            if SHOW_LOESS_FIG3
                for k = kk
                    if G(k).n_fit >= 10
                        plot(ax, G(k).appearance_sort, G(k).mu_loess, '-', 'Color', shade(colors(G(k).g,:), -0.35), ...
                             'LineWidth', 2);
                    end
                end
            end
            % individual y-axis limit for this panel
            mu_p = vertcat(G(kk).mu);
            [y_up_p, y_step_p] = nice_y_limit(max([mu_p; 0]) * 1.15);
            ht_fprintf(HT_FLOG, '  Figure 3, %s: y-axis 0-%g px/h (step %g)\n', c.label, y_up_p, y_step_p);
            style_panel(ax, x_lo, x_hi, y_up_p, y_step_p, fs_tick);
            title(ax, c.label, 'FontSize', fs_label, 'FontWeight','bold', 'Interpreter','none');
            subtitle(ax, sprintf('\\mu_{max}: %s     \\rho: %s', stars_na(c.p_mw_bh), stars_na(c.p_rho_bh)), ...
                     'FontSize', fs_tick - 1, 'Interpreter','tex');
            if pp + n_c3 > n_p3                % no panel below in this column
                xlabel(ax, 'Appearance time (h)', 'FontSize',fs_label, 'FontWeight','bold');
            end
            if c3 == 0
                ylabel(ax, 'Maximum growth rate (px/h)', 'FontSize',fs_label, 'FontWeight','bold');
            end
            legend(ax, h_leg, 'Location','northeast', 'Interpreter','none', ...
                   'FontSize', fs_tick - 1, 'Box','off');
            add_mirror_axes(ax, fs_tick);
        end
        note = sprintf(['\\mu_{max}: Mann–Whitney U test of maximum growth rate.  \\rho: comparison of the Spearman ' ...
                        '\\rho (appearance time vs \\mu_{max}) by Fisher z.  p_{BH}: Benjamini–Hochberg across the %d ' ...
                        'comparisons;  * p < 0.05,  ** p < 0.01,  *** p < 0.001,  ns = not significant'], n_p3);
        nax = axes('Parent', hfig, 'Position', [0.25/fig_w, 0.05/fig_h3, (fig_w - 0.5)/fig_w, (note_h - 0.05)/fig_h3], ...
                   'Visible','off', 'XLim',[0 1], 'YLim',[0 1]);
        text(nax, 0, 0.5, note, 'FontSize', 9, 'HorizontalAlignment','left', ...
             'VerticalAlignment','middle', 'Interpreter','tex');
        save_figure(hfig, p.out_dir, fig_stems{3}, grp_str, timestamp_log, fig_names{3}, HT_FLOG);
    end

    % ------------------------------------------------------------------
    %  Statistics table (xlsx)
    % ------------------------------------------------------------------
    xl_path = fullfile(p.out_dir, sprintf('AppearanceTime_vs_GrowthRate_stats_%s_%s.xlsx', grp_str, timestamp_log));
    try
        c_grp = cell(n_groups,1);  c_p = cell(n_groups,1);  c_sig = cell(n_groups,1);
        c_eq  = cell(n_groups,1);
        v_n = zeros(n_groups,1);  v_rho = nan(n_groups,1);  v_r2 = nan(n_groups,1);
        v_wr = nan(n_groups,1);   v_a = nan(n_groups,3);    v_med = nan(n_groups,1);
        for k = 1:n_groups
            st = ht.groups(G(k).g).growth;
            c_grp{k} = labels{G(k).g};
            v_n(k)   = G(k).n_fit;
            v_med(k) = median_or_nan(G(k).mu);
            v_rho(k) = st.spearman_rho;
            c_p{k}   = fmt_p(st.spearman_p, '%.4g');
            c_sig{k} = G(k).sig;
            v_r2(k)  = st.spearman_rho^2;
            v_wr(k)  = st.wt_spearman;
            v_a(k,:) = G(k).p_coef;
            c_eq{k}  = G(k).loess_eq;
        end
        T = table(c_grp, v_n, v_med, v_rho, c_p, c_sig, v_r2, v_wr, v_a(:,1), v_a(:,2), v_a(:,3), c_eq, ...
            'VariableNames', {'Group','n_fit','Median_mu_max','Spearman_rho','Spearman_p','Sig', ...
                              'Spearman_rho_squared','Weighted_Spearman_rho', ...
                              'LOESS_poly_a2','LOESS_poly_a1','LOESS_poly_a0','LOESS_poly_equation'});
        writetable(T, xl_path, 'Sheet', 'Stats');

        if ~isempty(C)
            n_c = numel(C);
            cr = cell(n_c, 16);
            for r = 1:n_c
                c = C(r);
                cr(r, :) = {c.label, c.n1, c.n2, c.med1, c.med2, c.med1 / c.med2, c.U, c.r_rb, ...
                            fmt_p(c.p_mw, '%.4g'), fmt_p(c.p_mw_bh, '%.4g'), stars_na(c.p_mw_bh), ...
                            c.rho1, c.rho2, fmt_p(c.p_rho, '%.4g'), fmt_p(c.p_rho_bh, '%.4g'), stars_na(c.p_rho_bh)};
            end
            writetable(cell2table(cr, 'VariableNames', {'Pair','n_1','n_2','Median_mu_max_1','Median_mu_max_2', ...
                       'Median_ratio','MW_U','Rank_biserial_r','MW_p','MW_p_BH','MW_Sig', ...
                       'Spearman_rho_1','Spearman_rho_2','Rho_diff_p','Rho_diff_p_BH','Rho_diff_Sig'}), ...
                       xl_path, 'Sheet', 'Comparisons');
        end

        readme = { ...
            'Column', 'Description'; ...
            'n_fit', 'Colonies with a successful Gompertz fit and a finite appearance time (points in every figure)'; ...
            'Median_mu_max', 'Median maximum growth rate of these colonies (px/h)'; ...
            'Spearman_rho / Spearman_p / Sig', 'Spearman correlation of appearance time vs mu_max per group (preprocess_pipeline_data STEP 6d); p as text, "<1e-300" if it underflows; Sig: *** p<0.001, ** p<0.01, * p<0.05, ns'; ...
            'Spearman_rho_squared', 'rho^2 (the squared Spearman rho, not an OLS R^2)'; ...
            'Weighted_Spearman_rho', 'Local-density-weighted Spearman rho (rho_w, STEP 6d)'; ...
            'LOESS_poly_a2/a1/a0 / equation', 'Degree-2 polynomial fitted to the LOESS curve (span 0.75): y = a2 x^2 + a1 x + a0; empty/NaN if n_fit < 10'; ...
            'Comparisons sheet', ['Pairs chosen by the user (first vs second group). MW_*: Mann–Whitney U test of mu_max (two-sided, normal approximation, ' ...
                                  'tie and continuity corrected); MW_U = U of the first group; Rank_biserial_r = 2U/(n_1 n_2) - 1 (>0: first group higher); ' ...
                                  'Median_ratio = first / second. Rho_diff_*: difference between the two Spearman rho, Fisher z test with Bonett–Wright variance (1 + rho^2/2)/(n - 3). ' ...
                                  'p_BH = Benjamini–Hochberg across the chosen pairs, separately per test; these are the stars in Figure 3. ' ...
                                  'Colonies pooled over positions and treated as independent.']; ...
        };
        writetable(cell2table(readme(2:end,:), 'VariableNames', readme(1,:)), xl_path, 'Sheet', 'README');
        ht_fprintf(HT_FLOG, '\nStatistics table saved: %s\n', xl_path);
    catch ME
        ht_fprintf(HT_FLOG, '\nWARNING: could not write xlsx (%s)\n', ME.message);
    end

    if HT_FLOG ~= -1
        fclose(HT_FLOG);
        fprintf('Log saved: %s\n', log_path);
    end
end


%% ================================================================
%  LOCAL: mann_whitney — two-sided Mann–Whitney U test (normal
%  approximation with tie and continuity correction), no toolbox.
%  U = U of x; r_rb = rank-biserial r = 2U/(n1 n2) - 1 (>0: x higher).
%% ================================================================
function [U, z, p, r_rb] = mann_whitney(x, y)
    x = x(isfinite(x));  y = y(isfinite(y));
    n1 = numel(x);  n2 = numel(y);
    U = NaN;  z = NaN;  p = NaN;  r_rb = NaN;
    if n1 < 2 || n2 < 2, return; end
    [rk, tie_sum] = tied_rank([x(:); y(:)]);
    N    = n1 + n2;
    U    = sum(rk(1:n1)) - n1*(n1 + 1)/2;
    r_rb = 2*U/(n1*n2) - 1;
    sig2 = n1*n2/12 * ((N + 1) - tie_sum/(N*(N - 1)));
    if sig2 <= 0, z = 0; p = 1; return; end
    d = U - n1*n2/2;
    z = (d - 0.5*sign(d)) / sqrt(sig2);
    p = min(1, erfc(abs(z)/sqrt(2)));
end


function [r, tie_sum] = tied_rank(v)
% Average ranks for ties; tie_sum = sum(t^3 - t) over tie groups.
    v = v(:);  N = numel(v);
    [s, o] = sort(v);
    r = zeros(N, 1);  tie_sum = 0;
    i = 1;
    while i <= N
        j = i;
        while j < N && s(j+1) == s(i), j = j + 1; end
        r(o(i:j)) = (i + j)/2;
        t = j - i + 1;
        tie_sum = tie_sum + (t^3 - t);
        i = j + 1;
    end
end


%% ================================================================
%  LOCAL: compare_spearman — two-sided test of rho1 = rho2 (independent
%  groups): Fisher z with the Bonett–Wright variance (1 + rho^2/2)/(n - 3).
%% ================================================================
function [z, p] = compare_spearman(r1, n1, r2, n2)
    z = NaN;  p = NaN;
    if ~all(isfinite([r1 r2])) || n1 < 4 || n2 < 4, return; end
    r1 = max(min(r1, 0.9999), -0.9999);
    r2 = max(min(r2, 0.9999), -0.9999);
    se = sqrt((1 + r1^2/2)/(n1 - 3) + (1 + r2^2/2)/(n2 - 3));
    z  = (atanh(r1) - atanh(r2)) / se;
    p  = min(1, erfc(abs(z)/sqrt(2)));
end


function q = bh_adjust(p)
    p = p(:)';  m = numel(p);
    q = nan(size(p));
    ok = isfinite(p);
    if ~any(ok), return; end
    pv = p(ok);  mm = numel(pv);
    [ps, o] = sort(pv);
    qs = ps .* mm ./ (1:mm);
    qs = min(1, fliplr(cummin(fliplr(qs))));
    qq = zeros(1, mm);  qq(o) = qs;
    q(ok) = qq;
    if m == 0, q = []; end
end


function s = stars_na(p)
    if isnan(p),      s = 'NA';
    elseif p < 0.001, s = '***';
    elseif p < 0.01,  s = '**';
    elseif p < 0.05,  s = '*';
    else,             s = 'ns';
    end
end


function m = median_or_nan(v)
    if isempty(v), m = NaN; else, m = median(v); end
end


%% ================================================================
%  LOCAL: style_panel — limits, ticks, dashed major + midpoint minor grid
%% ================================================================
function style_panel(ax, x_lo, x_hi, y_upper, y_step, fs_tick)
    xlim(ax, [x_lo, x_hi]);
    xticks(ax, x_lo:8:x_hi);
    ylim(ax, [0, y_upper]);
    yticks(ax, 0:y_step:y_upper);
    set(ax, 'XGrid','on','YGrid','on', 'XMinorGrid','on','YMinorGrid','on', ...
            'GridLineStyle','--', 'GridColor',[0.72 0.72 0.72],'GridAlpha',1.0, ...
            'MinorGridLineStyle','--', 'MinorGridColor',[0.87 0.87 0.87],'MinorGridAlpha',1.0);
    maj_x = xticks(ax);
    if numel(maj_x) >= 2
        ax.XAxis.MinorTickValues = (maj_x(1:end-1) + maj_x(2:end)) / 2;
    end
    maj_y = yticks(ax);
    if numel(maj_y) >= 2
        ax.YAxis.MinorTickValues = (maj_y(1:end-1) + maj_y(2:end)) / 2;
    end
    set(ax, 'FontSize',fs_tick, 'TickDir','in', 'Box','on', ...
            'XMinorTick','off', 'YMinorTick','off');
end


%% ================================================================
%  LOCAL: nice_y_limit — rounded upper limit and tick step for a y-axis
%  from 0 to at least y_max, with at most 6 tick intervals.
%% ================================================================
function [y_up, y_step] = nice_y_limit(y_max)
    if ~(y_max > 0), y_up = 1; y_step = 0.2; return; end
    mag   = 10^floor(log10(y_max));
    cands = [0.1 0.2 0.25 0.5 1 2 2.5 5 10] * mag;
    y_step = cands(end);
    for c = cands
        if ceil(y_max / c) <= 6, y_step = c; break; end
    end
    y_up = ceil(y_max / y_step) * y_step;
end


%% ================================================================
%  LOCAL: add_mirror_axes — ticks on the top + right frame
%% ================================================================
function add_mirror_axes(ax, fs_tick)
    ax2 = axes('Position', ax.Position, 'XAxisLocation','top', 'YAxisLocation','right', ...
               'Color','none', 'XColor',ax.XColor, 'YColor',ax.YColor, ...
               'FontSize',fs_tick, 'TickDir','in', 'Box','off', ...
               'XTick',ax.XTick, 'YTick',ax.YTick, 'XTickLabel',{}, 'YTickLabel',{}, ...
               'XLim',ax.XLim, 'YLim',ax.YLim);
    uistack(ax2, 'bottom');
    axes(ax);
end


%% ================================================================
%  LOCAL: save_figure — .fig / .pdf / .png
%% ================================================================
function save_figure(hfig, out_dir, stem, grp_str, ts, name, HT_FLOG)
    fig_path = fullfile(out_dir, sprintf('%s_%s_%s', stem, grp_str, ts));
    savefig(hfig, [fig_path '.fig']);
    exportgraphics(hfig, [fig_path '.pdf'], 'ContentType','vector');
    exportgraphics(hfig, [fig_path '.png'], 'Resolution',300);
    ht_fprintf(HT_FLOG, '\n%s saved:\n  %s.fig\n  %s.pdf\n  %s.png\n', ...
        name, fig_path, fig_path, fig_path);
    set(hfig, 'Visible', 'on');
    drawnow;
end


function s = sig_stars(p)
    if isnan(p),      s = 'ns';
    elseif p < 0.001, s = '***';
    elseif p < 0.01,  s = '**';
    elseif p < 0.05,  s = '*';
    else,             s = 'ns';
    end
end


function c = shade(c, f)
% f > 0 mixes colour c with white, f < 0 with black (|f| = fraction).
    if f >= 0, c = c + (1 - c) * f; else, c = c * (1 + f); end
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
%  p can underflow to exact 0.0 in double precision for large n — the
%  true value is just too small to represent — so that case is shown as
%  "<1e-300" rather than a misleading literal 0.
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
    ht_fprintf(HT_FLOG, '  plot_appearanceTime_vs_growthRate -- Group Selection\n');
    ht_fprintf(HT_FLOG, '============================================================\n');
    ht_fprintf(HT_FLOG, '  +-----+----------------+-------+----------+\n');
    ht_fprintf(HT_FLOG, '  | idx | label          | n_col | n_fitted |\n');
    ht_fprintf(HT_FLOG, '  +-----+----------------+-------+----------+\n');
    for k = 1:n_avail
        n_fit = sum(ht.groups(k).growth.fit_ok);
        ht_fprintf(HT_FLOG, '  | %3d | %-14s | %5d | %8d |\n', ...
            k, all_labels{k}, ht.groups(k).n_colonies, n_fit);
    end
    ht_fprintf(HT_FLOG, '  +-----+----------------+-------+----------+\n');
    ht_fprintf(HT_FLOG, '  Enter group labels (comma-separated), index numbers, or ALL:\n');

    while true
        raw = strtrim(input('  Your selection: ', 's'));
        if strcmpi(raw, 'all')
            sel_idx = 1:n_avail;
            ht_fprintf(HT_FLOG, '  -> All %d groups selected.\n', n_avail);
            return;
        end
        sel_idx = parse_group_tokens(raw, all_labels);
        if isempty(sel_idx)
            ht_fprintf(HT_FLOG, ['  WARNING: group(s) not recognised. Use labels (comma-separated if\n' ...
                                 '           they contain spaces, e.g. WT CTRL, 6S CTRL), index numbers or ALL.\n']);
            continue;
        end
        [~, ui]  = unique(sel_idx,'stable');
        sel_idx  = sel_idx(ui);
        ht_fprintf(HT_FLOG, '  -> Selected: %s\n', strjoin(all_labels(sel_idx), ', '));
        return;
    end
end


%% ====================================================================
%  LOCAL: ask_order — panel order of the selected groups. Numbers refer
%  to the list shown; groups not listed follow in their current order.
%  Enter keeps the current order.
%% ====================================================================
function sel = ask_order(HT_FLOG, sel, labels)
    n = numel(sel);
    if n < 2, return; end
    fprintf('\n  Display order — current order:\n');
    for k = 1:n, fprintf('  %3d | %s\n', k, labels{sel(k)}); end
    while true
        fprintf(['  New order as numbers from this list (e.g. 3 1 2; groups not listed follow\n' ...
                 '  in the current order) [Enter = keep]: ']);
        raw = strtrim(input('', 's'));
        if isempty(raw), break; end
        idx = parse_index_list(raw, n);
        if ~isempty(idx)
            sel = sel([idx, setdiff(1:n, idx, 'stable')]);
            break;
        end
        fprintf('  WARNING: enter numbers between 1 and %d. Try again.\n', n);
    end
    ht_fprintf(HT_FLOG, '  -> Order: %s\n', strjoin(labels(sel), ', '));
end


%% ====================================================================
%  LOCAL: ask_pairs — lists every pair of the selected groups (in display
%  order) and returns the chosen ones as [P x 2] indices into lbls.
%% ====================================================================
function P = ask_pairs(HT_FLOG, lbls)
    n = numel(lbls);
    P = zeros(0, 2);
    if n < 2
        ht_fprintf(HT_FLOG, '  (fewer than two groups: no pairwise comparisons)\n');
        return;
    end
    all_p = nchoosek(1:n, 2);
    pl = pair_labels(lbls, all_p);
    fprintf(['\n  Pairwise comparisons available (Mann–Whitney test of mu_max + comparison\n' ...
             '  of Spearman rho; each chosen pair also gets a panel in Figure 3):\n']);
    for r = 1:size(all_p, 1)
        fprintf('  %3d | %s\n', r, pl{r});
    end
    while true
        fprintf('  Comparisons — pair numbers (e.g. 1 4 6, ranges 1-3), ALL [Enter = none]: ');
        raw = strtrim(input('', 's'));
        if isempty(raw) || strcmpi(raw, 'none'), break; end
        if strcmpi(raw, 'all'), P = all_p; break; end
        idx = parse_index_list(raw, size(all_p, 1));
        if ~isempty(idx), P = all_p(sort(idx), :); break; end
        fprintf('  WARNING: enter pair numbers between 1 and %d. Try again.\n', size(all_p, 1));
    end
    if isempty(P)
        ht_fprintf(HT_FLOG, '  -> Comparisons: none\n');
    else
        ht_fprintf(HT_FLOG, '  -> Comparisons: %s\n', strjoin(pair_labels(lbls, P), ', '));
    end
end


function L = pair_labels(lbls, P)
    L = cell(1, size(P, 1));
    for r = 1:size(P, 1)
        L{r} = sprintf('%s vs %s', lbls{P(r, 1)}, lbls{P(r, 2)});
    end
end


%% ====================================================================
%  LOCAL: parse_index_list — "3 1 2", "1,4,6", "1-3 7" -> row vector of
%  unique integers in 1..nmax (input order kept); [] if anything invalid.
%% ====================================================================
function idx = parse_index_list(raw, nmax)
    idx = [];
    tok = regexp(strtrim(raw), '[,;\s]+', 'split');
    tok = tok(~cellfun(@isempty, tok));
    if isempty(tok), return; end
    out = [];
    for t = 1:numel(tok)
        rg = regexp(tok{t}, '^(\d+)-(\d+)$', 'tokens', 'once');
        if ~isempty(rg)
            a = str2double(rg{1});  b = str2double(rg{2});
            if a <= b, v = a:b; else, v = a:-1:b; end
        else
            v = str2double(tok{t});
            if ~isfinite(v) || v ~= floor(v), return; end
        end
        if any(v < 1 | v > nmax), return; end
        out = [out, v]; %#ok<AGROW>
    end
    [~, ui] = unique(out, 'stable');
    idx = out(sort(ui));
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


%% ====================================================================
%  LOCAL: parse_group_tokens
%  Group selection input: labels (comma-separated when labels contain
%  spaces, e.g. "WT CTRL, 6S CTRL"), labels separated by spaces, or
%  index numbers. Returns [] if anything is not recognised. The order of
%  the input is kept.
%% ====================================================================
function sel = parse_group_tokens(raw, all_labels)
    sel = [];
    norm = @(s) lower(regexprep(strtrim(char(s)), '\s+', ' '));
    lab  = cellfun(norm, all_labels, 'UniformOutput', false);
    tries = {regexp(raw, ',', 'split'), regexp(raw, '[,\s]+', 'split')};
    for tr = 1:numel(tries)
        tok = tries{tr};
        tok = tok(~cellfun(@(t) isempty(strtrim(t)), tok));
        s = zeros(1, numel(tok));
        ok = ~isempty(tok);
        for t = 1:numel(tok)
            m = find(strcmp(lab, norm(tok{t})), 1);
            if isempty(m)
                num = str2double(tok{t});
                if isfinite(num) && num == floor(num) && num >= 1 && num <= numel(all_labels)
                    m = num;
                end
            end
            if isempty(m), ok = false; break; end
            s(t) = m;
        end
        if ok
            [~, ui] = unique(s, 'stable');
            sel = s(sort(ui));
            return;
        end
    end
end

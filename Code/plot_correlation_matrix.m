%% HeteroTyper Pipeline for Bright Plates
% Produces two complementary statistical outputs:
%
%   (1) Pearson within-group correlation matrices
%       One heatmap per selected sample group showing how the selected
%       phenotypic features correlate with each other inside that group.
%
%   (2) Non-parametric between-group comparison (Kruskal-Wallis + Dunn)
%       For each selected feature, a Kruskal-Wallis test checks whether
%       any group differs. Dunn's post-hoc test with Bonferroni correction
%       then produces all pairwise adjusted p-values. Results are saved as:
%         • a p-value heatmap figure (.fig / .pdf / .png) per feature
%         • a single CSV table covering all features and all pairs
%
% Data are NOT assumed to be normally distributed — Kruskal-Wallis and
% Dunn's test are rank-based and valid for any continuous distribution.
%
% Requires preprocess_pipeline_data(data) to have been run first.
%
%% USAGE:
%   plot_correlation_matrix(ht);

function plot_correlation_matrix(ht)
    HT_FLOG = -1;

    % ------------------------------------------------------------------
    %  Open log file
    % ------------------------------------------------------------------
    timestamp_log = datestr(now, 'yyyymmdd_HHMMSS');
    log_path = fullfile(ht.params.out_dir, ...
                        sprintf('plot_correlation_matrix_%s.txt', timestamp_log));
    HT_FLOG = fopen(log_path, 'w');
    if HT_FLOG == -1
        warning('Could not open log file: %s', log_path);
        HT_FLOG = -1;
    end
    fprintf('Log file: %s\n', log_path);

    check_ht(ht);

    % ------------------------------------------------------------------
    %  Feature catalogue  (indices shared by both analysis sections)
    % ------------------------------------------------------------------
    all_feat_names  = {'Lag Time','Final Size','Area','Intensity', ...
                       'Mean Intensity','Int / Size','Perimeter', ...
                       'Circularity','Eccentricity','Solidity'};
    all_feat_fields = {'lag_time','size','area','intensity', ...
                       'mean_intensity','int_per_size','perimeter', ...
                       'circularity','eccentricity','solidity'};

    % ------------------------------------------------------------------
    %  User prompts  (shared for both outputs)
    % ------------------------------------------------------------------
    sel_groups = ask_group_selection(HT_FLOG, ht);
    sel_feats  = ask_feature_selection(HT_FLOG, all_feat_names);

    n_groups    = length(sel_groups);
    n_feats     = length(sel_feats);
    feat_names  = all_feat_names(sel_feats);
    feat_fields = all_feat_fields(sel_feats);

    ht_fprintf(HT_FLOG, '\n--- Analysis parameters ---\n');
    ht_fprintf(HT_FLOG, '  Groups   : %s\n', strjoin(ht.labels(sel_groups), ', '));
    ht_fprintf(HT_FLOG, '  Features : %s\n', strjoin(feat_names, ', '));
    ht_fprintf(HT_FLOG, '  Output   : %s\n', ht.params.out_dir);
    ht_fprintf(HT_FLOG, '---------------------------\n\n');

    % ==================================================================
    %  SECTION 1 — Within-group Pearson correlation matrices
    % ==================================================================
    ht_fprintf(HT_FLOG, '=== SECTION 1: Within-group Pearson correlation matrices ===\n\n');

    % Custom blue-white-red colormap
    blue   = [42  80  161] / 255;
    white  = [1   1   1  ];
    red    = [170 0   36 ] / 255;
    n_half = 128;
    custom_cmap = [ ...
        [linspace(blue(1),white(1),n_half)', linspace(blue(2),white(2),n_half)', linspace(blue(3),white(3),n_half)']; ...
        [linspace(white(1),red(1),n_half)',  linspace(white(2),red(2),n_half)',  linspace(white(3),red(3),n_half)']  ...
    ];

    n_cols    = min(n_groups, 2);
    n_rows    = ceil(n_groups / n_cols);
    cell_size = 0.50;
    panel_w   = n_feats * cell_size;
    panel_h   = n_feats * cell_size;
    gap_in    = 80 / 72;
    pad_in    = 60 / 72;
    fig_w     = 2*pad_in + n_cols*panel_w + (n_cols-1)*gap_in + 1.2;
    fig_h     = 2*pad_in + n_rows*panel_h + (n_rows-1)*gap_in + 0.6;

    hfig_corr = figure('Name','Correlation Matrices', ...
                       'Color','w', ...
                       'Units','inches', ...
                       'Position',[1 1 fig_w fig_h], ...
                       'PaperUnits','inches', ...
                       'PaperSize',[fig_w fig_h], ...
                       'PaperPosition',[0 0 fig_w fig_h]);

    tl = tiledlayout(hfig_corr, n_rows, n_cols, ...
                     'TileSpacing','loose', ...
                     'Padding','loose');

    fs_cell  = max(9,  round(11 - 0.3*n_feats));
    fs_label = max(11, round(13 - 0.2*n_feats));
    fs_title = 14;

    for k = 1:n_groups
        g   = sel_groups(k);
        grp = ht.groups(g);
        lbl = ht.labels{g};

        cols         = cellfun(@(f) grp.(f), feat_fields, 'UniformOutput', false);
        param_matrix = [cols{:}];
        param_matrix = param_matrix(all(isfinite(param_matrix), 2), :);

        nexttile(tl, k);

        if size(param_matrix, 1) < 3
            ht_fprintf(HT_FLOG, 'Group %s: insufficient data for correlation, skipping.\n', lbl);
            text(0.5, 0.5, sprintf('%s\n(no data)', lbl), ...
                 'Units','normalized', 'HorizontalAlignment','center');
            axis off;
            continue;
        end

        n_finite = size(param_matrix, 1);   % colonies with finite values on ALL selected features
        R  = corr(param_matrix, 'Type','Pearson');
        hm = heatmap(feat_names, feat_names, R, ...
                     'Colormap',        custom_cmap, ...
                     'ColorLimits',     [-1 1], ...
                     'CellLabelFormat', '%.2f', ...
                     'FontSize',        fs_cell);

        hm.Title    = sprintf('%s  -  Pearson Correlation', lbl);
        hm.FontSize = fs_label;
        try
            hm.NodeParent.Title.FontSize   = fs_title;
            hm.NodeParent.Title.FontWeight = 'bold';
        catch
        end

        % Report finite-row count (NaN-masked colonies excluded) and total
        % so the log makes clear why the two numbers may differ.
        ht_fprintf(HT_FLOG, 'Group %s: Pearson R matrix computed  (n=%d valid / %d total colonies).\n', ...
                   lbl, n_finite, grp.n_colonies);
    end

    % Save correlation matrix figure
    out_dir   = ht.params.out_dir;
    timestamp = datestr(now, 'yyyymmdd_HHMMSS');
    grp_str   = strjoin(ht.labels(sel_groups), '-');
    base_corr = sprintf('CorrelationMatrix_%s_%s', grp_str, timestamp);
    path_corr = fullfile(out_dir, base_corr);

    savefig(hfig_corr,  [path_corr '.fig']);
    exportgraphics(hfig_corr, [path_corr '.pdf'], 'ContentType','vector');
    exportgraphics(hfig_corr, [path_corr '.png'], 'Resolution',300);
    ht_fprintf(HT_FLOG, '\nCorrelation matrix figure saved:\n  %s.fig\n  %s.pdf\n  %s.png\n\n', ...
               path_corr, path_corr, path_corr);

    % ==================================================================
    %  SECTION 2 — Non-parametric between-group comparison
    %              Kruskal-Wallis omnibus  +  Dunn post-hoc (Bonferroni)
    % ==================================================================
    ht_fprintf(HT_FLOG, '=== SECTION 2: Non-parametric between-group comparison ===\n');
    ht_fprintf(HT_FLOG, '  Method  : Kruskal-Wallis omnibus + Dunn post-hoc\n');
    ht_fprintf(HT_FLOG, '  Correction: Bonferroni (per feature)\n\n');

    grp_labels = ht.labels(sel_groups);   % labels for selected groups only

    % Build all unique pairwise group combinations
    pair_idx  = nchoosek(1:n_groups, 2);   % [n_pairs x 2] indices into sel_groups
    n_pairs   = size(pair_idx, 1);
    pair_labels = cell(n_pairs, 1);
    for p_i = 1:n_pairs
        pair_labels{p_i} = sprintf('%s vs %s', ...
            grp_labels{pair_idx(p_i,1)}, grp_labels{pair_idx(p_i,2)});
    end

    % --- Storage for CSV output ---
    % Columns: Feature | KW_H | KW_p | Pair_1_padj | Pair_2_padj | ...
    csv_header = [{'Feature','KW_H_stat','KW_p_value'}, pair_labels'];
    csv_rows   = cell(n_feats, 3 + n_pairs);

    % Storage for combined figure (populated in stats loop, rendered after)
    pmat_all     = cell(n_feats, 1);
    log_pmat_all = cell(n_feats, 1);
    kw_H_all     = nan(n_feats, 1);
    kw_p_all     = nan(n_feats, 1);

    for fi = 1:n_feats
        fname  = feat_names{fi};
        ffield = feat_fields{fi};

        % Assemble per-group data vectors (finite values only)
        group_data = cell(n_groups, 1);
        for k = 1:n_groups
            v = ht.groups(sel_groups(k)).(ffield);
            group_data{k} = v(isfinite(v));
        end

        % --- Kruskal-Wallis omnibus test ---
        % Build concatenated data + group labels
        all_vals  = vertcat(group_data{:});
        all_grp   = arrayfun(@(k) repmat(k, numel(group_data{k}), 1), ...
                             (1:n_groups)', 'UniformOutput', false);
        all_grp   = vertcat(all_grp{:});

        if numel(all_vals) < 3 || n_groups < 2
            ht_fprintf(HT_FLOG, '[%s] Insufficient data — skipping.\n', fname);
            csv_rows(fi,:) = [{fname}, {NaN}, {NaN}, repmat({NaN}, 1, n_pairs)];
            continue;
        end

        % kruskalwallis returns: p-value, ANOVA table (cell array), stats struct.
        % The H-statistic (chi-squared approximation) lives in the table at
        % tbl{2,5} — it is NOT a field of the stats struct in any MATLAB version.
        [kw_p, kw_tbl, ~] = kruskalwallis(all_vals, all_grp, 'off');
        kw_H = kw_tbl{2,5};   % row 2 = "Groups" row, col 5 = Chi-sq / H value

        ht_fprintf(HT_FLOG, '[%s]  KW H=%.4f  p=%s\n', fname, kw_H, fmt_p(kw_p));

        % --- Dunn post-hoc with Bonferroni correction ---
        % Dunn's test: compare each pair using the pooled rank sum from KW.
        % The Bonferroni correction multiplies each raw p by the number of
        % pairs tested for this feature (m_pairs = n_groups*(n_groups-1)/2).
        % Adjusted p values are capped at 1.
        [all_ranks, ~] = tiedrank(all_vals);
        N = numel(all_vals);

        % Tie correction factor  C = 1 - sum(t^3 - t) / (N^3 - N)
        % where t is the size of each tied group.
        tie_correction = compute_tie_correction(all_ranks, N);

        padj_vals = nan(n_pairs, 1);
        for p_i = 1:n_pairs
            g1 = pair_idx(p_i, 1);
            g2 = pair_idx(p_i, 2);

            r1  = all_ranks(all_grp == g1);
            r2  = all_ranks(all_grp == g2);
            n1  = numel(r1);
            n2  = numel(r2);

            if n1 < 1 || n2 < 1
                continue;
            end

            % Dunn z-statistic
            mean_r  = (N + 1) / 2;
            se      = sqrt( tie_correction * N*(N+1)/12 * (1/n1 + 1/n2) );
            if se == 0
                padj_vals(p_i) = 1;
                continue;
            end
            z_stat  = (mean(r1) - mean(r2)) / se;
            raw_p   = 2 * (1 - normcdf(abs(z_stat)));

            % Bonferroni correction
            padj_vals(p_i) = min(1, raw_p * n_pairs);

            ht_fprintf(HT_FLOG, '    %s  z=%.3f  p_adj=%s%s\n', ...
                pair_labels{p_i}, z_stat, fmt_p(padj_vals(p_i)), ...
                sig_stars(padj_vals(p_i)));
        end

        % Store in CSV row
        csv_rows(fi,:) = [{fname}, {kw_H}, {kw_p}, num2cell(padj_vals')];

        % Build symmetric p-value matrix and -log10 transform — stored for
        % the combined figure drawn after the stats loop finishes.
        pmat_all{fi}     = nan(n_groups, n_groups);
        log_pmat_all{fi} = nan(n_groups, n_groups);
        kw_H_all(fi)     = kw_H;
        kw_p_all(fi)     = kw_p;
        for p_i = 1:n_pairs
            g1 = pair_idx(p_i,1);
            g2 = pair_idx(p_i,2);
            pmat_all{fi}(g1,g2) = padj_vals(p_i);
            pmat_all{fi}(g2,g1) = padj_vals(p_i);
        end
        % -log10 transform; clamp p away from 0 to avoid Inf
        pm_safe = max(pmat_all{fi}, 1e-300);
        lp      = -log10(pm_safe);
        lp(isnan(pmat_all{fi})) = NaN;           % restore diagonal as NaN
        log_pmat_all{fi} = lp;

        ht_fprintf(HT_FLOG, '\n');
    end  % end feature stats loop

    % ------------------------------------------------------------------
    %  Combined KW/Dunn figure — built with imagesc for reliable text.
    % ------------------------------------------------------------------

    clim_max = 4;   % -log10(0.0001); values above saturate to darkest red
    n_cm     = 256;
    kw_cmap  = [linspace(1, 0.6, n_cm)', ...
                linspace(1, 0.0, n_cm)', ...
                linspace(1, 0.0, n_cm)'];

    % Grid
    tl_cols = min(n_feats, 4);
    tl_rows = ceil(n_feats / tl_cols);

    % Panel geometry (all in inches) — kept tight
    cell_sz = 0.72;    % coloured cell size
    lbl_w   = 0.48;    % left margin  (y-axis labels)
    lbl_h   = 0.36;    % bottom margin (x-axis labels)
    ttl_h   = 0.58;    % top margin   (two-line title)
    gap_x   = 0.14;    % horizontal gap between panels (no per-panel colorbar)
    gap_y   = 0.08;    % vertical gap between panels
    pad     = 0.30;    % outer padding on left/bottom/top
    sup_h   = 0.44;    % super-title strip height
    leg_h   = 0.38;    % bottom legend strip height
    cb_w    = 0.16;    % single shared colorbar width
    cb_gap  = 0.22;    % gap between rightmost panel and colorbar
    cb_lbl  = 0.55;    % extra width for colorbar label text

    heat_w  = n_groups * cell_sz;
    heat_h  = n_groups * cell_sz;
    panel_w = lbl_w + heat_w;
    panel_h = ttl_h + heat_h + lbl_h;

    fig_w = pad + tl_cols*panel_w + (tl_cols-1)*gap_x + cb_gap + cb_w + cb_lbl + pad;
    fig_h = pad + sup_h + tl_rows*panel_h + (tl_rows-1)*gap_y + leg_h;

    hfig_kw = figure('Name','KW Dunn pairwise', ...
                     'Color','w', ...
                     'Units','inches', ...
                     'Position',[1 1 fig_w fig_h], ...
                     'PaperUnits','inches', ...
                     'PaperSize',[fig_w fig_h], ...
                     'PaperPosition',[0 0 fig_w fig_h]);

    colormap(hfig_kw, kw_cmap);

    % Super-title
    annotation(hfig_kw, 'textbox', ...
               [0, 1 - sup_h/fig_h, 1, sup_h/fig_h], ...
               'String', sprintf('Dunn p_{adj} (Bonferroni)  |  Groups: %s', ...
                                 strjoin(grp_labels, ', ')), ...
               'EdgeColor','none', 'FontSize',13, 'FontWeight','bold', ...
               'HorizontalAlignment','center', 'VerticalAlignment','middle', ...
               'FitBoxToText','off');

    fs_lbl   = 11;   % group name labels
    fs_stars = 15;   % stars / ns text inside cells
    fs_title =  9;   % per-panel title

    for fi = 1:n_feats
        row_i = ceil(fi / tl_cols);
        col_i = mod(fi-1, tl_cols) + 1;

        lp      = log_pmat_all{fi};
        pmat    = pmat_all{fi};
        kw_H    = kw_H_all(fi);
        kw_p    = kw_p_all(fi);
        lp_disp = min(lp, clim_max);

        lp_img             = lp_disp;
        lp_img(isnan(lp_img)) = -1;   % diagonal → below colour range → white

        % Cell-grid bottom-left in normalised figure coords
        x0 = (pad + (col_i-1)*(panel_w + gap_x) + lbl_w) / fig_w;
        y0 = (leg_h + (tl_rows - row_i)*(panel_h + gap_y) + lbl_h) / fig_h;
        w0 = heat_w / fig_w;
        h0 = heat_h / fig_h;

        ax = axes('Parent', hfig_kw, 'Units','normalized', ...
                  'Position', [x0, y0, w0, h0]);

        imagesc(ax, lp_img, [0, clim_max]);
        colormap(ax, kw_cmap);
        caxis(ax, [0, clim_max]);
        hold(ax, 'on');

        % Grey diagonal patches
        for d = 1:n_groups
            fill(ax, [d-0.5 d+0.5 d+0.5 d-0.5], ...
                     [d-0.5 d-0.5 d+0.5 d+0.5], ...
                 [0.88 0.88 0.88], 'EdgeColor','none');
        end

        % Grid lines
        for g = 0.5:1:n_groups+0.5
            plot(ax, [0.5 n_groups+0.5], [g g], 'k-', 'LineWidth',0.5);
            plot(ax, [g g], [0.5 n_groups+0.5], 'k-', 'LineWidth',0.5);
        end

        % Cell text: stars, 'ns', or '-' on diagonal
        for r = 1:n_groups
            for c = 1:n_groups
                if r == c
                    lbl = '-';
                    tc  = [0.5 0.5 0.5];
                else
                    stars = sig_stars(pmat(r,c));
                    if isempty(stars)
                        lbl = 'ns';
                    else
                        lbl = stars;
                    end
                    bg  = 0;
                    if ~isnan(lp_disp(r,c)) && lp_disp(r,c) > 0
                        bg = lp_disp(r,c) / clim_max;
                    end
                    tc = ternary_color(bg);
                end
                text(ax, c, r, lbl, ...
                     'HorizontalAlignment','center', ...
                     'VerticalAlignment',  'middle', ...
                     'FontSize',   fs_stars, ...
                     'FontWeight', 'bold', ...
                     'Color',      tc, ...
                     'Interpreter','none');
            end
        end

        set(ax, 'XLim',[0.5 n_groups+0.5], 'YLim',[0.5 n_groups+0.5], ...
                'YDir','reverse', 'XTick',[], 'YTick',[], ...
                'XColor','k', 'YColor','k', 'Box','on', 'LineWidth',0.8);

        % X-axis labels
        for c = 1:n_groups
            annotation(hfig_kw, 'textbox', ...
                [x0+(c-1)*cell_sz/fig_w, y0-lbl_h/fig_h, ...
                 cell_sz/fig_w, lbl_h/fig_h*0.92], ...
                'String',grp_labels{c}, 'EdgeColor','none', ...
                'FontSize',fs_lbl, 'FontWeight','bold', ...
                'HorizontalAlignment','center', 'VerticalAlignment','top', ...
                'FitBoxToText','off');
        end

        % Y-axis labels
        for r = 1:n_groups
            annotation(hfig_kw, 'textbox', ...
                [x0-lbl_w/fig_w, y0+(n_groups-r)*cell_sz/fig_h, ...
                 lbl_w/fig_w*0.95, cell_sz/fig_h], ...
                'String',grp_labels{r}, 'EdgeColor','none', ...
                'FontSize',fs_lbl, 'FontWeight','bold', ...
                'HorizontalAlignment','right', 'VerticalAlignment','middle', ...
                'FitBoxToText','off');
        end

        % Two-line panel title
        annotation(hfig_kw, 'textbox', ...
            [x0-lbl_w/fig_w, y0+h0, (heat_w+lbl_w)/fig_w, ttl_h/fig_h], ...
            'String', sprintf('%s\nKW H=%.1f, p=%s', feat_names{fi}, kw_H, fmt_p(kw_p)), ...
            'EdgeColor','none', 'FontSize',fs_title, 'FontWeight','bold', ...
            'HorizontalAlignment','center', 'VerticalAlignment','bottom', ...
            'FitBoxToText','off', 'Interpreter','none');
    end

    % ------------------------------------------------------------------
    %  Single shared colorbar — right edge of the figure, full height
    %  of the content area (top of last row to bottom of first row).
    % ------------------------------------------------------------------
    % Vertical extent: from bottom of bottom-row grid to top of top-row title
    cb_y_bot = (leg_h + lbl_h) / fig_h;
    cb_y_top = (leg_h + tl_rows*panel_h + (tl_rows-1)*gap_y - lbl_h*0.1) / fig_h;
    cb_x_pos = (pad + tl_cols*panel_w + (tl_cols-1)*gap_x + cb_gap) / fig_w;

    cb_ax = axes('Parent', hfig_kw, 'Units','normalized', ...
                 'Position', [cb_x_pos, cb_y_bot, cb_w/fig_w, cb_y_top - cb_y_bot]);
    colormap(cb_ax, kw_cmap);
    cb = colorbar(cb_ax, 'eastoutside');
    caxis(cb_ax, [0 clim_max]);
    cb.Ticks        = 0:1:clim_max;
    cb.TickLabels   = {'0','1','2','3','4 (≤0.0001)'};
    cb.FontSize     = 10;
    cb.Label.String = '-log_{10}(p_{adj})';
    cb.Label.FontSize   = 11;
    cb.Label.FontWeight = 'bold';
    set(cb_ax, 'Visible','off');

    % Bottom legend
    annotation(hfig_kw, 'textbox', [0.01 0.002 0.98 leg_h/fig_h*0.75], ...
               'String', ['ns = not significant   * p<0.05   ** p<0.01   *** p<0.001   ' ...
                          '(colour = -log_{10}(p_{adj}), saturates at 4)'], ...
               'EdgeColor','none', 'FontSize',10, 'FontWeight','bold', ...
               'HorizontalAlignment','center', 'VerticalAlignment','middle', ...
               'FitBoxToText','off');

    % Save combined figure
    base_kw  = sprintf('KWDunn_%s_%s', grp_str, timestamp);
    path_kw  = fullfile(out_dir, base_kw);
    savefig(hfig_kw,  [path_kw '.fig']);
    exportgraphics(hfig_kw, [path_kw '.pdf'], 'ContentType','vector');
    exportgraphics(hfig_kw, [path_kw '.png'], 'Resolution',300);
    ht_fprintf(HT_FLOG, '\nKW/Dunn combined figure saved:\n  %s.fig\n  %s.pdf\n  %s.png\n\n', ...
               path_kw, path_kw, path_kw);

    % ------------------------------------------------------------------
    %  Write combined CSV table
    % ------------------------------------------------------------------
    csv_path = fullfile(out_dir, ...
        sprintf('KWDunn_pairwise_%s_%s.csv', grp_str, timestamp));

    fcsv = fopen(csv_path, 'w');
    if fcsv == -1
        ht_fprintf(HT_FLOG, 'WARNING: Could not write CSV to %s\n', csv_path);
    else
        % Header row
        fprintf(fcsv, '%s\n', strjoin(csv_header, ','));
        % Data rows
        for fi = 1:n_feats
            row_parts = cell(1, 3 + n_pairs);
            row_parts{1} = csv_rows{fi,1};                           % Feature name
            row_parts{2} = sprintf('%.6g', csv_rows{fi,2});          % KW H
            row_parts{3} = fmt_p_csv(csv_rows{fi,3});                % KW p
            for p_i = 1:n_pairs
                val = csv_rows{fi, 3+p_i};
                if isnumeric(val) && ~isnan(val)
                    row_parts{3+p_i} = fmt_p_csv(val);
                else
                    row_parts{3+p_i} = 'NA';
                end
            end
            fprintf(fcsv, '%s\n', strjoin(row_parts, ','));
        end
        fclose(fcsv);
        ht_fprintf(HT_FLOG, 'Pairwise CSV saved: %s\n\n', csv_path);
    end

    ht_fprintf(HT_FLOG, '=== All outputs complete ===\n');
    ht_fprintf(HT_FLOG, '  Section 1 (Pearson matrices)     : %s  (.fig/.pdf/.png)\n', path_corr);
    ht_fprintf(HT_FLOG, '  Section 2 (KW + Dunn combined)   : %s  (.fig/.pdf/.png)\n', path_kw);
    ht_fprintf(HT_FLOG, '  Section 2 CSV summary             : %s\n', csv_path);

    % Close log
    if HT_FLOG ~= -1
        fclose(HT_FLOG);
        fprintf('Log saved: %s\n', log_path);
    end

end  % end main function


%% ================================================================
%  LOCAL FUNCTION: compute_tie_correction
%  Returns the Dunn/KW tie-correction factor:
%    C = 1 - sum(t_k^3 - t_k) / (N^3 - N)
%  where t_k is the size of the k-th group of tied ranks.
% ================================================================
function C = compute_tie_correction(ranks, N)
    % Count runs of equal values (ties)
    sorted_r = sort(ranks);
    tie_sum  = 0;
    i = 1;
    while i <= N
        j = i;
        while j < N && sorted_r(j+1) == sorted_r(j)
            j = j + 1;
        end
        t = j - i + 1;
        if t > 1
            tie_sum = tie_sum + (t^3 - t);
        end
        i = j + 1;
    end
    denom = N^3 - N;
    if denom == 0
        C = 1;
    else
        C = 1 - tie_sum / denom;
    end
end


%% ================================================================
%  LOCAL FUNCTION: fmt_p
%  Format p-value for display (e.g. '0.0032' or '<0.0001').
% ================================================================
function s = fmt_p(p)
    if isnan(p)
        s = 'NA';
    elseif p == 0
        % kruskalwallis/normcdf-based p can legitimately underflow to
        % exact 0.0 in double precision (huge n / huge effect) — the
        % true value is just too small to represent, not actually zero.
        s = '<1e-300';
    elseif p < 0.0001
        s = '<0.0001';
    else
        s = sprintf('%.4f', p);
    end
end


%% ================================================================
%  LOCAL FUNCTION: fmt_p_csv
%  Full-precision p-value text for the CSV export (%.6g instead of
%  fmt_p's 4-decimal/bucketed display). Exact 0.0 (double underflow)
%  is shown as "<1e-300" rather than a misleading literal 0.
% ================================================================
function s = fmt_p_csv(p)
    if isnan(p)
        s = 'NA';
    elseif p == 0
        s = '<1e-300';
    else
        s = sprintf('%.6g', p);
    end
end


%% ================================================================
%  LOCAL FUNCTION: sig_stars
%  Return significance star string for a p-value.
% ================================================================
function s = sig_stars(p)
    if isnan(p) || p >= 0.05
        s = '';
    elseif p < 0.001
        s = '***';
    elseif p < 0.01
        s = '**';
    else
        s = '*';
    end
end


%% ================================================================
%  LOCAL FUNCTION: ternary_color
%  Choose black or white text depending on background brightness.
%  bg is a scalar in [0,1] representing colormap position.
% ================================================================
function c = ternary_color(bg)
    if bg > 0.55
        c = [1 1 1];   % white text on dark background
    else
        c = [0 0 0];   % black text on light background
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
        tokens  = regexp(strtrim(raw), '[,\s]+', 'split');
        tokens  = tokens(~cellfun(@isempty, tokens));
        sel_idx = zeros(1, length(tokens));
        ok      = true;
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
%  LOCAL FUNCTION: ask_feature_selection
% ================================================================
function sel_idx = ask_feature_selection(HT_FLOG, feat_names)
    n = length(feat_names);
    ht_fprintf(HT_FLOG, '\nAvailable phenotypic features:\n');
    ht_fprintf(HT_FLOG, '  +-----+---------------------+\n');
    for k = 1:n
        ht_fprintf(HT_FLOG, '  | %3d | %-19s |\n', k, feat_names{k});
    end
    ht_fprintf(HT_FLOG, '  +-----+---------------------+\n');
    ht_fprintf(HT_FLOG, '  Enter feature numbers separated by commas  (e.g.  1,2,4,7)\n');
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
        ht_fprintf(HT_FLOG, '  -> Selected features: %s\n', strjoin(feat_names(sel_idx), ', '));
        return;
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
%  Writes to terminal and to the log file simultaneously.
% ================================================================
function ht_fprintf(HT_FLOG, varargin)
    fprintf(varargin{:});
    if ~isempty(HT_FLOG) && HT_FLOG ~= -1
        fprintf(HT_FLOG, varargin{:});
    end
end
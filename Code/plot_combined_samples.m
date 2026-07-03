%% HeteroTyper Pipeline for Bright Plates
% Plots per-group probability density histograms for user-selected
% phenotypic features, with Gini index and smart corner labels.
%
% Requires preprocess_pipeline_data(data) to have been run first.
%
%% USAGE:
%   plot_combined_samples(ht);

function plot_combined_samples(ht)

    % ------------------------------------------------------------------
    %  Validate input
    % ------------------------------------------------------------------
    if ~isstruct(ht) || ~isfield(ht, 'groups')
        error(['Input must be the ht struct produced by preprocess_pipeline_data.\n' ...
               'Run:  preprocess_pipeline_data(data)  first.']);
    end

    % ------------------------------------------------------------------
    %  Unpack everything already computed by the preprocessor
    % ------------------------------------------------------------------
    p      = ht.params;
    labels = ht.labels;
    colors = ht.colors;
    xb     = ht.xbins;

    incTime = p.incTime;
    max_lag = p.max_lag;
    out_dir = p.out_dir;

    % ------------------------------------------------------------------
    %  Open log file — mirrors everything printed to the terminal
    % ------------------------------------------------------------------
    timestamp_log = datestr(now, 'yyyymmdd_HHMMSS');
    log_path = fullfile(out_dir, ...
                        sprintf('plot_combined_samples_%s.txt', timestamp_log));
    flog = fopen(log_path, 'w');
    if flog == -1
        warning('Could not open log file: %s', log_path);
        flog = -1;
    end
    tlog('Log file: %s\n', log_path);

    function tlog(varargin)
        % Write to terminal AND log file simultaneously
        fprintf(varargin{:});
        if flog ~= -1
            fprintf(flog, varargin{:});
        end
    end
    % ------------------------------------------------------------------
    %  User inputs: which groups and which features to plot
    % ------------------------------------------------------------------
    [sel_groups, plot_feat_ix] = ask_selections(ht);

    tlog('\n--- Plot parameters ---\n');
    tlog('  Groups to plot    : %s\n',  strjoin(labels(sel_groups), ', '));
    tlog('  Features to plot  : [%s]\n', num2str(plot_feat_ix));
    tlog('  Output directory  : %s\n',   out_dir);
    tlog('-----------------------\n\n');

    % ------------------------------------------------------------------
    %  Figure layout — manual position grid so labels never clip
    % ------------------------------------------------------------------
    n_rows = length(sel_groups);
    n_cols = length(plot_feat_ix);

    panel_w  = 3.5;   panel_h  = 2.2;
    gap_x    = 0.55;  gap_y    = 0.65;
    margin_l = 1.00;  margin_r = 0.30;
    margin_b = 0.80;  margin_t = 0.30;
    ylbl_w   = 0.45;  % left space per inter-column gap for y-ticks + ylabel

    fig_w = margin_l + n_cols*panel_w + (n_cols-1)*(gap_x+ylbl_w) + margin_r;
    fig_h = margin_b + n_rows*panel_h + (n_rows-1)*gap_y + margin_t;

    hfig = figure('Name','Combined plots', ...
                  'Units','inches', ...
                  'Position',[1 1 fig_w fig_h], ...
                  'PaperUnits','inches', ...
                  'PaperSize',[fig_w fig_h], ...
                  'PaperPosition',[0 0 fig_w fig_h], ...
                  'Visible','off');

    ax_pos = zeros(n_rows, n_cols, 4);
    for r = 1:n_rows
        for c = 1:n_cols
            ax_pos(r,c,:) = [ ...
                (margin_l + (c-1)*(panel_w+gap_x+ylbl_w)) / fig_w, ...
                (margin_b + (n_rows-r)*(panel_h+gap_y)) / fig_h, ...
                panel_w / fig_w, ...
                panel_h / fig_h];
        end
    end

    fs_tick  = 12;
    fs_label = 13;
    fs_annot = 12;
    fs_gini  = 11;

    % ------------------------------------------------------------------
    %  Build x bin midpoint vectors (used for bar centres)
    % ------------------------------------------------------------------
    lag_step    = xb.lag_step;
    xp_lag      = xb.lag(1:end-1)         + lag_step/2;
    xp_size     = xb.size(1:end-1)        + xb.size_step/2;
    xp_area     = xb.area(1:end-1)        + xb.area_step/2;
    xp_peri     = xb.perimeter(1:end-1)   + xb.peri_step/2;
    xp_circ     = xb.circularity(1:end-1) + 0.005;
    xp_ecc      = xb.eccentricity(1:end-1)+ 0.005;
    xp_sol      = xb.solidity(1:end-1)    + 0.005;

    % Intensity, Mean Intensity, Int/Size: recompute with finer bins (~50)
    % rather than using the coarse precomputed vectors stored in ht.xbins.
    n_bins_int = 100;
    int_step   = round_to_sig(ht.global.int      / n_bins_int, 1);
    mint_step  = round_to_sig(ht.global.mean_int / n_bins_int, 2);
    ips_step   = round_to_sig(ht.global.int_size / n_bins_int, 1);

    int_upper  = ceil(ht.global.int      / int_step)  * int_step  + int_step;
    mint_upper = ceil(ht.global.mean_int / mint_step) * mint_step + mint_step;
    ips_upper  = ceil(ht.global.int_size / ips_step)  * ips_step  + ips_step;

    int_bins   = 0 : int_step  : int_upper;
    mint_bins  = 0 : mint_step : mint_upper;
    ips_bins   = 0 : ips_step  : ips_upper;

    xp_int     = int_bins(1:end-1)  + int_step/2;
    xp_mint    = mint_bins(1:end-1) + mint_step/2;
    xp_ips     = ips_bins(1:end-1)  + ips_step/2;

    tlog('Intensity bins (%d target): int_step=%.3g  mint_step=%.3g  ips_step=%.3g\n', ...
            n_bins_int, int_step, mint_step, ips_step);

    % ------------------------------------------------------------------
    %  Main plot loop — one row per group, one column per selected feature
    % ------------------------------------------------------------------
    for i = 1:n_rows
        g   = sel_groups(i);
        grp = ht.groups(g);

        % Print terminal range summary (median, IQR, Gini per feature)
        tlog('Group %s  (n=%d colonies):\n', labels{g}, grp.n_colonies);
        fns = {'lag_time','size','area','intensity','mean_intensity', ...
               'int_per_size','perimeter','circularity','eccentricity','solidity'};
        flab = {'Lag Time','Size','Area','Intensity','Mean Intensity', ...
                'Int/Size','Perimeter','Circularity','Eccentricity','Solidity'};
        tlog('  %-16s  %10s  %10s  %10s  %10s  %10s  %8s\n', ...
            'Feature', 'Min', 'Max', 'Median', 'IQR', 'Gini', 'n');
        tlog('  %s\n', repmat('-', 1, 82));
        for f = 1:length(fns)
            v = grp.(fns{f});
            if isempty(v), continue; end
            v_clean = v(~isnan(v));
            if isempty(v_clean), continue; end
            med_v  = median(v_clean);
            iqr_v  = iqr(v_clean);
            gini_v = grp.gini.(fns{f});
            tlog('  %-16s  %10.3g  %10.3g  %10.3g  %10.3g  %8.4f  %8d\n', ...
                flab{f}, min(v_clean), max(v_clean), med_v, iqr_v, gini_v, length(v_clean));
        end
        tlog('\n');

        % Full catalogue of all 10 panels
        % Columns: {norm_counts, xplot, xlims, xlabel_str, data_vec}
        % hist.* counts and xp_* midpoints are both derived from ht.xbins,
        % so they are guaranteed to be the same length.
        %
        % Intensity features use locally-recomputed bins (finer resolution);
        % their hist.* equivalents in ht use the same formula so lengths match.
        norm_int  = grp.hist.intensity;
        norm_mint = grp.hist.mean_intensity;
        norm_ips  = grp.hist.int_per_size;

        all_panels = { ...
            grp.hist.lag_time,     xp_lag,  [incTime, max_lag+lag_step], 'Lag time (h)',      grp.lag_time;       ...
            grp.hist.size,         xp_size, [0, xb.size_upper],          'Colony size (px)',  grp.size;           ...
            grp.hist.area,         xp_area, [0, xb.area_upper],          'Area (px)',         grp.area;           ...
            norm_int,              xp_int,  [0, int_upper],              'Intensity',         grp.intensity;      ...
            norm_mint,             xp_mint, [0, mint_upper],             'Mean intensity',    grp.mean_intensity; ...
            norm_ips,              xp_ips,  [0, ips_upper],              'Intensity / size',  grp.int_per_size;   ...
            grp.hist.perimeter,    xp_peri, [0, xb.peri_upper],          'Perimeter (px)',    grp.perimeter;      ...
            grp.hist.circularity,  xp_circ, [0, 1],                      'Circularity',       grp.circularity;    ...
            grp.hist.eccentricity, xp_ecc,  [0, 1],                      'Eccentricity',      grp.eccentricity;   ...
            grp.hist.solidity,     xp_sol,  [0, 1],                      'Solidity',          grp.solidity;       ...
        };

        for col_pos = 1:n_cols
            feat_id  = plot_feat_ix(col_pos);

            norm_p   = all_panels{feat_id, 1};
            xp       = all_panels{feat_id, 2};
            xlims    = all_panels{feat_id, 3};
            xlbl     = all_panels{feat_id, 4};
            data_vec = all_panels{feat_id, 5};

            ax = axes('Position', squeeze(ax_pos(i, col_pos, :))');

            % Bar + outline
            bar(xp, norm_p, 'FaceColor',colors(g,:),'EdgeColor','none','FaceAlpha',1);
            hold on;
            plot(xp, norm_p, 'Color','k','LineWidth',0.3);

            y_max = max(norm_p);
            if y_max == 0, y_max = 1; end
            axis([xlims(1) xlims(2) 0 y_max]);

            % Lag time: decade x-ticks
            if feat_id == 1
                xticks(incTime:10:max_lag);
            end

            % Median line (omit NaNs from censored/rejected colonies)
            med_val = median(data_vec(~isnan(data_vec)));
            if ~isempty(med_val) && ~isnan(med_val)
                line([med_val med_val], [0 y_max], ...
                     'Color','r','LineWidth',1.5);
            end

            % Gini index (pre-computed in ht)
            gini_val = grp.gini.(fns{feat_id});

            % Labels always top-left
            x_txt = xlims(1) + 0.02*(xlims(2)-xlims(1));

            text(x_txt, y_max*0.97, labels{g}, ...
                 'FontSize',fs_annot, 'FontWeight','bold', 'Color','k', ...
                 'HorizontalAlignment','left', 'VerticalAlignment','top', ...
                 'Clipping','off');
            text(x_txt, y_max*0.87, sprintf('%.3f', gini_val), ...
                 'FontSize',fs_gini, 'Color',colors(g,:), ...
                 'HorizontalAlignment','left', 'VerticalAlignment','top', ...
                 'Clipping','off');

            % x-label: bottom row only
            if i == n_rows
                xlabel(xlbl, 'FontSize',fs_label, 'FontWeight','bold');
            end

            % y-label: all columns
            ylabel('Probability density', 'FontSize',fs_label, 'FontWeight','bold');

            set(ax, 'FontSize',fs_tick, 'TickDir','in', 'Box','on', ...
                'XMinorTick','off', 'YMinorTick','off');
            % Mirror ticks on top and right without extra tick labels
            set(ax, 'XAxisLocation','bottom', 'YAxisLocation','left');
            % Draw mirrored tick marks on top and right axes
            ax2 = axes('Position', ax.Position, ...
                       'XAxisLocation','top', 'YAxisLocation','right', ...
                       'Color','none', 'XColor',ax.XColor, 'YColor',ax.YColor, ...
                       'FontSize',fs_tick, 'TickDir','in', 'Box','off', ...
                       'XTick',ax.XTick, 'YTick',ax.YTick, ...
                       'XTickLabel',{}, 'YTickLabel',{}, ...
                       'XLim',ax.XLim, 'YLim',ax.YLim);
            uistack(ax2, 'bottom');
            axes(ax);
        end
    end

    % Render the complete figure all at once now that all panels are drawn
    set(hfig, 'Visible', 'on');
    drawnow;

    % ------------------------------------------------------------------
    %  Biological Replicate Correlation Tables
    %  — Print to terminal and export as Excel workbook
    % ------------------------------------------------------------------
    if isfield(ht, 'biorep_corr') && ~isempty(ht.biorep_corr)
        % Check if STEP 6f was skipped entirely (no bio_rep_col was set)
        if isfield(ht.biorep_corr(1), 'skipped') && ht.biorep_corr(1).skipped
            tlog('\n============================================================\n');
            tlog('  BIOLOGICAL REPLICATE STATISTICS — SKIPPED\n');
            tlog('============================================================\n');
            tlog('  Reason: %s\n', ht.biorep_corr(1).skip_reason);
            tlog('\n  To generate biological replicate statistics:\n');
            tlog('  1. Re-run preprocess_pipeline_data(data)\n');
            tlog('  2. At parameter step 6, select the metadata column that\n');
            tlog('     identifies each biological replicate (e.g. "Set",\n');
            tlog('     "Strain", "BioRep") — do NOT press 0.\n');
            tlog('  3. Re-run plot_combined_samples(ht)\n');
            tlog('============================================================\n\n');
        else
        tlog('\n============================================================\n');
        tlog('  BIOLOGICAL REPLICATE CORRELATION & STATISTICAL ANALYSIS\n');
        tlog('============================================================\n');
        tlog('  Tests computed per replicate pair per morphology parameter:\n');
        tlog('    Pearson r     : Linear correlation of 20-quantile profiles.\n');
        tlog('                    r ~ 1 = strong linear agreement in shape.\n');
        tlog('    Pearson p     : Two-tailed p-value for Pearson r (df=18).\n');
        tlog('    Spearman rho  : Rank correlation of 20-quantile profiles.\n');
        tlog('                    rho ~ 1 = strong rank agreement.\n');
        tlog('    Spearman p    : Two-tailed p-value for Spearman rho (df=18).\n');
        tlog('    KS D          : Kolmogorov-Smirnov max-CDF distance [0-1].\n');
        tlog('                    D ~ 0 = near-identical distributions.\n');
        tlog('    KS p-value    : H0 = same distribution. p < 0.05 = sig. diff.\n');
        tlog('    JSD           : Jensen-Shannon Divergence [0-1].\n');
        tlog('                    JSD ~ 0 = same histogram shape.\n');
        tlog('    |dGini|       : |Gini_A - Gini_B| on raw colony values.\n');
        tlog('                    ~ 0 = replicates have similar heterogeneity.\n');
        tlog('    |Delta Median|: Absolute difference in medians.\n');
        tlog('    CV of medians : CV%% of per-replicate medians (group-level).\n');
        tlog('  Per-replicate quantities (not pairwise):\n');
        tlog('    Gini (values) : Classic Gini on raw colony values [0-1].\n');
        tlog('============================================================\n\n');

        % Prepare Excel output path
        timestamp_br = datestr(now, 'yyyymmdd_HHMMSS');
        xl_path = fullfile(out_dir, ...
            sprintf('biorep_correlation_stats_%s.xlsx', timestamp_br));

        % Track whether at least one group was written to Excel
        xl_written = false;

        for g = 1:length(ht.biorep_corr)
            bc = ht.biorep_corr(g);
            n_pairs = size(bc.spearman_rho, 1);
            n_feats = length(bc.features);

            tlog('--- Group: %s  (%d replicate(s), %d pair(s)) ---\n', ...
                bc.group_label, bc.n_reps, n_pairs);

            % Print bio-rep labels with plate counts
            rep_info = cell(bc.n_reps, 1);
            for r_idx = 1:bc.n_reps
                if isfield(bc, 'rep_n_plates') && ~isempty(bc.rep_n_plates)
                    rep_info{r_idx} = sprintf('%s (%d plate(s))', ...
                        bc.plate_labels{r_idx}, bc.rep_n_plates(r_idx));
                else
                    rep_info{r_idx} = bc.plate_labels{r_idx};
                end
            end
            tlog('  Bio-replicates: %s\n', strjoin(rep_info, ' | '));

            if n_pairs == 0
                tlog('  (Only 1 biological replicate — no pairwise comparison available.)\n');
                % Still print per-replicate Gini for the single replicate
                if isfield(bc, 'gini_values') && ~isempty(bc.gini_values)
                    tlog('  Per-replicate Gini (classic, raw colony values):\n');
                    tlog('  %-28s  %-20s  %12s\n', 'Replicate', 'Feature', 'Gini');
                    tlog('  %s\n', repmat('-', 1, 64));
                    for r_idx = 1:bc.n_reps
                        for fi = 1:n_feats
                            gv = bc.gini_values(r_idx, fi);
                            if ~isnan(gv)
                                tlog('  %-28s  %-20s  %12.4f\n', ...
                                    bc.plate_labels{r_idx}, bc.features{fi}, gv);
                            end
                        end
                    end
                end
                tlog('\n');
                continue;
            end

            % ---- Per-replicate Gini table ------------------------------
            if isfield(bc, 'gini_values') && ~isempty(bc.gini_values)
                tlog('  Per-bio-rep Gini (classic Gini on raw colony values):\n');
                tlog('  %-28s  %-20s  %12s\n', 'Replicate', 'Feature', 'Gini');
                tlog('  %s\n', repmat('-', 1, 64));
                for r_idx = 1:bc.n_reps
                    for fi = 1:n_feats
                        gv = bc.gini_values(r_idx, fi);
                        if ~isnan(gv)
                            tlog('  %-28s  %-20s  %12.4f\n', ...
                                bc.plate_labels{r_idx}, bc.features{fi}, gv);
                        end
                    end
                end
                tlog('\n');
            end

            % ---- Pairwise statistics table -----------------------------
            hdr = sprintf('  %-26s  %-18s  %8s %9s  %8s %9s  %6s %9s  %6s  %8s  %11s\n', ...
                'Pair', 'Feature', ...
                'Pearson_r', 'p', ...
                'Spear_rho', 'p', ...
                'KS_D', 'KS_p', 'JSD', '|dGini|', '|dMedian|');
            tlog('%s', hdr);
            tlog('  %s\n', repmat('-', 1, 116));

            for pi_p = 1:n_pairs
                for fi = 1:n_feats
                    pr_v   = bc.pearson_r(pi_p, fi);
                    pr_pv  = bc.pearson_p(pi_p, fi);
                    rho_v  = bc.spearman_rho(pi_p, fi);
                    sp_pv  = bc.spearman_p(pi_p, fi);
                    ks_d   = bc.ks_D(pi_p, fi);
                    ks_pv  = bc.ks_p(pi_p, fi);
                    jsd_v  = bc.jsd(pi_p, fi);
                    dg_v   = bc.delta_gini(pi_p, fi);
                    dmed   = bc.median_diff(pi_p, fi);

                    if isnan(pr_v), continue; end

                    % Significance star helper (shared for KS p) — a p that
                    % underflowed to exact 0.0 is still highly significant.
                    if ~isnan(ks_pv) && ks_pv < 0.001,      ks_sig = '***';
                    elseif ~isnan(ks_pv) && ks_pv < 0.01,   ks_sig = '** ';
                    elseif ~isnan(ks_pv) && ks_pv < 0.05,   ks_sig = '*  ';
                    else,                                     ks_sig = 'ns ';
                    end

                    tlog('  %-26s  %-18s  %+8.4f %9s  %+8.4f %9s  %6.4f %9s%s  %6.4f  %8.4f  %11.4g\n', ...
                        bc.pair_labels{pi_p}, bc.features{fi}, ...
                        pr_v, format_pval(pr_pv), rho_v, format_pval(sp_pv), ...
                        ks_d, format_pval(ks_pv), ks_sig, jsd_v, dg_v, dmed);
                end
                tlog('\n');
            end

            % CV of medians summary
            tlog('  CV of medians across all replicates:\n');
            for fi = 1:n_feats
                cv = bc.cv_medians(fi);
                if ~isnan(cv)
                    if cv <= 5
                        interp = 'excellent reproducibility';
                    elseif cv <= 15
                        interp = 'good reproducibility';
                    elseif cv <= 30
                        interp = 'moderate reproducibility';
                    else
                        interp = 'high variability';
                    end
                    tlog('    %-20s  CV = %6.2f%%  (%s)\n', ...
                        bc.features{fi}, cv, interp);
                end
            end
            tlog('\n');

            % ---- Excel sheet — pairwise stats for this group -----------
            sheet_name = regexprep(bc.group_label, '[\\\/\?\*\[\]:''\ ]', '_');
            sheet_name = sheet_name(1:min(31, length(sheet_name)));

            % Pairwise rows: one row per pair × feature
            % NOTE: row_pr_p / row_rho_p / row_ksp are stored as formatted
            % strings (not raw doubles) so a p-value that underflowed to
            % exact 0.0 in double precision is exported as "<1e-300"
            % rather than a misleading literal 0.
            row_pair   = {};  row_repA  = {};  row_repB  = {};  row_feat  = {};
            row_pr     = [];  row_pr_p  = {};
            row_rho    = [];  row_rho_p = {};
            row_ksD    = [];  row_ksp   = {};  row_kssig = {};
            row_jsd    = [];  row_dg    = [];  row_dmed  = [];

            for pi_p = 1:n_pairs
                tok = strsplit(bc.pair_labels{pi_p}, ' vs ');
                repA_str = tok{1};
                repB_str = ''; if length(tok) >= 2, repB_str = tok{2}; end

                for fi = 1:n_feats
                    pr_v   = bc.pearson_r(pi_p, fi);
                    pr_pv  = bc.pearson_p(pi_p, fi);
                    rho_v  = bc.spearman_rho(pi_p, fi);
                    sp_pv  = bc.spearman_p(pi_p, fi);
                    ks_d   = bc.ks_D(pi_p, fi);
                    ks_pv  = bc.ks_p(pi_p, fi);
                    jsd_v  = bc.jsd(pi_p, fi);
                    dg_v   = bc.delta_gini(pi_p, fi);
                    dmed   = bc.median_diff(pi_p, fi);

                    if ~isnan(ks_pv) && ks_pv < 0.001,      sig_xl = '***';
                    elseif ~isnan(ks_pv) && ks_pv < 0.01,   sig_xl = '**';
                    elseif ~isnan(ks_pv) && ks_pv < 0.05,   sig_xl = '*';
                    else,                                     sig_xl = 'ns';
                    end

                    row_pair{end+1,1}  = bc.pair_labels{pi_p};      %#ok<AGROW>
                    row_repA{end+1,1}  = repA_str;                   %#ok<AGROW>
                    row_repB{end+1,1}  = repB_str;                   %#ok<AGROW>
                    row_feat{end+1,1}  = bc.features{fi};            %#ok<AGROW>
                    row_pr(end+1,1)    = pr_v;                       %#ok<AGROW>
                    row_pr_p{end+1,1}  = format_pval(pr_pv);         %#ok<AGROW>
                    row_rho(end+1,1)   = rho_v;                      %#ok<AGROW>
                    row_rho_p{end+1,1} = format_pval(sp_pv);         %#ok<AGROW>
                    row_ksD(end+1,1)   = ks_d;                       %#ok<AGROW>
                    row_ksp{end+1,1}   = format_pval(ks_pv);         %#ok<AGROW>
                    row_kssig{end+1,1} = sig_xl;                     %#ok<AGROW>
                    row_jsd(end+1,1)   = jsd_v;                      %#ok<AGROW>
                    row_dg(end+1,1)    = dg_v;                       %#ok<AGROW>
                    row_dmed(end+1,1)  = dmed;                       %#ok<AGROW>
                end
            end

            % Append CV-of-medians rows at the bottom
            for fi = 1:n_feats
                row_pair{end+1,1}  = 'CV of medians (%)';
                row_repA{end+1,1}  = 'all replicates';
                row_repB{end+1,1}  = '';
                row_feat{end+1,1}  = bc.features{fi};
                row_pr(end+1,1)    = NaN;   row_pr_p{end+1,1}  = 'NA';
                row_rho(end+1,1)   = NaN;   row_rho_p{end+1,1} = 'NA';
                row_ksD(end+1,1)   = NaN;   row_ksp{end+1,1}   = 'NA';
                row_kssig{end+1,1} = '';
                row_jsd(end+1,1)   = NaN;   row_dg(end+1,1)    = NaN;
                row_dmed(end+1,1)  = bc.cv_medians(fi);
            end

            if ~isempty(row_feat)
                T_br = table(row_pair, row_repA, row_repB, row_feat, ...
                             row_pr, row_pr_p, row_rho, row_rho_p, ...
                             row_ksD, row_ksp, row_kssig, ...
                             row_jsd, row_dg, row_dmed, ...
                    'VariableNames', { ...
                        'Pair','RepA','RepB','Feature', ...
                        'Pearson_r','Pearson_p', ...
                        'Spearman_rho','Spearman_p', ...
                        'KS_D','KS_p','KS_sig', ...
                        'JSD','AbsDeltaGini', ...
                        'AbsDeltaMedian_or_CV_pct'});
                try
                    writetable(T_br, xl_path, 'Sheet', sheet_name);
                    xl_written = true;
                catch ME_xl
                    tlog('  WARNING: could not write Excel sheet "%s": %s\n', ...
                        sheet_name, ME_xl.message);
                end
            end

            % ---- Excel sheet — per-replicate Gini for this group -------
            if isfield(bc, 'gini_values') && ~isempty(bc.gini_values)
                gini_sheet = ['Gini_' sheet_name];
                gini_sheet = gini_sheet(1:min(31, length(gini_sheet)));

                g_rep  = {};  g_feat = {};  g_val = [];

                for r_idx = 1:bc.n_reps
                    for fi = 1:n_feats
                        gv = bc.gini_values(r_idx, fi);
                        g_rep{end+1,1}  = bc.plate_labels{r_idx};  %#ok<AGROW>
                        g_feat{end+1,1} = bc.features{fi};          %#ok<AGROW>
                        g_val(end+1,1)  = gv;                        %#ok<AGROW>
                    end
                end

                T_gini = table(g_rep, g_feat, g_val, ...
                    'VariableNames', {'Replicate','Feature','Gini'});
                try
                    writetable(T_gini, xl_path, 'Sheet', gini_sheet);
                    xl_written = true;
                catch ME_gini
                    tlog('  WARNING: could not write Gini sheet "%s": %s\n', ...
                        gini_sheet, ME_gini.message);
                end
            end
        end % groups loop

        % Write a README sheet explaining columns
        if xl_written
            readme_data = { ...
                'Column',                  'Description'; ...
                'Pair',                    'Label of the replicate pair being compared (RepA vs RepB). Rows labelled "CV of medians (%)" are group-level summaries, not pairwise.'; ...
                'RepA',                    'Label of replicate A (biological replicate identifier)'; ...
                'RepB',                    'Label of replicate B'; ...
                'Feature',                 'Morphology parameter name'; ...
                'Pearson_r',               'Pearson r on 20-quantile profiles of the two replicates. Parametric; sensitive to linear distributional shifts. r=+1: perfect linear agreement.'; ...
                'Pearson_p',               'Two-tailed p-value for Pearson r (df=18), as text. p<0.05: significantly correlated profiles. Shown as "<1e-300" instead of 0 when the true value underflows double precision.'; ...
                'Spearman_rho',            'Spearman rank correlation on same 20-quantile profiles. Non-parametric; robust to skewed distributions. rho=+1: perfect rank agreement.'; ...
                'Spearman_p',              'Two-tailed p-value for Spearman rho (df=18), as text (see Pearson_p note on "<1e-300").'; ...
                'KS_D',                    'Kolmogorov-Smirnov max-CDF distance [0,1]. D=0: identical distributions; D=1: completely non-overlapping.'; ...
                'KS_p',                    'Two-sample KS p-value, as text (see Pearson_p note on "<1e-300"). p<0.05: significant distributional difference between replicates.'; ...
                'KS_sig',                  'Significance stars for KS_p: ***p<0.001  **p<0.01  *p<0.05  ns=not significant'; ...
                'JSD',                     'Jensen-Shannon Divergence [0,1]. JSD=0: identical histogram shapes; JSD=1: maximally different.'; ...
                'AbsDeltaGini',            '|Gini_A - Gini_B| on raw colony values. Difference in within-replicate phenotypic heterogeneity. ~0 = both replicates have similar colony inequality.'; ...
                'AbsDeltaMedian_or_CV_pct','Pairwise rows: |median(A) - median(B)| in original units. CV rows: CV% of per-replicate medians across all replicates in the group.'; ...
                '', ''; ...
                'Gini_<group> sheet', ''; ...
                'Replicate',               'Biological replicate label'; ...
                'Feature',                 'Morphology parameter name'; ...
                'Gini',                    'Classic Gini coefficient on raw colony values [0,1]. 0=all colonies identical (no heterogeneity); 1=maximum inequality. Measures within-replicate phenotypic heterogeneity at the single-colony level.'; ...
                '', ''; ...
                'Interpretation guide', ''; ...
                'Pearson r > 0.95',        'Excellent linear agreement of distribution profiles'; ...
                'Pearson r 0.8-0.95',      'Good linear agreement'; ...
                'Pearson r < 0.8',         'Weak linear agreement — distributions differ in shape or scale'; ...
                'Spearman rho > 0.9',      'Excellent rank agreement between replicates'; ...
                'Spearman rho 0.7-0.9',    'Good agreement'; ...
                'Spearman rho < 0.7',      'Poor agreement — investigate plate-to-plate variability'; ...
                'KS D < 0.1',              'Distributions nearly identical'; ...
                'KS D > 0.2',              'Noticeable distributional shift between replicates'; ...
                'JSD < 0.05',              'Histogram shapes essentially identical'; ...
                'JSD > 0.2',               'Substantially different histogram shapes'; ...
                '|dGini| < 0.05',          'Very similar within-replicate heterogeneity between bio-reps'; ...
                '|dGini| > 0.15',          'One replicate has notably more colony-to-colony variation'; ...
                'CV%  < 5',                'Excellent reproducibility of central tendency across replicates'; ...
                'CV%  5-15',               'Good reproducibility'; ...
                'CV%  > 30',               'High variability — consider excluding outlier replicates'; ...
            };
            try
                readme_tbl = cell2table(readme_data(2:end,:), ...
                    'VariableNames', readme_data(1,:));
                writetable(readme_tbl, xl_path, 'Sheet', 'README');
                tlog('Biological replicate stats saved: %s\n', xl_path);
            catch ME_rm
                tlog('  WARNING: could not write README sheet: %s\n', ME_rm.message);
            end
        end  % if xl_written
        end  % if skipped ... else ... end
    else
        tlog('\n(No biorep_corr field found in ht — run preprocess_pipeline_data first.)\n');
    end

    % ------------------------------------------------------------------
    %  Distribution Statistics Table
    %  — Median, IQR, Gini, n for every group × feature
    %  — Written to a dedicated Excel workbook and logged to terminal
    % ------------------------------------------------------------------
    tlog('\n============================================================\n');
    tlog('  DISTRIBUTION STATISTICS  (Median | IQR | Gini | n)\n');
    tlog('============================================================\n');

    fns_ds  = {'lag_time','size','area','intensity','mean_intensity', ...
               'int_per_size','perimeter','circularity','eccentricity','solidity'};
    flab_ds = {'Lag Time','Size','Area','Intensity','Mean Intensity', ...
               'Int/Size','Perimeter','Circularity','Eccentricity','Solidity'};

    % Terminal header
    tlog('  %-12s  %-20s  %10s  %10s  %8s  %8s\n', ...
        'Group', 'Feature', 'Median', 'IQR', 'Gini', 'n');
    tlog('  %s\n', repmat('-', 1, 76));

    % Accumulate rows for Excel
    ds_group = {};  ds_feat = {};
    ds_med   = [];  ds_iqr  = [];  ds_gini = [];  ds_n = [];

    for gi = 1:length(sel_groups)
        g_idx = sel_groups(gi);
        grp_ds = ht.groups(g_idx);
        for fi = 1:length(fns_ds)
            fn  = fns_ds{fi};
            v   = grp_ds.(fn);
            if isempty(v), continue; end
            vc  = v(~isnan(v));
            if isempty(vc), continue; end
            med_v  = median(vc);
            iqr_v  = iqr(vc);
            gini_v = grp_ds.gini.(fn);
            n_v    = length(vc);
            tlog('  %-12s  %-20s  %10.4g  %10.4g  %8.4f  %8d\n', ...
                labels{g_idx}, flab_ds{fi}, med_v, iqr_v, gini_v, n_v);
            ds_group{end+1,1} = labels{g_idx};   %#ok<AGROW>
            ds_feat{end+1,1}  = flab_ds{fi};     %#ok<AGROW>
            ds_med(end+1,1)   = med_v;            %#ok<AGROW>
            ds_iqr(end+1,1)   = iqr_v;            %#ok<AGROW>
            ds_gini(end+1,1)  = gini_v;           %#ok<AGROW>
            ds_n(end+1,1)     = n_v;              %#ok<AGROW>
        end
        tlog('\n');
    end
    tlog('============================================================\n\n');

    % Write to Excel
    if ~isempty(ds_group)
        ts_ds   = datestr(now, 'yyyymmdd_HHMMSS');
        xl_ds   = fullfile(out_dir, sprintf('distribution_stats_%s.xlsx', ts_ds));
        T_ds    = table(ds_group, ds_feat, ds_med, ds_iqr, ds_gini, ds_n, ...
                    'VariableNames', {'Group','Feature','Median','IQR','Gini','n'});
        try
            writetable(T_ds, xl_ds, 'Sheet', 'Distribution_Stats');

            % README sheet
            readme_ds = { ...
                'Column',   'Description'; ...
                'Group',    'Sample group label (e.g. 3h, 7h, 24h, 48h)'; ...
                'Feature',  'Phenotypic feature name'; ...
                'Median',   'Median value across all valid colonies in the group (NaN-excluded)'; ...
                'IQR',      'Interquartile range (Q75 - Q25) — robust spread measure for non-normal distributions'; ...
                'Gini',     'Gini coefficient [0,1]: 0 = all colonies identical, 1 = maximum heterogeneity'; ...
                'n',        'Number of valid (non-NaN) colonies contributing to this row'; ...
            };
            readme_ds_tbl = cell2table(readme_ds(2:end,:), 'VariableNames', readme_ds(1,:));
            writetable(readme_ds_tbl, xl_ds, 'Sheet', 'README');
            tlog('Distribution statistics saved: %s\n', xl_ds);
        catch ME_ds
            tlog('  WARNING: could not write distribution stats Excel file: %s\n', ME_ds.message);
        end
    end

    % ------------------------------------------------------------------
    %  Pairwise Group Hypothesis Tests
    %  — precomputed in preprocess_pipeline_data (ht.pairwise); this
    %  section only filters to the user-selected groups/features and
    %  prints + exports. No statistics are calculated here.
    % ------------------------------------------------------------------
    tlog('\n============================================================\n');
    tlog('  PAIRWISE GROUP HYPOTHESIS TESTS\n');
    tlog('============================================================\n');
    tlog('  Test choice per pair per feature (Lilliefors normality test on\n');
    tlog('  each group''s raw per-colony values, computed in\n');
    tlog('  preprocess_pipeline_data):\n');
    tlog('    both groups normal  -> Welch two-sample t-test\n');
    tlog('    otherwise           -> Wilcoxon rank-sum test (Mann-Whitney U)\n');
    tlog('  p       = uncorrected p-value.\n');
    tlog('  p_FDR   = Benjamini-Hochberg FDR-corrected p, per feature, across\n');
    tlog('            ALL group pairs tested for that feature (the full set\n');
    tlog('            computed in preprocess_pipeline_data — not just the\n');
    tlog('            pairs/groups shown here if you selected a subset).\n');
    tlog('============================================================\n\n');

    if ~isfield(ht, 'pairwise') || ~isfield(ht.pairwise, 'p')
        tlog('(No hypothesis-test fields in ht.pairwise — re-run preprocess_pipeline_data.)\n\n');
    else
        pw = ht.pairwise;
        has_fdr = isfield(pw, 'p_fdr');

        % pw_pval / pw_pfdr are formatted strings (not raw doubles) so a
        % p-value that underflowed to exact 0.0 exports as "<1e-300"
        % rather than a misleading literal 0.
        pw_pair = {};  pw_feat_c = {};  pw_test = {};
        pw_pval = {};  pw_pfdr = {};  pw_stat_v = [];  pw_sig = {};  pw_sigfdr = {};
        pw_nA   = [];  pw_nB     = [];
        pw_medA = [];  pw_medB   = [];  pw_medDiff = [];

        hdr = sprintf('  %-24s  %-18s  %-17s  %10s %4s  %10s %4s  %6s  %6s  %10s  %10s\n', ...
            'Pair', 'Feature', 'Test', 'p', 'sig', 'p_FDR', 'sig', 'n_A', 'n_B', 'Median_A', 'Median_B');
        tlog('%s', hdr);
        tlog('  %s\n', repmat('-', 1, 134));

        for pi_p = 1:length(pw.pair_labels)
            ga = pw.pair_idx(pi_p, 1);
            gb = pw.pair_idx(pi_p, 2);
            if ~ismember(ga, sel_groups) || ~ismember(gb, sel_groups)
                continue;
            end
            for fi = plot_feat_ix
                pval = pw.p(pi_p, fi);
                if isnan(pval), continue; end
                tname = pw.test_name{pi_p, fi};

                if has_fdr, pfdr = pw.p_fdr(pi_p, fi); else, pfdr = NaN; end

                % Stars are computed from the raw numeric p/p_fdr (a p that
                % underflowed to exact 0.0 is still maximally significant);
                % only the DISPLAYED number is swapped for "<1e-300".
                sig    = sig_stars(pval);
                sigfdr = sig_stars(pfdr);

                tlog('  %-24s  %-18s  %-17s  %10s %4s  %10s %4s  %6d  %6d  %10.4g  %10.4g\n', ...
                    pw.pair_labels{pi_p}, pw.feat_names{fi}, tname, ...
                    format_pval(pval), sig, format_pval(pfdr), sigfdr, ...
                    pw.n_A(pi_p,fi), pw.n_B(pi_p,fi), pw.median_A(pi_p,fi), pw.median_B(pi_p,fi));

                pw_pair{end+1,1}    = pw.pair_labels{pi_p};        %#ok<AGROW>
                pw_feat_c{end+1,1}  = pw.feat_names{fi};            %#ok<AGROW>
                pw_test{end+1,1}    = tname;                        %#ok<AGROW>
                pw_pval{end+1,1}    = format_pval(pval);             %#ok<AGROW>
                pw_pfdr{end+1,1}    = format_pval(pfdr);              %#ok<AGROW>
                pw_stat_v(end+1,1)  = pw.stat(pi_p,fi);                %#ok<AGROW>
                pw_sig{end+1,1}     = strtrim(sig);                     %#ok<AGROW>
                pw_sigfdr{end+1,1}  = strtrim(sigfdr);                   %#ok<AGROW>
                pw_nA(end+1,1)      = pw.n_A(pi_p,fi);                    %#ok<AGROW>
                pw_nB(end+1,1)      = pw.n_B(pi_p,fi);                     %#ok<AGROW>
                pw_medA(end+1,1)    = pw.median_A(pi_p,fi);                 %#ok<AGROW>
                pw_medB(end+1,1)    = pw.median_B(pi_p,fi);                  %#ok<AGROW>
                pw_medDiff(end+1,1) = pw.median_diff(pi_p,fi);                %#ok<AGROW>
            end
        end
        tlog('\n============================================================\n\n');

        if ~isempty(pw_pair)
            ts_pw = datestr(now, 'yyyymmdd_HHMMSS');
            xl_pw = fullfile(out_dir, sprintf('pairwise_group_tests_%s.xlsx', ts_pw));
            T_pw = table(pw_pair, pw_feat_c, pw_test, pw_pval, pw_sig, pw_pfdr, pw_sigfdr, ...
                         pw_stat_v, pw_nA, pw_nB, pw_medA, pw_medB, pw_medDiff, ...
                'VariableNames', {'Pair','Feature','Test','p','Sig','p_FDR','Sig_FDR', ...
                                  'Statistic','n_A','n_B','Median_A','Median_B','AbsMedianDiff'});
            try
                writetable(T_pw, xl_pw, 'Sheet', 'Pairwise_Tests');

                readme_pw = { ...
                    'Column',        'Description'; ...
                    'Pair',          'Sample-group pair being compared (A vs B)'; ...
                    'Feature',       'Phenotypic feature name'; ...
                    'Test',          '"t-test (Welch)" if both groups'' raw values pass a Lilliefors normality test, else "Wilcoxon rank-sum" (Mann-Whitney U)'; ...
                    'p',             'Uncorrected p-value of the chosen test, as text. Shown as "<1e-300" instead of 0 when the true value underflows double precision (huge n / huge effect size).'; ...
                    'Sig',           'Stars for uncorrected p: *** p<0.001, ** p<0.01, * p<0.05, ns = not significant'; ...
                    'p_FDR',         'Benjamini-Hochberg FDR-corrected p-value, as text (see p column note on "<1e-300"), corrected per feature across ALL group pairs tested for that feature in preprocess_pipeline_data (not only the rows shown/exported here if groups/features were subset-selected)'; ...
                    'Sig_FDR',       'Stars for p_FDR, same thresholds as Sig'; ...
                    'Statistic',     't-statistic (Welch t-test) or z/rank-sum statistic (Wilcoxon)'; ...
                    'n_A / n_B',     'Number of valid (non-NaN) colonies used from group A / B'; ...
                    'Median_A / B',  'Median of each group''s raw (unbinned) values'; ...
                    'AbsMedianDiff', '|Median_A - Median_B|'; ...
                };
                readme_pw_tbl = cell2table(readme_pw(2:end,:), 'VariableNames', readme_pw(1,:));
                writetable(readme_pw_tbl, xl_pw, 'Sheet', 'README');
                tlog('Pairwise group hypothesis tests saved: %s\n\n', xl_pw);
            catch ME_pw
                tlog('  WARNING: could not write pairwise tests Excel file: %s\n\n', ME_pw.message);
            end
        end
    end

    % ------------------------------------------------------------------
    %  Save figure (.fig, .pdf, .png)
    % ------------------------------------------------------------------
    feat_str  = strrep(num2str(plot_feat_ix), '  ', '-');
    timestamp = datestr(now, 'yyyymmdd_HHMMSS');
    fig_name  = sprintf('HeteroTyper_features_%s_%s', feat_str, timestamp);
    fig_path  = fullfile(out_dir, fig_name);

    savefig(hfig, [fig_path '.fig']);
    exportgraphics(hfig, [fig_path '.pdf'], 'ContentType','vector');
    exportgraphics(hfig, [fig_path '.png'], 'Resolution',300);

    tlog('Figure saved:\n  %s.fig\n  %s.pdf\n  %s.png\n', ...
            fig_path, fig_path, fig_path);

    % Close log file
    if flog ~= -1
        fclose(flog);
        fprintf('Log saved: %s\n', log_path);
    end

end  % end main function


%% ================================================================
%  LOCAL FUNCTION: sig_stars
%  Significance stars for a p-value. Returns 'NA ' for NaN (e.g. a
%  comparison never run, or no FDR field available).
% ================================================================
function s = sig_stars(pval)
    if isnan(pval),      s = 'NA ';
    elseif pval < 0.001, s = '***';
    elseif pval < 0.01,  s = '** ';
    elseif pval < 0.05,  s = '*  ';
    else,                s = 'ns ';
    end
end


%% ================================================================
%  LOCAL FUNCTION: format_pval
%  Text display for a p-value. A p-value computed from ranksum/ttest2/
%  corr/kstest2 etc. can legitimately underflow to exact 0.0 in double
%  precision when n is large and the effect is strong (the true value
%  is just too small to represent, not actually zero) — display that
%  case as "<1e-300" instead of a misleading literal 0.
% ================================================================
function s = format_pval(p)
    if isnan(p)
        s = 'NA';
    elseif p == 0
        s = '<1e-300';
    else
        s = sprintf('%.4g', p);
    end
end


%% ================================================================
%  LOCAL FUNCTION: local_normhist
%  Normalised histogram (probability density) from a data vector
%  and bin-edge vector.  Mirrors the preprocessor's logic exactly.
% ================================================================
function nc = local_normhist(vec, bin_edges)
    counts = histcounts(vec, bin_edges);
    total  = sum(counts);
    if total > 0
        nc = counts / total;
    else
        nc = counts;
    end
end


%% ================================================================
%  LOCAL FUNCTION: round_to_sig
%  Round x to n significant figures (for bin step sizing).
%  Duplicated from preprocess_pipeline_data to keep this file
%  self-contained — no dependency on the preprocessor being on path.
% ================================================================
function y = round_to_sig(x, n)
    if x == 0, y = 1; return; end
    d = ceil(log10(abs(x)));
    p = n - d;
    y = round(x * 10^p) / 10^p;
    if y <= 0, y = 10^(-p); end
end



%% ================================================================
%  LOCAL FUNCTION: ask_selections
%  Shows groups and features together in one screen, asks for both
%  selections before returning.
% ================================================================
function [sel_groups, plot_feat_ix] = ask_selections(ht)

    all_labels = ht.labels;
    n_avail    = length(all_labels);

    % Print the full menu in one shot
    fprintf('\n============================================================\n');
    fprintf('  plot_combined_samples — Selection\n');
    fprintf('============================================================\n');

    fprintf('\n  GROUPS\n');
    fprintf('  +-----+--------+-------------+\n');
    fprintf('  | idx | label  | n colonies  |\n');
    fprintf('  +-----+--------+-------------+\n');
    for k = 1:n_avail
        fprintf('  | %3d | %-6s | %11d |\n', k, all_labels{k}, ht.groups(k).n_colonies);
    end
    fprintf('  +-----+--------+-------------+\n');

    fprintf('\n  FEATURES\n');
    fprintf('  +----+--------------------------+\n');
    fprintf('  |  1 | Lag time                 |\n');
    fprintf('  |  2 | Colony size              |\n');
    fprintf('  |  3 | Area                     |\n');
    fprintf('  |  4 | Intensity                |\n');
    fprintf('  |  5 | Mean intensity           |\n');
    fprintf('  |  6 | Intensity / size         |\n');
    fprintf('  |  7 | Perimeter                |\n');
    fprintf('  |  8 | Circularity              |\n');
    fprintf('  |  9 | Eccentricity             |\n');
    fprintf('  | 10 | Solidity                 |\n');
    fprintf('  +----+--------------------------+\n');
    fprintf('\n============================================================\n');

    % --- Group selection ---
    while true
        fprintf('  Groups  — enter labels e.g. [ 3h 24h ]  or  ALL: ');
        raw_g = strtrim(input('', 's'));
        if strcmpi(raw_g, 'all')
            sel_groups = 1:n_avail;
            fprintf('  -> All %d groups selected.\n', n_avail);
            break;
        end
        tokens = regexp(raw_g, '[,\s]+', 'split');
        tokens = tokens(~cellfun(@isempty, tokens));
        sel_groups = zeros(1, length(tokens));
        ok = true;
        for t = 1:length(tokens)
            match = find(strcmpi(all_labels, tokens{t}));
            if isempty(match)
                fprintf('  WARNING  "%s" not recognised. Try again.\n', tokens{t});
                ok = false; break;
            end
            sel_groups(t) = match(1);
        end
        if ~ok, continue; end
        [~, ui]    = unique(sel_groups, 'stable');
        sel_groups = sel_groups(ui);
        fprintf('  -> Groups: %s\n', strjoin(all_labels(sel_groups), ', '));
        break;
    end

    % --- Feature selection ---
    while true
        fprintf('  Features — enter numbers e.g. [ 1 2 3 ]  or  ALL: ');
        raw_f  = strtrim(input('', 's'));
        if strcmpi(raw_f, 'all')
            plot_feat_ix = 1:10;
            fprintf('  -> All 10 features selected.\n');
            return;
        end
        tokens = regexp(raw_f, '[,\s]+', 'split');
        nums   = str2double(tokens);
        if any(isnan(nums)) || any(nums < 1) || any(nums > 10) || any(nums ~= floor(nums))
            fprintf('  WARNING  Enter integers 1-10 or ALL.\n');
            continue;
        end
        [~, ui]      = unique(nums, 'stable');
        plot_feat_ix = nums(ui);
        fprintf('  -> Features: [%s]\n', num2str(plot_feat_ix));
        return;
    end
end
% ================================================================
%  LOCAL FUNCTION: ask_feature_selection
%  Only user prompt still needed — everything else is in ht.
% ================================================================
function plot_feat_ix = ask_feature_selection()

    fprintf('\nSelect phenotypic features to visualise\n');
    fprintf('   +----+--------------------------+\n');
    fprintf('   |  1 | Lag time                 |\n');
    fprintf('   |  2 | Colony size              |\n');
    fprintf('   |  3 | Area                     |\n');
    fprintf('   |  4 | Intensity                |\n');
    fprintf('   |  5 | Mean intensity           |\n');
    fprintf('   |  6 | Intensity / size         |\n');
    fprintf('   |  7 | Perimeter                |\n');
    fprintf('   |  8 | Circularity              |\n');
    fprintf('   |  9 | Eccentricity             |\n');
    fprintf('   | 10 | Solidity                 |\n');
    fprintf('   +----+--------------------------+\n');
    fprintf('   Enter numbers separated by commas (e.g.  1,2,3) or namespace (e.g. 1 2 3)\n');
    fprintf('   or type  ALL  to include all 10 features.\n');

    while true
        raw = strtrim(input('   Your selection: ', 's'));

        if strcmpi(raw, 'all')
            plot_feat_ix = 1:10;
            fprintf('   -> Plotting all 10 features.\n');
            return;
        end

        tokens = regexp(strtrim(raw), '[,\s]+', 'split');
        nums   = str2double(tokens);

        if any(isnan(nums)) || any(nums < 1) || any(nums > 10) || any(nums ~= floor(nums))
            fprintf('   WARNING  Invalid input. Enter integers 1-10 separated by commas, or ALL.\n');
            continue;
        end

        [~, ui]      = unique(nums, 'stable');
        plot_feat_ix = nums(ui);
        fprintf('   -> Plotting features: [%s]\n', num2str(plot_feat_ix));
        return;
    end

end
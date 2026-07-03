%% HeteroTyper Pipeline
%  plot_combined_samples_biological_replicates
%  Identical panel layout to plot_combined_samples, but produces ONE
%  FIGURE PER BIOLOGICAL REPLICATE instead of one figure for all groups.
%  Within each figure the rows are sample groups and columns are the
%  selected phenotypic features, showing only the colonies that belong
%  to that biological replicate. Each replicate's sample groups are also
%  compared pairwise (e.g. Rep1: 3h vs 7h, 3h vs 24h, ...) using the
%  normality-aware hypothesis test (Welch t-test / Wilcoxon rank-sum)
%  computed in preprocess_pipeline_data (STEP 6i) — NOT the same
%  comparison as ht.biorep_corr, which compares replicates against each
%  other WITHIN one sample group instead.
%
%  All pooling, histograms, and hypothesis tests are precomputed in
%  preprocess_pipeline_data(data) -> ht.biorep_group_tests. This script
%  only renders and prints; it performs no statistical calculations.
%
%%  USAGE
%  -----
%    plot_combined_samples_biological_replicates(ht)
%
%  The function asks the user to choose:
%    (a) which sample groups to include
%    (b) which features to plot
%  (The biological replicate column itself was already chosen when
%  preprocess_pipeline_data(data) ran — STEP 6i uses ht.params.bio_rep_col.)
%
%  Outputs per replicate (saved to ht.params.out_dir):
%    biorep_<label>_features_<ids>_<ts>.fig / .pdf / .png
%  Plus:
%    plot_biorep_samples_<ts>.txt                 - combined log
%    biorep_group_pairwise_tests_<ts>.xlsx        - one sheet per replicate

function plot_combined_samples_biological_replicates(ht)

    %% ----------------------------------------------------------------
    %  0.  Validate
    %% ----------------------------------------------------------------
    if ~isstruct(ht) || ~isfield(ht, 'groups')
        error(['Input must be the ht struct from preprocess_pipeline_data.\n' ...
               'Run:  preprocess_pipeline_data(data)  first.']);
    end
    if ~isfield(ht, 'biorep_group_tests') || isempty(ht.biorep_group_tests)
        error(['ht.biorep_group_tests is missing.\n' ...
               'Re-run preprocess_pipeline_data(data) to populate it (STEP 6i).']);
    end
    if isfield(ht.biorep_group_tests(1), 'skipped') && ht.biorep_group_tests(1).skipped
        error(['STEP 6i was skipped during preprocessing (%s).\n' ...
               'Re-run preprocess_pipeline_data(data) and select a biological\n' ...
               'replicate column at parameter step 6 (do NOT press 0).'], ...
               ht.biorep_group_tests(1).skip_reason);
    end

    p       = ht.params;
    labels  = ht.labels;
    colors  = ht.colors;
    xb      = ht.xbins;
    out_dir = p.out_dir;
    incTime = p.incTime;
    max_lag = p.max_lag;
    bgt     = ht.biorep_group_tests;   % precomputed per-replicate data + tests (STEP 6i)
    n_reps  = length(bgt);

    %% ----------------------------------------------------------------
    %  1.  Log file
    %% ----------------------------------------------------------------
    ts_log   = datestr(now, 'yyyymmdd_HHMMSS');
    log_path = fullfile(out_dir, sprintf('plot_biorep_samples_%s.txt', ts_log));
    flog = fopen(log_path, 'w');
    if flog == -1
        warning('Could not open log file: %s', log_path);
        flog = -1;
    end

    function tlog(varargin)
        fprintf(varargin{:});
        if flog ~= -1
            fprintf(flog, varargin{:});
        end
    end

    tlog('plot_combined_samples_biological_replicates\n');
    tlog('Run: %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
    tlog('Output directory: %s\n', out_dir);
    tlog('Biological replicate column: "%s"  (fixed at preprocessing time — STEP 6i)\n\n', ...
        p.bio_rep_col);

    %% ----------------------------------------------------------------
    %  2.  Select groups and features
    %% ----------------------------------------------------------------
    [sel_groups, plot_feat_ix] = ask_selections(ht);

    tlog('Groups   : %s\n',  strjoin(labels(sel_groups), ', '));
    tlog('Features : [%s]\n\n', num2str(plot_feat_ix));
    tlog('Biological replicates found (%d): %s\n\n', ...
        n_reps, strjoin({bgt.rep_label}, ', '));

    %% ----------------------------------------------------------------
    %  3.  Bin midpoints / axis limits — rendering only. The actual
    %      histogram counts were already computed in STEP 6i using the
    %      same bin edges as ht.groups(g).hist.*.
    %% ----------------------------------------------------------------
    lag_step = xb.lag_step;
    xp_lag   = xb.lag(1:end-1)          + lag_step/2;
    xp_size  = xb.size(1:end-1)         + xb.size_step/2;
    xp_area  = xb.area(1:end-1)         + xb.area_step/2;
    xp_peri  = xb.perimeter(1:end-1)    + xb.peri_step/2;
    xp_circ  = xb.circularity(1:end-1)  + 0.005;
    xp_ecc   = xb.eccentricity(1:end-1) + 0.005;
    xp_sol   = xb.solidity(1:end-1)     + 0.005;

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

    % Panel catalogue: {xp midpoints, x-axis limits, x-label, feature column 1-10}
    % Feature column indexing matches ht.pairwise.feat_names / bgt(r).features order.
    % Row index == feature column (1-10), matching ht.pairwise.feat_names /
    % bgt(r).features order, so panel_catalog{feat_id, ...} lines up with
    % bgt(r).hist{a_idx, feat_id} etc. directly.
    panel_catalog = { ...
        xp_lag,  [incTime, max_lag+lag_step], 'Lag time (h)';      ...
        xp_size, [0, xb.size_upper],          'Colony size (px)';  ...
        xp_area, [0, xb.area_upper],          'Area (px)';         ...
        xp_int,  [0, int_upper],              'Intensity';         ...
        xp_mint, [0, mint_upper],             'Mean intensity';    ...
        xp_ips,  [0, ips_upper],              'Intensity / size';  ...
        xp_peri, [0, xb.peri_upper],          'Perimeter (px)';    ...
        xp_circ, [0, 1],                      'Circularity';       ...
        xp_ecc,  [0, 1],                      'Eccentricity';      ...
        xp_sol,  [0, 1],                      'Solidity';          ...
    };

    %% ----------------------------------------------------------------
    %  4.  Figure layout (identical to plot_combined_samples)
    %% ----------------------------------------------------------------
    n_rows = length(sel_groups);
    n_cols = length(plot_feat_ix);

    panel_w  = 3.5;   panel_h  = 2.2;
    gap_x    = 0.55;  gap_y    = 0.65;
    margin_l = 1.00;  margin_r = 0.30;
    margin_b = 0.80;  margin_t = 0.55;   % extra top margin for sgtitle
    ylbl_w   = 0.45;

    fig_w = margin_l + n_cols*panel_w + (n_cols-1)*(gap_x+ylbl_w) + margin_r;
    fig_h = margin_b + n_rows*panel_h + (n_rows-1)*gap_y + margin_t;

    ax_pos = zeros(n_rows, n_cols, 4);
    for row = 1:n_rows
        for col = 1:n_cols
            ax_pos(row,col,:) = [ ...
                (margin_l + (col-1)*(panel_w+gap_x+ylbl_w)) / fig_w, ...
                (margin_b + (n_rows-row)*(panel_h+gap_y))   / fig_h, ...
                panel_w / fig_w, ...
                panel_h / fig_h];
        end
    end

    fs_tick  = 12;
    fs_label = 13;
    fs_annot = 12;
    fs_gini  = 11;

    %% ----------------------------------------------------------------
    %  5.  Pairwise hypothesis-test table for one replicate — printed
    %      and appended to the shared Excel workbook. Nested function so
    %      it can share this function's workspace (bgt, labels, xl_path, ...).
    %% ----------------------------------------------------------------
    xl_path    = fullfile(out_dir, sprintf('biorep_group_pairwise_tests_%s.xlsx', ts_log));
    xl_written = false;

    function write_pairwise_table(r)
        rep_lbl       = bgt(r).rep_label;
        pair_labels_r = bgt(r).pair_labels;

        tlog('  Pairwise group hypothesis tests (Rep %s):\n', rep_lbl);
        if isempty(pair_labels_r)
            tlog('    (fewer than 2 sample groups present in this replicate — no comparison.)\n\n');
            return;
        end

        has_fdr = isfield(bgt(r), 'p_fdr') && ~isempty(bgt(r).p_fdr);

        hdr = sprintf('    %-24s  %-18s  %-17s  %10s %4s  %10s %4s  %6s  %6s  %10s  %10s\n', ...
            'Pair', 'Feature', 'Test', 'p', 'sig', 'p_FDR', 'sig', 'n_A', 'n_B', 'Median_A', 'Median_B');
        tlog('%s', hdr);
        tlog('    %s\n', repmat('-', 1, 132));
        tlog('    (p_FDR: Benjamini-Hochberg corrected per feature across this replicate''s own pairs only.)\n');

        % row_p / row_pfdr are formatted strings (not raw doubles) so a
        % p-value that underflowed to exact 0.0 exports as "<1e-300"
        % rather than a misleading literal 0.
        row_pair = {};  row_feat = {};  row_test = {};  row_sig = {};  row_sigfdr = {};
        row_p = {};  row_pfdr = {};  row_stat = [];  row_nA = [];  row_nB = [];
        row_medA = [];  row_medB = [];  row_medDiff = [];

        for pi_p = 1:length(pair_labels_r)
            tok = strsplit(pair_labels_r{pi_p}, ' vs ');
            if length(tok) < 2, continue; end
            gA = find(strcmp(labels, tok{1}), 1);
            gB = find(strcmp(labels, tok{2}), 1);
            if isempty(gA) || isempty(gB), continue; end
            if ~ismember(gA, sel_groups) || ~ismember(gB, sel_groups), continue; end

            for fi = plot_feat_ix
                pval = bgt(r).p(pi_p, fi);
                if isnan(pval), continue; end
                tname = bgt(r).test_name{pi_p, fi};

                if has_fdr, pfdr = bgt(r).p_fdr(pi_p, fi); else, pfdr = NaN; end

                % Stars use the raw numeric p/p_fdr; only the DISPLAYED
                % number is swapped for "<1e-300" on exact underflow.
                sig    = sig_stars(pval);
                sigfdr = sig_stars(pfdr);

                tlog('    %-24s  %-18s  %-17s  %10s %4s  %10s %4s  %6d  %6d  %10.4g  %10.4g\n', ...
                    pair_labels_r{pi_p}, bgt(r).features{fi}, tname, ...
                    format_pval(pval), sig, format_pval(pfdr), sigfdr, ...
                    bgt(r).n_A(pi_p,fi), bgt(r).n_B(pi_p,fi), ...
                    bgt(r).median_A(pi_p,fi), bgt(r).median_B(pi_p,fi));

                row_pair{end+1,1}    = pair_labels_r{pi_p};          %#ok<AGROW>
                row_feat{end+1,1}    = bgt(r).features{fi};           %#ok<AGROW>
                row_test{end+1,1}    = tname;                         %#ok<AGROW>
                row_p{end+1,1}       = format_pval(pval);              %#ok<AGROW>
                row_sig{end+1,1}     = strtrim(sig);                    %#ok<AGROW>
                row_pfdr{end+1,1}    = format_pval(pfdr);                %#ok<AGROW>
                row_sigfdr{end+1,1}  = strtrim(sigfdr);                  %#ok<AGROW>
                row_stat(end+1,1)    = bgt(r).stat(pi_p,fi);              %#ok<AGROW>
                row_nA(end+1,1)      = bgt(r).n_A(pi_p,fi);                %#ok<AGROW>
                row_nB(end+1,1)      = bgt(r).n_B(pi_p,fi);                 %#ok<AGROW>
                row_medA(end+1,1)    = bgt(r).median_A(pi_p,fi);             %#ok<AGROW>
                row_medB(end+1,1)    = bgt(r).median_B(pi_p,fi);              %#ok<AGROW>
                row_medDiff(end+1,1) = bgt(r).median_diff(pi_p,fi);            %#ok<AGROW>
            end
        end
        tlog('\n');

        if isempty(row_pair)
            return;
        end

        T_r = table(row_pair, row_feat, row_test, row_p, row_sig, row_pfdr, row_sigfdr, ...
                    row_stat, row_nA, row_nB, row_medA, row_medB, row_medDiff, ...
            'VariableNames', {'Pair','Feature','Test','p','Sig','p_FDR','Sig_FDR', ...
                              'Statistic','n_A','n_B','Median_A','Median_B','AbsMedianDiff'});

        sheet_name = matlab.lang.makeValidName(rep_lbl);
        sheet_name = sheet_name(1:min(31, length(sheet_name)));
        try
            writetable(T_r, xl_path, 'Sheet', sheet_name);
            xl_written = true;
        catch ME_r
            tlog('    WARNING: could not write Excel sheet "%s": %s\n', sheet_name, ME_r.message);
        end
    end

    %% ----------------------------------------------------------------
    %  6.  Draw one figure per biological replicate
    %% ----------------------------------------------------------------
    for r = 1:n_reps
        rep_lbl = bgt(r).rep_label;
        rep_fn  = matlab.lang.makeValidName(rep_lbl);

        tlog('--- Replicate: %s ---\n', rep_lbl);
        tlog('  Sample groups present: %s\n', strjoin(bgt(r).groups_present, ', '));

        hfig = figure('Name', sprintf('Biorep: %s', rep_lbl), ...
                      'Color','w', 'Units','inches', ...
                      'Position',[1 1 fig_w fig_h], ...
                      'PaperUnits','inches', ...
                      'PaperSize',[fig_w fig_h], ...
                      'PaperPosition',[0 0 fig_w fig_h], ...
                      'Visible','off');

        for gi = 1:n_rows
            g       = sel_groups(gi);
            grp_col = colors(g, :);
            grp_lbl = labels{g};
            a_idx   = find(strcmp(bgt(r).groups_present, grp_lbl), 1);

            for ci = 1:n_cols
                feat_id = plot_feat_ix(ci);
                xp      = panel_catalog{feat_id, 1};
                xlims   = panel_catalog{feat_id, 2};
                xlbl    = panel_catalog{feat_id, 3};

                ax = axes('Position', squeeze(ax_pos(gi, ci, :))'); %#ok<LAXES>

                norm_p = [];
                if ~isempty(a_idx)
                    norm_p = bgt(r).hist{a_idx, feat_id};
                end

                if isempty(norm_p) || all(norm_p == 0)
                    text(0.5, 0.5, 'no data', 'Units','normalized', ...
                         'HorizontalAlignment','center', ...
                         'Color',[0.55 0.55 0.55], 'FontSize', fs_tick);
                    set(ax, 'XLim', xlims, 'YLim', [0 1], ...
                            'FontSize', fs_tick, 'TickDir','in', 'Box','on');
                    if gi == n_rows
                        xlabel(xlbl, 'FontSize', fs_label, 'FontWeight','bold');
                    end
                    ylabel('Probability density', 'FontSize', fs_label, 'FontWeight','bold');
                    continue;
                end

                % Histogram — precomputed normalised counts (STEP 6i)
                bar(xp, norm_p, 'FaceColor',grp_col, 'EdgeColor','none', 'FaceAlpha',1);
                hold on;
                plot(xp, norm_p, 'Color','k', 'LineWidth',0.3);

                y_max = max(norm_p);
                if y_max == 0, y_max = 1; end
                axis([xlims(1) xlims(2) 0 y_max]);

                % Lag time: decade x-ticks
                if feat_id == 1
                    xticks(incTime:10:max_lag);
                end

                % Median line — precomputed
                med_val = bgt(r).median(a_idx, feat_id);
                if ~isnan(med_val)
                    line([med_val med_val], [0 y_max], 'Color','r', 'LineWidth',1.5);
                end

                % Group label + colony count + Gini — precomputed
                gini_val = bgt(r).gini(a_idx, feat_id);
                n_val    = bgt(r).n(a_idx, feat_id);
                x_txt    = xlims(1) + 0.02*(xlims(2)-xlims(1));

                text(x_txt, y_max*0.97, sprintf('%s  (n=%d)', grp_lbl, n_val), ...
                     'FontSize',fs_annot, 'FontWeight','bold', 'Color','k', ...
                     'HorizontalAlignment','left', 'VerticalAlignment','top', ...
                     'Clipping','off');
                text(x_txt, y_max*0.84, sprintf('Gini=%.3f', gini_val), ...
                     'FontSize',fs_gini, 'Color',grp_col, ...
                     'HorizontalAlignment','left', 'VerticalAlignment','top', ...
                     'Clipping','off');

                % x-label on bottom row only
                if gi == n_rows
                    xlabel(xlbl, 'FontSize', fs_label, 'FontWeight','bold');
                end
                ylabel('Probability density', 'FontSize', fs_label, 'FontWeight','bold');

                set(ax, 'FontSize',fs_tick, 'TickDir','in', 'Box','on', ...
                    'XMinorTick','off', 'YMinorTick','off');

                % Mirror tick marks on top and right (cosmetic, no labels)
                ax2 = axes('Position', ax.Position, ...
                           'XAxisLocation','top', 'YAxisLocation','right', ...
                           'Color','none', 'XColor',ax.XColor, 'YColor',ax.YColor, ...
                           'FontSize',fs_tick, 'TickDir','in', 'Box','off', ...
                           'XTick',ax.XTick, 'YTick',ax.YTick, ...
                           'XTickLabel',{}, 'YTickLabel',{}, ...
                           'XLim',ax.XLim, 'YLim',ax.YLim);
                uistack(ax2, 'bottom');
                axes(ax); %#ok<LAXES>

                tlog('  %-12s  %-20s  n=%d  median=%.4g  Gini=%.4f\n', ...
                    grp_lbl, bgt(r).features{feat_id}, n_val, med_val, gini_val);
            end
        end

        % Figure title showing replicate name
        sgtitle(hfig, sprintf('Biological Replicate: %s', rep_lbl), ...
                'FontSize', 15, 'FontWeight','bold');

        set(hfig, 'Visible','on');
        drawnow;

        % Save
        feat_str  = strrep(num2str(plot_feat_ix), '  ', '-');
        ts_fig    = datestr(now, 'yyyymmdd_HHMMSS');
        fig_stem  = sprintf('biorep_%s_features_%s_%s', rep_fn, feat_str, ts_fig);
        fig_path  = fullfile(out_dir, fig_stem);

        savefig(hfig, [fig_path '.fig']);
        exportgraphics(hfig, [fig_path '.pdf'], 'ContentType','vector');
        exportgraphics(hfig, [fig_path '.png'], 'Resolution',300);

        tlog('Saved: %s  (.fig/.pdf/.png)\n\n', fig_path);

        % Pairwise group hypothesis-test table for this replicate
        write_pairwise_table(r);
    end

    %% ----------------------------------------------------------------
    %  7.  README sheet for the pairwise-tests workbook
    %% ----------------------------------------------------------------
    if xl_written
        readme_pw = { ...
            'Column',        'Description'; ...
            'Pair',          'Sample-group pair compared WITHIN this replicate (A vs B)'; ...
            'Feature',       'Phenotypic feature name'; ...
            'Test',          '"t-test (Welch)" if both groups'' raw values pass a Lilliefors normality test, else "Wilcoxon rank-sum" (Mann-Whitney U)'; ...
            'p',             'Uncorrected p-value of the chosen test, as text. Shown as "<1e-300" instead of 0 when the true value underflows double precision (huge n / huge effect size).'; ...
            'Sig',           'Stars for uncorrected p: *** p<0.001, ** p<0.01, * p<0.05, ns = not significant'; ...
            'p_FDR',         'Benjamini-Hochberg FDR-corrected p-value, as text (see p column note on "<1e-300"), corrected per feature across ALL pairs tested WITHIN this replicate (this replicate''s correction family only — not pooled with other replicates)'; ...
            'Sig_FDR',       'Stars for p_FDR, same thresholds as Sig'; ...
            'Statistic',     't-statistic (Welch t-test) or z/rank-sum statistic (Wilcoxon)'; ...
            'n_A / n_B',     'Number of valid (non-NaN) colonies used from group A / B, this replicate only'; ...
            'Median_A / B',  'Median of each group''s raw (unbinned) values, this replicate only'; ...
            'AbsMedianDiff', '|Median_A - Median_B|'; ...
            '', ''; ...
            'Sheet per replicate', 'Each sheet holds the pairwise tests for one biological replicate (e.g. Rep1: 3h vs 7h, 3h vs 24h, 7h vs 24h, ...), independent of every other replicate.'; ...
        };
        try
            readme_pw_tbl = cell2table(readme_pw(2:end,:), 'VariableNames', readme_pw(1,:));
            writetable(readme_pw_tbl, xl_path, 'Sheet', 'README');
            tlog('Per-replicate pairwise hypothesis tests saved: %s\n', xl_path);
        catch ME_rm
            tlog('  WARNING: could not write README sheet: %s\n', ME_rm.message);
        end
    end

    %% ----------------------------------------------------------------
    %  8.  Close log
    %% ----------------------------------------------------------------
    tlog('Done.\n');
    if flog ~= -1
        fclose(flog);
        fprintf('Log saved: %s\n', log_path);
    end

end  % end main function


%% ====================================================================
%  LOCAL: ask_selections  (same as plot_combined_samples)
%% ====================================================================
function [sel_groups, plot_feat_ix] = ask_selections(ht)

    all_labels = ht.labels;
    n_avail    = length(all_labels);

    fprintf('\n============================================================\n');
    fprintf('  plot_combined_samples_biological_replicates — Selection\n');
    fprintf('============================================================\n');

    fprintf('\n  GROUPS\n');
    fprintf('  +-----+------------+-------------+\n');
    fprintf('  | idx | label      | n colonies  |\n');
    fprintf('  +-----+------------+-------------+\n');
    for k = 1:n_avail
        fprintf('  | %3d | %-10s | %11d |\n', k, all_labels{k}, ht.groups(k).n_colonies);
    end
    fprintf('  +-----+------------+-------------+\n');

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

    % Group selection
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
                fprintf('  WARNING: "%s" not recognised. Try again.\n', tokens{t});
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

    % Feature selection
    while true
        fprintf('  Features — enter numbers e.g. [ 1 2 3 ]  or  ALL: ');
        raw_f = strtrim(input('', 's'));
        if strcmpi(raw_f, 'all')
            plot_feat_ix = 1:10;
            fprintf('  -> All 10 features selected.\n');
            return;
        end
        tokens = regexp(raw_f, '[,\s]+', 'split');
        nums   = str2double(tokens);
        if any(isnan(nums)) || any(nums < 1) || any(nums > 10) || any(nums ~= floor(nums))
            fprintf('  WARNING: enter integers 1-10 or ALL.\n');
            continue;
        end
        [~, ui]      = unique(nums, 'stable');
        plot_feat_ix = nums(ui);
        fprintf('  -> Features: [%s]\n', num2str(plot_feat_ix));
        return;
    end
end


%% ====================================================================
%  LOCAL: round_to_sig
%% ====================================================================
function y = round_to_sig(x, n)
    if x == 0, y = 1; return; end
    d = ceil(log10(abs(x)));
    p = n - d;
    y = round(x * 10^p) / 10^p;
    if y <= 0, y = 10^(-p); end
end


%% ====================================================================
%  LOCAL: sig_stars
%  Significance stars for a p-value. Returns 'NA ' for NaN (e.g. no
%  FDR field available, or comparison never run).
%% ====================================================================
function s = sig_stars(pval)
    if isnan(pval),      s = 'NA ';
    elseif pval < 0.001, s = '***';
    elseif pval < 0.01,  s = '** ';
    elseif pval < 0.05,  s = '*  ';
    else,                s = 'ns ';
    end
end


%% ====================================================================
%  LOCAL: format_pval
%  Text display for a p-value. A p-value from ranksum/ttest2 can
%  legitimately underflow to exact 0.0 in double precision when n is
%  large and the effect is strong (the true value is just too small to
%  represent) — display that case as "<1e-300" instead of a misleading
%  literal 0.
%% ====================================================================
function s = format_pval(p)
    if isnan(p)
        s = 'NA';
    elseif p == 0
        s = '<1e-300';
    else
        s = sprintf('%.4g', p);
    end
end

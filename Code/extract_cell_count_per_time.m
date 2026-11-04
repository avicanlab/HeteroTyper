%% HeteroTyper Pipeline
% Two independent, mutually-exclusive analyses. The script asks once, at
% the very start, which one to run — only that one executes.
%
%  STEP 1 — Colony appearance counts per timepoint
%    For each selected sample group, counts how many colonies newly reach
%    a pixel-size threshold at each imaging timepoint — i.e. the number
%    of colonies whose smoothed size timecourse FIRST crosses the
%    threshold at that timepoint. A colony counted at one timepoint is
%    never counted again later, so summing all timepoints reproduces the
%    same valid-colony count used by plot_combined_samples.
%
%  STEP 2 — EARLY vs LATE appearance timepoints: growth rate, doubling
%           time, and size growth over +5h
%    Within ONE selected sample group, pick how many "EARLY" and how many
%    "LATE" imaging timepoints to compare (not fixed to 2+2 — e.g. 1 vs 3
%    is fine), then that many timepoint indices for each phase (into that
%    group's time grid). Colonies are assigned to whichever selected
%    timepoint they first reached the size threshold at (same rule as
%    STEP 1). For the SAME colonies in each timepoint bin, two figures
%    are produced, each a jittered scatter across all selected timepoints
%    (x = timepoint, colour = EARLY/LATE):
%      Figure 1 — Growth rate  (ht.groups(g).growth.mu_max)
%      Figure 2 — Doubling time (ht.groups(g).doublingTime)
%    plus a 3rd figure/analysis:
%      Figure 3 — Colony size growth over +5h: for the SAME colonies in
%        each bin, size at their own appearance timepoint next to size
%        ~5h later (paired scatter with a connecting line per colony),
%        plus a main slope fit (px/h) — one for the pooled EARLY bins'
%        before/+5h points, one for the pooled LATE bins' — with its
%        equation, Pearson r, R^2 and p, both logged and drawn on the plot.
%    plus a 4th figure/analysis:
%      Figures 4, 5 and 6 — Colony size at THREE timepoints per bin: Figure 4 =
%        X-5h, X, X+5h; Figure 5 = X-3h, X, X+5h; Figure 6 = X-2h, X, X+2h
%        (same colonies; only colonies valid at all three; same analysis,
%        different window lengths), with TWO
%        main slopes per cohort (EARLY, LATE): pre-appearance X-L -> X and
%        post-appearance X -> X+A, each with equation, R^2, r and p.
%        EARLY vs LATE within a window = independent slopes (compare_slopes,
%        side-panel bars); pre vs post within a cohort = the same colonies,
%        so each colony's own pre/post slope is compared with the paired
%        test (Lilliefors on differences -> paired t-test / signrank).
%      Figure 7 — Pre-appearance windows only: X-3h -> X and X-2h -> X, each
%        analysed on its own (no post window) and shown as two subpanels of
%        ONE figure with a shared y-axis. Per window and cohort: pooled OLS
%        slope (equation, R^2, r, p), EARLY vs LATE slope comparison
%        (compare_slopes) and EARLY vs LATE per-colony slopes (Lilliefors-
%        gated Welch/rank-sum + KS). Each window's subpanel has its own
%        steepness bar panel beside it (EARLY vs LATE slope +/- SE, with the
%        compare_slopes p), bar panels sharing one y-axis.
%
%    Statistics (per feature, Figures 1-2, and the +5h size GAIN in
%    Figure 3 pooled EARLY vs LATE): every pairwise timepoint comparison,
%    plus one planned contrast (pooled EARLY vs pooled LATE). The test
%    choice is decided PER PAIR from a Lilliefors normality test on each
%    of that pair's two samples: both normal -> Welch's two-sample t-test
%    (ttest2, 'Vartype','unequal'); otherwise Wilcoxon rank-sum (ranksum).
%    A two-sample Kolmogorov-Smirnov test (kstest2) is additionally
%    reported for whole-distribution shape, alongside the location test —
%    same convention as ht.pairwise in preprocess_pipeline_data. BH-FDR
%    corrects the pairwise p-values (and separately the KS p-values) per
%    feature, across the pairwise comparisons run for that feature.
%    See STATISTICS_NOTES.md decision log (2026-09-29 / 2026-09-30).
%
%    Figure 3's before/+5h comparison is a MATCHED-PAIRS design (the same
%    colonies measured twice), so Welch/rank-sum would ignore the pairing
%    and is NOT used there. Instead: Lilliefors on the paired differences
%    -> paired t-test (ttest), or Wilcoxon signed-rank (signrank)
%    otherwise. The +5h size GAIN compared between pooled EARLY and pooled
%    LATE (different colonies in each cohort) IS an independent-samples
%    comparison, so that one reuses the Welch/rank-sum + KS convention.
%    The main slope per cohort is a separate, additional analysis: an
%    ordinary-least-squares fit (polyfit, degree 1) of colony size against
%    elapsed time since that colony's own appearance timepoint (0 or the
%    bin's actual +5h offset), pooled across every colony in every EARLY
%    (or every LATE) bin — same Pearson r / R^2 / slope-equation
%    convention as plot_QC_full_dataset.m's OLS section.
%
%    Colonies within a timepoint bin can come from more than one plate
%    (biological/technical replicate) pooled into the same sample group.
%    This is a potential source of clustering/batch bias that a simple
%    two-sample test does not account for, so the per-plate composition
%    of every bin is also reported — check it before over-interpreting a
%    significant p-value driven by one plate.
%
% Requires preprocess_pipeline_data(data) to have been run first.
% STEP 2 additionally requires ht.groups(g).growth (STEP 6d).
%
%% USAGE:
%   extract_cell_count_per_time(ht);

function extract_cell_count_per_time(ht)

    % ------------------------------------------------------------------
    %  Validate input
    % ------------------------------------------------------------------
    if ~isstruct(ht) || ~isfield(ht, 'groups')
        error(['Input must be the ht struct produced by preprocess_pipeline_data.\n' ...
               'Run:  preprocess_pipeline_data(data)  first.']);
    end

    p       = ht.params;
    out_dir = p.out_dir;

    if isfield(p, 'size_threshold')
        default_threshold = p.size_threshold;
    else
        default_threshold = 100;
    end

    % ------------------------------------------------------------------
    %  Ask which step to run — mutually exclusive
    % ------------------------------------------------------------------
    step_choice = ask_step_choice();

    % ------------------------------------------------------------------
    %  Open log file — mirrors everything printed to the terminal
    % ------------------------------------------------------------------
    timestamp_log = datestr(now, 'yyyymmdd_HHMMSS');
    log_path = fullfile(out_dir, ...
        sprintf('extract_cell_count_per_time_Step%d_%s.txt', step_choice, timestamp_log));
    flog = fopen(log_path, 'w');
    if flog == -1
        warning('Could not open log file: %s', log_path);
        flog = -1;
    end
    tlog('Log file: %s\n', log_path);
    tlog('Step selected: %d\n\n', step_choice);

    function tlog(varargin)
        % Write to terminal AND log file simultaneously
        fprintf(varargin{:});
        if flog ~= -1
            fprintf(flog, varargin{:});
        end
    end

    tlog_fn = @tlog;

    if step_choice == 1
        run_step1(ht, default_threshold, tlog_fn);
    else
        run_step2(ht, default_threshold, tlog_fn);
    end

    tlog('Done.\n');

    % Close log file
    if flog ~= -1
        fclose(flog);
        fprintf('Log saved: %s\n', log_path);
    end

end  % end main function


%% ================================================================
%  STEP 1 — Colony appearance counts per timepoint
% ================================================================
function run_step1(ht, default_threshold, tlog)

    labels  = ht.labels;
    out_dir = ht.params.out_dir;

    [sel_groups, threshold] = ask_selections_step1(ht, default_threshold);

    tlog('\n--- Parameters ---\n');
    tlog('  Groups to analyze      : %s\n', strjoin(labels(sel_groups), ', '));
    tlog('  Appearance threshold   : %g px\n', threshold);
    tlog('  Preprocessing filter   : %g px (ht.params.size_threshold)\n', default_threshold);
    if threshold < default_threshold
        tlog(['  WARNING: the appearance threshold is LOWER than the size filter\n' ...
              '           already applied in preprocess_pipeline_data. Colonies whose\n' ...
              '           final size was below %g px were NaN-masked across their\n' ...
              '           ENTIRE timecourse there and cannot be recovered here.\n' ...
              '           Re-run preprocess_pipeline_data with a lower cfg.size_threshold\n' ...
              '           if those colonies are needed.\n'], default_threshold);
    end
    tlog('  Output directory       : %s\n', out_dir);
    tlog('-------------------------\n\n');

    group_results = struct('label', {}, 'time_h', {}, 'counts', {}, 'n_valid', {});

    for gi = 1:length(sel_groups)
        g   = sel_groups(gi);
        grp = ht.groups(g);
        lbl = labels{g};

        tlog('============================================================\n');
        tlog('  Group: %s\n', lbl);
        tlog('============================================================\n');

        tc  = grp.size_timecourse;   % [T x N]
        t_h = grp.time_h;            % [T x 1]

        if isempty(tc) || isempty(t_h)
            tlog('  No size_timecourse data available for this group — skipped.\n\n');
            continue;
        end

        [T, N] = size(tc);
        if N == 0
            tlog('  No colonies in this group — skipped.\n\n');
            continue;
        end

        % Colonies counted elsewhere in the pipeline (plot_combined_samples,
        % etc.) as "valid" are exactly the non-NaN entries of grp.size —
        % rejected (final size < filter) and right-censored colonies are
        % already NaN there.
        n_valid = sum(~isnan(grp.size));

        first_idx = compute_first_appearance(tc, threshold);

        counts = zeros(T, 1);
        for t = 1:T
            counts(t) = sum(first_idx == t);
        end

        n_counted  = sum(~isnan(first_idx));
        n_no_cross = n_valid - n_counted;

        tlog('  Colonies total (incl. rejected/censored)              : %d\n', grp.n_colonies);
        tlog('  Colonies valid (same filtering as plot_combined_samples): %d\n', n_valid);
        tlog('  Colonies counted (reached %g px in the timecourse)     : %d\n', threshold, n_counted);
        if n_no_cross > 0
            tlog(['  NOTE: %d valid colonies never reached %g px in their smoothed\n' ...
                  '        timecourse (threshold raised above the preprocessing filter,\n' ...
                  '        or smoothing/final-Area mismatch) — excluded from the counts below.\n'], ...
                  n_no_cross, threshold);
        elseif n_no_cross < 0
            tlog('  NOTE: %d colonies reached threshold despite being NaN in grp.size (unexpected).\n', -n_no_cross);
        end
        tlog('\n');

        tlog('  %10s  %8s  %10s  %10s\n', 'Time (h)', 'New', 'Cumulative', '% of valid');
        tlog('  %s\n', repmat('-', 1, 46));
        cum = 0;
        for t = 1:T
            cum = cum + counts(t);
            if n_valid > 0
                pct = 100 * cum / n_valid;
            else
                pct = 0;
            end
            tlog('  %10.2f  %8d  %10d  %9.1f%%\n', t_h(t), counts(t), cum, pct);
        end
        tlog('\n');

        group_results(end+1) = struct( ...
            'label',   lbl, ...
            'time_h',  t_h, ...
            'counts',  counts, ...
            'n_valid', n_valid); %#ok<AGROW>
    end

    % ------------------------------------------------------------------
    %  Combined table across groups, only if they share the same time grid
    % ------------------------------------------------------------------
    if length(group_results) > 1
        same_grid = true;
        ref_t = group_results(1).time_h;
        for k = 2:length(group_results)
            gt = group_results(k).time_h;
            if length(gt) ~= length(ref_t) || any(abs(gt(:) - ref_t(:)) > 1e-9)
                same_grid = false;
                break;
            end
        end

        if same_grid
            tlog('============================================================\n');
            tlog('  COMBINED — new colonies per timepoint, all selected groups\n');
            tlog('============================================================\n');

            hdr = sprintf('  %10s', 'Time (h)');
            for k = 1:length(group_results)
                hdr = [hdr sprintf('  %10s', group_results(k).label)]; %#ok<AGROW>
            end
            tlog('%s\n', hdr);
            tlog('  %s\n', repmat('-', 1, 12 + 12*length(group_results)));

            for t = 1:length(ref_t)
                row = sprintf('  %10.2f', ref_t(t));
                for k = 1:length(group_results)
                    row = [row sprintf('  %10d', group_results(k).counts(t))]; %#ok<AGROW>
                end
                tlog('%s\n', row);
            end
            tlog('\n');
        end
    end
end


%% ================================================================
%  STEP 2 — EARLY vs LATE: growth rate, final size, doubling time,
%           and size growth over +5h
% ================================================================
function run_step2(ht, default_threshold, tlog)

    labels  = ht.labels;
    out_dir = ht.params.out_dir;

    if ~isfield(ht.groups(1), 'growth')
        error(['Step 2 requires ht.groups(g).growth (max Gompertz growth rate fits).\n' ...
               'Re-run preprocess_pipeline_data (STEP 6d) first.']);
    end

    % --- Pick exactly one sample group ---
    g   = ask_single_group(ht);
    grp = ht.groups(g);
    lbl = labels{g};

    tc  = grp.size_timecourse;
    t_h = grp.time_h;
    if isempty(tc) || isempty(t_h)
        error('Selected group "%s" has no size_timecourse data.', lbl);
    end
    [T, ~] = size(tc);

    % --- Appearance threshold (same concept as Step 1) ---
    threshold = ask_threshold(default_threshold);
    first_idx = compute_first_appearance(tc, threshold);

    % --- Show per-timepoint appearance counts to help pick indices ---
    fprintf('\n  Group %s — colonies newly appearing per timepoint (threshold=%g px):\n', lbl, threshold);
    fprintf('  +-----+----------+----------+\n');
    fprintf('  | idx | time (h) | n_appear |\n');
    fprintf('  +-----+----------+----------+\n');
    for t = 1:T
        fprintf('  | %3d | %8.2f | %8d |\n', t, t_h(t), sum(first_idx == t));
    end
    fprintf('  +-----+----------+----------+\n');

    % --- Pick EARLY + LATE timepoint indices (user chooses how many of each) ---
    [early_idx, late_idx] = ask_timepoint_selection(T);

    n_early = length(early_idx);
    n_late  = length(late_idx);
    n_bins  = n_early + n_late;

    sel_tp      = [early_idx, late_idx];   % 1 x n_bins, EARLY-1..EARLY-n_early, LATE-1..LATE-n_late
    phase_names = [arrayfun(@(i) sprintf('EARLY-%d', i), 1:n_early, 'UniformOutput', false), ...
                   arrayfun(@(i) sprintf('LATE-%d',  i), 1:n_late,  'UniformOutput', false)];
    is_early    = [true(1, n_early), false(1, n_late)];

    EARLY_COLOR = [42  80  161] / 255;   % blue  (matches plot_correlation_matrix)
    LATE_COLOR  = [170 0   36 ] / 255;   % red   (matches plot_correlation_matrix)

    tlog('\n--- Parameters ---\n');
    tlog('  Sample group           : %s\n', lbl);
    tlog('  Appearance threshold   : %g px\n', threshold);
    tlog('  Preprocessing filter   : %g px (ht.params.size_threshold)\n', default_threshold);
    if threshold < default_threshold
        tlog(['  WARNING: threshold is lower than the preprocessing size filter — colonies\n' ...
              '           already NaN-masked upstream cannot be recovered here.\n']);
    end
    tlog('  EARLY / LATE counts    : %d EARLY, %d LATE timepoint(s)\n', n_early, n_late);
    for k = 1:n_bins
        tlog('  %-8s : timepoint idx %2d  (t = %.2f h)\n', phase_names{k}, sel_tp(k), t_h(sel_tp(k)));
    end
    tlog('  Output directory       : %s\n', out_dir);
    tlog('-------------------------\n\n');

    tlog(['  NOTE ON GROWTH RATE / DOUBLING TIME: neither is recomputed by this\n' ...
          '  script. Each is reused as-is from the value already computed once per\n' ...
          '  colony in preprocess_pipeline_data (mu_max: Gompertz fit over that\n' ...
          '  colony''s FULL size timecourse, STEP 6d; doublingTime: log-linear fit\n' ...
          '  over the post-appearance timecourse) — the same values\n' ...
          '  plot_appearanceTime_vs_growthRate.m / plot_combined_samples would show for that\n' ...
          '  colony. Binning colonies by appearance timepoint here only chooses WHICH\n' ...
          '  colonies to compare; it does not restrict either feature to the window\n' ...
          '  between the selected timepoints, so Figures 1-2 compare "colonies that\n' ...
          '  appeared early/late" by their overall (whole-trajectory) phenotype, not a\n' ...
          '  quantity measured only in that window. Figure 3 (size growth over +5h) is\n' ...
          '  the exception — that one DOES use only the size measured at the two\n' ...
          '  specific +5h-apart timepoints.\n\n']);

    % ------------------------------------------------------------------
    %  Extract colonies + both features for each of the n_bins timepoint
    %  bins. colony_idx is the full "appeared at this timepoint" set (same
    %  rule as Step 1); mu/dt keep NaN where that specific feature is
    %  unavailable for a colony (e.g. Gompertz fit failed, or too few
    %  post-appearance points for a doubling-time fit) — every downstream test
    %  drops its own NaNs, so one shared colony_idx per bin is enough.
    % ------------------------------------------------------------------
    bins = struct('label', {}, 'time_h', {}, 'tp_idx', {}, 'is_early', {}, ...
                   'colony_idx', {}, 'n_appear', {}, 'mu', {}, 'dt', {}, ...
                   'n_plates', {}, 'plate_ids', {}, 'plate_counts', {});

    for k = 1:n_bins
        tp      = sel_tp(k);
        col_idx = find(first_idx == tp);
        n_appear = length(col_idx);

        mu_vals = nan(n_appear, 1);
        dt_vals = nan(n_appear, 1);
        if n_appear > 0
            mu_vals = grp.growth.mu_max(col_idx);
            dt_vals = grp.doublingTime(col_idx);
        end

        [n_plates, plate_ids, plate_counts] = plate_breakdown(grp, col_idx);

        bins(k).label        = phase_names{k};
        bins(k).time_h       = t_h(tp);
        bins(k).tp_idx       = tp;
        bins(k).is_early     = is_early(k);
        bins(k).colony_idx   = col_idx;
        bins(k).n_appear     = n_appear;
        bins(k).mu           = mu_vals(:);
        bins(k).dt           = dt_vals(:);
        bins(k).n_plates     = n_plates;
        bins(k).plate_ids    = plate_ids;
        bins(k).plate_counts = plate_counts;
    end

    % ------------------------------------------------------------------
    %  Log: per-bin summary + plate-composition bias check
    % ------------------------------------------------------------------
    tlog('============================================================\n');
    tlog('  BIN SUMMARY  (colonies whose FIRST appearance = this timepoint)\n');
    tlog('============================================================\n');
    tlog('  %-8s  %8s  %8s  %6s  %6s  %8s  %s\n', ...
        'Bin', 'Time(h)', 'n_appear', 'n_mu', 'n_dt', 'n_plates', 'plate composition (plate:n)');
    tlog('  %s\n', repmat('-', 1, 92));
    for k = 1:n_bins
        n_mu = sum(~isnan(bins(k).mu));
        n_dt = sum(~isnan(bins(k).dt));
        comp_str = plate_comp_string(bins(k).plate_ids, bins(k).plate_counts);
        tlog('  %-8s  %8.2f  %8d  %6d  %6d  %8d  %s\n', ...
            bins(k).label, bins(k).time_h, bins(k).n_appear, n_mu, n_dt, bins(k).n_plates, comp_str);
    end
    tlog('\n');
    tlog(['  CAVEAT: colonies within a bin can come from more than one plate\n' ...
          '  (biological/technical replicate) pooled into this sample group. The\n' ...
          '  two-sample tests below treat every colony as an independent draw and\n' ...
          '  do NOT model plate-of-origin. If "n_plates" is 1 for a bin, or counts\n' ...
          '  are dominated by a single plate, treat any significant p-value with\n' ...
          '  caution — it may reflect a plate/batch effect rather than a genuine\n' ...
          '  EARLY-vs-LATE (or timepoint-vs-timepoint) difference.\n\n']);

    for k = 1:n_bins
        n_mu = sum(~isnan(bins(k).mu));
        n_dt = sum(~isnan(bins(k).dt));
        if min([n_mu, n_dt]) < 3
            tlog('  WARNING: %s has n_mu=%d, n_dt=%d — statistics for whichever is <3 will be unreliable.\n', ...
                bins(k).label, n_mu, n_dt);
        end
    end
    tlog('\n');

    % ------------------------------------------------------------------
    %  Figures 1-2: growth rate, doubling time — same pairwise-comparison
    %  + pooled-contrast statistics for each, Lilliefors-gated per pair
    %  (Welch t-test / Wilcoxon rank-sum) + KS.
    % ------------------------------------------------------------------
    feat_defs = struct( ...
        'field', {'mu',                    'dt'}, ...
        'label', {'Growth rate (mu_max)',  'Doubling time'}, ...
        'unit',  {'\mu_{max} (px/h)',      'Doubling time (h)'}, ...
        'tag',   {'GrowthRate',            'DoublingTime'});

    fig_paths = {};
    for fi = 1:numel(feat_defs)
        fd = feat_defs(fi);
        tlog('============================================================\n');
        tlog('  STATISTICS — %s\n', fd.label);
        tlog('============================================================\n');
        [pw_res, res_pooled] = run_bin_stats(tlog, bins, is_early, fd.field, fd.label);
        fp = make_bin_comparison_figure(out_dir, lbl, threshold, bins, is_early, fd.field, ...
                 fd.unit, fd.label, fd.tag, pw_res, res_pooled, EARLY_COLOR, LATE_COLOR);
        fig_paths{end+1} = fp;
        tlog('Figure saved:\n  %s.fig\n  %s.pdf\n  %s.png\n\n', fp, fp, fp);
    end

    % ------------------------------------------------------------------
    %  Figure 3: colony size at each bin's own appearance timepoint vs
    %  ~5h later (SAME colonies — matched pairs), plus the independent-
    %  samples comparison of the size GAIN between pooled EARLY and
    %  pooled LATE cohorts, plus a main slope (OLS fit) per cohort.
    % ------------------------------------------------------------------
    tlog('============================================================\n');
    tlog('  STATISTICS — Colony size growth over +5h\n');
    tlog('============================================================\n');
    tlog(['  Same colonies, measured at their own bin timepoint and again ~5h later.\n' ...
          '  This is a MATCHED-PAIRS comparison (not independent samples), so Welch/\n' ...
          '  rank-sum + KS (used above) do not apply here: a paired test is used\n' ...
          '  instead — Lilliefors on the paired differences -> paired t-test (ttest),\n' ...
          '  or Wilcoxon signed-rank (signrank) otherwise. The +5h size GAIN compared\n' ...
          '  between the pooled EARLY and pooled LATE cohorts (different colonies) IS\n' ...
          '  an independent-samples comparison, so that one reuses Welch/rank-sum + KS.\n\n']);

    step = median(diff(t_h));
    for k = 1:n_bins
        target = bins(k).time_h + 5;
        [idx5, ok5] = nearest_timepoint(t_h, target, step);
        bins(k).tp5_ok = ok5;
        if ok5
            bins(k).tp5_idx     = idx5;
            bins(k).time_h5     = t_h(idx5);
            bins(k).size_before = tc(bins(k).tp_idx, bins(k).colony_idx)';
            bins(k).size_after  = tc(idx5,           bins(k).colony_idx)';
        else
            bins(k).tp5_idx     = NaN;
            bins(k).time_h5     = NaN;
            bins(k).size_before = [];
            bins(k).size_after  = [];
        end
    end

    if ~any([bins.tp5_ok])
        tlog('  WARNING: no bin has a usable +5h timepoint within the imaging window — Figure 3 skipped.\n\n');
    else
        tlog('  %-8s  %8s  %8s  %-20s  %10s %4s\n', 'Bin', 't(h)', 't+5h(h)', 'Test', 'p', 'sig');
        tlog('  %s\n', repmat('-', 1, 68));
        for k = 1:n_bins
            if ~bins(k).tp5_ok
                tlog('  %-8s  %8.2f  %8s  (no +5h timepoint available)\n', bins(k).label, bins(k).time_h, '--');
                continue;
            end
            res_p = compare_paired(bins(k).size_before, bins(k).size_after);
            bins(k).res_pair = res_p;
            tlog('  %-8s  %8.2f  %8.2f  %-20s  %10s %4s\n', bins(k).label, bins(k).time_h, bins(k).time_h5, ...
                blank_if_empty(res_p.test_name), fmt_p(res_p.p), sig_stars(res_p.p));
        end
        tlog('\n');

        d_early = [];
        d_late  = [];
        rel_t_early = []; size_early_all = [];
        rel_t_late  = []; size_late_all  = [];
        for k = 1:n_bins
            if ~bins(k).tp5_ok, continue; end
            valid = ~isnan(bins(k).size_before) & ~isnan(bins(k).size_after);
            before_v = bins(k).size_before(valid);
            after_v  = bins(k).size_after(valid);
            d        = after_v - before_v;
            dt_k     = bins(k).time_h5 - bins(k).time_h;
            n_v      = sum(valid);
            if bins(k).is_early
                d_early        = [d_early; d]; 
                rel_t_early     = [rel_t_early; zeros(n_v,1); repmat(dt_k, n_v, 1)]; 
                size_early_all  = [size_early_all; before_v; after_v]; 
            else
                d_late          = [d_late; d]; 
                rel_t_late      = [rel_t_late; zeros(n_v,1); repmat(dt_k, n_v, 1)]; 
                size_late_all   = [size_late_all; before_v; after_v]; 
            end
        end
        res_growth_cmp = two_sample_test(d_early, d_late);
        tlog('  --- Size GAIN over +5h: pooled EARLY vs pooled LATE (independent cohorts) ---\n');
        tlog('  %-18s  %10s %4s  %8s  %10s %4s  %6s  %6s  %10s  %10s\n', ...
            'Test', 'p', 'sig', 'KS_D', 'KS_p', 'sig', 'n_A', 'n_B', 'Median_A', 'Median_B');
        tlog('  %s\n', repmat('-', 1, 100));
        tlog('  %-18s  %10s %4s  %8.4f  %10s %4s  %6d  %6d  %10.4g  %10.4g\n', ...
            blank_if_empty(res_growth_cmp.test_name), fmt_p(res_growth_cmp.p), sig_stars(res_growth_cmp.p), ...
            res_growth_cmp.ks_D, fmt_p(res_growth_cmp.ks_p), sig_stars(res_growth_cmp.ks_p), ...
            res_growth_cmp.n_a, res_growth_cmp.n_b, res_growth_cmp.median_a, res_growth_cmp.median_b);
        tlog('\n');

        % ---- Main slope: OLS fit of size vs elapsed time since the
        % colony's own appearance timepoint (0, or the bin's actual +5h
        % offset), pooled across every colony in every EARLY (or every
        % LATE) bin — same Pearson r / R^2 / slope-equation convention as
        % plot_QC_full_dataset.m's OLS section. With n_early or n_late > 1
        % this pools across bins on an elapsed-time basis (not absolute
        % clock time), so different bins' own +5h offsets can differ
        % slightly (grid snapping) without biasing the fit.
        fit_early = fit_growth_rate(rel_t_early, size_early_all);
        fit_late  = fit_growth_rate(rel_t_late,  size_late_all);
        tlog('  --- Main slope: colony size vs elapsed time since appearance (OLS) ---\n');
        log_growth_fit(tlog, 'EARLY', fit_early);
        log_growth_fit(tlog, 'LATE',  fit_late);

        % ---- Slope comparison: is EARLY reliably steeper/shallower than
        % LATE, or is the numeric difference within the noise of each fit?
        slope_cmp = compare_slopes(fit_early, fit_late);
        tlog('  --- Slope comparison: EARLY vs LATE steepness ---\n');
        if isnan(slope_cmp.t)
            tlog('    Insufficient data in one or both fits for a slope comparison.\n\n');
        else
            tlog('    Delta slope (EARLY - LATE) = %.4g px/h\n', slope_cmp.diff);
            tlog('    SE(delta)                  = %.4g px/h\n', slope_cmp.se);
            tlog('    t(%d) = %.4f,  p = %s %s\n\n', ...
                slope_cmp.df, slope_cmp.t, fmt_p(slope_cmp.p), strtrim(sig_stars(slope_cmp.p)));
        end

        fp = make_growth5h_figure(out_dir, lbl, bins, res_growth_cmp, fit_early, fit_late, slope_cmp, EARLY_COLOR, LATE_COLOR);
        fig_paths{end+1} = fp;
        tlog('Figure saved:\n  %s.fig\n  %s.pdf\n  %s.png\n\n', fp, fp, fp);
    end

    % ------------------------------------------------------------------
    %  Figures 4, 5 and 6: three timepoints per bin — X-L, X (appearance),
    %  X+A — for the SAME colonies, with TWO slopes per growth category
    %  (EARLY / LATE): pre-appearance (X-L -> X) and post-appearance
    %  (X -> X+A). Figure 4: L=5h, A=5h. Figure 5: L=3h, A=5h.
    %  Figure 6: L=2h, A=2h (a short, symmetric window around X).
    % ------------------------------------------------------------------
    win_cfgs = [5 5; 3 5; 2 2];     % [lookback L, lookahead A] in hours, one row per figure
    for cfg_i = 1:size(win_cfgs, 1)
        LOOKBACK_H  = win_cfgs(cfg_i, 1);
        LOOKAHEAD_H = win_cfgs(cfg_i, 2);
        pre_lbl  = sprintf('X-%dh', LOOKBACK_H);
        post_lbl = sprintf('X+%dh', LOOKAHEAD_H);
        fig_no   = 3 + cfg_i;

        tlog('============================================================\n');
        tlog('  STATISTICS — Figure %d: colony size %s -> X -> %s (two slopes per cohort)\n', ...
            fig_no, pre_lbl, post_lbl);
        tlog('============================================================\n');
        tlog(['  Same colonies at three timepoints: %dh BEFORE their appearance timepoint X\n' ...
              '  (colonies are still below the size threshold then), at X, and %dh after.\n' ...
              '  Only colonies with a valid size at ALL THREE timepoints are used, so the\n' ...
              '  pre and post slopes describe the same colony set. Per cohort (EARLY, LATE),\n' ...
              '  pooled across its bins: one OLS slope for %s -> X and one for X -> %s\n' ...
              '  (same fit convention as above). EARLY vs LATE within one window = different\n' ...
              '  colonies -> independent-slopes test (compare_slopes). Pre vs post WITHIN a\n' ...
              '  cohort = the same colonies -> NOT independent, so compare_slopes would be\n' ...
              '  wrong there; instead each colony''s own pre and post slope (px/h) are\n' ...
              '  compared with the paired test (Lilliefors on the differences -> paired\n' ...
              '  t-test or Wilcoxon signed-rank).\n\n'], LOOKBACK_H, LOOKAHEAD_H, pre_lbl, post_lbl);

        for k = 1:n_bins
            bins(k).tp3_ok       = false;
            bins(k).time_h_pre   = NaN;
            bins(k).time_h_post  = NaN;
            bins(k).size_pre     = [];
            bins(k).size_x       = [];
            bins(k).size_post    = [];
            [idx_pre,  ok_pre]  = nearest_timepoint(t_h, bins(k).time_h - LOOKBACK_H,  step);
            [idx_post, ok_post] = nearest_timepoint(t_h, bins(k).time_h + LOOKAHEAD_H, step);
            if ok_pre && ok_post && idx_pre < bins(k).tp_idx && idx_post > bins(k).tp_idx
                bins(k).tp3_ok      = true;
                bins(k).time_h_pre  = t_h(idx_pre);
                bins(k).time_h_post = t_h(idx_post);
                bins(k).size_pre    = tc(idx_pre,         bins(k).colony_idx)';
                bins(k).size_x      = tc(bins(k).tp_idx,  bins(k).colony_idx)';
                bins(k).size_post   = tc(idx_post,        bins(k).colony_idx)';
            end
        end

        for k = 1:n_bins
            if ~bins(k).tp3_ok
                tlog('  %-8s  (no %s and/or %s timepoint around t=%.2f h within the imaging window — excluded)\n', ...
                    bins(k).label, pre_lbl, post_lbl, bins(k).time_h);
            else
                tlog('  %-8s  t: %.2f h -> %.2f h -> %.2f h\n', bins(k).label, ...
                    bins(k).time_h_pre, bins(k).time_h, bins(k).time_h_post);
            end
        end
        tlog('\n');

        if ~any([bins.tp3_ok])
            tlog('  WARNING: no bin has both a %s and a %s timepoint inside the imaging window — Figure %d skipped.\n\n', ...
                pre_lbl, post_lbl, fig_no);
            continue;
        end

        P_E = pool_three_point_cohort(bins, true);
        P_L = pool_three_point_cohort(bins, false);

        S = struct();
        S.lookback       = LOOKBACK_H;
        S.lookahead      = LOOKAHEAD_H;
        S.pre_lbl        = pre_lbl;
        S.post_lbl       = post_lbl;
        S.fit_pre_early  = fit_growth_rate(P_E.rt_pre,  P_E.y_pre);
        S.fit_post_early = fit_growth_rate(P_E.rt_post, P_E.y_post);
        S.fit_pre_late   = fit_growth_rate(P_L.rt_pre,  P_L.y_pre);
        S.fit_post_late  = fit_growth_rate(P_L.rt_post, P_L.y_post);

        tlog('  --- Main slopes (OLS, size vs elapsed time since the start of each window) ---\n');
        log_growth_fit(tlog, sprintf('EARLY  %s -> X ', pre_lbl),   S.fit_pre_early);
        log_growth_fit(tlog, sprintf('EARLY  X -> %s ', post_lbl),  S.fit_post_early);
        log_growth_fit(tlog, sprintf('LATE   %s -> X ', pre_lbl),   S.fit_pre_late);
        log_growth_fit(tlog, sprintf('LATE   X -> %s ', post_lbl),  S.fit_post_late);

        S.cmp_pre  = compare_slopes(S.fit_pre_early,  S.fit_pre_late);
        S.cmp_post = compare_slopes(S.fit_post_early, S.fit_post_late);
        tlog('  --- EARLY vs LATE steepness, per window (independent colonies) ---\n');
        log_slope_comparison(tlog, sprintf('%s -> X ', pre_lbl),  S.cmp_pre);
        log_slope_comparison(tlog, sprintf('X -> %s ', post_lbl), S.cmp_post);

        S.acc_early = compare_paired(P_E.s_pre, P_E.s_post);
        S.acc_late  = compare_paired(P_L.s_pre, P_L.s_post);
        tlog('  --- Pre vs post slope WITHIN each cohort (same colonies, per-colony slopes, paired) ---\n');
        log_paired_slopes(tlog, 'EARLY', S.acc_early, pre_lbl, post_lbl);
        log_paired_slopes(tlog, 'LATE',  S.acc_late,  pre_lbl, post_lbl);

        fp = make_growth3pt_figure(out_dir, lbl, bins, S, EARLY_COLOR, LATE_COLOR);
        fig_paths{end+1} = fp;
        tlog('Figure %d saved:\n  %s.fig\n  %s.pdf\n  %s.png\n\n', fig_no, fp, fp, fp);
    end

    % ------------------------------------------------------------------
    %  Figure 7: pre-appearance windows ONLY — X-3h -> X and X-2h -> X —
    %  each analysed individually (no post-appearance window) and shown as
    %  two subpanels of ONE figure with a shared y-axis.
    % ------------------------------------------------------------------
    tlog('============================================================\n');
    tlog('  STATISTICS — Figure 7: pre-appearance windows X-3h -> X and X-2h -> X\n');
    tlog('============================================================\n');
    tlog(['  Each window is analysed on its own (subpanels of one figure, shared y-axis).\n' ...
          '  Per cohort: one pooled OLS slope (colonies valid at both timepoints), EARLY\n' ...
          '  vs LATE via compare_slopes (independent colonies), plus EARLY vs LATE on the\n' ...
          '  per-colony slopes (Lilliefors-gated Welch/rank-sum + KS).\n' ...
          '  CAVEAT: colonies are binned by the timepoint at which they FIRST reach the\n' ...
          '  size threshold, so their size at X is >= threshold by construction. These\n' ...
          '  pre-window slopes therefore partly reflect that selection, not a free\n' ...
          '  measurement of growth rate — compare cohorts/windows, do not read the\n' ...
          '  absolute px/h as an unbiased rate.\n\n']);

    pre_panels = struct('L', {}, 'pre_lbl', {}, 'pb', {}, 'fit_e', {}, 'fit_l', {}, ...
                        'cmp', {}, 'cs', {}, 'n_e', {}, 'n_l', {});
    for L_pre = [3 2]
        Pn = compute_pre_window_panel(bins, tc, t_h, step, L_pre);
        pre_panels(end+1) = Pn; %#ok<AGROW>
        tlog('  --- Window %s -> X  (n colonies: EARLY %d, LATE %d) ---\n', Pn.pre_lbl, Pn.n_e, Pn.n_l);
        log_growth_fit(tlog, sprintf('EARLY  %s -> X ', Pn.pre_lbl), Pn.fit_e);
        log_growth_fit(tlog, sprintf('LATE   %s -> X ', Pn.pre_lbl), Pn.fit_l);
        log_slope_comparison(tlog, sprintf('%s -> X, EARLY vs LATE slope ', Pn.pre_lbl), Pn.cmp);
        if isempty(Pn.cs.test_name)
            tlog('  Per-colony slope, EARLY vs LATE: insufficient data.\n\n');
        else
            tlog('  Per-colony slope, EARLY vs LATE: %s p=%s %s; KS D=%.4f p=%s %s; median %.4g vs %.4g px/h (n=%d, %d)\n\n', ...
                Pn.cs.test_name, fmt_p(Pn.cs.p), strtrim(sig_stars(Pn.cs.p)), ...
                Pn.cs.ks_D, fmt_p(Pn.cs.ks_p), strtrim(sig_stars(Pn.cs.ks_p)), ...
                Pn.cs.median_a, Pn.cs.median_b, Pn.cs.n_a, Pn.cs.n_b);
        end
    end

    if any(arrayfun(@(P) any([P.pb.ok]), pre_panels))
        fp = make_pre_window_figure(out_dir, lbl, pre_panels, EARLY_COLOR, LATE_COLOR);
        fig_paths{end+1} = fp;
        tlog('Figure 7 saved:\n  %s.fig\n  %s.pdf\n  %s.png\n\n', fp, fp, fp);
    else
        tlog('  WARNING: no bin has an X-3h or X-2h timepoint inside the imaging window — Figure 7 skipped.\n\n');
    end
end


%% ================================================================
%  LOCAL FUNCTION: compute_pre_window_panel
%  Pre-appearance window ONLY (X-L -> X) for every bin, pooled per cohort
%  (EARLY / LATE): one OLS slope each (fit_growth_rate), EARLY vs LATE
%  slope comparison (compare_slopes), and the EARLY vs LATE comparison of
%  the per-colony slopes (two_sample_test: Lilliefors-gated Welch/rank-sum
%  + KS). Only colonies with a valid size at both X-L and X are used.
% ================================================================
function P = compute_pre_window_panel(bins, tc, t_h, step, L)
    P.L       = L;
    P.pre_lbl = sprintf('X-%dh', L);
    n = numel(bins);
    pb = struct('label', cell(1, n), 'is_early', false, 'ok', false, ...
                'time_pre', NaN, 'time_x', NaN, 'y_pre', [], 'y_x', []);
    for k = 1:n
        pb(k).label    = bins(k).label;
        pb(k).is_early = bins(k).is_early;
        [idx_pre, ok_pre] = nearest_timepoint(t_h, bins(k).time_h - L, step);
        if ok_pre && idx_pre < bins(k).tp_idx
            pb(k).ok       = true;
            pb(k).time_pre = t_h(idx_pre);
            pb(k).time_x   = bins(k).time_h;
            pb(k).y_pre    = tc(idx_pre,        bins(k).colony_idx)';
            pb(k).y_x      = tc(bins(k).tp_idx, bins(k).colony_idx)';
        end
    end
    P.pb = pb;

    Pe = pool_pre_cohort(pb, true);
    Pl = pool_pre_cohort(pb, false);
    P.fit_e = fit_growth_rate(Pe.rt, Pe.y);
    P.fit_l = fit_growth_rate(Pl.rt, Pl.y);
    P.cmp   = compare_slopes(P.fit_e, P.fit_l);
    P.cs    = two_sample_test(Pe.s, Pl.s);
    P.n_e   = Pe.n;
    P.n_l   = Pl.n;
end

function Q = pool_pre_cohort(pb, want_early)
    Q.rt = []; Q.y = []; Q.s = []; Q.n = 0;
    for k = 1:numel(pb)
        if ~pb(k).ok || pb(k).is_early ~= want_early, continue; end
        pre = pb(k).y_pre(:);
        xs  = pb(k).y_x(:);
        v   = ~isnan(pre) & ~isnan(xs);
        if ~any(v), continue; end
        pre = pre(v); xs = xs(v);
        n_v = sum(v);
        dt  = pb(k).time_x - pb(k).time_pre;
        Q.rt = [Q.rt; zeros(n_v,1); repmat(dt, n_v, 1)]; %#ok<AGROW>
        Q.y  = [Q.y;  pre; xs];                          %#ok<AGROW>
        Q.s  = [Q.s;  (xs - pre) / dt];                  %#ok<AGROW>
        Q.n  = Q.n + n_v;
    end
end


%% ================================================================
%  LOCAL FUNCTION: make_pre_window_figure
%  ONE figure, one subpanel per pre-appearance window (e.g. X-3h -> X and
%  X-2h -> X), SHARED y-axis across the subpanels so the windows are
%  directly comparable. In each: EARLY cluster and LATE cluster, a thin
%  line per colony from its size at X-L to its size at X, a dashed pooled
%  OLS slope line per cohort (equation / R^2 / r / p in the corners), and
%  brackets for the EARLY-vs-LATE slope comparison (compare_slopes) and
%  per-colony-slope comparison (two_sample_test).
% ================================================================
function fig_path = make_pre_window_figure(out_dir, lbl, panels, EARLY_COLOR, LATE_COLOR)

    n_pan = numel(panels);

    % Shared y-range over every panel's data
    all_vals = [];
    n_max_bins = 1;
    for pi = 1:n_pan
        for k = 1:numel(panels(pi).pb)
            if ~panels(pi).pb(k).ok, continue; end
            all_vals = [all_vals; panels(pi).pb(k).y_pre(:); panels(pi).pb(k).y_x(:)]; %#ok<AGROW>
        end
        n_max_bins = max(n_max_bins, sum([panels(pi).pb.ok]));
    end
    all_vals = all_vals(~isnan(all_vals));
    if isempty(all_vals)
        y_min = 0; y_max = 1;
    else
        y_min = min(all_vals); y_max = max(all_vals);
    end
    if y_max <= y_min, y_max = y_min + 1; end
    yr = y_max - y_min;

    fig_w = n_pan * max(6.2, 1.4*n_max_bins + 4.8);
    fig_h = 5.4;
    hfig = figure('Name', sprintf('Size growth, pre-appearance windows — %s', lbl), ...
                  'Color', 'w', 'Units', 'inches', ...
                  'Position', [1 1 fig_w fig_h], ...
                  'PaperUnits', 'inches', ...
                  'PaperPosition', [0 0 fig_w fig_h], ...
                  'PaperSize', [fig_w fig_h], 'Visible', 'off');

    % Each window gets one block: scatter panel (left) + its bar panel (right)
    block_w = 0.92 / n_pan;
    jitter = 0.10;
    for pi = 1:n_pan
        P  = panels(pi);
        ax = axes('Position', [0.06 + (pi-1)*block_w, 0.28, 0.66*block_w, 0.56]);
        hold(ax, 'on');

        ok     = find([P.pb.ok]);
        n_ok   = numel(ok);
        if n_ok == 0
            text(ax, 0.5, 0.5, sprintf('No bin has a usable %s timepoint', P.pre_lbl), ...
                'Units', 'normalized', 'HorizontalAlignment', 'center', 'FontSize', 10);
            axis(ax, 'off');
            continue;
        end
        is_e   = arrayfun(@(k) P.pb(k).is_early, ok);
        n_e_ok = sum(is_e);
        n_l_ok = n_ok - n_e_ok;

        x_pos = zeros(n_ok, 2);
        cursor = 1;
        last_is_early = [];
        for m = 1:n_ok
            if ~isempty(last_is_early) && last_is_early ~= is_e(m)
                cursor = cursor + 1.0;
            end
            x_pos(m, :) = [cursor, cursor + 0.7];
            cursor = cursor + 0.7 + 0.5;
            last_is_early = is_e(m);
        end

        xtl = cell(1, 2*n_ok);
        for m = 1:n_ok
            k = ok(m);
            if is_e(m), col = EARLY_COLOR; else, col = LATE_COLOR; end
            line_col = 0.55*col + 0.45*[1 1 1];
            pre = P.pb(k).y_pre(:);
            xs  = P.pb(k).y_x(:);
            v   = ~isnan(pre) & ~isnan(xs);
            pre = pre(v); xs = xs(v);
            n_k = numel(pre);
            if n_k > 0
                xa = x_pos(m,1) + (rand(n_k,1) - 0.5) * 2*jitter;
                xb = x_pos(m,2) + (rand(n_k,1) - 0.5) * 2*jitter;
                plot(ax, [xa'; xb'], [pre'; xs'], 'Color', line_col, 'LineWidth', 0.6);
                scatter(ax, xa, pre, 20, col, 'filled', 'MarkerFaceAlpha', 0.6, 'MarkerEdgeColor', 'none');
                scatter(ax, xb, xs,  20, col, 'filled', 'MarkerFaceAlpha', 0.6, 'MarkerEdgeColor', 'none');
                plot(ax, [x_pos(m,1)-0.18 x_pos(m,1)+0.18], median(pre)*[1 1], 'Color', [0.1 0.1 0.1], 'LineWidth', 2);
                plot(ax, [x_pos(m,2)-0.18 x_pos(m,2)+0.18], median(xs)*[1 1],  'Color', [0.1 0.1 0.1], 'LineWidth', 2);
            end
            multi = (is_e(m) && n_e_ok > 1) || (~is_e(m) && n_l_ok > 1);
            if multi
                xtl{2*m-1} = sprintf('%s\n%.1fh', P.pb(k).label, P.pb(k).time_pre);
            else
                xtl{2*m-1} = sprintf('%.1fh', P.pb(k).time_pre);
            end
            xtl{2*m} = sprintf('%.1fh', P.pb(k).time_x);
        end

        % Pooled slope line per cohort (dashed), from mean "pre" slot to mean "X" slot
        fits = {P.fit_e, P.fit_l};
        cols = {EARLY_COLOR, LATE_COLOR};
        for c = 1:2
            rows = (is_e == (c == 1));
            if ~any(rows) || isnan(fits{c}.slope), continue; end
            dts = arrayfun(@(k) P.pb(k).time_x - P.pb(k).time_pre, ok(rows));
            plot(ax, [mean(x_pos(rows,1)) mean(x_pos(rows,2))], ...
                [fits{c}.intercept, fits{c}.intercept + fits{c}.slope * mean(dts)], ...
                '--', 'Color', cols{c}, 'LineWidth', 3);
        end

        xlim(ax, [0.3, cursor - 0.3]);
        ylim(ax, [y_min - 0.06*yr, y_max + 0.80*yr]);
        set(ax, 'XTick', reshape(x_pos', 1, []), 'XTickLabel', xtl, ...
            'FontSize', 9, 'TickDir', 'in', 'Box', 'on');
        ylabel(ax, 'Colony size (px)', 'FontSize', 12, 'FontWeight', 'bold');
        title(ax, sprintf('%s -> X', P.pre_lbl), 'FontSize', 13, 'FontWeight', 'bold');

        % Brackets: EARLY vs LATE (slope fit, then per-colony slopes)
        if n_e_ok > 0 && n_l_ok > 0
            e_mid = mean(reshape(x_pos(is_e, :), 1, []));
            l_mid = mean(reshape(x_pos(~is_e, :), 1, []));
            if ~isnan(P.cmp.p)
                txt = sprintf('slope EARLY vs LATE: p=%s %s', fmt_p(P.cmp.p), strtrim(sig_stars(P.cmp.p)));
                draw_bracket(ax, e_mid, l_mid, y_max + 0.06*yr, 0.015*yr, txt, 8, 1.2);
            end
            if ~isempty(P.cs.test_name)
                txt = sprintf('per-colony slope, %s: p=%s %s', P.cs.test_name, fmt_p(P.cs.p), strtrim(sig_stars(P.cs.p)));
                draw_bracket(ax, e_mid, l_mid, y_max + 0.20*yr, 0.015*yr, txt, 8, 1.0);
            end
        end

        text(ax, 0.02, 0.97, growth_fit_annotation_str('EARLY', P.fit_e), 'Units', 'normalized', ...
            'FontSize', 8, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', ...
            'Color', EARLY_COLOR, 'BackgroundColor', 'w', 'Margin', 3);
        text(ax, 0.98, 0.97, growth_fit_annotation_str('LATE', P.fit_l), 'Units', 'normalized', ...
            'FontSize', 8, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'right', ...
            'Color', LATE_COLOR, 'BackgroundColor', 'w', 'Margin', 3);

        % EARLY / LATE heading under this panel's tick labels
        drawnow;
        xl_final = xlim(ax);
        pos      = ax.Position;
        hd_w = 0.12; hd_h = 0.05; hd_y = max(0.01, pos(2) - 0.10);
        if n_e_ok > 0
            frac = (mean(reshape(x_pos(is_e, :), 1, [])) - xl_final(1)) / diff(xl_final);
            annotation(hfig, 'textbox', [pos(1)+frac*pos(3)-hd_w/2, hd_y, hd_w, hd_h], ...
                'String', 'EARLY', 'EdgeColor', 'none', 'HorizontalAlignment', 'center', ...
                'FontWeight', 'bold', 'FontSize', 12, 'Color', EARLY_COLOR);
        end
        if n_l_ok > 0
            frac = (mean(reshape(x_pos(~is_e, :), 1, [])) - xl_final(1)) / diff(xl_final);
            annotation(hfig, 'textbox', [pos(1)+frac*pos(3)-hd_w/2, hd_y, hd_w, hd_h], ...
                'String', 'LATE', 'EdgeColor', 'none', 'HorizontalAlignment', 'center', ...
                'FontWeight', 'bold', 'FontSize', 12, 'Color', LATE_COLOR);
        end
    end

    % ------------------------------------------------------------------
    %  Bar panels (one per window, to the right of its scatter panel):
    %  EARLY vs LATE pooled slope +/- SE, compare_slopes p bracketed
    %  between the two bars; shared y-axis across the bar panels.
    % ------------------------------------------------------------------
    bar_tops = []; bar_bot = 0;
    for pi = 1:n_pan
        for f = {panels(pi).fit_e, panels(pi).fit_l}
            fi = f{1};
            if isnan(fi.slope), continue; end
            se0 = fi.se_slope;
            if isnan(se0), se0 = 0; end
            bar_tops(end+1) = fi.slope + se0; %#ok<AGROW>
            bar_bot = min(bar_bot, fi.slope);
        end
    end
    if isempty(bar_tops) || max(bar_tops) <= 0
        bar_top = 1;
    else
        bar_top = max(bar_tops);
    end
    bar_base = bar_top * 1.08;

    bcols = {EARLY_COLOR, LATE_COLOR};
    for pi = 1:n_pan
        P   = panels(pi);
        axb = axes('Position', [0.06 + (pi-1)*block_w + 0.74*block_w, 0.28, 0.21*block_w, 0.56]);
        hold(axb, 'on');
        fits_pair = {P.fit_e, P.fit_l};
        for c = 1:2
            fi = fits_pair{c};
            if isnan(fi.slope), continue; end
            bar(axb, c, fi.slope, 0.6, 'FaceColor', bcols{c}, 'EdgeColor', 'none');
            if ~isnan(fi.se_slope)
                errorbar(axb, c, fi.slope, fi.se_slope, 'Color', 'k', 'LineWidth', 1.2, ...
                    'CapSize', 8, 'LineStyle', 'none');
            end
        end
        xlim(axb, [0.3 2.7]);
        ylim(axb, [bar_bot * 1.2, bar_base + 0.30*bar_top]);
        set(axb, 'XTick', [1 2], 'XTickLabel', {'EARLY', 'LATE'}, ...
            'FontSize', 9, 'TickDir', 'in', 'Box', 'on');
        ylabel(axb, 'Slope (px/h)', 'FontSize', 10, 'FontWeight', 'bold');
        title(axb, sprintf('%s -> X', P.pre_lbl), 'FontSize', 10, 'FontWeight', 'bold');
        grid(axb, 'on');
        if ~isnan(P.cmp.p)
            txt_b = sprintf('p=%s %s', fmt_p(P.cmp.p), strtrim(sig_stars(P.cmp.p)));
            draw_bracket(axb, 1, 2, bar_base, 0.02*bar_top, txt_b, 8, 1.0);
        end
    end

    annotation(hfig, 'textbox', [0 0.93 1 0.06], ...
        'String', sprintf('Colony size before appearance — Group %s (shared y-axes)', lbl), ...
        'EdgeColor', 'none', 'HorizontalAlignment', 'center', ...
        'FontWeight', 'bold', 'FontSize', 13);

    set(hfig, 'Visible', 'on');
    drawnow;

    timestamp_fig = datestr(now, 'yyyymmdd_HHMMSS');
    fig_name = sprintf('SizeGrowth_PreWindows_EarlyLate_%s_%s', lbl, timestamp_fig);
    fig_path = fullfile(out_dir, fig_name);
    savefig(hfig, [fig_path '.fig']);
    exportgraphics(hfig, [fig_path '.pdf'], 'ContentType', 'vector');
    exportgraphics(hfig, [fig_path '.png'], 'Resolution', 300);
end


%% ================================================================
%  LOCAL FUNCTION: compute_first_appearance
%  First row index per colony (column) where the smoothed size
%  timecourse reaches >= threshold. NaN if it never does (rejected /
%  censored colonies are already all-NaN, so they never satisfy the
%  comparison and are automatically excluded).
% ================================================================
function first_idx = compute_first_appearance(tc, threshold)
    N = size(tc, 2);
    first_idx = nan(N, 1);
    for ci = 1:N
        idx = find(tc(:, ci) >= threshold, 1, 'first');
        if ~isempty(idx)
            first_idx(ci) = idx;
        end
    end
end


%% ================================================================
%  LOCAL FUNCTION: nearest_timepoint
%  Nearest grid index to `target`; ok=false if target is beyond the
%  imaging window or no grid point is within one sampling step of it.
% ================================================================
function [idx, ok] = nearest_timepoint(t_h, target, step)
    if target > max(t_h) + 0.5*step || target < min(t_h) - 0.5*step
        idx = NaN; ok = false;
        return;
    end
    [dmin, idx] = min(abs(t_h - target));
    ok = dmin <= step;
end


%% ================================================================
%  LOCAL FUNCTION: fit_growth_rate
%  Ordinary-least-squares fit of size (y) vs elapsed time (x), plus
%  Pearson r / R^2 / slope-equation — same convention as
%  plot_QC_full_dataset.m's OLS section (polyfit degree 1, R^2 from
%  Pearson r, "y = mx + c" equation string). Returns NaN/'n/a' fields
%  when there are fewer than 3 points or all x values are identical
%  (a degree-1 fit is then undefined).
% ================================================================
function fit = fit_growth_rate(x, y)
    x = x(:); y = y(:);
    valid = ~isnan(x) & ~isnan(y);
    x = x(valid); y = y(valid);

    fit.n         = length(x);
    fit.slope     = NaN;
    fit.intercept = NaN;
    fit.r         = NaN;
    fit.p         = NaN;
    fit.R2        = NaN;
    fit.eq_str    = 'n/a';
    fit.se_slope  = NaN;   % standard error of the slope (for compare_slopes)
    fit.dof       = NaN;   % residual degrees of freedom (n - 2)

    if fit.n < 3 || numel(unique(x)) < 2
        return;
    end

    [r, pv] = corr(x, y, 'Type', 'Pearson');
    fit.r  = r;
    fit.p  = pv;
    fit.R2 = r^2;

    coeff         = polyfit(x, y, 1);
    fit.slope     = coeff(1);
    fit.intercept = coeff(2);
    if fit.intercept >= 0
        fit.eq_str = sprintf('y = %.3fx + %.1f', fit.slope, fit.intercept);
    else
        fit.eq_str = sprintf('y = %.3fx - %.1f', fit.slope, abs(fit.intercept));
    end

    % Standard error of the slope (standard OLS formula), needed to test
    % whether two independently-fitted slopes differ (see compare_slopes).
    fit.dof = fit.n - 2;
    if fit.dof > 0
        y_pred = polyval(coeff, x);
        s_e    = sqrt(sum((y - y_pred).^2) / fit.dof);
        ss_x   = sum((x - mean(x)).^2);
        if ss_x > 0
            fit.se_slope = s_e / sqrt(ss_x);
        end
    end
end


%% ================================================================
%  LOCAL FUNCTION: log_growth_fit
% ================================================================
function log_growth_fit(tlog, cohort_label, fit)
    tlog('  %s (pooled across all colonies in that cohort''s bins):\n', cohort_label);
    if isnan(fit.slope)
        tlog('    Insufficient data for a linear fit (n=%d points).\n\n', fit.n);
        return;
    end
    tlog('    n points       = %d\n', fit.n);
    tlog('    Slope          = %.4g px/h\n', fit.slope);
    tlog('    Intercept      = %.4g px\n', fit.intercept);
    tlog('    Pearson r      = %.4f,  p = %s\n', fit.r, fmt_p(fit.p));
    tlog('    R^2            = %.4f\n', fit.R2);
    tlog('    Equation       : %s\n\n', fit.eq_str);
end


%% ================================================================
%  LOCAL FUNCTION: compare_slopes
%  Tests whether two independently-fitted regression slopes differ:
%  t = (slope_a - slope_b) / sqrt(SE_a^2 + SE_b^2), df = dof_a + dof_b
%  (standard test for a difference between two independent regression
%  coefficients, using each fit's own OLS standard error). This answers
%  "is EARLY reliably steeper/shallower than LATE", which the two slope
%  numbers alone do not — a numerically bigger slope can still be
%  statistically indistinguishable from the other given the noise in
%  each fit. NaN fields if either fit was itself insufficient (see
%  fit_growth_rate) or the combined standard error is exactly zero.
% ================================================================
function res = compare_slopes(fit_a, fit_b)
    res.diff = NaN;
    res.se   = NaN;
    res.t    = NaN;
    res.df   = NaN;
    res.p    = NaN;

    if isnan(fit_a.slope) || isnan(fit_b.slope) || isnan(fit_a.se_slope) || isnan(fit_b.se_slope)
        return;
    end

    res.diff = fit_a.slope - fit_b.slope;
    res.se   = sqrt(fit_a.se_slope^2 + fit_b.se_slope^2);
    res.df   = fit_a.dof + fit_b.dof;

    if res.se == 0 || res.df <= 0
        return;
    end
    res.t = res.diff / res.se;
    res.p = 2 * (1 - tcdf(abs(res.t), res.df));
end


%% ================================================================
%  LOCAL FUNCTION: run_bin_stats
%  All pairwise timepoint-bin comparisons for one feature (BH-FDR
%  corrected, p and KS p as separate families), plus the planned
%  pooled-EARLY-vs-pooled-LATE contrast. Logs a table for each and
%  returns the pairwise results + pooled result for the figure to reuse.
% ================================================================
function [pw_res, res_pooled] = run_bin_stats(tlog, bins, is_early, feat_field, feat_label)

    n_bins  = length(bins);
    pairs   = nchoosek(1:n_bins, 2);
    n_pairs = size(pairs, 1);

    pw_p   = nan(n_pairs, 1);
    pw_ksp = nan(n_pairs, 1);
    pw_res = cell(n_pairs, 1);
    for pi = 1:n_pairs
        a = pairs(pi, 1); b = pairs(pi, 2);
        res = two_sample_test(bins(a).(feat_field), bins(b).(feat_field));
        pw_res{pi} = res;
        pw_p(pi)   = res.p;
        pw_ksp(pi) = res.ks_p;
    end
    pw_pfdr    = bh_fdr(pw_p);
    pw_ks_pfdr = bh_fdr(pw_ksp);

    tlog('  --- %s: all pairwise timepoint comparisons ---\n', feat_label);
    tlog('  (test decided per pair: Lilliefors on both samples -> Welch t-test if\n');
    tlog('   both normal, else Wilcoxon rank-sum; + two-sample KS test)\n');
    tlog('  %-18s  %-18s  %10s %4s  %10s %4s  %8s  %10s %4s  %6s  %6s  %10s  %10s\n', ...
        'Pair', 'Test', 'p', 'sig', 'p_FDR', 'sig', 'KS_D', 'KS_p', 'sig', 'n_A', 'n_B', 'Median_A', 'Median_B');
    tlog('  %s\n', repmat('-', 1, 140));
    for pi = 1:n_pairs
        a = pairs(pi, 1); b = pairs(pi, 2);
        res = pw_res{pi};
        pair_lbl = sprintf('%s vs %s', bins(a).label, bins(b).label);
        tlog('  %-18s  %-18s  %10s %4s  %10s %4s  %8.4f  %10s %4s  %6d  %6d  %10.4g  %10.4g\n', ...
            pair_lbl, blank_if_empty(res.test_name), fmt_p(res.p), sig_stars(res.p), ...
            fmt_p(pw_pfdr(pi)), sig_stars(pw_pfdr(pi)), res.ks_D, fmt_p(res.ks_p), sig_stars(res.ks_p), ...
            res.n_a, res.n_b, res.median_a, res.median_b);
    end
    tlog('\n');

    early_vals = vertcat(bins(is_early).(feat_field));
    late_vals  = vertcat(bins(~is_early).(feat_field));
    res_pooled = two_sample_test(early_vals, late_vals);

    tlog('  --- %s: planned contrast, pooled EARLY vs pooled LATE ---\n', feat_label);
    tlog('  %-18s  %10s %4s  %8s  %10s %4s  %6s  %6s  %10s  %10s\n', ...
        'Test', 'p', 'sig', 'KS_D', 'KS_p', 'sig', 'n_A', 'n_B', 'Median_A', 'Median_B');
    tlog('  %s\n', repmat('-', 1, 100));
    tlog('  %-18s  %10s %4s  %8.4f  %10s %4s  %6d  %6d  %10.4g  %10.4g\n', ...
        blank_if_empty(res_pooled.test_name), fmt_p(res_pooled.p), sig_stars(res_pooled.p), ...
        res_pooled.ks_D, fmt_p(res_pooled.ks_p), sig_stars(res_pooled.ks_p), ...
        res_pooled.n_a, res_pooled.n_b, res_pooled.median_a, res_pooled.median_b);
    tlog('\n');
end


%% ================================================================
%  LOCAL FUNCTION: draw_bracket
%  Significance bracket: a horizontal line from x1 to x2 at height y
%  with short downward end-ticks, labelled with p-value + stars text
%  centred just above the line.
% ================================================================
function draw_bracket(ax, x1, x2, y, tick_h, txt, font_sz, line_w)
    line(ax, [x1 x1 x2 x2], [y-tick_h, y, y, y-tick_h], ...
         'Color', 'k', 'LineWidth', line_w);
    text(ax, (x1+x2)/2, y, txt, ...
         'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
         'FontSize', font_sz);
end


%% ================================================================
%  LOCAL FUNCTION: make_bin_comparison_figure
%  Jittered scatter of one feature across all selected timepoint bins
%  (x = timepoint, colour = EARLY/LATE), median bar per bin, and
%  significance brackets for every pairwise comparison + the pooled
%  EARLY-vs-LATE contrast. Saves .fig/.pdf/.png.
% ================================================================
function fig_path = make_bin_comparison_figure(out_dir, lbl, threshold, bins, is_early, feat_field, ...
        y_unit_label, feat_title, fig_tag, pw_res, res_pooled, EARLY_COLOR, LATE_COLOR)

    n_bins  = length(bins);
    pairs   = nchoosek(1:n_bins, 2);
    n_pairs = size(pairs, 1);

    AX_BOTTOM      = 0.24;   % axes bottom, normalized fig units — leaves
                             % room below for tick labels + EARLY/LATE row
    GROUP_LABEL_DY = 0.10;   % how far below the axes bottom the EARLY/LATE
                             % heading sits (normalized fig units)

    all_vals = vertcat(bins.(feat_field));
    all_vals = all_vals(~isnan(all_vals));
    y_data_top = 1;
    if ~isempty(all_vals) && max(all_vals) > 0
        y_data_top = max(all_vals) * 1.08;
    end
    y_data_bottom = 0;
    if ~isempty(all_vals) && min(all_vals) < 0
        y_data_bottom = min(all_vals) * 1.08;
    end

    n_levels    = n_pairs + 1;                 % all pairwise + 1 pooled
    yspan       = y_data_top - y_data_bottom;
    if yspan <= 0, yspan = 1; end
    step        = 0.085 * yspan;
    tick_h      = 0.018 * yspan;
    y_top_total = y_data_top + (n_levels + 1) * step;

    fig_w = 7.0; fig_h = 5.0 + 0.28*n_levels;
    hfig = figure('Name', sprintf('%s: EARLY vs LATE', feat_title), ...
                  'Color', 'w', ...
                  'Units', 'inches', ...
                  'Position', [1 1 fig_w fig_h], ...
                  'PaperUnits', 'inches', ...
                  'PaperSize', [fig_w fig_h], ...
                  'PaperPosition', [0 0 fig_w fig_h], ...
                  'Visible', 'off');
    ax = axes('Position', [0.13 AX_BOTTOM 0.80 0.98-AX_BOTTOM-0.06]);
    ax.PositionConstraint = 'innerposition';   % keep this exact box fixed —
                                                % titles/labels must not shrink
                                                % it, or the EARLY/LATE heading
                                                % below would drift out of place
    hold(ax, 'on');

    tick_t = [bins.time_h];
    if numel(unique(tick_t)) > 1
        min_gap = min(diff(sort(unique(tick_t))));
    else
        min_gap = max(1, tick_t(1)*0.1);
    end
    jitter_w = 0.30 * min_gap;
    x_pad    = max(2.5*jitter_w, 0.06*range(tick_t));
    if x_pad == 0, x_pad = 1; end
    xlim(ax, [min(tick_t)-x_pad, max(tick_t)+x_pad]);

    h_early = gobjects(0); h_late = gobjects(0);
    for k = 1:n_bins
        vals_k = bins(k).(feat_field);
        vals_k = vals_k(~isnan(vals_k));
        n_k = length(vals_k);
        if n_k == 0, continue; end
        x_k = bins(k).time_h + (rand(n_k, 1) - 0.5) * jitter_w;
        if bins(k).is_early
            col = EARLY_COLOR;
        else
            col = LATE_COLOR;
        end
        hpt = scatter(ax, x_k, vals_k, 26, col, 'filled', 'MarkerFaceAlpha', 0.75);
        if bins(k).is_early && isempty(h_early), h_early = hpt; end
        if ~bins(k).is_early && isempty(h_late), h_late = hpt; end

        % Median marker
        med_k = median(vals_k);
        line(ax, [bins(k).time_h - jitter_w*0.65, bins(k).time_h + jitter_w*0.65], ...
             [med_k med_k], 'Color', 'k', 'LineWidth', 2);
    end

    ylim(ax, [y_data_bottom, y_top_total]);

    % --- Significance brackets: all pairwise (narrowest first, stacked
    %     low to high) then the pooled EARLY-vs-LATE contrast on top ---
    bracket_x1 = zeros(n_pairs, 1); bracket_x2 = zeros(n_pairs, 1);
    for pi = 1:n_pairs
        a = pairs(pi, 1); b = pairs(pi, 2);
        bracket_x1(pi) = min(bins(a).time_h, bins(b).time_h);
        bracket_x2(pi) = max(bins(a).time_h, bins(b).time_h);
    end
    [~, span_order] = sort(bracket_x2 - bracket_x1);
    for lvl = 1:n_pairs
        pi  = span_order(lvl);
        res = pw_res{pi};
        txt = sprintf('p=%s %s', fmt_p(res.p), strtrim(sig_stars(res.p)));
        y_lvl = y_data_top + lvl*step;
        draw_bracket(ax, bracket_x1(pi), bracket_x2(pi), y_lvl, tick_h, txt, 8.5, 1.0);
    end
    pooled_x1  = mean(tick_t(is_early));
    pooled_x2  = mean(tick_t(~is_early));
    pooled_txt = sprintf('EARLY vs LATE: p=%s %s', fmt_p(res_pooled.p), strtrim(sig_stars(res_pooled.p)));
    y_pooled   = y_data_top + (n_pairs+1)*step;
    draw_bracket(ax, pooled_x1, pooled_x2, y_pooled, tick_h*1.5, pooled_txt, 10, 1.8);

    % --- Tick labels: real time (h) + n colonies (this feature), one line ---
    xticks(ax, sort(tick_t));
    [~, ord] = sort(tick_t);
    xtl = cell(1, n_bins);
    for kk = 1:n_bins
        k = ord(kk);
        vals_k = bins(k).(feat_field);
        n_valid_k = sum(~isnan(vals_k));
        xtl{kk} = sprintf('%.1f h (n=%d)', bins(k).time_h, n_valid_k);
    end
    xticklabels(ax, xtl);

    ylabel(ax, y_unit_label, 'FontSize', 12, 'FontWeight', 'bold', 'Interpreter', 'tex');
    title(ax, sprintf('%s — Group %s — threshold %g px', feat_title, lbl, threshold), 'FontSize', 13);
    set(ax, 'Box', 'on', 'TickDir', 'in', 'FontSize', 11);
    grid(ax, 'on');

    leg_handles = gobjects(0);
    leg_labels  = {};
    if ~isempty(h_early), leg_handles(end+1) = h_early; leg_labels{end+1} = 'EARLY'; end
    if ~isempty(h_late),  leg_handles(end+1) = h_late;  leg_labels{end+1} = 'LATE';  end
    if ~isempty(leg_handles)
        legend(ax, leg_handles, leg_labels, 'Location', 'northeast', 'Box', 'off');
    end

    % --- EARLY / LATE group heading, spanning its two ticks, below the
    %     tick-label row (figure-normalized coords so it doesn't depend
    %     on the data y-scale) ---
    drawnow;
    xl_final = xlim(ax);
    frac_e = (pooled_x1 - xl_final(1)) / diff(xl_final);
    frac_l = (pooled_x2 - xl_final(1)) / diff(xl_final);
    pos    = ax.Position;   % [left bottom width height], normalized
    lbl_w  = 0.20; lbl_h = 0.055;
    fig_y  = max(0.01, pos(2) - GROUP_LABEL_DY);
    annotation(hfig, 'textbox', [pos(1)+frac_e*pos(3)-lbl_w/2, fig_y, lbl_w, lbl_h], ...
        'String', 'EARLY', 'EdgeColor', 'none', 'HorizontalAlignment', 'center', ...
        'FontWeight', 'bold', 'FontSize', 12, 'Color', EARLY_COLOR);
    annotation(hfig, 'textbox', [pos(1)+frac_l*pos(3)-lbl_w/2, fig_y, lbl_w, lbl_h], ...
        'String', 'LATE', 'EdgeColor', 'none', 'HorizontalAlignment', 'center', ...
        'FontWeight', 'bold', 'FontSize', 12, 'Color', LATE_COLOR);

    set(hfig, 'Visible', 'on');
    drawnow;

    timestamp_fig = datestr(now, 'yyyymmdd_HHMMSS');
    fig_name = sprintf('%s_EarlyLate_%s_%s', fig_tag, lbl, timestamp_fig);
    fig_path = fullfile(out_dir, fig_name);
    savefig(hfig, [fig_path '.fig']);
    exportgraphics(hfig, [fig_path '.pdf'], 'ContentType', 'vector');
    exportgraphics(hfig, [fig_path '.png'], 'Resolution', 300);
end


%% ================================================================
%  LOCAL FUNCTION: make_growth5h_figure
%  Paired scatter: for every bin with a usable +5h timepoint, the SAME
%  colonies' size at their own appearance timepoint next to their size
%  ~5h later, with a connecting line per colony. EARLY bins are grouped
%  together, LATE bins grouped together, with a gap between the clusters.
%  A thick dashed "main slope" line per cohort overlays the OLS fit
%  (fit_early / fit_late, from fit_growth_rate) evaluated at that
%  cohort's own mean before/after x-positions, with its equation, R^2
%  and Pearson p annotated in a fixed corner of the axes. A side panel
%  (steepness bar chart, ax2) compares the two slopes directly, with the
%  slope_cmp t-test on their difference (compare_slopes).
% ================================================================
function fig_path = make_growth5h_figure(out_dir, lbl, bins, res_growth_cmp, fit_early, fit_late, slope_cmp, EARLY_COLOR, LATE_COLOR)

    ok_bins = find([bins.tp5_ok]);
    n_ok    = length(ok_bins);
    is_early_ok = arrayfun(@(k) bins(k).is_early, ok_bins);
    n_early_ok  = sum(is_early_ok);
    n_late_ok   = n_ok - n_early_ok;

    % x-positions: each ok bin gets 2 slots (before, after); a gap opens
    % between the EARLY cluster and the LATE cluster.
    gap    = 1.0;
    x_pos  = zeros(n_ok, 2);
    cursor = 1;
    last_is_early = [];
    for m = 1:n_ok
        k = ok_bins(m);
        if ~isempty(last_is_early) && last_is_early ~= bins(k).is_early
            cursor = cursor + gap;
        end
        x_pos(m, :) = [cursor, cursor + 0.7];
        cursor = cursor + 0.7 + 0.5;
        last_is_early = bins(k).is_early;
    end

    % Extra width reserved for the slope-steepness side panel (ax2).
    side_w_in = 2.2;
    fig_w = max(6.4, 1.3*n_ok + 2) + side_w_in;
    hfig = figure('Name', sprintf('Size growth over +5h — %s', lbl), ...
                  'Color', 'w', 'Units', 'inches', ...
                  'Position', [1 1 fig_w 5.2], ...
                  'PaperUnits', 'inches', ...
                  'PaperPosition', [0 0 fig_w 5.2], ...
                  'PaperSize', [fig_w 5.2], 'Visible', 'off');
    ax = axes('Position', [0.09 0.30 0.62 0.56]);
    hold(ax, 'on');

    all_vals = [];
    for m = 1:n_ok
        k = ok_bins(m);
        all_vals = [all_vals; bins(k).size_before(:); bins(k).size_after(:)]; %#ok<AGROW>
    end
    all_vals = all_vals(~isnan(all_vals));
    if isempty(all_vals)
        y_min = 0; y_max = 1;
    else
        y_min = min(all_vals); y_max = max(all_vals);
    end
    if y_max <= y_min, y_max = y_min + 1; end
    yr = y_max - y_min;

    jitter = 0.10;
    xtl = cell(1, 2*n_ok);
    for m = 1:n_ok
        k = ok_bins(m);
        if bins(k).is_early
            col = EARLY_COLOR; line_col = 0.55*EARLY_COLOR + 0.45*[1 1 1];
        else
            col = LATE_COLOR;  line_col = 0.55*LATE_COLOR  + 0.45*[1 1 1];
        end
        before = bins(k).size_before(:);
        after  = bins(k).size_after(:);
        n_k    = length(before);
        xb = x_pos(m,1) + (rand(n_k,1) - 0.5) * 2*jitter;
        xa = x_pos(m,2) + (rand(n_k,1) - 0.5) * 2*jitter;

        for i = 1:n_k
            if isnan(before(i)) || isnan(after(i)), continue; end
            plot(ax, [xb(i) xa(i)], [before(i) after(i)], 'Color', line_col, 'LineWidth',0.6);
        end
        valid_b = ~isnan(before); valid_a = ~isnan(after);
        scatter(ax, xb(valid_b), before(valid_b), 20, col, 'filled', ...
            'MarkerFaceAlpha',0.6, 'MarkerEdgeColor','none');
        scatter(ax, xa(valid_a), after(valid_a), 20, col, 'filled', ...
            'MarkerFaceAlpha',0.6, 'MarkerEdgeColor','none');

        if any(valid_b)
            med_b = median(before(valid_b));
            plot(ax, [x_pos(m,1)-0.18 x_pos(m,1)+0.18], [med_b med_b], 'Color',[0.1 0.1 0.1], 'LineWidth',2);
        end
        if any(valid_a)
            med_a = median(after(valid_a));
            plot(ax, [x_pos(m,2)-0.18 x_pos(m,2)+0.18], [med_a med_a], 'Color',[0.1 0.1 0.1], 'LineWidth',2);
        end

        % Paired-test bracket for this bin (before vs after, same colonies)
        if isfield(bins(k), 'res_pair') && ~isempty(bins(k).res_pair.test_name)
            rp   = bins(k).res_pair;
            ybk  = y_max + 0.05*yr;
            txt  = sprintf('%s: p=%s %s', rp.test_name, fmt_p(rp.p), strtrim(sig_stars(rp.p)));
            draw_bracket(ax, x_pos(m,1), x_pos(m,2), ybk, 0.015*yr, txt, 8, 1.0);
        end

        % Primary tick label is just the real time (h) — "X" / "X+5" for
        % EARLY, "Y" / "Y+5" for LATE — matching Figures 1-2's convention
        % of showing real time on the axis. The EARLY/LATE heading itself
        % is drawn below the axis (see after title(), below), not baked
        % into this tick label. If more than one bin shares a phase, the
        % "before" tick also gets that bin's own label so they can still
        % be told apart (still redundant with the BIN SUMMARY log table).
        if (bins(k).is_early && n_early_ok > 1) || (~bins(k).is_early && n_late_ok > 1)
            xtl{2*m-1} = sprintf('%s\n%.1fh', bins(k).label, bins(k).time_h);
        else
            xtl{2*m-1} = sprintf('%.1fh', bins(k).time_h);
        end
        xtl{2*m} = sprintf('%.1fh', bins(k).time_h5);
    end

    xlim(ax, [0.3, cursor - 0.3]);
    ylim(ax, [y_min - 0.06*yr, y_max + 0.30*yr]);
    set(ax, 'XTick', reshape(x_pos', 1, []), 'XTickLabel', xtl, ...
        'FontSize', 9, 'TickDir', 'in', 'Box', 'on');
    ylabel(ax, 'Colony size (px)', 'FontSize', 13, 'FontWeight', 'bold');
    title(ax, sprintf('Colony size: appearance vs +5h — Group %s', lbl), ...
        'FontSize', 13, 'FontWeight', 'bold');

    % --- EARLY / LATE group heading, spanning its own cluster of ticks,
    %     below the tick-label row (figure-normalized coords) — same
    %     convention as make_bin_comparison_figure's heading. ---
    drawnow;
    xl_final = xlim(ax);
    pos      = ax.Position;
    hd_w = 0.16; hd_h = 0.05; hd_dy = 0.10;
    hd_y = max(0.01, pos(2) - hd_dy);
    if n_early_ok > 0
        early_span = x_pos(is_early_ok, :);
        frac_e = (mean(early_span(:)) - xl_final(1)) / diff(xl_final);
        annotation(hfig, 'textbox', [pos(1)+frac_e*pos(3)-hd_w/2, hd_y, hd_w, hd_h], ...
            'String', 'EARLY', 'EdgeColor', 'none', 'HorizontalAlignment', 'center', ...
            'FontWeight', 'bold', 'FontSize', 12, 'Color', EARLY_COLOR);
    end
    if n_late_ok > 0
        late_span = x_pos(~is_early_ok, :);
        frac_l = (mean(late_span(:)) - xl_final(1)) / diff(xl_final);
        annotation(hfig, 'textbox', [pos(1)+frac_l*pos(3)-hd_w/2, hd_y, hd_w, hd_h], ...
            'String', 'LATE', 'EdgeColor', 'none', 'HorizontalAlignment', 'center', ...
            'FontWeight', 'bold', 'FontSize', 12, 'Color', LATE_COLOR);
    end

    % Independent-samples bracket: pooled EARLY size-gain vs pooled LATE
    early_xs = []; late_xs = [];
    for m = 1:n_ok
        k = ok_bins(m);
        if bins(k).is_early
            early_xs = [early_xs, x_pos(m,:)]; %#ok<AGROW>
        else
            late_xs = [late_xs, x_pos(m,:)]; %#ok<AGROW>
        end
    end
    if ~isempty(early_xs) && ~isempty(late_xs)
        early_mid = mean(early_xs);
        late_mid  = mean(late_xs);
        yb2  = y_max + 0.18*yr;
        txt2 = sprintf('Size GAIN, EARLY vs LATE — %s: p=%s %s', blank_if_empty(res_growth_cmp.test_name), ...
            fmt_p(res_growth_cmp.p), strtrim(sig_stars(res_growth_cmp.p)));
        draw_bracket(ax, early_mid, late_mid, yb2, 0.02*yr, txt2, 9, 1.5);
    end

    % --- Main slope: pooled OLS fit line per cohort, drawn from that
    %     cohort's mean "before" x-position to its mean "after" x-position
    %     (categorical layout, so real elapsed time only enters via the
    %     fit itself, evaluated at rel_t=0 and rel_t=mean dt for the fit
    %     line's two endpoints) ---
    draw_growth_fit_line(ax, x_pos, ok_bins, bins, true,  fit_early, EARLY_COLOR);
    draw_growth_fit_line(ax, x_pos, ok_bins, bins, false, fit_late,  LATE_COLOR);

    ann_early = growth_fit_annotation_str('EARLY', fit_early);
    ann_late  = growth_fit_annotation_str('LATE',  fit_late);
    text(ax, 0.02, 0.97, ann_early, 'Units','normalized', 'FontSize',8, ...
        'VerticalAlignment','top', 'HorizontalAlignment','left', ...
        'Color', EARLY_COLOR, 'BackgroundColor','w', 'Margin',3);
    text(ax, 0.98, 0.97, ann_late, 'Units','normalized', 'FontSize',8, ...
        'VerticalAlignment','top', 'HorizontalAlignment','right', ...
        'Color', LATE_COLOR, 'BackgroundColor','w', 'Margin',3);

    % ------------------------------------------------------------------
    %  Side panel: is the EARLY slope steeper/shallower than the LATE
    %  slope? A 2-bar comparison (slope +/- its own SE) with the
    %  slope_cmp t-test (compare_slopes) bracketed between them — reading
    %  the two "main slope" numbers off the main panel only says which is
    %  numerically bigger, not whether that difference is more than noise.
    % ------------------------------------------------------------------
    ax2 = axes('Position', [0.78 0.30 0.19 0.56]);
    hold(ax2, 'on');

    slopes = [fit_early.slope, fit_late.slope];
    ses    = [fit_early.se_slope, fit_late.se_slope];
    cols   = {EARLY_COLOR, LATE_COLOR};
    for i = 1:2
        if isnan(slopes(i)), continue; end
        bar(ax2, i, slopes(i), 0.6, 'FaceColor', cols{i}, 'EdgeColor', 'none');
        if ~isnan(ses(i))
            errorbar(ax2, i, slopes(i), ses(i), 'Color', 'k', 'LineWidth', 1.2, ...
                'CapSize', 8, 'LineStyle', 'none');
        end
    end
    xlim(ax2, [0.3 2.7]);
    set(ax2, 'XTick', [1 2], 'XTickLabel', {'EARLY', 'LATE'}, ...
        'FontSize', 9, 'TickDir', 'in', 'Box', 'on');
    ylabel(ax2, 'Slope (px/h)', 'FontSize', 10, 'FontWeight', 'bold');
    title(ax2, 'Steepness', 'FontSize', 10, 'FontWeight', 'bold');
    grid(ax2, 'on');

    finite_slopes = slopes(~isnan(slopes));
    if ~isempty(finite_slopes) && ~isnan(slope_cmp.p)
        yl2   = ylim(ax2);
        y_bar = max([finite_slopes(:); 0]) + 0.12*range(yl2);
        txt_s = sprintf('p=%s %s', fmt_p(slope_cmp.p), strtrim(sig_stars(slope_cmp.p)));
        draw_bracket(ax2, 1, 2, y_bar, 0.03*range(yl2), txt_s, 8, 1.0);
        ylim(ax2, [min(0, yl2(1)), max(yl2(2), y_bar + 0.15*range(yl2))]);
    end

    set(hfig, 'Visible', 'on');
    drawnow;

    timestamp_fig = datestr(now, 'yyyymmdd_HHMMSS');
    fig_name = sprintf('SizeGrowth5h_EarlyLate_%s_%s', lbl, timestamp_fig);
    fig_path = fullfile(out_dir, fig_name);
    savefig(hfig, [fig_path '.fig']);
    exportgraphics(hfig, [fig_path '.pdf'], 'ContentType', 'vector');
    exportgraphics(hfig, [fig_path '.png'], 'Resolution', 300);
end


%% ================================================================
%  LOCAL FUNCTION: draw_growth_fit_line
%  Draws the pooled OLS fit (fit_growth_rate) for one cohort (EARLY or
%  LATE) as a thick dashed line from that cohort's mean "before" x
%  position to its mean "after" x position — the two endpoints' y-values
%  come from evaluating the fit at rel_t=0 and rel_t=mean(actual +5h
%  offset), NOT from the raw medians, since this line represents the
%  population-level trend, not a per-bin summary.
% ================================================================
function draw_growth_fit_line(ax, x_pos, ok_bins, bins, want_early, fit, color)
    if isnan(fit.slope), return; end

    before_xs = []; after_xs = []; dts = [];
    for m = 1:length(ok_bins)
        k = ok_bins(m);
        if bins(k).is_early ~= want_early, continue; end
        before_xs(end+1) = x_pos(m, 1); %#ok<AGROW>
        after_xs(end+1)  = x_pos(m, 2); %#ok<AGROW>
        dts(end+1)       = bins(k).time_h5 - bins(k).time_h; %#ok<AGROW>
    end
    if isempty(before_xs), return; end

    mean_bx = mean(before_xs);
    mean_ax = mean(after_xs);
    mean_dt = mean(dts);

    y0 = fit.intercept;
    y1 = fit.intercept + fit.slope * mean_dt;
    plot(ax, [mean_bx mean_ax], [y0 y1], '--', 'Color', color, 'LineWidth', 3);
end


%% ================================================================
%  LOCAL FUNCTION: growth_fit_annotation_str
% ================================================================
function s = growth_fit_annotation_str(label, fit)
    if isnan(fit.slope)
        s = sprintf('%s: insufficient data for a fit (n=%d)', label, fit.n);
        return;
    end
    s = sprintf('%s slope fit (n=%d)\n%s\nR^2 = %.3f\nr = %.3f, p = %s', ...
        label, fit.n, fit.eq_str, fit.R2, fit.r, fmt_p(fit.p));
end


%% ================================================================
%  LOCAL FUNCTION: log_slope_comparison / log_paired_slopes
% ================================================================
function log_slope_comparison(tlog, window_label, cmp)
    tlog('  %s:\n', window_label);
    if isnan(cmp.t)
        tlog('    Insufficient data in one or both fits for a slope comparison.\n\n');
        return;
    end
    tlog('    Delta slope (EARLY - LATE) = %.4g px/h\n', cmp.diff);
    tlog('    SE(delta)                  = %.4g px/h\n', cmp.se);
    tlog('    t(%d) = %.4f,  p = %s %s\n\n', cmp.df, cmp.t, fmt_p(cmp.p), strtrim(sig_stars(cmp.p)));
end

function log_paired_slopes(tlog, label, res, pre_lbl, post_lbl)
    tlog('  %s:\n', label);
    if isempty(res.test_name)
        tlog('    Insufficient paired data (n=%d colonies; need >=3).\n\n', res.n);
        return;
    end
    if res.normal_diff, nrm = 'normal'; else, nrm = 'not normal'; end
    tlog('    n colonies              = %d\n', res.n);
    tlog('    Median slope %s -> X = %.4g px/h\n', pre_lbl, res.median_before);
    tlog('    Median slope X -> %s = %.4g px/h\n', post_lbl, res.median_after);
    tlog('    Median (post - pre)     = %.4g px/h\n', res.median_diff);
    tlog('    Normality of differences (Lilliefors): %s\n', nrm);
    tlog('    %s: p=%s %s\n\n', res.test_name, fmt_p(res.p), strtrim(sig_stars(res.p)));
end


%% ================================================================
%  LOCAL FUNCTION: pool_three_point_cohort
%  Pools every EARLY (or LATE) bin's colonies that have a valid size at
%  ALL THREE timepoints (X-L, X, X+A). For each window, x = elapsed
%  time since the start of that window (0 at the window start, the bin's
%  actual window length at the end) and y = size, ready for
%  fit_growth_rate; s_pre / s_post are each colony's own slope (px/h)
%  per window, for the paired pre-vs-post comparison.
% ================================================================
function P = pool_three_point_cohort(bins, want_early)
    P.rt_pre = []; P.y_pre = []; P.rt_post = []; P.y_post = [];
    P.s_pre  = []; P.s_post = []; P.n = 0;
    for k = 1:numel(bins)
        if ~bins(k).tp3_ok || bins(k).is_early ~= want_early, continue; end
        pre = bins(k).size_pre(:);
        bef = bins(k).size_x(:);
        aft = bins(k).size_post(:);
        v   = ~isnan(pre) & ~isnan(bef) & ~isnan(aft);
        if ~any(v), continue; end
        pre = pre(v); bef = bef(v); aft = aft(v);
        n_v     = sum(v);
        dt_pre  = bins(k).time_h      - bins(k).time_h_pre;
        dt_post = bins(k).time_h_post - bins(k).time_h;

        P.rt_pre  = [P.rt_pre;  zeros(n_v,1); repmat(dt_pre,  n_v, 1)];
        P.y_pre   = [P.y_pre;   pre; bef];                              
        P.rt_post = [P.rt_post; zeros(n_v,1); repmat(dt_post, n_v, 1)]; 
        P.y_post  = [P.y_post;  bef; aft];                              
        P.s_pre   = [P.s_pre;   (bef - pre) / dt_pre];                  
        P.s_post  = [P.s_post;  (aft - bef) / dt_post];                
        P.n       = P.n + n_v;
    end
end


%% ================================================================
%  LOCAL FUNCTION: make_growth3pt_figure
%  Paired scatter at THREE timepoints per bin (X-L, X, X+A; labels from
%  S.pre_lbl / S.post_lbl) for the
%  same colonies, EARLY cluster then LATE cluster, shared y-axis, one
%  polyline per colony. Per cohort two main-slope lines are overlaid
%  (dotted = X-L -> X, dashed = X -> X+A) with equations/R^2 in the
%  corners, and a side panel compares the two cohorts' slopes within
%  each window (compare_slopes). The per-cohort paired pre-vs-post
%  per-colony-slope test is bracketed above each cohort.
% ================================================================
function fig_path = make_growth3pt_figure(out_dir, lbl, bins, S, EARLY_COLOR, LATE_COLOR)

    ok_bins     = find([bins.tp3_ok]);
    n_ok        = numel(ok_bins);
    is_early_ok = arrayfun(@(k) bins(k).is_early, ok_bins);
    n_early_ok  = sum(is_early_ok);
    n_late_ok   = n_ok - n_early_ok;

    gap    = 1.0;
    x_pos  = zeros(n_ok, 3);
    cursor = 1;
    last_is_early = [];
    for m = 1:n_ok
        k = ok_bins(m);
        if ~isempty(last_is_early) && last_is_early ~= bins(k).is_early
            cursor = cursor + gap;
        end
        x_pos(m, :) = cursor + [0 0.7 1.4];
        cursor = cursor + 1.4 + 0.5;
        last_is_early = bins(k).is_early;
    end

    side_w_in = 3.2;
    fig_w = max(7.2, 1.8*n_ok + 2.5) + side_w_in;
    fig_h = 5.4;
    hfig = figure('Name', sprintf('Size growth %s, X, %s — %s', S.pre_lbl, S.post_lbl, lbl), ...
                  'Color', 'w', 'Units', 'inches', ...
                  'Position', [1 1 fig_w fig_h], ...
                  'PaperUnits', 'inches', ...
                  'PaperPosition', [0 0 fig_w fig_h], ...
                  'PaperSize', [fig_w fig_h], 'Visible', 'off');
    ax = axes('Position', [0.08 0.30 0.62 0.56]);
    hold(ax, 'on');

    all_vals = [];
    for m = 1:n_ok
        k = ok_bins(m);
        all_vals = [all_vals; bins(k).size_pre(:); bins(k).size_x(:); bins(k).size_post(:)];
    end
    all_vals = all_vals(~isnan(all_vals));
    if isempty(all_vals)
        y_min = 0; y_max = 1;
    else
        y_min = min(all_vals); y_max = max(all_vals);
    end
    if y_max <= y_min, y_max = y_min + 1; end
    yr = y_max - y_min;

    jitter = 0.10;
    xtl    = cell(1, 3*n_ok);
    for m = 1:n_ok
        k = ok_bins(m);
        if bins(k).is_early
            col = EARLY_COLOR;
        else
            col = LATE_COLOR;
        end
        line_col = 0.55*col + 0.45*[1 1 1];

        pre = bins(k).size_pre(:);
        bef = bins(k).size_x(:);
        aft = bins(k).size_post(:);
        v   = ~isnan(pre) & ~isnan(bef) & ~isnan(aft);
        pre = pre(v); bef = bef(v); aft = aft(v);
        n_k = numel(pre);
        if n_k > 0
            xp = x_pos(m,1) + (rand(n_k,1) - 0.5) * 2*jitter;
            xb = x_pos(m,2) + (rand(n_k,1) - 0.5) * 2*jitter;
            xa = x_pos(m,3) + (rand(n_k,1) - 0.5) * 2*jitter;
            plot(ax, [xp'; xb'; xa'], [pre'; bef'; aft'], 'Color', line_col, 'LineWidth', 0.6);
            scatter(ax, xp, pre, 20, col, 'filled', 'MarkerFaceAlpha', 0.6, 'MarkerEdgeColor', 'none');
            scatter(ax, xb, bef, 20, col, 'filled', 'MarkerFaceAlpha', 0.6, 'MarkerEdgeColor', 'none');
            scatter(ax, xa, aft, 20, col, 'filled', 'MarkerFaceAlpha', 0.6, 'MarkerEdgeColor', 'none');
            meds = [median(pre) median(bef) median(aft)];
            for j = 1:3
                plot(ax, [x_pos(m,j)-0.18 x_pos(m,j)+0.18], [meds(j) meds(j)], ...
                    'Color', [0.1 0.1 0.1], 'LineWidth', 2);
            end
        end

        if (bins(k).is_early && n_early_ok > 1) || (~bins(k).is_early && n_late_ok > 1)
            xtl{3*m-2} = sprintf('%s\n%.1fh', bins(k).label, bins(k).time_h_pre);
        else
            xtl{3*m-2} = sprintf('%.1fh', bins(k).time_h_pre);
        end
        xtl{3*m-1} = sprintf('%.1fh', bins(k).time_h);
        xtl{3*m}   = sprintf('%.1fh', bins(k).time_h_post);
    end

    draw_cohort_window_lines(ax, x_pos, ok_bins, bins, true,  S.fit_pre_early, S.fit_post_early, EARLY_COLOR);
    draw_cohort_window_lines(ax, x_pos, ok_bins, bins, false, S.fit_pre_late,  S.fit_post_late,  LATE_COLOR);

    xlim(ax, [0.3, cursor - 0.3]);
    ylim(ax, [y_min - 0.06*yr, y_max + 0.95*yr]);
    set(ax, 'XTick', reshape(x_pos', 1, []), 'XTickLabel', xtl, ...
        'FontSize', 9, 'TickDir', 'in', 'Box', 'on');
    ylabel(ax, 'Colony size (px)', 'FontSize', 13, 'FontWeight', 'bold');
    title(ax, sprintf('Colony size: %s, X, %s — Group %s', S.pre_lbl, S.post_lbl, lbl), ...
        'FontSize', 13, 'FontWeight', 'bold');

    % Paired pre-vs-post slope test, bracketed from the middle of the
    % pre-window to the middle of the post-window of each cohort.
    cohorts = {true, S.acc_early; false, S.acc_late};
    for c = 1:2
        rows = (is_early_ok == cohorts{c,1});
        res  = cohorts{c,2};
        if ~any(rows) || isempty(res.test_name), continue; end
        x1 = mean(mean(x_pos(rows, 1:2)));
        x2 = mean(mean(x_pos(rows, 2:3)));
        txt = sprintf('pre vs post slope (paired): p=%s %s', fmt_p(res.p), strtrim(sig_stars(res.p)));
        draw_bracket(ax, x1, x2, y_max + 0.06*yr, 0.015*yr, txt, 8, 1.2);
    end

    text(ax, 0.02, 0.97, three_pt_annotation_str('EARLY', S.fit_pre_early, S.fit_post_early, S.pre_lbl, S.post_lbl), ...
        'Units', 'normalized', 'FontSize', 8, 'VerticalAlignment', 'top', ...
        'HorizontalAlignment', 'left', 'Color', EARLY_COLOR, 'BackgroundColor', 'w', 'Margin', 3);
    text(ax, 0.98, 0.97, three_pt_annotation_str('LATE', S.fit_pre_late, S.fit_post_late, S.pre_lbl, S.post_lbl), ...
        'Units', 'normalized', 'FontSize', 8, 'VerticalAlignment', 'top', ...
        'HorizontalAlignment', 'right', 'Color', LATE_COLOR, 'BackgroundColor', 'w', 'Margin', 3);

    % EARLY / LATE heading below the tick labels (same convention as above)
    drawnow;
    xl_final = xlim(ax);
    pos      = ax.Position;
    hd_w = 0.16; hd_h = 0.05; hd_y = max(0.01, pos(2) - 0.10);
    if n_early_ok > 0
        span   = x_pos(is_early_ok, :);
        frac_e = (mean(span(:)) - xl_final(1)) / diff(xl_final);
        annotation(hfig, 'textbox', [pos(1)+frac_e*pos(3)-hd_w/2, hd_y, hd_w, hd_h], ...
            'String', 'EARLY', 'EdgeColor', 'none', 'HorizontalAlignment', 'center', ...
            'FontWeight', 'bold', 'FontSize', 12, 'Color', EARLY_COLOR);
    end
    if n_late_ok > 0
        span   = x_pos(~is_early_ok, :);
        frac_l = (mean(span(:)) - xl_final(1)) / diff(xl_final);
        annotation(hfig, 'textbox', [pos(1)+frac_l*pos(3)-hd_w/2, hd_y, hd_w, hd_h], ...
            'String', 'LATE', 'EdgeColor', 'none', 'HorizontalAlignment', 'center', ...
            'FontWeight', 'bold', 'FontSize', 12, 'Color', LATE_COLOR);
    end

    % ------------------------------------------------------------------
    %  Side panel: both slopes of each cohort as bars — per cohort a light
    %  bar (X-L -> X) and a dark bar (X -> X+A), slope +/- SE. Brackets:
    %   level 1  pre vs post WITHIN each cohort (paired per-colony slopes)
    %   level 2  EARLY vs LATE, X-L -> X window (compare_slopes)
    %   level 3  EARLY vs LATE, X -> X+A window (compare_slopes)
    % ------------------------------------------------------------------
    ax2 = axes('Position', [0.77 0.30 0.20 0.56]);
    hold(ax2, 'on');
    sl = [S.fit_pre_early.slope,    S.fit_post_early.slope; ...
          S.fit_pre_late.slope,     S.fit_post_late.slope];
    se = [S.fit_pre_early.se_slope, S.fit_post_early.se_slope; ...
          S.fit_pre_late.se_slope,  S.fit_post_late.se_slope];
    se(isnan(se)) = 0;
    cols  = {EARLY_COLOR, LATE_COLOR};
    off   = 0.2;
    bar_x = zeros(2, 2);
    for c = 1:2
        for w = 1:2
            xx = c + (2*w - 3) * off;       % w=1 (pre) left, w=2 (post) right
            bar_x(c, w) = xx;
            if isnan(sl(c, w)), continue; end
            if w == 1
                fc = 0.5*cols{c} + 0.5*[1 1 1];
            else
                fc = cols{c};
            end
            bar(ax2, xx, sl(c, w), 0.34, 'FaceColor', fc, 'EdgeColor', cols{c}, 'LineWidth', 1);
            errorbar(ax2, xx, sl(c, w), se(c, w), 'Color', 'k', 'LineWidth', 1.2, ...
                'CapSize', 6, 'LineStyle', 'none');
        end
    end
    xlim(ax2, [0.4 2.6]);
    set(ax2, 'XTick', [1 2], 'XTickLabel', {'EARLY', 'LATE'}, ...
        'FontSize', 9, 'TickDir', 'in', 'Box', 'on');
    ylabel(ax2, 'Slope (px/h)', 'FontSize', 10, 'FontWeight', 'bold');
    xlabel(ax2, sprintf('light: %s->X    dark: X->%s', S.pre_lbl, S.post_lbl), 'FontSize', 8);
    title(ax2, 'Steepness', 'FontSize', 10, 'FontWeight', 'bold');
    grid(ax2, 'on');

    fin = ~isnan(sl);
    if any(fin(:))
        top = max(sl(fin) + se(fin));
        bot = min([0; sl(fin)]);
        if top <= 0, top = 1; end
        stp  = 0.10 * top;
        base = top + 0.06*top;
        y_top2 = base + 3*stp + 0.05*top;
        ylim(ax2, [bot * 1.2, y_top2]);

        accs = {S.acc_early, S.acc_late};
        for c = 1:2
            if isempty(accs{c}.test_name), continue; end
            txt_a = sprintf('pre vs post: p=%s %s', fmt_p(accs{c}.p), strtrim(sig_stars(accs{c}.p)));
            draw_bracket(ax2, bar_x(c,1), bar_x(c,2), base, 0.015*top, txt_a, 7, 1.0);
        end
        if ~isnan(S.cmp_pre.p)
            txt_p = sprintf('pre: p=%s %s', fmt_p(S.cmp_pre.p), strtrim(sig_stars(S.cmp_pre.p)));
            draw_bracket(ax2, bar_x(1,1), bar_x(2,1), base + stp, 0.015*top, txt_p, 7, 1.0);
        end
        if ~isnan(S.cmp_post.p)
            txt_q = sprintf('post: p=%s %s', fmt_p(S.cmp_post.p), strtrim(sig_stars(S.cmp_post.p)));
            draw_bracket(ax2, bar_x(1,2), bar_x(2,2), base + 2*stp, 0.015*top, txt_q, 7, 1.0);
        end
    end

    set(hfig, 'Visible', 'on');
    drawnow;

    timestamp_fig = datestr(now, 'yyyymmdd_HHMMSS');
    fig_name = sprintf('SizeGrowth3pt_Xminus%dh_Xplus%dh_EarlyLate_%s_%s', S.lookback, S.lookahead, lbl, timestamp_fig);
    fig_path = fullfile(out_dir, fig_name);
    savefig(hfig, [fig_path '.fig']);
    exportgraphics(hfig, [fig_path '.pdf'], 'ContentType', 'vector');
    exportgraphics(hfig, [fig_path '.png'], 'Resolution', 300);
end


%% ================================================================
%  LOCAL FUNCTION: draw_cohort_window_lines
%  For one cohort, draws the pre-window fit (dotted, x-slot 1 -> 2) and
%  post-window fit (dashed, x-slot 2 -> 3). Each line runs between the
%  cohort's mean slot positions; its y-values come from evaluating the
%  fit at 0 and at the cohort's mean window length.
% ================================================================
function draw_cohort_window_lines(ax, x_pos, ok_bins, bins, want_early, fit_pre, fit_post, color)
    xs1 = []; xs2 = []; xs3 = []; dt_pre = []; dt_post = [];
    for m = 1:numel(ok_bins)
        k = ok_bins(m);
        if bins(k).is_early ~= want_early, continue; end
        xs1(end+1)    = x_pos(m, 1);                       
        xs2(end+1)    = x_pos(m, 2);                       
        xs3(end+1)    = x_pos(m, 3);                        
        dt_pre(end+1)  = bins(k).time_h  - bins(k).time_h_pre; 
        dt_post(end+1) = bins(k).time_h_post - bins(k).time_h;
    end
    if isempty(xs1), return; end
    if ~isnan(fit_pre.slope)
        plot(ax, [mean(xs1) mean(xs2)], ...
            [fit_pre.intercept, fit_pre.intercept + fit_pre.slope * mean(dt_pre)], ...
            ':', 'Color', color, 'LineWidth', 3.5);
    end
    if ~isnan(fit_post.slope)
        plot(ax, [mean(xs2) mean(xs3)], ...
            [fit_post.intercept, fit_post.intercept + fit_post.slope * mean(dt_post)], ...
            '--', 'Color', color, 'LineWidth', 3);
    end
end


%% ================================================================
%  LOCAL FUNCTION: three_pt_annotation_str
% ================================================================
function s = three_pt_annotation_str(label, fit_pre, fit_post, pre_lbl, post_lbl)
    s = sprintf('%s\n%s\n%s', label, ...
        window_fit_str(sprintf('%s->X  (dotted)', pre_lbl), fit_pre), ...
        window_fit_str(sprintf('X->%s (dashed)', post_lbl), fit_post));
end

function s = window_fit_str(name, fit)
    if isnan(fit.slope)
        s = sprintf('%s: insufficient data (n=%d)', name, fit.n);
        return;
    end
    s = sprintf('%s n=%d\n  %s\n  R^2=%.3f, p=%s', name, fit.n, fit.eq_str, fit.R2, fmt_p(fit.p));
end


%% ================================================================
%  LOCAL FUNCTION: plate_breakdown
%  Maps a list of within-group colony indices to their plate of origin
%  using grp.plate_seg_start/plate_seg_end/plate_indices, and returns
%  the distinct plate count plus a per-plate colony tally.
% ================================================================
function [n_plates, plate_ids, plate_counts] = plate_breakdown(grp, colony_idx)
    if isempty(colony_idx) || ~isfield(grp, 'plate_seg_start') || isempty(grp.plate_seg_start)
        n_plates = 0; plate_ids = []; plate_counts = [];
        return;
    end
    seg_start = grp.plate_seg_start;
    seg_end   = grp.plate_seg_end;
    plates    = grp.plate_indices;

    plate_of = nan(length(colony_idx), 1);
    for m = 1:length(colony_idx)
        ci = colony_idx(m);
        k = find(ci >= seg_start & ci <= seg_end, 1);
        if ~isempty(k)
            plate_of(m) = plates(k);
        end
    end
    plate_ids    = unique(plate_of(~isnan(plate_of)));
    plate_counts = arrayfun(@(x) sum(plate_of == x), plate_ids);
    n_plates     = length(plate_ids);
end


function s = plate_comp_string(plate_ids, plate_counts)
    if isempty(plate_ids)
        s = '(none)';
        return;
    end
    parts = cell(1, length(plate_ids));
    for i = 1:length(plate_ids)
        parts{i} = sprintf('%d:%d', plate_ids(i), plate_counts(i));
    end
    s = strjoin(parts, ', ');
end


function s = blank_if_empty(s)
    if isempty(s), s = 'n/a'; end
end


%% ================================================================
%  LOCAL FUNCTION: two_sample_test
%  Independent-samples comparison — mirrors preprocess_pipeline_data's
%  ht.pairwise test selection (Lilliefors, decided PER PAIR from these
%  two specific samples -> Welch ttest2 if both normal, else ranksum),
%  plus a two-sample KS test (kstest2) for whole-distribution shape
%  (STATISTICS_NOTES.md decision log, 2026-09-29 / 2026-09-30). Returns
%  NaN/'' fields if either sample has fewer than 3 valid values.
% ================================================================
function res = two_sample_test(va, vb)
    va = va(:); va = va(~isnan(va));
    vb = vb(:); vb = vb(~isnan(vb));

    res.n_a       = length(va);
    res.n_b       = length(vb);
    res.median_a  = NaN;
    res.median_b  = NaN;
    res.p         = NaN;
    res.stat      = NaN;
    res.test_name = '';
    res.normal_a  = false;
    res.normal_b  = false;
    res.ks_D      = NaN;
    res.ks_p      = NaN;

    if res.n_a < 3 || res.n_b < 3
        return;
    end

    res.median_a = median(va);
    res.median_b = median(vb);
    res.normal_a = local_is_normal(va);
    res.normal_b = local_is_normal(vb);

    if res.normal_a && res.normal_b
        [~, pv, ~, stats] = ttest2(va, vb, 'Vartype', 'unequal');
        res.p         = pv;
        res.stat      = stats.tstat;
        res.test_name = 't-test (Welch)';
    else
        [pv, ~, stats] = ranksum(va, vb);
        res.p = pv;
        if isfield(stats, 'zval')
            res.stat = stats.zval;
        else
            res.stat = stats.ranksum;
        end
        res.test_name = 'Wilcoxon rank-sum';
    end

    try
        [~, ks_pv] = kstest2(va, vb);
        res.ks_p = ks_pv;
        all_x = sort([va; vb]);
        cdf_a = arrayfun(@(x) mean(va <= x), all_x);
        cdf_b = arrayfun(@(x) mean(vb <= x), all_x);
        res.ks_D = max(abs(cdf_a - cdf_b));
    catch
        % leave ks_D / ks_p as NaN
    end
end


%% ================================================================
%  LOCAL FUNCTION: compare_paired
%  Matched-pairs comparison (SAME colonies measured twice — size at an
%  appearance timepoint vs size 5h later). Paired analogue of
%  two_sample_test: Lilliefors on the paired differences -> paired
%  ttest, or Wilcoxon signrank otherwise. No KS test here — kstest2
%  compares independent samples and would discard the pairing.
% ================================================================
function res = compare_paired(before, after)
    before = before(:);
    after  = after(:);
    valid  = ~isnan(before) & ~isnan(after);
    before = before(valid);
    after  = after(valid);
    d      = after - before;

    res.n             = length(d);
    res.median_before = NaN;
    res.median_after  = NaN;
    res.median_diff   = NaN;
    res.p             = NaN;
    res.stat          = NaN;
    res.test_name     = '';
    res.normal_diff   = false;

    if res.n < 3
        return;
    end

    res.median_before = median(before);
    res.median_after  = median(after);
    res.median_diff   = median(d);
    res.normal_diff   = local_is_normal(d);

    if res.normal_diff
        [~, pv, ~, st] = ttest(before, after);
        res.p         = pv;
        res.stat      = st.tstat;
        res.test_name = 't-test (paired)';
    else
        [pv, ~, st] = signrank(before, after);
        res.p = pv;
        if isfield(st, 'zval')
            res.stat = st.zval;
        else
            res.stat = st.signedrank;
        end
        res.test_name = 'Wilcoxon signed-rank';
    end
end


%% ================================================================
%  LOCAL FUNCTION: local_is_normal
%  Lilliefors normality test. h==1 rejects normality at alpha=0.05;
%  treated conservatively as "not normal" on error or too little data.
% ================================================================
function tf = local_is_normal(x)
    tf = false;
    x = x(~isnan(x));
    if length(x) < 4 || std(x) == 0
        return;
    end
    try
        h  = lillietest(x);
        tf = (h == 0);
    catch
        tf = false;
    end
end


%% ================================================================
%  LOCAL FUNCTION: bh_fdr
%  Benjamini-Hochberg FDR correction (duplicated from
%  preprocess_pipeline_data). NaN entries are ignored and returned NaN.
% ================================================================
function q = bh_fdr(pvals)
    q = nan(size(pvals));
    valid = ~isnan(pvals);
    m = sum(valid);
    if m == 0
        return;
    end

    idx         = find(valid);
    [ps, order] = sort(pvals(idx));
    ranks       = (1:m)';
    q_sorted    = ps(:) .* m ./ ranks;

    for i = m-1:-1:1
        q_sorted(i) = min(q_sorted(i), q_sorted(i+1));
    end
    q_sorted = min(q_sorted, 1);

    q(idx(order)) = q_sorted;
end


%% ================================================================
%  LOCAL FUNCTION: sig_stars / fmt_p
% ================================================================
function s = sig_stars(pval)
    if isnan(pval),      s = 'NA ';
    elseif pval < 0.001, s = '***';
    elseif pval < 0.01,  s = '** ';
    elseif pval < 0.05,  s = '*  ';
    else,                s = 'ns ';
    end
end

function s = fmt_p(p)
    if isnan(p)
        s = 'NA';
    elseif p == 0
        s = '<1e-300';
    else
        s = sprintf('%.4g', p);
    end
end


%% ================================================================
%  LOCAL FUNCTION: ask_step_choice
% ================================================================
function step = ask_step_choice()
    fprintf('\n============================================================\n');
    fprintf('  extract_cell_count_per_time\n');
    fprintf('============================================================\n');
    fprintf('  STEP 1 : Colony appearance counts per timepoint\n');
    fprintf('           How many colonies newly reach the size threshold at\n');
    fprintf('           each imaging timepoint, per sample group.\n\n');
    fprintf('  STEP 2 : EARLY vs LATE appearance timepoints\n');
    fprintf('           Pick how many EARLY and LATE timepoints in one sample\n');
    fprintf('           group; compare growth rate, final size, doubling time,\n');
    fprintf('           and size growth over +5h for colonies that first\n');
    fprintf('           appeared at each of those timepoints.\n');
    fprintf('============================================================\n');
    while true
        raw = strtrim(input('  Run which step? [1/2]: ', 's'));
        if strcmpi(raw, '1'), step = 1; return; end
        if strcmpi(raw, '2'), step = 2; return; end
        fprintf('  WARNING  Enter 1 or 2.\n');
    end
end


%% ================================================================
%  LOCAL FUNCTION: ask_threshold
%  Shared by Step 1 and Step 2.
% ================================================================
function threshold = ask_threshold(default_threshold)
    fprintf('\n  COLONY SIZE THRESHOLD (appear control)\n');
    fprintf('  A colony "appears" the first time its smoothed size timecourse\n');
    fprintf('  reaches this many pixels. Default matches the size_threshold\n');
    fprintf('  used in preprocess_pipeline_data (%g px).\n', default_threshold);
    while true
        raw_t = strtrim(input(sprintf( ...
            '  Threshold in px [default = %g, press Enter to keep]: ', default_threshold), 's'));
        if isempty(raw_t)
            threshold = default_threshold;
            fprintf('  -> Using default threshold: %g px\n', threshold);
            return;
        end
        val = str2double(raw_t);
        if isnan(val) || val <= 0
            fprintf('  WARNING  Enter a positive number, or press Enter for default.\n');
            continue;
        end
        threshold = val;
        fprintf('  -> Using threshold: %g px\n', threshold);
        return;
    end
end


%% ================================================================
%  LOCAL FUNCTION: ask_selections_step1
%  Groups (multi-select, or ALL) + appearance threshold.
% ================================================================
function [sel_groups, threshold] = ask_selections_step1(ht, default_threshold)

    all_labels = ht.labels;
    n_avail    = length(all_labels);

    fprintf('\n============================================================\n');
    fprintf('  STEP 1 — Selection\n');
    fprintf('============================================================\n');

    fprintf('\n  GROUPS\n');
    fprintf('  +-----+--------+-------------+\n');
    fprintf('  | idx | label  | n colonies  |\n');
    fprintf('  +-----+--------+-------------+\n');
    for k = 1:n_avail
        fprintf('  | %3d | %-6s | %11d |\n', k, all_labels{k}, ht.groups(k).n_colonies);
    end
    fprintf('  +-----+--------+-------------+\n');
    fprintf('\n============================================================\n');

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

    threshold = ask_threshold(default_threshold);
end


%% ================================================================
%  LOCAL FUNCTION: ask_single_group
%  Step 2 needs exactly ONE sample group (its own time grid + growth fits).
% ================================================================
function g = ask_single_group(ht)
    all_labels = ht.labels;
    n_avail    = length(all_labels);

    fprintf('\n============================================================\n');
    fprintf('  STEP 2 — Sample group selection (choose exactly ONE)\n');
    fprintf('============================================================\n');
    fprintf('  +-----+--------+-------+----------+\n');
    fprintf('  | idx | label  | n_col | n_fitted |\n');
    fprintf('  +-----+--------+-------+----------+\n');
    for k = 1:n_avail
        if isfield(ht.groups(k), 'growth') && isfield(ht.groups(k).growth, 'fit_ok')
            n_fit = sum(ht.groups(k).growth.fit_ok);
        else
            n_fit = 0;
        end
        fprintf('  | %3d | %-6s | %5d | %8d |\n', k, all_labels{k}, ht.groups(k).n_colonies, n_fit);
    end
    fprintf('  +-----+--------+-------+----------+\n');

    while true
        raw = strtrim(input('  Enter ONE group label: ', 's'));
        tokens = regexp(raw, '[,\s]+', 'split');
        tokens = tokens(~cellfun(@isempty, tokens));
        if length(tokens) ~= 1
            fprintf('  WARNING  Enter exactly one group label.\n');
            continue;
        end
        match = find(strcmpi(all_labels, tokens{1}));
        if isempty(match)
            fprintf('  WARNING  "%s" not recognised. Try again.\n', tokens{1});
            continue;
        end
        g = match(1);
        fprintf('  -> Selected group: %s\n', all_labels{g});
        return;
    end
end


%% ================================================================
%  LOCAL FUNCTION: ask_timepoint_selection
%  User-chosen NUMBER of EARLY and LATE timepoints (not fixed to 2+2),
%  then that many distinct timepoint indices into the group's time grid
%  for each phase. More timepoints per phase means more pairwise
%  comparisons (and a busier plot) — the user is warned, not blocked.
% ================================================================
function [early_idx, late_idx] = ask_timepoint_selection(T)
    fprintf('\n  TIMEPOINT SELECTION\n');
    fprintf('  Choose how many EARLY and how many LATE timepoints to compare\n');
    fprintf('  (need not be equal, e.g. 1 EARLY vs 3 LATE is fine).\n');

    n_early = ask_phase_count('EARLY', T);
    n_late  = ask_phase_count('LATE',  T, n_early);

    n_total_pairs = nchoosek(n_early + n_late, 2);
    if n_total_pairs > 10
        fprintf(['  NOTE: %d + %d timepoints -> %d pairwise comparisons will be drawn\n' ...
                  '        on each plot as stacked brackets — it will look busy.\n'], ...
                  n_early, n_late, n_total_pairs);
    end

    while true
        fprintf('  EARLY timepoint indices (%d needed, e.g. 2 4): ', n_early);
        raw_e = strtrim(input('', 's'));
        e_idx = parse_n_indices(raw_e, n_early, T);
        if isempty(e_idx), continue; end

        fprintf('  LATE timepoint indices  (%d needed, e.g. 10 14): ', n_late);
        raw_l = strtrim(input('', 's'));
        l_idx = parse_n_indices(raw_l, n_late, T);
        if isempty(l_idx), continue; end

        if any(ismember(e_idx, l_idx))
            fprintf('  WARNING  EARLY and LATE indices must be distinct (no overlap). Try again.\n');
            continue;
        end

        early_idx = sort(e_idx);
        late_idx  = sort(l_idx);

        if max(early_idx) >= min(late_idx)
            fprintf(['  WARNING  EARLY timepoints are not strictly earlier than LATE timepoints.\n' ...
                      '           Proceeding anyway — verify this is intended.\n']);
        end
        fprintf('  -> EARLY idx: [%s]   LATE idx: [%s]\n', ...
            num2str(early_idx), num2str(late_idx));
        return;
    end
end

function n = ask_phase_count(phase_name, T, other_n)
    if nargin < 3, other_n = 0; end
    while true
        raw = strtrim(input(sprintf('  How many %s timepoints? [default = 2]: ', phase_name), 's'));
        if isempty(raw)
            n = 2;
        else
            n = str2double(raw);
        end
        if isnan(n) || n < 1 || n ~= floor(n)
            fprintf('  WARNING  Enter a positive integer.\n');
            continue;
        end
        if n + other_n > T
            fprintf('  WARNING  %d %s + %d already chosen exceeds the %d available timepoints.\n', ...
                n, phase_name, other_n, T);
            continue;
        end
        return;
    end
end

function idx = parse_n_indices(raw, n, T)
    idx = [];
    tokens = regexp(raw, '[,\s]+', 'split');
    tokens = tokens(~cellfun(@isempty, tokens));
    if length(tokens) ~= n
        fprintf('  WARNING  Enter exactly %d integer(s).\n', n);
        return;
    end
    nums = str2double(tokens);
    if any(isnan(nums)) || any(nums < 1) || any(nums > T) || any(nums ~= floor(nums))
        fprintf('  WARNING  Enter integers between 1 and %d.\n', T);
        return;
    end
    if length(unique(nums)) ~= n
        fprintf('  WARNING  All %d indices must be different.\n', n);
        return;
    end
    idx = nums;
end

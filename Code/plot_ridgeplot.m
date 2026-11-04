%% HeteroTyper Pipeline
%  plot_ridgeplot  —  ridgeline + interval plots of colony features for any
%  set of sample groups (replaces plot_appearanceTime_ridgeplot,
%  plot_appearanceTime_ridgeplot_6s_stress and plot_appearanceTime_ridgeplot_6srna_growth).
%
%  One figure per selected feature:
%    Grey density per sample group with, on its baseline, bars in the group
%    colour covering the central 95 % (2.5–97.5th percentile, lightest),
%    80 % (10–90th) and 50 % of colonies (25–75th = IQR, darkest); black
%    dot = group median. Vertical dashed line = median of all selected
%    groups pooled. n, median, IQR, MAD_med and Gini_counts are written
%    next to every ridge; Brown–Forsythe results below the axes.
%
%  INTERACTIVE CHOICES (in this order)
%    1. Sample groups           labels / index numbers / ALL
%    2. Display order           top to bottom, as numbers of the selected
%                               groups (Enter = keep selection order)
%    3. Features                appearance time, size, area, intensity, ... — every
%                               feature field present in ht.groups (Enter =
%                               appearance time only)
%    4. Statistical comparisons all possible pairs of the selected groups
%                               are listed and numbered; enter the pair
%                               numbers to test (e.g. 1 4 6, ranges 1-3,
%                               ALL; Enter = none). The same pairs are
%                               tested for every selected feature.
%    5. Appearance threshold    for Gini_counts (px)
%
%  DATA
%    Feature value = ht.groups(g).<feature>, finite values only (appearance time:
%    ht.groups(g).appearance_time, non-NaN, as in plot_combined_samples_newGini).
%    x-axis: appearance time from incTime to max_appearance (5 h ticks); other features
%    from the 0.5th to the 99.5th percentile of all selected groups pooled
%    (5 % padding), so extreme outliers can lie off-axis — they are still
%    included in every statistic.
%
%  GINI (Gini_counts)
%    Appearance time : Gini of the appearance-count table (colonies newly reaching
%               the threshold at each imaging timepoint, not cumulative).
%    Others   : Gini of the counts of the appearance-counted colonies over
%               GINI_NBINS equal-width bins spanning the pooled range of the
%               selected groups (shared edges, so groups are comparable).
%               Edit gini_edges() to use the pipeline's histogram edges
%               instead, if those should be matched exactly.
%
%  BROWN–FORSYTHE TEST (spread, robust to skew)
%    z = |x - group median| per colony; one-way ANOVA F-test on z.
%    Reported for all selected groups together and for the chosen pairs,
%    Benjamini–Hochberg adjusted across the chosen pairs (per feature), with
%    the MAD_med ratio and IQR ratio as effect sizes. Tests ABSOLUTE spread
%    in the units of the feature: independent of a shift in the median.
%    Colonies of all positions are pooled and treated as independent.
%
%  Densities are Gaussian kernel density estimates (Silverman bandwidth;
%  for appearance time never narrower than the imaging interval, because appearance
%  times sit on the imaging grid). Percentiles use the same definition as
%  MATLAB prctile. No toolbox function is needed.
%
%  USAGE
%    plot_ridgeplot(ht)
%
%  Outputs (ht.params.out_dir):
%    Ridgeplot_<feature>_<ts>.fig/.pdf/.png   (one per feature)
%    plot_ridgeplot_<ts>.txt                  (log)
%    ridgeplot_stats_<ts>.xlsx                (Stats, BrownForsythe, README)

function plot_ridgeplot(ht)

    GINI_NBINS = 30;      % bins for Gini_counts of features other than appearance time

    if ~isstruct(ht) || ~isfield(ht, 'groups')
        error(['Input must be the ht struct produced by preprocess_pipeline_data.\n' ...
               'Run:  preprocess_pipeline_data(data)  first.']);
    end

    p       = ht.params;
    labels  = ht.labels;
    colors  = ht.colors;
    out_dir = p.out_dir;
    incTime = p.incTime;
    max_appearance = p.max_appearance;
    if isfield(p, 'img_int'),        img_int = p.img_int;          else, img_int = 0.5; end
    if isfield(p, 'size_threshold'), size_filter = p.size_threshold; else, size_filter = 100; end

    % ------------------------------------------------------------------
    %  Log
    % ------------------------------------------------------------------
    ts       = datestr(now, 'yyyymmdd_HHMMSS');
    log_path = fullfile(out_dir, sprintf('plot_ridgeplot_%s.txt', ts));
    xl_path  = fullfile(out_dir, sprintf('ridgeplot_stats_%s.xlsx', ts));
    flog = fopen(log_path, 'w');
    if flog == -1, warning('Could not open log file: %s', log_path); end

    function tlog(varargin)
        fprintf(varargin{:});
        if flog ~= -1, fprintf(flog, varargin{:}); end
    end

    % ------------------------------------------------------------------
    %  Interactive choices
    % ------------------------------------------------------------------
    sel       = ask_groups(ht, 'plot_ridgeplot');
    sel       = ask_order(sel, labels);
    feats     = ask_features(ht);
    pairs     = ask_pairs(labels(sel), 'Brown–Forsythe test of spread');   % indices into sel
    threshold = ask_threshold(size_filter);
    n_g       = numel(sel);

    tlog('\nplot_ridgeplot\n');
    tlog('  Groups (top to bottom): %s\n', strjoin(labels(sel), ', '));
    tlog('  Features   : %s\n', strjoin({feats.name}, ', '));
    if isempty(pairs)
        tlog('  Comparisons: none (all-group test only)\n');
    else
        tlog('  Comparisons: %s\n', strjoin(pair_labels(labels(sel), pairs), ', '));
    end
    tlog('  Gini: Gini_counts, appearance threshold %g px; %d bins for features other than appearance time\n', ...
         threshold, GINI_NBINS);
    tlog('  Output: %s\n\n', out_dir);

    fs_tick = 11;  fs_lbl = 12;  fs_grp = 12;  fs_note = 9;

    stat_rows = cell(0, 15);
    bf_rows   = cell(0, 11);

    % ==================================================================
    %  One figure per feature
    % ==================================================================
    for fi = 1:numel(feats)
        F  = feats(fi);
        us = unit_suffix(F.unit);

        % ---------------- values ----------------
        vals = cell(1, n_g);
        for k = 1:n_g
            v = double(ht.groups(sel(k)).(F.field));
            v = v(:);
            vals{k} = v(isfinite(v));
        end
        pooled = vertcat(vals{:});
        if isempty(pooled)
            tlog('  %s: no finite values in the selected groups — skipped.\n\n', F.name);
            continue;
        end

        % ---------------- x axis ----------------
        if strcmp(F.field, 'appearance_time')
            x_step  = 5;                                 % ticks: incTime, +5 h, ...
            x_end   = incTime + ceil((max_appearance - incTime) / x_step) * x_step;
            x_lim   = [incTime, x_end];
            x_ticks = incTime : x_step : x_end;
            min_bw  = img_int;
        else
            [x_lim, x_ticks] = feature_axis(pooled);
            min_bw = diff(x_lim) / 1000;
        end
        x_grid = linspace(x_lim(1), x_lim(2), 512);

        ng_edges = cell(1, 10);
        if F.gi > 1
            ng_edges{F.gi} = gini_edges(pooled, GINI_NBINS);
        end

        % ---------------- per-group statistics ----------------
        S = struct('lbl', {}, 'col', {}, 'v', {}, 'n', {}, 'med', {}, 'q', {}, ...
                   'dens', {}, 'gini', {}, 'n_gini', {}, 'mad', {});
        for k = 1:n_g
            g  = sel(k);
            v  = vals{k};
            ng = newgini_colony_set(ht.groups(g), [], threshold, ng_edges);

            S(k).lbl    = labels{g};
            S(k).col    = colors(g, :);
            S(k).v      = v;
            S(k).n      = numel(v);
            S(k).med    = median_or_nan(v);
            S(k).q      = qtile(v, [2.5 10 25 75 90 97.5]);
            S(k).dens   = kde_gauss(v, x_grid, min_bw);
            S(k).gini   = ng.counts(F.gi);
            S(k).n_gini = ng.n(F.gi);
        end
        med_all = median_or_nan(pooled);

        tlog('====================================================================\n');
        tlog('  FEATURE: %s%s\n', F.name, unit_paren(F.unit));
        tlog('====================================================================\n');
        tlog('  %-14s  %6s  %10s  %10s  %10s  %10s  %10s  %10s  %10s  %10s  %11s\n', ...
            'Group', 'n', 'Median', 'Q25', 'Q75', 'IQR', 'P10', 'P90', 'P2.5', 'P97.5', 'Gini_counts');
        tlog('  %s\n', repmat('-', 1, 128));
        for k = 1:n_g
            q = S(k).q;
            tlog('  %-14s  %6d  %10.4g  %10.4g  %10.4g  %10.4g  %10.4g  %10.4g  %10.4g  %10.4g  %11.4f\n', ...
                S(k).lbl, S(k).n, S(k).med, q(3), q(4), q(4)-q(3), q(2), q(5), q(1), q(6), S(k).gini);
        end
        tlog('  Pooled median of all selected groups: %.4g%s\n\n', med_all, us);

        % ---------------- Brown–Forsythe ----------------
        BF = brown_forsythe_all({S.v}, {S.lbl}, pairs);
        for k = 1:n_g, S(k).mad = BF.mad(k); end

        tlog('  BROWN–FORSYTHE TEST — spread of %s\n', lower(F.name));
        tlog('  (z = |x - group median|; ANOVA on z. MAD_med = mean z per group%s.)\n', unit_comma(F.unit));
        tlog('  %-14s  %6s  %10s\n', 'Group', 'n', 'MAD_med');
        for k = 1:n_g
            tlog('  %-14s  %6d  %10.4g\n', S(k).lbl, S(k).n, S(k).mad);
        end
        tlog('\n  All groups : F(%d, %d) = %.3f, p = %s %s\n\n', ...
            BF.all.df1, BF.all.df2, BF.all.F, fmt_p(BF.all.p), stars(BF.all.p));
        if ~isempty(BF.pair)
            tlog('  %-32s  %9s  %8s  %10s  %10s %4s  %11s  %9s\n', ...
                'Pair', 'F', 'df2', 'p', 'p_BH', 'sig', 'MAD ratio', 'IQR ratio');
            tlog('  %s\n', repmat('-', 1, 108));
            for i = 1:numel(BF.pair)
                b = BF.pair(i);
                tlog('  %-32s  %9.3f  %8d  %10s  %10s %4s  %11.3f  %9.3f\n', ...
                    b.label, b.F, b.df2, fmt_p(b.p), fmt_p(b.p_bh), stars(b.p_bh), ...
                    b.mad_ratio, iqr_ratio(S(b.i), S(b.j)));
            end
            tlog(['  MAD / IQR ratio = first group / second group (>1: first group more spread).\n' ...
                  '  p_BH: Benjamini–Hochberg across the chosen pairs (within this feature). Colonies\n' ...
                  '  pooled over positions and treated as independent; large n makes small differences significant.\n\n']);
        end
        bf_lines = bf_text_lines(BF, S, F.name);

        % ---------------- figure: grey ridge + 95/80/50 % bars + median dot ----------------
        n_bf   = numel(bf_lines);
        bf_hgt = 0.35 + 0.19*n_bf;              % room below the x-axis for the BF block (in)
        ridge_h2 = 0.75;   bar_h = 0.16;
        left_w = 1.45;  info_w = 1.7;  right_m = 0.3;
        fig_w2 = 10.5;  fig_h2 = 1.1 + 1.0*n_g + 1.6 + bf_hgt;   % top: legend; bottom: BF block
        ax_b2  = 0.75 + bf_hgt;
        hf = figure('Name', sprintf('Ridgeline + intervals — %s', F.name), 'Color','w', ...
                    'Units','inches', 'Position',[1.5 1 fig_w2 fig_h2], 'PaperUnits','inches', ...
                    'PaperSize',[fig_w2 fig_h2], 'PaperPosition',[0 0 fig_w2 fig_h2], ...
                    'Visible','off');
        ax_w = fig_w2 - left_w - info_w - right_m;
        ax_h = fig_h2 - ax_b2 - 1.75;
        ax = axes('Position', [left_w/fig_w2, ax_b2/fig_h2, ax_w/fig_w2, ax_h/fig_h2]);
        hold(ax, 'on');

        base  = n_g - (1:n_g);
        y_top = max(base) + ridge_h2 + 0.35;
        if isfinite(med_all)
            plot(ax, [med_all med_all], [-0.45, y_top - 0.1], '--', 'Color', [0.25 0.25 0.25], 'LineWidth', 1.4);
            text(ax, med_all + 0.01*diff(x_lim), y_top - 0.1, ...
                 sprintf('Median, all groups (%s%s)', fmt_v(med_all, F), us), ...
                 'FontSize', fs_note+1, 'VerticalAlignment','top', 'HorizontalAlignment','left');
        end

        for k = 1:n_g
            d = S(k).dens;
            b = base(k);
            if ~isempty(d) && max(d) > 0
                y = b + ridge_h2 * d / max(d);
                keep = d > 1e-3 * max(d);
                xs = x_grid(keep);  ys = y(keep);
                patch(ax, [xs, fliplr(xs)], [ys, b*ones(size(ys))], [0.83 0.83 0.83], ...
                      'EdgeColor', 'none', 'FaceAlpha', 0.9);
            end
            q = S(k).q;
            if all(isfinite(q))
                shades = [shade(S(k).col, 0.70); shade(S(k).col, 0.40); S(k).col];   % 95, 80, 50 %
                lohi   = [q(1) q(6); q(2) q(5); q(3) q(4)];
                for r = 1:3
                    patch(ax, [lohi(r,1) lohi(r,2) lohi(r,2) lohi(r,1)], ...
                              [b - bar_h/2, b - bar_h/2, b + bar_h/2, b + bar_h/2], ...
                              shades(r,:), 'EdgeColor', 'none');
                end
                plot(ax, S(k).med, b, 'o', 'MarkerSize', 7, ...
                     'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
            end
            % info column right of the axes: n, median [IQR], MAD_med, Gini_counts
            text(ax, x_lim(2) + 0.02*diff(x_lim), b + 0.05, ...
                 sprintf('n = %d\nmedian %s%s  [IQR %s]\nMAD_{med} = %s%s\nGini = %.3f', ...
                         S(k).n, fmt_v(S(k).med, F), us, fmt_v(q(4) - q(3), F), ...
                         fmt_v(S(k).mad, F), us, S(k).gini), ...
                 'FontSize', fs_note, 'VerticalAlignment','middle', 'Clipping','off');
        end

        set(ax, 'YTick', fliplr(base), 'YTickLabel', fliplr({S.lbl}), ...
                'TickLabelInterpreter','none', 'FontSize', fs_tick, 'Box','off', ...
                'TickDir','out', 'XLim', x_lim, 'XTick', x_ticks, ...
                'YLim', [-0.45, y_top], 'XGrid','on', 'GridColor',[0.85 0.85 0.85], ...
                'GridAlpha', 1, 'YColor', 'none');
        for k = 1:n_g     % group names in group colour (YColor none hides tick labels)
            text(ax, x_lim(1) - 0.02*diff(x_lim), base(k), S(k).lbl, ...
                 'HorizontalAlignment','right', 'VerticalAlignment','middle', ...
                 'FontSize', fs_grp, 'FontWeight','bold', 'Color', S(k).col, ...
                 'Interpreter','none', 'Clipping','off');
        end
        xlabel(ax, [F.name unit_paren(F.unit)], 'FontSize', fs_lbl, 'FontWeight','bold');

        % Legend (top, neutral grey so it applies to every group colour)
        lg_w = 2.6;  lg_h = 1.15;      % legend graphic; its text runs to the right (~4 in)
        lax = axes('Position', [left_w/fig_w2, (fig_h2 - 0.2 - lg_h)/fig_h2, ...
                                lg_w/fig_w2, lg_h/fig_h2], 'XLim',[0 10], 'YLim',[0 1.6], ...
                   'Visible','off');
        hold(lax, 'on');
        xd = linspace(0.6, 9.4, 200);
        dd = exp(-0.5*((xd-4.3)/1.1).^2) + 0.35*exp(-0.5*((xd-6.6)/1.0).^2);
        patch(lax, [xd fliplr(xd)], [0.55 + 0.9*dd/max(dd), 0.55*ones(size(xd))], [0.83 0.83 0.83], 'EdgeColor','none');
        gcol = [0.35 0.35 0.35];
        lq = [1.2 8.9; 2.2 7.6; 3.3 5.9];   lmed = 4.6;
        shg = [shade(gcol,0.70); shade(gcol,0.40); gcol];
        for r = 1:3
            patch(lax, [lq(r,1) lq(r,2) lq(r,2) lq(r,1)], [0.47 0.47 0.63 0.63], shg(r,:), 'EdgeColor','none');
        end
        plot(lax, lmed, 0.55, 'o', 'MarkerSize', 6, 'MarkerFaceColor','k', 'MarkerEdgeColor','k');
        text(lax, 10.1, 1.25, sprintf('Distribution of %s', lower(F.name)), 'FontSize', fs_note, 'Clipping','off');
        text(lax, 10.1, 0.95, 'Median (dot)', 'FontSize', fs_note, 'Clipping','off');
        text(lax, 10.1, 0.65, 'Bars: central 50 % (IQR, darkest), 80 % and 95 % (lightest) of colonies', ...
             'FontSize', fs_note, 'Clipping','off');
        text(lax, 10.1, 0.35, 'Bars are drawn in shades of each sample-group colour', 'FontSize', fs_note, ...
             'Clipping','off');

        draw_text_block(hf, [0.25, 0.1, fig_w2 - 0.5, bf_hgt - 0.05], fig_w2, fig_h2, bf_lines, fs_note);
        stem = sprintf('Ridgeplot_%s_%s', F.field, ts);
        save_fig(hf, fullfile(out_dir, stem));
        tlog('Figure saved: %s (.fig/.pdf/.png)\n\n', stem);

        % ---------------- rows for the xlsx ----------------
        for k = 1:n_g
            q = S(k).q;
            stat_rows(end+1, :) = {F.name, F.unit, S(k).lbl, S(k).n, S(k).med, q(3), q(4), q(4)-q(3), ...
                                   q(2), q(5), q(1), q(6), S(k).mad, S(k).n_gini, S(k).gini}; %#ok<AGROW>
        end
        bf_rows(end+1, :) = {F.name, 'All groups', strjoin({S.lbl}, ' | '), BF.all.F, BF.all.df1, ...
                             BF.all.df2, fmt_p(BF.all.p), fmt_p(BF.all.p), strtrim(stars(BF.all.p)), NaN, NaN}; %#ok<AGROW>
        for i = 1:numel(BF.pair)
            b = BF.pair(i);
            bf_rows(end+1, :) = {F.name, 'Pair', b.label, b.F, b.df1, b.df2, fmt_p(b.p), fmt_p(b.p_bh), ...
                                 strtrim(stars(b.p_bh)), b.mad_ratio, iqr_ratio(S(b.i), S(b.j))}; %#ok<AGROW>
        end
    end

    % ------------------------------------------------------------------
    %  Excel
    % ------------------------------------------------------------------
    if ~isempty(stat_rows)
        try
            writetable(cell2table(stat_rows, 'VariableNames', {'Feature','Unit','Group','n','Median', ...
                       'Q25','Q75','IQR','P10','P90','P2_5','P97_5','MAD_med','n_Gini','Gini_counts'}), ...
                       xl_path, 'Sheet', 'Stats');
            writetable(cell2table(bf_rows, 'VariableNames', {'Feature','Comparison','Groups','F','df1','df2', ...
                       'p','p_BH','Sig','MAD_ratio','IQR_ratio'}), xl_path, 'Sheet', 'BrownForsythe');
            readme = { ...
                'Column', 'Description'; ...
                'Feature / Unit', 'Colony feature (ht.groups(g).<field>) and its unit; all values of a row are in this unit'; ...
                'n', 'Colonies with a finite value of the feature'; ...
                'Median / Q25 / Q75 / IQR', 'Median, quartiles and interquartile range; the darkest bar in the figure spans Q25-Q75'; ...
                'P10 / P90', '10th / 90th percentile: the 80 % bar in the figure'; ...
                'P2_5 / P97_5', '2.5th / 97.5th percentile: the 95 % bar in the figure'; ...
                'MAD_med', 'Mean absolute deviation from the group median: the per-group quantity compared by the Brown–Forsythe test'; ...
                'BrownForsythe sheet', ['Brown–Forsythe test of equal spread (ANOVA on |x - group median|). All groups together and the pairs chosen by the user; ' ...
                                        'p_BH = Benjamini–Hochberg adjusted across the chosen pairs within each feature (= p for the all-group test). ' ...
                                        'MAD_ratio / IQR_ratio = first / second group of the pair (>1: first more spread). ' ...
                                        'Colonies pooled over positions and treated as independent.']; ...
                'n_Gini / Gini_counts', sprintf(['Appearance time: Gini of the Count column of the appearance table (colonies newly reaching %g px at each imaging timepoint, ' ...
                                         'not cumulative; zero rows included). Other features: Gini of the counts of the same colonies over %d equal-width bins ' ...
                                         'spanning the pooled range of the selected groups (zero bins included).'], threshold, GINI_NBINS); ...
            };
            writetable(cell2table(readme(2:end,:), 'VariableNames', readme(1,:)), xl_path, 'Sheet', 'README');
            tlog('Statistics saved: %s\n', xl_path);
        catch ME
            tlog('WARNING: could not write xlsx (%s)\n', ME.message);
        end
    end

    tlog('Done.\n');
    if flog ~= -1
        fclose(flog);
        fprintf('Log saved: %s\n', log_path);
    end
end


%% ====================================================================
%  LOCAL: feature_table — every feature the script can plot.
%  gi = column of the feature in newgini_colony_set (keep this order).
%  Edit name / unit / fmt here if the units in your dataset differ.
%% ====================================================================
function T = feature_table()
    T = struct( ...
        'field', {'appearance_time','size','area','intensity','mean_intensity', ...
                  'int_per_size','perimeter','circularity','eccentricity','solidity'}, ...
        'name',  {'Appearance time','Size','Area','Intensity','Mean intensity', ...
                  'Intensity per size','Perimeter','Circularity','Eccentricity','Solidity'}, ...
        'unit',  {'h','px','px','a.u.','a.u.','a.u./px','px','','',''}, ...
        'fmt',   {'%.1f','%.3g','%.3g','%.3g','%.3g','%.3g','%.3g','%.3g','%.3g','%.3g'}, ...
        'gi',    num2cell(1:10));
end


%% ====================================================================
%  LOCAL: ask_features — numbered list of the features present in
%  ht.groups; numbers / names / ALL; Enter = appearance time.
%% ====================================================================
function feats = ask_features(ht)
    T  = feature_table();
    ok = false(1, numel(T));
    for i = 1:numel(T)
        ok(i) = isfield(ht.groups, T(i).field);
    end
    T = T(ok);
    if isempty(T), error('plot_ridgeplot: none of the feature fields is present in ht.groups.'); end

    fprintf('\n  Features available:\n');
    for i = 1:numel(T)
        fprintf('  %3d | %-20s %s\n', i, T(i).name, unit_paren(T(i).unit));
    end
    def = find(strcmp({T.field}, 'appearance_time'), 1);
    if isempty(def), def = 1; end
    while true
        fprintf('  Features — numbers or names (comma-separated), or ALL [Enter = %s]: ', T(def).name);
        raw = strtrim(input('', 's'));
        if isempty(raw),         idx = def;           break; end
        if strcmpi(raw, 'all'),  idx = 1:numel(T);    break; end
        idx = parse_index_list(raw, numel(T));
        if isempty(idx)          % try names
            nrm = @(s) lower(regexprep(char(s), '[\s_\-]', ''));
            tok = strtrim(regexp(raw, ',', 'split'));
            idx = zeros(1, numel(tok));
            for t = 1:numel(tok)
                m = find(strcmp(cellfun(nrm, {T.name}, 'UniformOutput', false), nrm(tok{t})) | ...
                         strcmp(cellfun(nrm, {T.field}, 'UniformOutput', false), nrm(tok{t})), 1);
                if isempty(m), idx = []; break; end
                idx(t) = m;
            end
            if ~isempty(idx), [~, ui] = unique(idx, 'stable'); idx = idx(sort(ui)); end
        end
        if ~isempty(idx), break; end
        fprintf('  WARNING: feature(s) not recognised. Try again.\n');
    end
    feats = T(idx);
    fprintf('  -> Features: %s\n', strjoin({feats.name}, ', '));
end


%% ====================================================================
%  LOCAL: ask_order — display order of the selected groups (top to
%  bottom). Numbers refer to the list shown; groups not listed follow in
%  their current order. Enter keeps the current order.
%% ====================================================================
function sel = ask_order(sel, labels)
    n = numel(sel);
    if n < 2, return; end
    fprintf('\n  Display order (top to bottom) — current order:\n');
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
    fprintf('  -> Order: %s\n', strjoin(labels(sel), ', '));
end


%% ====================================================================
%  LOCAL: ask_pairs — lists every pair of the selected groups (in display
%  order) and returns the chosen ones as [P x 2] indices into lbls.
%% ====================================================================
function P = ask_pairs(lbls, test_name)
    n = numel(lbls);
    P = zeros(0, 2);
    if n < 2
        fprintf('  (fewer than two groups: no pairwise comparisons)\n');
        return;
    end
    all_p = nchoosek(1:n, 2);
    fprintf('\n  Pairwise comparisons available (%s):\n', test_name);
    pl = pair_labels(lbls, all_p);
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
        fprintf('  -> Comparisons: none\n');
    else
        fprintf('  -> Comparisons: %s\n', strjoin(pair_labels(lbls, P), ', '));
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


%% ====================================================================
%  LOCAL: feature_axis — x limits and ticks for a feature other than appearance
%  time: 0.5–99.5th percentile of the pooled values, 5 % padding, rounded
%  to a 1/2/2.5/5 step (never below 0 for non-negative features).
%% ====================================================================
function [lim, ticks] = feature_axis(v)
    lo = qtile(v, 0.5);  hi = qtile(v, 99.5);
    if ~(hi > lo), lo = min(v);  hi = max(v); end
    if ~(hi > lo), lo = lo - 0.5;  hi = hi + 0.5; end
    pad = 0.05 * (hi - lo);
    lo  = lo - pad;  hi = hi + pad;
    if min(v) >= 0, lo = max(lo, 0); end
    step  = nice_step((hi - lo) / 6);
    lo    = floor(lo / step) * step;
    hi    = ceil(hi / step) * step;
    lim   = [lo, hi];
    ticks = lo : step : hi;
end


function s = nice_step(raw)
    mag = 10^floor(log10(raw));
    m   = raw / mag;
    if     m <= 1,   s = 1;
    elseif m <= 2,   s = 2;
    elseif m <= 2.5, s = 2.5;
    elseif m <= 5,   s = 5;
    else,            s = 10;
    end
    s = s * mag;
end


%% ====================================================================
%  LOCAL: gini_edges — shared bin edges for Gini_counts of a feature
%  other than appearance time (equal width over the pooled range).
%% ====================================================================
function e = gini_edges(v, nb)
    lo = min(v);  hi = max(v);
    if ~(hi > lo), lo = lo - 0.5;  hi = hi + 0.5; end
    e = linspace(lo, hi, nb + 1);
end


function s = unit_suffix(u)
    if isempty(u), s = ''; else, s = [' ' u]; end
end

function s = unit_paren(u)
    if isempty(u), s = ''; else, s = [' (' u ')']; end
end

function s = unit_comma(u)
    if isempty(u), s = ''; else, s = [', in ' u]; end
end

function s = fmt_v(x, F)
    if ~isfinite(x), s = 'NA'; else, s = sprintf(F.fmt, x); end
end


%% ====================================================================
%  LOCAL: save_fig
%% ====================================================================
function save_fig(h, path_stem)
    set(h, 'Visible', 'on');
    drawnow;
    savefig(h, [path_stem '.fig']);
    exportgraphics(h, [path_stem '.pdf'], 'ContentType','vector');
    exportgraphics(h, [path_stem '.png'], 'Resolution',300);
end


%% ====================================================================
%  LOCAL: brown_forsythe_all
%  Brown–Forsythe test (median-centred Levene): z = |x - median of its
%  group|, one-way ANOVA on z. All groups together + the chosen pairs,
%  with Benjamini–Hochberg adjustment across the pairs.
%% ====================================================================
function BF = brown_forsythe_all(vals, lbls, pairs)
% pairs : [P x 2] indices into vals of the pairs to test
    k = numel(vals);
    BF.mad = nan(1, k);
    for i = 1:k
        v = vals{i}(:);
        if ~isempty(v), BF.mad(i) = mean(abs(v - median(v))); end
    end
    [BF.all.F, BF.all.df1, BF.all.df2, BF.all.p] = bf_core(vals);

    BF.pair = struct('i', {}, 'j', {}, 'label', {}, 'F', {}, 'df1', {}, 'df2', {}, ...
                     'p', {}, 'p_bh', {}, 'mad_ratio', {});
    for r = 1:size(pairs, 1)
        i = pairs(r, 1);  j = pairs(r, 2);
        [F, d1, d2, pv] = bf_core(vals([i j]));
        BF.pair(end+1) = struct('i', i, 'j', j, ...
            'label', sprintf('%s vs %s', lbls{i}, lbls{j}), 'F', F, 'df1', d1, ...
            'df2', d2, 'p', pv, 'p_bh', NaN, 'mad_ratio', BF.mad(i) / BF.mad(j)); %#ok<AGROW>
    end
    if ~isempty(BF.pair)
        padj = bh_adjust([BF.pair.p]);
        for i = 1:numel(BF.pair), BF.pair(i).p_bh = padj(i); end
    end
end


function [F, df1, df2, p] = bf_core(vals)
    keep = cellfun(@(v) numel(v) >= 2, vals);
    vals = vals(keep);
    k = numel(vals);
    F = NaN;  p = NaN;  df1 = k - 1;  df2 = 0;
    if k < 2, return; end
    z = cell(1, k);
    for i = 1:k
        v = vals{i}(:);
        z{i} = abs(v - median(v));
    end
    n    = cellfun(@numel, z);
    N    = sum(n);
    df2  = N - k;
    zbar = cellfun(@mean, z);
    zall = mean(vertcat(z{:}));
    ssb  = sum(n .* (zbar - zall).^2);
    ssw  = sum(cellfun(@(zi) sum((zi - mean(zi)).^2), z));
    if ssw <= 0, return; end
    F = (ssb / df1) / (ssw / df2);
    p = f_upper(F, df1, df2);
end


function p = f_upper(F, d1, d2)
% Upper tail of the F distribution via the regularised incomplete beta
% function (base MATLAB, no toolbox): P(F' > F).
    if ~isfinite(F) || F <= 0, p = 1; return; end
    p = betainc(d2 / (d2 + d1 * F), d2 / 2, d1 / 2);
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


function r = iqr_ratio(a, b)
    ia = a.q(4) - a.q(3);  ib = b.q(4) - b.q(3);
    if ib > 0, r = ia / ib; else, r = NaN; end
end


function s = fmt_p(p)
    if isnan(p),      s = 'NA';
    elseif p == 0,    s = '<1e-300';
    elseif p < 1e-4,  s = sprintf('%.2e', p);
    else,             s = sprintf('%.4f', p);
    end
end


function s = stars(p)
    if isnan(p),      s = 'NA ';
    elseif p < 0.001, s = '***';
    elseif p < 0.01,  s = '** ';
    elseif p < 0.05,  s = '*  ';
    else,             s = 'ns ';
    end
end


%% ====================================================================
%  LOCAL: bf_text_lines — Brown–Forsythe summary lines for the figure
%% ====================================================================
function L = bf_text_lines(BF, S, fname)
    L = {sprintf(['Brown–Forsythe test of %s spread (ANOVA on |x - group median|):  ' ...
                  'all groups F(%d, %d) = %.2f, p = %s (%s)'], ...
                 lower(fname), BF.all.df1, BF.all.df2, BF.all.F, fmt_p(BF.all.p), strtrim(stars(BF.all.p)))};
    for i = 1:numel(BF.pair)
        b = BF.pair(i);
        L{end+1} = sprintf('   %s:  F(1, %d) = %.2f,  p_{BH} = %s (%s),  MAD ratio %.2f,  IQR ratio %.2f', ...
            strrep(b.label, '_', '\_'), b.df2, b.F, fmt_p(b.p_bh), strtrim(stars(b.p_bh)), b.mad_ratio, ...
            iqr_ratio(S(b.i), S(b.j))); %#ok<AGROW>
    end
    if numel(BF.pair) > 1
        L{end+1} = sprintf('   p_{BH}: Benjamini–Hochberg across the %d chosen comparisons', numel(BF.pair));
    end
end


function draw_text_block(hfig, pos_in, fig_w, fig_h, lines, fs)
% Text lines in an invisible axes at pos_in = [left bottom width height] (in).
    ax = axes('Parent', hfig, 'Position', [pos_in(1)/fig_w, pos_in(2)/fig_h, ...
              pos_in(3)/fig_w, pos_in(4)/fig_h], 'Visible','off', 'XLim',[0 1], 'YLim',[0 1]);
    n = numel(lines);
    for i = 1:n
        text(ax, 0, 1 - (i - 0.5)/max(n, 1), lines{i}, 'FontSize', fs, ...
             'HorizontalAlignment','left', 'VerticalAlignment','middle', 'Interpreter','tex');
    end
end


%% ====================================================================
%  LOCAL: kde_gauss — Gaussian kernel density on grid xg.
%  Bandwidth: Silverman's rule, but at least min_bw (appearance time: the
%  imaging interval, because appearance times are discrete on the imaging grid).
%% ====================================================================
function d = kde_gauss(v, xg, min_bw)
    v = v(isfinite(v));
    n = numel(v);
    if n < 2, d = zeros(size(xg)); return; end
    s  = std(v);
    iq = diff(qtile(v, [25 75]));
    sig = min([s, iq/1.34]);
    if ~(sig > 0), sig = max(s, iq/1.34); end
    bw = 0.9 * sig * n^(-1/5);
    if ~(bw > 0) || bw < min_bw, bw = min_bw; end
    d = zeros(size(xg));
    for i = 1:n
        d = d + exp(-0.5 * ((xg - v(i)) / bw).^2);
    end
    d = d / (n * bw * sqrt(2*pi));
end


%% ====================================================================
%  LOCAL: qtile — percentiles, same definition as MATLAB prctile
%% ====================================================================
function q = qtile(v, pcts)
    v = sort(v(isfinite(v)));
    n = numel(v);
    q = nan(size(pcts));
    if n == 0, return; end
    if n == 1, q(:) = v; return; end
    pos = n * pcts / 100 + 0.5;          % prctile: value i sits at 100*(i-0.5)/n
    pos = min(max(pos, 1), n);
    lo  = floor(pos);  hi = ceil(pos);
    q   = v(lo)' + (pos - lo) .* (v(hi)' - v(lo)');
end


function c = shade(c, f)
% Mix colour c with white by fraction f (0 = c, 1 = white).
    c = c + (1 - c) * f;
end


function m = median_or_nan(v)
    if isempty(v), m = NaN; else, m = median(v); end
end


%% ====================================================================
%  LOCAL: ask_groups
%% ====================================================================
function sel = ask_groups(ht, script_name)
    all_labels = ht.labels;
    fprintf('\n============================================================\n');
    fprintf('  %s — Group selection\n', script_name);
    fprintf('============================================================\n');
    for k = 1:numel(all_labels)
        fprintf('  %3d | %-14s | n colonies %d\n', k, all_labels{k}, ht.groups(k).n_colonies);
    end
    while true
        fprintf('  Groups — labels (comma-separated) / index numbers, or ALL: ');
        raw = strtrim(input('', 's'));
        if strcmpi(raw, 'all')
            sel = 1:numel(all_labels);
            break;
        end
        sel = parse_group_tokens(raw, all_labels);
        if ~isempty(sel), break; end
        fprintf('  WARNING: group(s) not recognised. Try again.\n');
    end
    fprintf('  -> Groups: %s\n', strjoin(all_labels(sel), ', '));
end


%% ================================================================
%  NEW GINI — LOCAL FUNCTIONS  (Gini_counts)
%
%  Appearance time  : appearance table. Rows = every imaging timepoint;
%              Count = colonies whose smoothed size timecourse FIRST
%              reaches >= threshold px at that timepoint (newly
%              appearing, not cumulative; each colony counted once).
%  Other     : the colonies counted in the appearance table, binned on
%  features    the given bin edges; Count = colonies per bin.
%  Gini_counts = Gini of the Count column itself (zero rows included):
%              how unevenly colonies are spread over the timepoints / bins.
% ================================================================
function ng = newgini_colony_set(grp, col_idx, threshold, edges)
% grp       : one ht.groups(g)
% col_idx   : colony indices (columns of size_timecourse) to use; [] = all
% threshold : appearance threshold (px)
% edges     : {1 x 10} bin edges for features 2-10 (edges{1} unused)
    fields = {'appearance_time','size','area','intensity','mean_intensity', ...
              'int_per_size','perimeter','circularity','eccentricity','solidity'};
    ng = struct('ok', false, 'counts', nan(1,10), ...
                'n', zeros(1,10), 'time_h', [], 'appearance_counts', [], ...
                'n_counted', 0, 'n_appearance', 0, 'n_appearance_nocnt', 0, 'n_cnt_noappearance', 0);

    if ~isfield(grp, 'size_timecourse') || isempty(grp.size_timecourse) || ...
       ~isfield(grp, 'time_h') || isempty(grp.time_h)
        return;
    end
    tc  = grp.size_timecourse;     % [T x N]
    t_h = grp.time_h(:);           % [T x 1]
    [T, N] = size(tc);
    if isempty(col_idx), col_idx = (1:N)'; end
    col_idx = col_idx(:);
    col_idx = col_idx(col_idx >= 1 & col_idx <= N);

    % ---------------- appearance (appearance time) table ----------------
    first_idx = compute_first_appearance(tc(:, col_idx), threshold);
    counted   = ~isnan(first_idx);
    if any(counted)
        counts = accumarray(first_idx(counted), 1, [T 1]);
    else
        counts = zeros(T, 1);
    end

    appearance = grp.appearance_time(:);
    if numel(appearance) == N
        appearance = appearance(col_idx);
        ng.n_appearance       = sum(~isnan(appearance));
        ng.n_appearance_nocnt = sum(~counted & ~isnan(appearance));
        ng.n_cnt_noappearance = sum(counted & isnan(appearance));
    end

    ng.ok         = true;
    ng.time_h     = t_h;
    ng.appearance_counts = counts;
    ng.n_counted  = sum(counts);
    ng.n(1)       = ng.n_counted;
    if ng.n(1) > 0
        ng.counts(1) = compute_gini_plain(counts);   % Gini of the Count column
    end

    % ---------------- other colony features ----------------
    for fi = 2:10
        if ~isfield(grp, fields{fi}) || numel(edges) < fi || numel(edges{fi}) < 2
            continue;
        end
        v = grp.(fields{fi});
        v = v(:);
        if numel(v) ~= N, continue; end
        v    = v(col_idx);
        x    = v(counted & ~isnan(v));
        e    = edges{fi}(:)';
        fcnt = histcounts(x, e)';
        ng.n(fi) = sum(fcnt);
        if ng.n(fi) > 0
            ng.counts(fi) = compute_gini_plain(fcnt);   % Gini of the Count column
        end
    end
end


%% ================================================================
%  LOCAL FUNCTION: compute_first_appearance
%  First row index per colony where the smoothed size timecourse
%  reaches >= threshold. NaN if never (size-rejected / censored
%  colonies are all-NaN in size_timecourse and are excluded).
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
%  LOCAL FUNCTION: compute_gini_plain
%  Classic Gini (0 = perfectly equal, 1 = maximally unequal), no origin.
%  Same as compute_gini in compute_gini_from_appearance_counts.m.
% ================================================================
function g = compute_gini_plain(x)
    x = x(isfinite(x) & x >= 0);
    if isempty(x) || sum(x) == 0
        g = 0;
        return;
    end
    x = sort(x(:));
    n = length(x);
    g = (2 * sum((1:n)' .* x) / (n * sum(x))) - (n+1)/n;
    g = max(0, min(1, g));
end


%% ================================================================
%  LOCAL FUNCTION: ask_threshold
%  Appearance threshold (px); Enter keeps the default.
% ================================================================
function thr = ask_threshold(default_thr)
    while true
        fprintf('  Appearance threshold in px for Gini_counts [Enter = %g]: ', default_thr);
        raw = strtrim(input('', 's'));
        if isempty(raw)
            thr = default_thr;
            break;
        end
        v = str2double(raw);
        if isfinite(v) && v > 0
            thr = v;
            break;
        end
        fprintf('  WARNING: enter a positive number or press Enter.\n');
    end
    fprintf('  -> Appearance threshold: %g px\n', thr);
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

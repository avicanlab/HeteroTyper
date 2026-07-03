%% HeteroTyper Pipeline
%  plot_morphology_colonies_Gini  —  Violin + scatter of lag time AND final
%  colony size per plate, grouped by sample time.  For each metric the Gini
%  coefficient is annotated on every plate tick, and the correlation between
%  colony count and Gini is shown in the subplot title.
%
%  Four figures are produced:
%    Figure 1 — Lag time,      Spearman corr(N, Gini) in title
%    Figure 2 — Colony size,   Spearman corr(N, Gini) in title
%    Figure 3 — Lag time,      Pearson  corr(N, Gini) in title
%    Figure 4 — Colony size,   Pearson  corr(N, Gini) in title
%
%%  USAGE
%  -----
%    stats = plot_morphology_colonies_Gini(ht, your_data_variable)
%
%  Example:
%    stats = plot_morphology_colonies_Gini(ht, data_26C_050626)
%
%  On first call the function asks for bio/tech replicate column labels.
%  To re-select: ht = rmfield(ht, 'plate_labels');  then re-run.
%
%  DEPENDENCY
%  ----------
%  Per-plate colony size is read from ht.per_plate_size (populated by
%  preprocess_pipeline_data).  If that field is absent the function falls
%  back to extracting size on the fly from raw_data.

function stats = plot_morphology_colonies_Gini(ht, raw_data)

    %% ----------------------------------------------------------------
    %  0.  Validate inputs
    %% ----------------------------------------------------------------
    if nargin < 1 || ~isstruct(ht) || ~isfield(ht, 'groups')
        error('First argument must be the ht struct from preprocess_pipeline_data.');
    end
    if nargin < 2 || ~isstruct(raw_data) || ~isfield(raw_data, 'processed')
        error(['Second argument must be your raw data struct.\n' ...
               'Usage: stats = plot_morphology_colonies_Gini(ht, data_26C_050626)']);
    end

    p         = ht.params;
    labels    = ht.labels;
    n_groups  = length(labels);
    incTime   = p.incTime;
    max_lag   = p.max_lag;
    out_dir   = p.out_dir;
    nr_plates = ht.qc.nr_plates;

    %% ----------------------------------------------------------------
    %  1.  Plate labels
    %% ----------------------------------------------------------------
    need_labels = true;
    if isfield(ht, 'plate_labels') && ~isempty(ht.plate_labels)
        first_lbl = '';
        for ii = 1:length(ht.plate_labels)
            if ~isempty(ht.plate_labels{ii})
                first_lbl = ht.plate_labels{ii};
                break;
            end
        end
        is_default = ~isempty(regexp(first_lbl, '^(Plate|P)\d+$', 'once'));
        if ~is_default
            need_labels = false;
            fprintf('[plot_morphology_colonies_Gini] Using cached ht.plate_labels (e.g. "%s").\n', first_lbl);
            fprintf('  To re-select: ht = rmfield(ht,''plate_labels''); then re-run.\n');
        else
            fprintf('[plot_morphology_colonies_Gini] Cached labels are defaults — rebuilding from metadata.\n');
        end
    end

    if need_labels
        if isfield(ht, 'metadata') && istable(ht.metadata) && height(ht.metadata) > 0
            meta_tbl = ht.metadata;
        elseif isfield(raw_data, 'metadata') && isfield(raw_data.metadata, 'original') && ...
               istable(raw_data.metadata.original) && height(raw_data.metadata.original) > 0
            meta_tbl = raw_data.metadata.original;
        else
            meta_tbl = table();
        end
        ht.plate_labels = build_plate_labels_interactive(meta_tbl, nr_plates);
        assignin('base', 'ht', ht);
    end

    plate_labels = ht.plate_labels;

    %% ----------------------------------------------------------------
    %  2.  Per-plate lag-time cache
    %% ----------------------------------------------------------------
    if ~isfield(ht, 'per_plate') || isempty(fieldnames(ht.per_plate))
        ht = build_per_plate_lag_cache(ht, raw_data, incTime, max_lag);
        assignin('base', 'ht', ht);
    end

    %% ----------------------------------------------------------------
    %  3.  Per-plate colony-size cache
    %% ----------------------------------------------------------------
    if ~isfield(ht, 'per_plate_size') || isempty(ht.per_plate_size) || ...
       all(cellfun(@isempty, ht.per_plate_size))
        fprintf('[plot_morphology_colonies_Gini] ht.per_plate_size not found — building from raw_data.\n');
        fprintf('  (Re-run preprocess_pipeline_data to cache this automatically.)\n');
        ht = build_per_plate_size_cache(ht, raw_data);
        assignin('base', 'ht', ht);
    else
        fprintf('[plot_morphology_colonies_Gini] Using ht.per_plate_size from preprocess_pipeline_data.\n');
    end

    %% ----------------------------------------------------------------
    %  4.  Pre-compute per-plate data for all groups
    %      pd(g) contains everything needed to draw any figure variant.
    %% ----------------------------------------------------------------
    fprintf('\n[plot_morphology_colonies_Gini] Computing per-plate Gini values...\n');

    pd(n_groups) = struct();   % plate_data array

    for g = 1:n_groups
        plate_ids = ht.groups(g).plate_indices;
        n_plates  = length(plate_ids);

        tick_labels = cell(n_plates, 1);
        for pi = 1:n_plates
            pidx = plate_ids(pi);
            if pidx <= length(plate_labels) && ~isempty(plate_labels{pidx})
                tick_labels{pi} = plate_labels{pidx};
            else
                tick_labels{pi} = sprintf('P%d', pidx);
            end
        end

        lag_vals  = cell(n_plates, 1);   % per-plate lag vectors
        size_vals = cell(n_plates, 1);   % per-plate size vectors
        N_col     = nan(n_plates, 1);    % colony count (same for lag and size)
        gini_lag  = nan(n_plates, 1);
        gini_size = nan(n_plates, 1);

        for pi = 1:n_plates
            pidx = plate_ids(pi);

            % --- Lag time ---
            fld_lag = sprintf('plate_%d', pidx);
            if isfield(ht.per_plate, fld_lag) && ~isempty(ht.per_plate.(fld_lag))
                lag = ht.per_plate.(fld_lag);
                lag_vals{pi} = lag;
                N_col(pi)   = numel(lag);
                if numel(lag) >= 3
                    gini_lag(pi) = gini_coeff(lag);
                end
            else
                N_col(pi) = 0;
            end

            % --- Colony size ---
            if pidx <= length(ht.per_plate_size) && ~isempty(ht.per_plate_size{pidx})
                sz = ht.per_plate_size{pidx}(:);
                sz = sz(isfinite(sz));
                size_vals{pi} = sz;
                if numel(sz) >= 3
                    gini_size(pi) = gini_coeff(sz);
                end
            end
        end

        pd(g).label       = labels{g};
        pd(g).fn          = matlab.lang.makeValidName(labels{g});
        pd(g).plate_ids   = plate_ids;
        pd(g).n_plates    = n_plates;
        pd(g).tick_labels = tick_labels;
        pd(g).lag_vals    = lag_vals;
        pd(g).size_vals   = size_vals;
        pd(g).N_col       = N_col;
        pd(g).gini_lag    = gini_lag;
        pd(g).gini_size   = gini_size;

        % Export per-group Gini tables to workspace
        T_lag = table((1:n_plates)', tick_labels(:), N_col, gini_lag, ...
                      'VariableNames',{'PlatePos','Label','N_colonies','Gini_lag'});
        assignin('base', ['GiniTable_lag_'  pd(g).fn], T_lag);

        T_size = table((1:n_plates)', tick_labels(:), N_col, gini_size, ...
                       'VariableNames',{'PlatePos','Label','N_colonies','Gini_size'});
        assignin('base', ['GiniTable_size_' pd(g).fn], T_size);
    end

    %% ----------------------------------------------------------------
    %  5.  Y-axis ranges
    %% ----------------------------------------------------------------
    y_lo_lag   = incTime;
    y_hi_lag   = max_lag + 6;
    y_gini_lag = max_lag + 2.0;

    all_sizes_pooled = [];
    for g = 1:n_groups
        for pi = 1:pd(g).n_plates
            if ~isempty(pd(g).size_vals{pi})
                all_sizes_pooled = [all_sizes_pooled; pd(g).size_vals{pi}]; %#ok<AGROW>
            end
        end
    end
    if isempty(all_sizes_pooled)
        y_hi_size = 5000;
    else
        y_hi_size = prctile(all_sizes_pooled, 99) * 1.25;
    end
    y_lo_size   = 0;
    y_gini_size = y_hi_size * 0.92;

    %% ----------------------------------------------------------------
    %  6.  Compute correlation tables (Spearman + Pearson, lag + size)
    %% ----------------------------------------------------------------
    [CorrTable_lag_sp,   pool_lag_sp]  = compute_corr_table(pd, n_groups, 'lag',  'Spearman');
    [CorrTable_lag_pe,   pool_lag_pe]  = compute_corr_table(pd, n_groups, 'lag',  'Pearson');
    [CorrTable_size_sp,  pool_size_sp] = compute_corr_table(pd, n_groups, 'size', 'Spearman');
    [CorrTable_size_pe,  pool_size_pe] = compute_corr_table(pd, n_groups, 'size', 'Pearson');

    % Print and export all four tables
    print_corr_table('Lag time',     'Spearman', CorrTable_lag_sp,  pool_lag_sp);
    print_corr_table('Lag time',     'Pearson',  CorrTable_lag_pe,  pool_lag_pe);
    print_corr_table('Colony size',  'Spearman', CorrTable_size_sp, pool_size_sp);
    print_corr_table('Colony size',  'Pearson',  CorrTable_size_pe, pool_size_pe);

    assignin('base', 'GiniCorr_NvsGini_lag_Spearman',  CorrTable_lag_sp);
    assignin('base', 'GiniCorr_NvsGini_lag_Pearson',   CorrTable_lag_pe);
    assignin('base', 'GiniCorr_NvsGini_size_Spearman', CorrTable_size_sp);
    assignin('base', 'GiniCorr_NvsGini_size_Pearson',  CorrTable_size_pe);

    fprintf('Correlation tables saved to workspace:\n');
    fprintf('  GiniCorr_NvsGini_lag_Spearman\n');
    fprintf('  GiniCorr_NvsGini_lag_Pearson\n');
    fprintf('  GiniCorr_NvsGini_size_Spearman\n');
    fprintf('  GiniCorr_NvsGini_size_Pearson\n\n');

    %% ----------------------------------------------------------------
    %  7.  Build stats struct
    %% ----------------------------------------------------------------
    stats = struct();
    for g = 1:n_groups
        fn  = pd(g).fn;
        all_lag  = vertcat_cell(pd(g).lag_vals);
        all_size = vertcat_cell(pd(g).size_vals);
        stats.(fn).label       = pd(g).label;
        stats.(fn).N           = sum(pd(g).N_col(~isnan(pd(g).N_col)));
        stats.(fn).lag         = all_lag;
        stats.(fn).median_lag  = safe_median(all_lag);
        stats.(fn).iqr_lag     = safe_iqr(all_lag);
        stats.(fn).size        = all_size;
        stats.(fn).median_size = safe_median(all_size);
        stats.(fn).iqr_size    = safe_iqr(all_size);
    end

    %% ----------------------------------------------------------------
    %  8.  Draw four figures
    %% ----------------------------------------------------------------
    ts = datestr(now, 'yyyymmdd_HHMMSS');

    cfg_lag_sp  = make_cfg('lag',  'Spearman', y_lo_lag,  y_hi_lag,  y_gini_lag,  CorrTable_lag_sp);
    cfg_lag_pe  = make_cfg('lag',  'Pearson',  y_lo_lag,  y_hi_lag,  y_gini_lag,  CorrTable_lag_pe);
    cfg_size_sp = make_cfg('size', 'Spearman', y_lo_size, y_hi_size, y_gini_size, CorrTable_size_sp);
    cfg_size_pe = make_cfg('size', 'Pearson',  y_lo_size, y_hi_size, y_gini_size, CorrTable_size_pe);

    hfig1 = draw_figure(pd, n_groups, ht.colors, cfg_lag_sp,  labels, n_groups);
    hfig2 = draw_figure(pd, n_groups, ht.colors, cfg_size_sp, labels, n_groups);
    hfig3 = draw_figure(pd, n_groups, ht.colors, cfg_lag_pe,  labels, n_groups);
    hfig4 = draw_figure(pd, n_groups, ht.colors, cfg_size_pe, labels, n_groups);

    %% ----------------------------------------------------------------
    %  9.  Save all four figures
    %% ----------------------------------------------------------------
    fnames = { ...
        sprintf('morphology_Gini_lagTime_Spearman_%s',  ts), ...
        sprintf('morphology_Gini_colonySize_Spearman_%s', ts), ...
        sprintf('morphology_Gini_lagTime_Pearson_%s',   ts), ...
        sprintf('morphology_Gini_colonySize_Pearson_%s',  ts)};
    hfigs = {hfig1, hfig2, hfig3, hfig4};

    for k = 1:4
        fp = fullfile(out_dir, fnames{k});
        save_fig(hfigs{k}, [fp '.fig'], 'fig');
        save_fig(hfigs{k}, [fp '.png'], 'png');
        save_fig(hfigs{k}, [fp '.svg'], 'svg');
        save_fig(hfigs{k}, [fp '.pdf'], 'pdf');
    end

    %% ----------------------------------------------------------------
    %  10.  Write log file
    %% ----------------------------------------------------------------
    log_path = fullfile(out_dir, sprintf('morphology_Gini_log_%s.txt', ts));
    fid = fopen(log_path, 'w');
    if fid ~= -1
        fprintf(fid, '================================================\n');
        fprintf(fid, '  plot_morphology_colonies_Gini  —  Run log\n');
        fprintf(fid, '  %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
        fprintf(fid, '================================================\n\n');

        fprintf(fid, 'Output directory: %s\n\n', out_dir);

        fprintf(fid, 'Files saved:\n');
        for k = 1:4
            for ext = {'.fig','.png','.svg','.pdf'}
                fprintf(fid, '  %s%s\n', fnames{k}, ext{1});
            end
        end
        fprintf(fid, '\n');

        % Per-group summary
        fprintf(fid, '------------------------------------------------\n');
        fprintf(fid, '  Per-group summary\n');
        fprintf(fid, '------------------------------------------------\n');
        fprintf(fid, '  %-12s  %6s  %12s  %10s  %12s  %10s\n', ...
            'Group', 'N', 'Med lag(h)', 'IQR lag', 'Med size(px)', 'IQR size');
        for g = 1:n_groups
            fn = pd(g).fn;
            fprintf(fid, '  %-12s  %6d  %12.2f  %10.2f  %12.1f  %10.1f\n', ...
                pd(g).label, stats.(fn).N, ...
                stats.(fn).median_lag,  stats.(fn).iqr_lag, ...
                stats.(fn).median_size, stats.(fn).iqr_size);
        end
        fprintf(fid, '\n');

        % Correlation tables
        corr_configs = { ...
            'Lag time',    'Spearman', CorrTable_lag_sp;  ...
            'Lag time',    'Pearson',  CorrTable_lag_pe;  ...
            'Colony size', 'Spearman', CorrTable_size_sp; ...
            'Colony size', 'Pearson',  CorrTable_size_pe};

        for ci = 1:size(corr_configs, 1)
            met  = corr_configs{ci,1};
            ctyp = corr_configs{ci,2};
            T    = corr_configs{ci,3};
            if strcmp(ctyp,'Spearman'), rlbl = 'rho'; else, rlbl = 'r'; end
            fprintf(fid, '------------------------------------------------\n');
            fprintf(fid, '  %s corr: N_colonies vs Gini (%s)\n', ctyp, met);
            fprintf(fid, '------------------------------------------------\n');
            fprintf(fid, '  %-12s  %8s  %8s  %8s\n', 'Group','N plates', rlbl, 'p-value');
            for k = 1:height(T)
                if isnan(T.r(k))
                    fprintf(fid, '  %-12s  %8d  %8s  %8s\n', T.Group{k}, T.N_plates(k), 'NaN','NaN');
                else
                    fprintf(fid, '  %-12s  %8d  %8.3f  %8s  %s\n', ...
                        T.Group{k}, T.N_plates(k), T.r(k), fmt_p(T.P_value(k)), pval2star(T.P_value(k)));
                end
            end
            fprintf(fid, '\n');
        end

        % Per-plate Gini values
        fprintf(fid, '------------------------------------------------\n');
        fprintf(fid, '  Per-plate Gini coefficients\n');
        fprintf(fid, '------------------------------------------------\n');
        for g = 1:n_groups
            fprintf(fid, '  Group: %s\n', pd(g).label);
            fprintf(fid, '    %-20s  %6s  %8s  %8s\n', 'Plate','N','Gini_lag','Gini_sz');
            for pi = 1:pd(g).n_plates
                fprintf(fid, '    %-20s  %6d  %8.3f  %8.3f\n', ...
                    pd(g).tick_labels{pi}, pd(g).N_col(pi), ...
                    pd(g).gini_lag(pi), pd(g).gini_size(pi));
            end
            fprintf(fid, '\n');
        end

        fclose(fid);
        fprintf('Log saved  ->  %s\n', log_path);
    else
        warning('Could not write log file to: %s', log_path);
    end

    fprintf('\nDone — outputs saved to: %s\n', out_dir);

end


%% ====================================================================
%  LOCAL: make_cfg
%  Packages the draw parameters for one figure variant into a struct.
%% ====================================================================
function cfg = make_cfg(metric, corr_type, y_lo, y_hi, y_gini, corr_table)
    cfg.metric      = metric;      % 'lag' | 'size'
    cfg.corr_type   = corr_type;   % 'Spearman' | 'Pearson'
    cfg.y_lo        = y_lo;
    cfg.y_hi        = y_hi;
    cfg.y_gini      = y_gini;
    cfg.corr_table  = corr_table;  % table with Group / Spearman_r (or Pearson_r) / P_value
    if strcmp(metric, 'lag')
        cfg.ylabel_str = 'Lag time (h)';
        cfg.fig_name   = sprintf('Lag time — colonies & Gini (%s)', corr_type);
    else
        cfg.ylabel_str = 'Final colony size (px)';
        cfg.fig_name   = sprintf('Colony size — colonies & Gini (%s)', corr_type);
    end
end


%% ====================================================================
%  LOCAL: draw_figure
%  Draws one violin+Gini figure for a given metric + correlation type.
%% ====================================================================
function hfig = draw_figure(pd, n_groups, group_colors, cfg, labels, n_grp)

    fig_w = 14;   fig_h = 11;

    if n_grp <= 4
        n_rows_fig = ceil(n_grp / 2);
        n_cols_fig = min(2, n_grp);
    else
        n_rows_fig = n_grp;
        n_cols_fig = 1;
    end

    hfig = figure('Name', cfg.fig_name, 'Color','w', ...
                  'Units','inches','Position',[1 1 fig_w fig_h], ...
                  'PaperUnits','inches','PaperSize',[fig_w fig_h], ...
                  'PaperPosition',[0 0 fig_w fig_h]);

    use_lag = strcmp(cfg.metric, 'lag');

    for g = 1:n_grp
        ax = subplot(n_rows_fig, n_cols_fig, g);
        hold(ax, 'on');

        n_plates    = pd(g).n_plates;
        tick_labels = pd(g).tick_labels;
        violin_col  = group_colors(g, :);
        dot_col     = [0.15 0.15 0.15];

        % ---- Per-group y-axis range ----------------------------------------
        % Collect all values in this group to set the y range independently
        % for each subplot, so each group uses its own data scale.
        all_grp_vals = [];
        for pi_tmp = 1:n_plates
            if use_lag
                v_tmp = pd(g).lag_vals{pi_tmp};
            else
                v_tmp = pd(g).size_vals{pi_tmp};
            end
            if ~isempty(v_tmp)
                all_grp_vals = [all_grp_vals; v_tmp(:)]; %#ok<AGROW>
            end
        end

        if isempty(all_grp_vals)
            y_lo_grp = cfg.y_lo;
            y_hi_grp = cfg.y_hi;
        else
            y_lo_grp = cfg.y_lo;                         % keep fixed lower bound
            y_hi_grp = max(all_grp_vals);                % true data max this group
        end
        % All Gini labels in this group start at the same y: group max + 5%.
        % Extra headroom (20%) keeps the rotated labels inside the axes.
        data_span  = y_hi_grp - y_lo_grp;
        y_lbl_grp  = y_hi_grp + 0.05 * data_span;   % shared label baseline
        label_gap  = 0.20 * data_span;
        y_ax_hi    = y_hi_grp + label_gap;
        % -----------------------------------------------------------------------

        for pi = 1:n_plates
            if use_lag
                vals = pd(g).lag_vals{pi};
            else
                vals = pd(g).size_vals{pi};
            end

            if isempty(vals), continue; end

            n_col = numel(vals);

            % Violin
            if n_col > 5
                try
                    if use_lag
                        % Lag time: bounded support keeps the violin within
                        % the valid [incTime, max_lag] range.
                        bw = max(0.4, std(vals) * 0.18);
                        [f, yi] = ksdensity(vals, 'Bandwidth', bw, ...
                                            'Support', [cfg.y_lo, cfg.y_hi]);
                    else
                        % Colony size: evaluate only over the data's own range
                        % (+/- 3 bw) to avoid near-zero density at evaluation
                        % points far from the data mass.
                        bw = max(10, std(vals) * 0.18);
                        yi = linspace(max(cfg.y_lo, min(vals) - 3*bw), ...
                                      min(y_ax_hi,  max(vals) + 3*bw), 200)';
                        f  = ksdensity(vals, yi, 'Bandwidth', bw);
                    end
                    f = max(f, 0);              % clamp negatives from ksdensity
                    % Smooth with a moving average (~5% of grid length) to
                    % fill bimodal pinch points before the floor is applied.
                    n_sm = max(3, round(length(f) * 0.05));
                    sm_win = ones(n_sm, 1) / n_sm;
                    f = conv(f, sm_win, 'same');
                    f = max(f, 0);
                    % 4% floor ensures the polygon is always visibly wide
                    % — too small to distort shape, large enough to close gaps.
                    f = max(f, max(f) * 0.04);
                    if max(f) > eps
                        f = f ./ max(f) * 0.38;
                        patch(ax, [pi-f, fliplr(pi+f)], [yi, fliplr(yi)], ...
                              violin_col, 'EdgeColor','none','FaceAlpha',0.45);
                    end
                catch
                end
            end

            jx = pi + 0.06 * randn(n_col, 1);
            scatter(ax, jx, vals, 7, dot_col, 'filled', 'MarkerFaceAlpha', 0.40);
            line(ax, [pi-0.25 pi+0.25], [median(vals) median(vals)], ...
                 'Color', violin_col*0.6, 'LineWidth', 2.0);

            % Gini label — placed just above this plate's own maximum so it
            % never overlaps the data points and is consistent per plate.
            if n_col >= 3
                if use_lag
                    G = pd(g).gini_lag(pi);
                else
                    G = pd(g).gini_size(pi);
                end
                if ~isnan(G)
                    text(ax, pi, y_lbl_grp, sprintf('G=%.2f N=%d', G, n_col), ...
                         'Color',[0.05 0.20 0.70], ...
                         'HorizontalAlignment','left', ...
                         'VerticalAlignment','bottom', ...
                         'Rotation', 30, ...
                         'FontSize', 9, 'FontWeight', 'bold');
                end
            end
        end

        % Correlation result for this group from the pre-computed table
        grp_label = pd(g).label;
        row_match = strcmp(cfg.corr_table.Group, grp_label);
        if any(row_match)
            rho_g  = cfg.corr_table.r(row_match);
            pval_g = cfg.corr_table.P_value(row_match);
        else
            rho_g = NaN;  pval_g = NaN;
        end

        sig_str = pval2star(pval_g);
        if isnan(rho_g)
            corr_line = 'insufficient data';
        else
            if strcmp(cfg.corr_type, 'Spearman')
                corr_line = sprintf('Spearman r_s=%.3f, p=%s (%s)', rho_g, fmt_p(pval_g, '%.3g'), sig_str);
            else
                corr_line = sprintf('Pearson r=%.3f, p=%s (%s)', rho_g, fmt_p(pval_g, '%.3g'), sig_str);
            end
        end

        ylabel(ax, cfg.ylabel_str, 'FontSize', 14, 'FontWeight', 'bold');
        title(ax, {grp_label, corr_line}, 'FontSize', 14, 'FontWeight', 'bold');
        xlim(ax, [0.5, n_plates+0.5]);
        ylim(ax, [y_lo_grp, y_ax_hi]);
        set(ax, 'XTick', 1:n_plates, 'XTickLabel', tick_labels, ...
                'TickLabelInterpreter','none','FontSize', 11, ...
                'Box','off','GridAlpha',0.3);
        xtickangle(ax, 35);
        grid(ax, 'on');

        if g == n_grp || (n_cols_fig==2 && g >= n_grp-1)
            xlabel(ax, 'Biological replicate', 'FontSize', 14, 'FontWeight', 'bold');
        end
    end
end


%% ====================================================================
%  LOCAL: compute_corr_table
%  Computes correlation (Spearman or Pearson) between N_colonies and
%  Gini for each group + pooled, returns a table and pooled-row struct.
%% ====================================================================
function [T, pool] = compute_corr_table(pd, n_groups, metric, corr_type)

    use_lag = strcmp(metric, 'lag');

    groups_cell = cell(n_groups+1, 1);
    n_plates_vec = zeros(n_groups+1, 1);
    r_vec        = nan(n_groups+1, 1);
    p_vec        = nan(n_groups+1, 1);

    pool_N    = [];
    pool_Gini = [];

    for g = 1:n_groups
        N_col = pd(g).N_col;
        if use_lag
            gini_vec = pd(g).gini_lag;
        else
            gini_vec = pd(g).gini_size;
        end

        ok = ~isnan(N_col) & ~isnan(gini_vec);

        groups_cell{g}   = pd(g).label;
        n_plates_vec(g)  = nnz(ok);

        if nnz(ok) >= 3
            [r_vec(g), p_vec(g)] = corr(N_col(ok), gini_vec(ok), 'Type', corr_type);
        end

        pool_N    = [pool_N;    N_col(ok)];    %#ok<AGROW>
        pool_Gini = [pool_Gini; gini_vec(ok)]; %#ok<AGROW>
    end

    % Pooled row
    groups_cell{n_groups+1}  = 'POOLED';
    n_plates_vec(n_groups+1) = length(pool_N);
    if length(pool_N) >= 3
        [r_vec(n_groups+1), p_vec(n_groups+1)] = corr(pool_N, pool_Gini, 'Type', corr_type);
    end

    T = table(groups_cell, n_plates_vec, r_vec, p_vec, ...
              'VariableNames', {'Group','N_plates','r','P_value'});

    pool.N    = pool_N;
    pool.Gini = pool_Gini;
    pool.r    = r_vec(n_groups+1);
    pool.p    = p_vec(n_groups+1);
end


%% ====================================================================
%  LOCAL: print_corr_table
%% ====================================================================
function print_corr_table(metric_name, corr_type, T, pool)
    if strcmp(corr_type, 'Spearman')
        r_label = 'rho';
    else
        r_label = 'r';
    end

    fprintf('\n');
    fprintf('========================================================\n');
    fprintf('  %s corr: Colony count vs Gini (%s)\n', corr_type, metric_name);
    fprintf('========================================================\n');
    fprintf('  %-10s  %6s  %10s  %8s\n', 'Group', 'n plates', r_label, 'p-value');
    fprintf('  %s\n', repmat('-', 1, 44));

    for k = 1:height(T)
        if isnan(T.r(k))
            fprintf('  %-10s  %6d  %10s  %8s\n', T.Group{k}, T.N_plates(k), 'NaN', 'NaN');
        else
            sig = pval2star(T.P_value(k));
            fprintf('  %-10s  %6d  %10.3f  %8s  %s\n', ...
                T.Group{k}, T.N_plates(k), T.r(k), fmt_p(T.P_value(k)), sig);
        end
    end
    fprintf('========================================================\n\n');
end


%% ====================================================================
%  LOCAL: build_plate_labels_interactive
%% ====================================================================
function plate_labels = build_plate_labels_interactive(meta_tbl, nr_plates)

    plate_labels = cell(nr_plates, 1);

    if isempty(meta_tbl) || ~istable(meta_tbl) || height(meta_tbl) == 0 || ...
       width(meta_tbl) == 0
        for i = 1:nr_plates
            plate_labels{i} = sprintf('Plate%d', i);
        end
        return;
    end

    col_names = meta_tbl.Properties.VariableNames;
    n_cols    = length(col_names);

    fprintf('\n============================================================\n');
    fprintf('  X-axis label column selection\n');
    fprintf('  ht.metadata has %d rows x %d columns.\n', height(meta_tbl), n_cols);
    fprintf('============================================================\n\n');
    fprintf('  %-4s  %-28s  %s\n', 'No.', 'Column name', 'Unique values (preview)');
    fprintf('  %s\n', repmat('-',1,70));

    for k = 1:n_cols
        try
            raw_vals = meta_tbl.(col_names{k});
            if iscell(raw_vals)
                uv = unique(raw_vals(~cellfun(@isempty, raw_vals)));
                preview = strjoin(uv(1:min(5,end)), ' | ');
            elseif isnumeric(raw_vals) || islogical(raw_vals)
                uv = unique(raw_vals(~isnan(double(raw_vals))));
                preview = strtrim(num2str(uv(1:min(5,end))'));
            else
                uv = unique(raw_vals);
                preview = char(uv(1));
            end
        catch
            preview = '(preview unavailable)';
        end
        fprintf('  [%2d]  %-28s  %s\n', k, col_names{k}, preview);
    end

    fprintf('\n');

    bio_col = '';
    fprintf('  Step A — BIOLOGICAL replicate column (enter 0 to use plate indices)\n\n');
    while true
        raw = strtrim(input('  Biological replicate column number: ', 's'));
        idx = str2double(raw);
        if ~isnan(idx) && idx == 0
            fprintf('  -> Skipped. Using plate indices.\n\n');  break;
        end
        if ~isnan(idx) && floor(idx)==idx && idx>=1 && idx<=n_cols
            bio_col = col_names{idx};
            fprintf('  -> Biological replicate: "%s"\n\n', bio_col);  break;
        end
        fprintf('  WARNING: enter a number between 0 and %d.\n', n_cols);
    end

    tech_col = '';
    fprintf('  Step B — TECHNICAL replicate column (enter 0 to skip)\n\n');
    while true
        raw = strtrim(input('  Technical replicate column number (0 to skip): ', 's'));
        idx = str2double(raw);
        if ~isnan(idx) && idx == 0
            fprintf('  -> Skipped.\n\n');  break;
        end
        if ~isnan(idx) && floor(idx)==idx && idx>=1 && idx<=n_cols
            tech_col = col_names{idx};
            fprintf('  -> Technical replicate: "%s"\n\n', tech_col);  break;
        end
        fprintf('  WARNING: enter a number between 0 and %d.\n', n_cols);
    end

    dil_col = '';
    fprintf('  Step C — DILUTION column (enter 0 to skip)\n\n');
    while true
        raw = strtrim(input('  Dilution column number (0 to skip): ', 's'));
        idx = str2double(raw);
        if ~isnan(idx) && idx == 0
            fprintf('  -> Skipped.\n\n');  break;
        end
        if ~isnan(idx) && floor(idx)==idx && idx>=1 && idx<=n_cols
            dil_col = col_names{idx};
            fprintf('  -> Dilution: "%s"\n\n', dil_col);  break;
        end
        fprintf('  WARNING: enter a number between 0 and %d.\n', n_cols);
    end

    n_meta = height(meta_tbl);
    for i = 1:nr_plates
        if i > n_meta
            plate_labels{i} = sprintf('P%d', i);  continue;
        end
        if isempty(bio_col)
            lbl = sprintf('P%d', i);
        else
            lbl = val2str(meta_tbl.(bio_col)(i));
        end
        if ~isempty(tech_col)
            lbl = [lbl '_' val2str(meta_tbl.(tech_col)(i))];
        end
        if ~isempty(dil_col)
            lbl = [lbl '_' val2str(meta_tbl.(dil_col)(i))];
        end
        plate_labels{i} = lbl;
    end

    fprintf('  Labels assigned — first 8 plates:\n');
    for i = 1:min(8, nr_plates)
        fprintf('    Plate %2d  ->  "%s"\n', i, plate_labels{i});
    end
    if nr_plates > 8
        fprintf('    ... (%d more)\n', nr_plates-8);
    end
    fprintf('\n');
end


%% ====================================================================
%  LOCAL: build_per_plate_lag_cache
%% ====================================================================
function ht = build_per_plate_lag_cache(ht, raw_data, incTime, max_lag)

    nr_plates = ht.qc.nr_plates;
    min_col   = ht.params.min_col;
    max_col   = ht.params.max_col;
    cens_tol  = 0.25;
    ht.per_plate = struct();
    n_ok = 0;  n_bad = 0;

    fprintf('\n[build_per_plate_lag_cache]  %d plates...\n', nr_plates);

    for i = 1:nr_plates
        fld = sprintf('plate_%d', i);

        if i > length(raw_data.processed)
            ht.per_plate.(fld) = [];  n_bad = n_bad+1;  continue;
        end

        proc = raw_data.processed{i};

        if ~isfield(proc,'growth_quant') || ~proc.growth_quant
            ht.per_plate.(fld) = [];  n_bad = n_bad+1;  continue;
        end

        if ~isfield(proc,'colonies') || ~isfield(proc.colonies,'new') || ...
           ~isfield(proc.colonies.new,'lag_time')
            ht.per_plate.(fld) = [];  n_bad = n_bad+1;  continue;
        end

        lag_raw = proc.colonies.new.lag_time(:) + incTime;
        lag_raw(lag_raw >= max_lag - cens_tol) = NaN;
        n_valid = sum(~isnan(lag_raw));

        if n_valid < min_col || n_valid > max_col
            ht.per_plate.(fld) = [];  n_bad = n_bad+1;  continue;
        end

        ht.per_plate.(fld) = lag_raw(~isnan(lag_raw));
        n_ok = n_ok + 1;
    end

    fprintf('[build_per_plate_lag_cache]  Done: %d OK, %d skipped\n\n', n_ok, n_bad);
    assignin('base', 'ht', ht);
end


%% ====================================================================
%  LOCAL: build_per_plate_size_cache
%  Fallback: extracts per-plate final colony size from raw_data when
%  ht.per_plate_size is absent.
%% ====================================================================
function ht = build_per_plate_size_cache(ht, raw_data)

    nr_plates      = ht.qc.nr_plates;
    size_thr       = ht.params.size_threshold;
    per_plate_size = cell(nr_plates, 1);
    n_ok = 0;  n_bad = 0;

    fprintf('\n[build_per_plate_size_cache]  %d plates...\n', nr_plates);

    for i = 1:nr_plates
        if i > length(raw_data.processed)
            n_bad = n_bad + 1;  continue;
        end

        proc = raw_data.processed{i};

        if ~isfield(proc,'growth_quant') || ~proc.growth_quant
            n_bad = n_bad + 1;  continue;
        end

        if ~isfield(proc,'colonies') || ~isfield(proc.colonies,'new') || ...
           ~isfield(proc.colonies.new,'lag_time')
            n_bad = n_bad + 1;  continue;
        end

        colonies = proc.colonies;
        c_new    = colonies.new;
        n_col    = length(c_new.lag_time);
        if n_col == 0,  n_bad = n_bad + 1;  continue;  end

        rp_clean = [];
        if isfield(colonies, 'region_props_clean')
            rp_clean = colonies.region_props_clean;
        elseif isfield(colonies, 'region_props')
            rp_all = colonies.region_props;
            if isfield(colonies, 'flag_colony_ok') && ...
               length(colonies.flag_colony_ok) == height(rp_all)
                cidx = find(colonies.flag_colony_ok);
                if length(cidx) == n_col
                    rp_clean = rp_all(cidx, :);
                end
            end
        end

        if isempty(rp_clean) || height(rp_clean) ~= n_col
            n_bad = n_bad + 1;  continue;
        end

        sz = rp_clean.Area;
        sz(sz < size_thr) = NaN;

        % Also exclude colonies with undetectable lag time (censored).
        lt_raw_s = proc.colonies.new.lag_time(:) + ht.params.incTime;
        cens_tol_s = 0.5 * ht.params.lag_step;
        is_cens_s = ~isnan(lt_raw_s) & (lt_raw_s >= ht.params.max_lag - cens_tol_s);
        sz(is_cens_s) = NaN;

        per_plate_size{i} = sz(~isnan(sz));
        n_ok = n_ok + 1;
    end

    ht.per_plate_size = per_plate_size;
    fprintf('[build_per_plate_size_cache]  Done: %d OK, %d skipped\n\n', n_ok, n_bad);
    assignin('base', 'ht', ht);
end


%% ====================================================================
%  LOCAL helpers
%% ====================================================================
function s = pval2star(p)
    if isnan(p),       s = '';
    elseif p < 0.001,  s = '***';
    elseif p < 0.01,   s = '**';
    elseif p < 0.05,   s = '*';
    else,              s = 'ns';
    end
end

% Text display for a p-value. Pearson/Spearman p (via corr()) can
% legitimately underflow to exact 0.0 in double precision for large n /
% strong correlation — the true value is just too small to represent,
% not actually zero — shown as "<1e-300" instead. Only affects display;
% pval2star(...) above must keep receiving the raw numeric p.
function s = fmt_p(p, fmt)
    if nargin < 2, fmt = '%.4f'; end
    if isnan(p)
        s = 'NaN';
    elseif p == 0
        s = '<1e-300';
    else
        s = sprintf(fmt, p);
    end
end

function G = gini_coeff(x)
    x = x(isfinite(x) & x >= 0);
    n = numel(x);
    if n < 2 || sum(x) == 0, G = NaN; return; end
    x = sort(x(:));
    i = (1:n)';
    G = (2*sum(i.*x)) / (n*sum(x)) - (n+1)/n;
    G = max(0, min(1, G));
end

function v = vertcat_cell(c)
    v = [];
    for k = 1:length(c)
        if ~isempty(c{k}), v = [v; c{k}(:)]; end %#ok<AGROW>
    end
end

function m = safe_median(v)
    v = v(isfinite(v));
    if isempty(v), m = NaN; else, m = median(v); end
end

function q = safe_iqr(v)
    v = v(isfinite(v));
    if isempty(v), q = NaN; else, q = iqr(v); end
end

function s = val2str(v)
    if iscell(v),          s = strtrim(char(v{1}));
    elseif isnumeric(v),   s = strtrim(num2str(v));
    elseif ischar(v),      s = strtrim(v);
    elseif isstring(v),    s = strtrim(char(v));
    else,                  s = '?';
    end
end

function save_fig(hfig, fpath, fmt)
    try
        switch fmt
            case 'fig', savefig(hfig, fpath);
            case 'png', print(hfig, fpath, '-dpng', '-r300');
            case 'svg', print(hfig, fpath, '-dsvg');
            case 'pdf', print(hfig, fpath, '-dpdf', '-bestfit');
        end
        fprintf('Saved .%s  ->  %s\n', fmt, fpath);
    catch ME
        warning('Could not save .%s: %s', fmt, ME.message);
    end
end

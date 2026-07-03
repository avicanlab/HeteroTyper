%% HeteroTyper Pipeline
%  plot_lagTime_colonies_Gini  —  Violin + scatter of lag time per plate,
%  grouped by sample time. X-axis shows bio+tech replicate labels read
%  directly from ht.metadata (the full metadata table stored in ht).
%
%%  USAGE
%  -----
%    stats = plot_lagTime_colonies_Gini(ht, your_data_variable)
%
%  Example:
%    stats = plot_lagTime_colonies_Gini(ht, data_26C_050626)
%
%  On first call the function asks:
%    A) Which metadata column = biological replicate  (e.g. Set)
%    B) Which metadata column = technical replicate   (e.g. Replicate)
%  Labels are combined as "Set1_R1" on the x-axis and cached in ht so
%  subsequent calls skip the prompt.
%
%  To re-select columns:
%    ht = rmfield(ht, 'plate_labels');  then re-run.

function stats = plot_lagTime_colonies_Gini(ht, raw_data)

    %% ----------------------------------------------------------------
    %  0.  Validate inputs
    %% ----------------------------------------------------------------
    if nargin < 1 || ~isstruct(ht) || ~isfield(ht, 'groups')
        error('First argument must be the ht struct from preprocess_pipeline_data.');
    end
    if nargin < 2 || ~isstruct(raw_data) || ~isfield(raw_data, 'processed')
        error(['Second argument must be your raw data struct.\n' ...
               'Usage: stats = plot_lagTime_colonies_Gini(ht, data_26C_050626)']);
    end

    p        = ht.params;
    labels   = ht.labels;
    n_groups = length(labels);
    incTime  = p.incTime;
    max_lag  = p.max_lag;
    out_dir  = p.out_dir;
    nr_plates = ht.qc.nr_plates;

    %% ----------------------------------------------------------------
    %  1.  Plate labels — built interactively from ht.metadata each time
    %      unless already cached with real (non-default) values.
    %      Cache is invalidated if all labels look like "Plate<N>" or "P<N>".
    %% ----------------------------------------------------------------

    % Check if cached labels exist AND contain real metadata values
    need_labels = true;
    if isfield(ht, 'plate_labels') && ~isempty(ht.plate_labels)
        % Inspect first non-empty label — if it matches Plate<N> or P<N> pattern,
        % the cache is stale (defaults were stored, not real metadata)
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
            fprintf('[plot_lagTime_colonies_Gini] Using cached ht.plate_labels (e.g. "%s").\n', first_lbl);
            fprintf('  To re-select: ht = rmfield(ht,''plate_labels''); then re-run.\n');
        else
            fprintf('[plot_lagTime_colonies_Gini] Cached labels are defaults ("%s") — rebuilding from metadata.\n', first_lbl);
        end
    end

    if need_labels
        % Load metadata: prefer ht.metadata (patched preprocess), else raw_data
        if isfield(ht, 'metadata') && istable(ht.metadata) && height(ht.metadata) > 0
            meta_tbl = ht.metadata;
            fprintf('[plot_lagTime_colonies_Gini] Reading metadata from ht.metadata (%d rows, %d cols).\n', ...
                height(meta_tbl), width(meta_tbl));
        elseif isfield(raw_data, 'metadata') && isfield(raw_data.metadata, 'original') && ...
               istable(raw_data.metadata.original) && height(raw_data.metadata.original) > 0
            meta_tbl = raw_data.metadata.original;
            fprintf('[plot_lagTime_colonies_Gini] Reading metadata from raw_data.metadata.original (%d rows, %d cols).\n', ...
                height(meta_tbl), width(meta_tbl));
        else
            meta_tbl = table();
            fprintf('[plot_lagTime_colonies_Gini] WARNING: no metadata found — will use plate indices.\n');
        end

        ht.plate_labels = build_plate_labels_interactive(meta_tbl, nr_plates);
        assignin('base', 'ht', ht);
    end

    plate_labels = ht.plate_labels;   % {nr_plates x 1} cell of strings

    %% ----------------------------------------------------------------
    %  2.  Per-plate lag-time cache
    %% ----------------------------------------------------------------
    if ~isfield(ht, 'per_plate') || isempty(fieldnames(ht.per_plate))
        ht = build_per_plate_cache(ht, raw_data, incTime, max_lag);
        assignin('base', 'ht', ht);
    end

    %% ----------------------------------------------------------------
    %  3.  Figure layout
    %% ----------------------------------------------------------------
    fig_w = 14;   fig_h = 11;

    hfig = figure('Name','Lag time — colonies & Gini','Color','w', ...
                  'Units','inches','Position',[1 1 fig_w fig_h], ...
                  'PaperUnits','inches','PaperSize',[fig_w fig_h], ...
                  'PaperPosition',[0 0 fig_w fig_h]);

    if n_groups <= 4
        n_rows_fig = ceil(n_groups / 2);
        n_cols_fig = min(2, n_groups);
    else
        n_rows_fig = n_groups;
        n_cols_fig = 1;
    end

    stats        = struct();
    group_colors = ht.colors;
    y_lo         = incTime;
    y_hi         = max_lag + 6;
    y_gini       = max_lag + 2.0;   % fixed y for all Gini labels

    % Accumulators for pooled (across all groups) correlation
    pool_N    = [];   % colony counts across ALL plates
    pool_Gini = [];   % Gini values  across ALL plates
    pool_grp  = {};   % group label  for each plate (for reference)

    % Per-group correlation results table (filled inside loop)
    corr_group  = {};
    corr_rho    = [];
    corr_pval   = [];
    corr_nplates = [];

    %% ----------------------------------------------------------------
    %  4.  Plot loop — one subplot per time group
    %% ----------------------------------------------------------------
    for g = 1:n_groups

        ax = subplot(n_rows_fig, n_cols_fig, g);
        hold(ax, 'on');

        grp_label = labels{g};
        fn        = matlab.lang.makeValidName(grp_label);

        plate_ids = ht.groups(g).plate_indices;
        n_plates  = length(plate_ids);

        if n_plates == 0
            title(ax, grp_label);
            continue;
        end

        % ---- X-axis label: bio+tech rep from plate_labels ----
        tick_labels = cell(n_plates, 1);
        for pi = 1:n_plates
            pidx = plate_ids(pi);
            if pidx <= length(plate_labels) && ~isempty(plate_labels{pidx})
                tick_labels{pi} = plate_labels{pidx};
            else
                tick_labels{pi} = sprintf('P%d', pidx);
            end
        end

        % ---- Per-plate rendering ----
        all_lag    = [];
        all_pos    = [];
        plate_N    = nan(n_plates, 1);
        plate_Gini = nan(n_plates, 1);

        violin_col = group_colors(g, :);
        dot_col    = [0.15 0.15 0.15];

        for pi = 1:n_plates

            pidx = plate_ids(pi);
            fld  = sprintf('plate_%d', pidx);

            if isfield(ht.per_plate, fld)
                lag = ht.per_plate.(fld);
            else
                lag = [];
            end

            if isempty(lag)
                plate_N(pi) = 0;
                fprintf('  %s | %s (plate %d): no data\n', grp_label, tick_labels{pi}, pidx);
                continue;
            end

            n_col       = numel(lag);
            plate_N(pi) = n_col;

            % violin
            if n_col > 5
                try
                    bw = max(0.4, std(lag) * 0.18);
                    [f, yi] = ksdensity(lag, 'Bandwidth', bw, 'Support', [y_lo y_hi]);
                    f = f ./ max(f) * 0.38;
                    patch(ax, [pi-f, fliplr(pi+f)], [yi, fliplr(yi)], ...
                          violin_col, 'EdgeColor','none','FaceAlpha',0.45);
                catch err_v
                    fprintf('  violin failed %s plate %d: %s\n', grp_label, pidx, err_v.message);
                end
            end

            % scatter
            jx = pi + 0.06 * randn(n_col, 1);
            scatter(ax, jx, lag, 7, dot_col, 'filled', 'MarkerFaceAlpha', 0.40);

            % median line
            line(ax, [pi-0.25 pi+0.25], [median(lag) median(lag)], ...
                 'Color', violin_col*0.6, 'LineWidth', 2.0);

            % Gini label — fixed y, 30 deg, one per plate tick
            if n_col >= 3
                G              = gini_coeff(lag);
                plate_Gini(pi) = G;

                fprintf('%s | %s: Gini=%.3f  N=%d\n', grp_label, tick_labels{pi}, G, n_col);

                text(ax, pi, y_gini, sprintf('G=%.2f N=%d', G, n_col), ...
                     'Color',[0.05 0.20 0.70], ...
                     'HorizontalAlignment','left', ...
                     'VerticalAlignment','bottom', ...
                     'Rotation', 30, ...
                     'FontSize', 9, 'FontWeight', 'bold');
            end

            all_lag = [all_lag; lag];              %#ok<AGROW>
            all_pos = [all_pos; pi*ones(n_col,1)]; %#ok<AGROW>

        end % plate loop

        % Gini table
        T = table((1:n_plates)', tick_labels(:), plate_N, plate_Gini, ...
                  'VariableNames',{'PlatePos','Label','N_colonies','Gini_lag'});
        assignin('base', ['GiniTable_' fn], T);

        % Spearman corr(N, Gini) — per group
        ok = ~isnan(plate_N) & ~isnan(plate_Gini);
        if nnz(ok) >= 3
            [rho, pval] = corr(plate_N(ok), plate_Gini(ok), 'Type','Spearman');
            corr_group{end+1}   = grp_label;  %#ok<AGROW>
            corr_rho(end+1)     = rho;         %#ok<AGROW>
            corr_pval(end+1)    = pval;        %#ok<AGROW>
            corr_nplates(end+1) = nnz(ok);    %#ok<AGROW>
        else
            rho  = NaN;
            pval = NaN;
            corr_group{end+1}   = grp_label;  %#ok<AGROW>
            corr_rho(end+1)     = NaN;         %#ok<AGROW>
            corr_pval(end+1)    = NaN;         %#ok<AGROW>
            corr_nplates(end+1) = nnz(ok);    %#ok<AGROW>
        end

        % Accumulate for pooled correlation
        pool_N    = [pool_N;    plate_N(ok)];                      %#ok<AGROW>
        pool_Gini = [pool_Gini; plate_Gini(ok)];                   %#ok<AGROW>
        pool_grp  = [pool_grp;  repmat({grp_label}, nnz(ok), 1)]; %#ok<AGROW>

        % Stats
        stats.(fn).label  = grp_label;
        stats.(fn).lag    = all_lag;
        stats.(fn).plate  = all_pos;
        stats.(fn).N      = numel(all_lag);
        stats.(fn).median = median(all_lag);
        stats.(fn).iqr    = iqr(all_lag);

        fprintf('%s:  N=%d  |  median=%.2f h  |  IQR=%.2f h\n', ...
            grp_label, stats.(fn).N, stats.(fn).median, stats.(fn).iqr);

        % Axes — title includes Spearman rho and p-value
        ylabel(ax, 'Lag time (h)', 'FontSize', 14, 'FontWeight', 'bold');

        % Build significance star string
        if isnan(pval)
            sig_str = '';
        elseif pval < 0.001
            sig_str = '***';
        elseif pval < 0.01
            sig_str = '**';
        elseif pval < 0.05
            sig_str = '*';
        else
            sig_str = 'ns';
        end

        if isnan(rho)
            corr_line = 'Spearman: insufficient data';
        else
            corr_line = sprintf('Spearman r_s=%.3f, p=%s (%s)', rho, fmt_p(pval), sig_str);
        end

        title(ax, {grp_label, corr_line}, 'FontSize', 14, 'FontWeight', 'bold');
        xlim(ax, [0.5, n_plates+0.5]);
        ylim(ax, [y_lo, y_hi + 6.0]);
        set(ax, 'XTick', 1:n_plates, 'XTickLabel', tick_labels, ...
                'TickLabelInterpreter','none','FontSize', 11, ...
                'Box','off','GridAlpha',0.3);
        xtickangle(ax, 35);
        grid(ax, 'on');

        if g == n_groups || (n_cols_fig==2 && g >= n_groups-1)
            xlabel(ax, 'Biological replicate', 'FontSize', 14, 'FontWeight', 'bold');
        end

    end % group loop

    %% ----------------------------------------------------------------
    %  5.  Correlation summary — per group + pooled
    %% ----------------------------------------------------------------

    % --- Pooled Spearman (all groups combined) ---
    if length(pool_N) >= 3
        [rho_pool, pval_pool] = corr(pool_N, pool_Gini, 'Type','Spearman');
    else
        rho_pool  = NaN;
        pval_pool = NaN;
    end

    % --- Print formatted table ---
    fprintf('\n');
    fprintf('========================================================\n');
    fprintf('  Spearman correlation: Colony count vs Gini coefficient\n');
    fprintf('  (per time group, then pooled across all groups)\n');
    fprintf('========================================================\n');
    fprintf('  %-10s  %6s  %10s  %8s  %s\n', ...
        'Group', 'n plates', 'rho', 'p-value', 'Interpretation');
    fprintf('  %s\n', repmat('-', 1, 60));

    for gi = 1:length(corr_group)
        if isnan(corr_rho(gi))
            interp = 'insufficient data';
            fprintf('  %-10s  %6d  %10s  %8s  %s\n', ...
                corr_group{gi}, corr_nplates(gi), 'NaN', 'NaN', interp);
        else
            if corr_pval(gi) < 0.001,     sig = '***';
            elseif corr_pval(gi) < 0.01,  sig = '**';
            elseif corr_pval(gi) < 0.05,  sig = '*';
            else,                          sig = 'ns';
            end
            interp = sprintf('%s (p%s0.05)', sig, char('<' * (corr_pval(gi)<0.05) + '>' * (corr_pval(gi)>=0.05)));
            fprintf('  %-10s  %6d  %10.3f  %8s  %s\n', ...
                corr_group{gi}, corr_nplates(gi), corr_rho(gi), fmt_p(corr_pval(gi), '%.4f'), interp);
        end
    end

    fprintf('  %s\n', repmat('-', 1, 60));
    if ~isnan(rho_pool)
        if pval_pool < 0.001,     sig_p = '***';
        elseif pval_pool < 0.01,  sig_p = '**';
        elseif pval_pool < 0.05,  sig_p = '*';
        else,                     sig_p = 'ns';
        end
        interp_p = sprintf('%s (p%s0.05)', sig_p, char('<' * (pval_pool<0.05) + '>' * (pval_pool>=0.05)));
        fprintf('  %-10s  %6d  %10.3f  %8s  %s\n', ...
            'POOLED', length(pool_N), rho_pool, fmt_p(pval_pool, '%.4f'), interp_p);
    else
        fprintf('  %-10s  %6d  %10s  %8s  insufficient data\n', ...
            'POOLED', length(pool_N), 'NaN', 'NaN');
    end
    fprintf('========================================================\n\n');

    % --- Export correlation table to workspace ---
    % P_value is stored as formatted text (not a raw double) so a
    % Spearman p that underflowed to exact 0.0 in double precision
    % displays as "<1e-300" instead of a misleading literal 0 when the
    % table is viewed/printed.
    all_pval_disp = arrayfun(@(x) fmt_p(x, '%.4g'), [corr_pval(:); pval_pool], 'UniformOutput', false);
    CorrTable = table( ...
        [corr_group(:);  {'POOLED'}], ...
        [corr_nplates(:); length(pool_N)], ...
        [corr_rho(:);    rho_pool], ...
        all_pval_disp, ...
        'VariableNames', {'Group','N_plates','Spearman_rho','P_value'});
    assignin('base', 'GiniCorr_NvsGini', CorrTable);
    fprintf('  Correlation table saved to workspace as: GiniCorr_NvsGini\n\n');

    %% ----------------------------------------------------------------
    %  6.  Save outputs
    %% ----------------------------------------------------------------
    ts = datestr(now, 'yyyymmdd_HHMMSS');
    bn = sprintf('lagTime_colonies_Gini_%s', ts);

    save_fig(hfig, fullfile(out_dir, [bn '.fig']), 'fig');
    save_fig(hfig, fullfile(out_dir, [bn '.png']), 'png');
    save_fig(hfig, fullfile(out_dir, [bn '.svg']), 'svg');

    fprintf('\nDone — outputs saved to: %s\n', out_dir);

end


%% ====================================================================
%  LOCAL: build_plate_labels_interactive
%  Shows all metadata columns with value previews, asks user to pick
%  bio-rep column (required) and tech-rep column (optional), then builds
%  a {nr_plates x 1} cell of combined label strings.
%% ====================================================================
function plate_labels = build_plate_labels_interactive(meta_tbl, nr_plates)

    plate_labels = cell(nr_plates, 1);

    % Default fallback if no metadata available
    if isempty(meta_tbl) || ~istable(meta_tbl) || height(meta_tbl) == 0 || ...
       width(meta_tbl) == 0
        fprintf('  No metadata available — using plate index numbers.\n');
        for i = 1:nr_plates
            plate_labels{i} = sprintf('Plate%d', i);
        end
        return;
    end

    col_names = meta_tbl.Properties.VariableNames;
    n_cols    = length(col_names);

    fprintf('\n');
    fprintf('============================================================\n');
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

    % --- Biological replicate column (required) ---
    fprintf('  Step A — BIOLOGICAL replicate column\n');
    fprintf('  (groups different samples/strains/conditions on the x-axis,\n');
    fprintf('   e.g. "Set" -> Set1, Set2 ...)\n');
    fprintf('  Enter 0 to use plate index numbers only.\n\n');

    bio_col = '';
    while true
        raw = strtrim(input('  Biological replicate column number: ', 's'));
        idx = str2double(raw);
        if ~isnan(idx) && idx == 0
            fprintf('  -> Skipped. Using plate indices.\n\n');
            break;
        end
        if ~isnan(idx) && floor(idx)==idx && idx>=1 && idx<=n_cols
            bio_col = col_names{idx};
            fprintf('  -> Biological replicate: "%s"\n\n', bio_col);
            break;
        end
        fprintf('  WARNING: enter a number between 0 and %d.\n', n_cols);
    end

    % --- Technical replicate column (optional) ---
    fprintf('  Step B — TECHNICAL replicate column\n');
    fprintf('  (distinguishes replicate plates of the same sample,\n');
    fprintf('   e.g. "Replicate" -> R1, R2, R3 ...)\n');
    fprintf('  Enter 0 to skip.\n\n');

    tech_col = '';
    while true
        raw = strtrim(input('  Technical replicate column number (0 to skip): ', 's'));
        idx = str2double(raw);
        if ~isnan(idx) && idx == 0
            fprintf('  -> Skipped.\n\n');
            break;
        end
        if ~isnan(idx) && floor(idx)==idx && idx>=1 && idx<=n_cols
            tech_col = col_names{idx};
            fprintf('  -> Technical replicate: "%s"\n\n', tech_col);
            break;
        end
        fprintf('  WARNING: enter a number between 0 and %d.\n', n_cols);
    end

    % --- Dilution column (optional) ---
    fprintf('  Step C — DILUTION column\n');
    fprintf('  (e.g. "Dilution" -> -4(4), -5 ...)\n');
    fprintf('  Enter 0 to skip.\n\n');

    dil_col = '';
    while true
        raw = strtrim(input('  Dilution column number (0 to skip): ', 's'));
        idx = str2double(raw);
        if ~isnan(idx) && idx == 0
            fprintf('  -> Skipped.\n\n');
            break;
        end
        if ~isnan(idx) && floor(idx)==idx && idx>=1 && idx<=n_cols
            dil_col = col_names{idx};
            fprintf('  -> Dilution: "%s"\n\n', dil_col);
            break;
        end
        fprintf('  WARNING: enter a number between 0 and %d.\n', n_cols);
    end

    % --- Build label per plate: BioRep_TechRep_Dilution ---
    n_meta = height(meta_tbl);

    for i = 1:nr_plates
        if i > n_meta
            plate_labels{i} = sprintf('P%d', i);
            continue;
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

    % Confirm preview
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
%  LOCAL: val2str  — convert one table cell to a plain trimmed string
%% ====================================================================
function s = val2str(v)
    if iscell(v),          s = strtrim(char(v{1}));
    elseif isnumeric(v),   s = strtrim(num2str(v));
    elseif ischar(v),      s = strtrim(v);
    elseif isstring(v),    s = strtrim(char(v));
    else,                  s = '?';
    end
end


%% ====================================================================
%  LOCAL: build_per_plate_cache
%% ====================================================================
function ht = build_per_plate_cache(ht, raw_data, incTime, max_lag)

    nr_plates = ht.qc.nr_plates;
    min_col   = ht.params.min_col;
    max_col   = ht.params.max_col;
    cens_tol  = 0.25;
    ht.per_plate = struct();
    n_ok = 0;  n_bad = 0;

    fprintf('\n[build_per_plate_cache]  %d plates...\n', nr_plates);

    for i = 1:nr_plates
        fld = sprintf('plate_%d', i);

        if i > length(raw_data.processed)
            ht.per_plate.(fld) = [];  n_bad = n_bad+1;
            fprintf('  Plate %3d: out of range\n', i);  continue;
        end

        proc = raw_data.processed{i};

        if ~isfield(proc,'growth_quant') || ~proc.growth_quant
            ht.per_plate.(fld) = [];  n_bad = n_bad+1;
            fprintf('  Plate %3d: growth_quant=false\n', i);  continue;
        end

        if ~isfield(proc,'colonies') || ~isfield(proc.colonies,'new') || ...
           ~isfield(proc.colonies.new,'lag_time')
            ht.per_plate.(fld) = [];  n_bad = n_bad+1;
            fprintf('  Plate %3d: lag_time missing\n', i);  continue;
        end

        lag_raw = proc.colonies.new.lag_time(:) + incTime;
        lag_raw(lag_raw >= max_lag - cens_tol) = NaN;
        n_valid = sum(~isnan(lag_raw));

        if n_valid < min_col || n_valid > max_col
            ht.per_plate.(fld) = [];  n_bad = n_bad+1;
            fprintf('  Plate %3d: n=%d outside [%d,%d]\n', i, n_valid, min_col, max_col);
            continue;
        end

        lag = lag_raw(~isnan(lag_raw));
        ht.per_plate.(fld) = lag;
        n_ok = n_ok+1;
        fprintf('  Plate %3d: OK  n=%d  median=%.1f h\n', i, numel(lag), median(lag));
    end

    fprintf('[build_per_plate_cache]  Done: %d OK, %d skipped\n\n', n_ok, n_bad);
    assignin('base', 'ht', ht);
end


%% ====================================================================
%  LOCAL: fmt_p
%  Text display for a p-value. Spearman p (via corr()) can legitimately
%  underflow to exact 0.0 in double precision for large n / strong
%  correlation — the true value is just too small to represent, not
%  actually zero — shown as "<1e-300" instead.
%% ====================================================================
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


%% ====================================================================
%  LOCAL: gini_coeff
%% ====================================================================
function G = gini_coeff(x)
    x = x(isfinite(x) & x >= 0);
    n = numel(x);
    if n < 2 || sum(x) == 0, G = NaN; return; end
    x = sort(x(:));
    i = (1:n)';
    G = (2*sum(i.*x)) / (n*sum(x)) - (n+1)/n;
    G = max(0, min(1, G));
end


%% ====================================================================
%  LOCAL: save_fig
%% ====================================================================
function save_fig(hfig, fpath, fmt)
    try
        switch fmt
            case 'fig', savefig(hfig, fpath);
            case 'png', print(hfig, fpath, '-dpng', '-r300');
            case 'svg', print(hfig, fpath, '-dsvg');
        end
        fprintf('Saved .%s  ->  %s\n', fmt, fpath);
    catch ME
        warning('Could not save .%s: %s', fmt, ME.message);
    end
end
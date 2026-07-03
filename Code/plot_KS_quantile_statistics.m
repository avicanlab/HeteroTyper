%% HeteroTyper Pipeline for Bright Plates
% Plots pairwise Kolmogorov-Smirnov (D) statistics and 90th-percentile
% quantile differences for selected sample groups and features.
%
% All statistics are pre-computed by preprocess_pipeline_data and
% stored in ht.pairwise — this script only selects and plots.
%
% Output:
%   - Figure 1 : KS statistic D  (x-axis fixed 0–1, step 0.5)
%   - Figure 2 : Quantile difference Q90  (x-axis auto-scaled per feature)
%   - .txt log file mirroring all terminal output
%   - .xlsx file with all numeric values used in the figures
%
% Requires preprocess_pipeline_data(data) to have been run first.
%
%% USAGE:
%   plot_KS_quantile_statistics(ht)

function plot_KS_quantile_statistics(ht)

    % ------------------------------------------------------------------
    %  Validate input
    % ------------------------------------------------------------------
    if ~isstruct(ht) || ~isfield(ht,'groups')
        error(['Input must be the ht struct produced by preprocess_pipeline_data.\n' ...
               'Run:  preprocess_pipeline_data(data)  first.']);
    end
    if ~isfield(ht,'pairwise')
        error(['ht.pairwise not found.\n' ...
               'Please re-run preprocess_pipeline_data to generate pairwise statistics.']);
    end

    p        = ht.params;
    out_dir  = p.out_dir;
    pw       = ht.pairwise;   % shorthand

    % ------------------------------------------------------------------
    %  Open log file
    % ------------------------------------------------------------------
    timestamp_log = datestr(now,'yyyymmdd_HHMMSS');
    log_path = fullfile(out_dir, ...
        sprintf('plot_KS_quantile_statistics_%s.txt', timestamp_log));
    HT_FLOG = fopen(log_path,'w');
    if HT_FLOG == -1
        warning('Could not open log file: %s', log_path);
        HT_FLOG = -1;
    end
    ht_fprintf(HT_FLOG, 'Log file: %s\n', log_path);

    % ------------------------------------------------------------------
    %  User selections — groups and features shown together, chosen once
    % ------------------------------------------------------------------
    [sel_groups, sel_feats] = ask_selections(HT_FLOG, ht);

    n_groups   = length(sel_groups);
    n_feats    = length(sel_feats);
    feat_names = pw.feat_names(sel_feats);

    ht_fprintf(HT_FLOG, '\n--- KS & Quantile plot parameters ---\n');
    ht_fprintf(HT_FLOG, '  Groups   : %s\n', strjoin(ht.labels(sel_groups),', '));
    ht_fprintf(HT_FLOG, '  Features : %s\n', strjoin(feat_names,', '));
    ht_fprintf(HT_FLOG, '  Output   : %s\n', out_dir);
    ht_fprintf(HT_FLOG, '-------------------------------------\n\n');

    % ------------------------------------------------------------------
    %  Filter pairwise data to only the selected groups
    %  (pw contains ALL pairs; we keep only pairs where both groups
    %   are in sel_groups)
    % ------------------------------------------------------------------
    sel_labels = ht.labels(sel_groups);

    keep = false(size(pw.pair_idx,1), 1);
    for pi = 1:size(pw.pair_idx,1)
        a_in = ismember(pw.pair_idx(pi,1), sel_groups);
        b_in = ismember(pw.pair_idx(pi,2), sel_groups);
        keep(pi) = a_in && b_in;
    end

    pair_labels = pw.pair_labels(keep);
    KS_D        = pw.KS_D(keep, sel_feats);
    Q90         = pw.Q90(keep, sel_feats);
    n_pairs     = sum(keep);

    if n_pairs == 0
        error('No pairs found for the selected groups. Select at least 2 groups.');
    end

    % ------------------------------------------------------------------
    %  Print values to log
    % ------------------------------------------------------------------
    ht_fprintf(HT_FLOG, '%-28s  %-22s  %8s  %12s\n', ...
               'Pair','Feature','KS D','Q90 diff');
    ht_fprintf(HT_FLOG, '%s\n', repmat('-',1,76));
    for f = 1:n_feats
        for pi = 1:n_pairs
            if isnan(KS_D(pi,f)), continue; end
            ht_fprintf(HT_FLOG, '%-28s  %-22s  %8.4f  %12.4f\n', ...
                pair_labels{pi}, feat_names{f}, KS_D(pi,f), Q90(pi,f));
        end
    end
    ht_fprintf(HT_FLOG, '\n');

    % ------------------------------------------------------------------
    %  Figure layout
    % ------------------------------------------------------------------
    panel_w  = 2.4;   panel_h  = max(1.8, 0.32*n_pairs + 0.9);
    gap_x    = 0.50;   % wider gap prevents x-tick labels colliding between panels

    % margin_l scales with longest pair label to prevent y-tick text overflow
    max_label_len = max(cellfun(@length, pair_labels));
    margin_l = max(2.60, 0.11 * max_label_len);
    margin_r = 0.35;
    % margin_b: room for x-tick labels (FontSize 11 ~0.18") + x-title (FontSize 12 ~0.22") + gap
    margin_b = 0.75;
    % margin_t: room for panel title (FontSize 11 ~0.20") + gap above
    margin_t = 0.55;

    fig_w = margin_l + n_feats*panel_w + (n_feats-1)*gap_x + margin_r;
    fig_h = margin_b + panel_h + margin_t;

    % One uniform grey shade per feature column, evenly spaced light→dark.
    gray_levels  = linspace(0.80, 0.20, n_feats);
    gray_palette = repmat(gray_levels', 1, 3);   % [n_feats x 3] RGB

    % ------------------------------------------------------------------
    %  Figure 1 — KS statistic D  (x fixed 0..1, tick at 0.5)
    % ------------------------------------------------------------------
    hfig_ks = figure('Name','KS Statistics (D)', ...
                     'Color','w','Units','inches', ...
                     'Position',[1 5 fig_w fig_h], ...
                     'PaperUnits','inches','PaperSize',[fig_w fig_h], ...
                     'PaperPosition',[0 0 fig_w fig_h], ...
                     'Visible','off');

    for f = 1:n_feats
        left   = (margin_l + (f-1)*(panel_w+gap_x)) / fig_w;
        bottom = margin_b / fig_h;
        ax = axes('Position',[left, bottom, panel_w/fig_w, panel_h/fig_h], ...
                  'Parent',hfig_ks);
        hold(ax,'on');

        for pi = 1:n_pairs
            if isnan(KS_D(pi,f)), continue; end
            barh(ax, pi, KS_D(pi,f), 0.65, ...
                 'FaceColor',gray_palette(f,:), 'EdgeColor','none');
        end

        xlim(ax,[0 1]);
        xticks(ax,[0 0.5 1]);
        xline(ax, 0.5,'--','Color',[0.7 0.7 0.7],'LineWidth',0.8);
        ylim(ax,[0.4 n_pairs+0.6]);
        set(ax,'YTick',1:n_pairs,'YTickLabel',{},'YDir','reverse', ...
               'FontSize',11,'TickDir','in','Box','on');

        if f == 1
            set(ax,'YTickLabel',pair_labels);
            ax.YAxis.FontSize = 10;
        end

        title(ax, feat_names{f},'FontSize',11,'FontWeight','bold','Interpreter','none');
        ax.TitleFontSizeMultiplier = 1;  % prevent MATLAB from scaling title
    end

    % Fixed 0.30-inch band at the bottom for the x-axis label
    ann_h_ks = 0.30 / fig_h;
    annotation(hfig_ks,'textbox',[0, 0, 1, ann_h_ks], ...
               'String','Kolmogorov-Smirnov (D)', ...
               'HorizontalAlignment','center','VerticalAlignment','middle', ...
               'FontSize',12,'FontWeight','bold','EdgeColor','none');

    set(hfig_ks,'Visible','on'); drawnow;

    % ------------------------------------------------------------------
    %  Figure 2 — Quantile difference Q90  (x auto-scaled per feature)
    % ------------------------------------------------------------------
    hfig_q = figure('Name','Quantile Difference (Q90)', ...
                    'Color','w','Units','inches', ...
                    'Position',[1 0.5 fig_w fig_h], ...
                    'PaperUnits','inches','PaperSize',[fig_w fig_h], ...
                    'PaperPosition',[0 0 fig_w fig_h], ...
                    'Visible','off');

    for f = 1:n_feats
        left   = (margin_l + (f-1)*(panel_w+gap_x)) / fig_w;
        bottom = margin_b / fig_h;
        ax = axes('Position',[left, bottom, panel_w/fig_w, panel_h/fig_h], ...
                  'Parent',hfig_q); %#ok<LAXES>
        hold(ax,'on');

        vals = Q90(:,f);
        finite_vals = vals(isfinite(vals));

        for pi = 1:n_pairs
            if isnan(vals(pi)), continue; end
            barh(ax, pi, vals(pi), 0.65, ...
                 'FaceColor',gray_palette(f,:), 'EdgeColor','none');
        end

        % Symmetric x-axis scaled to data
        if ~isempty(finite_vals)
            abs_max = max(abs(finite_vals)) * 1.15;
            if abs_max == 0, abs_max = 1; end
            [~, tick_step] = nice_axis(abs_max);
            tick_step_int = max(1, round(tick_step));
            abs_max_int   = ceil(abs_max / tick_step_int) * tick_step_int;
            tick_vals     = -abs_max_int : tick_step_int : abs_max_int;
            xlim(ax, [-abs_max_int abs_max_int]);
            xticks(ax, tick_vals);
            xticklabels(ax, arrayfun(@(v) sprintf('%d', round(v)), tick_vals, ...
                            'UniformOutput', false));
        end
        xline(ax,0,'-','Color',[0.4 0.4 0.4],'LineWidth',0.8);

        ylim(ax,[0.4 n_pairs+0.6]);
        set(ax,'YTick',1:n_pairs,'YTickLabel',{},'YDir','reverse', ...
               'FontSize',11,'TickDir','in','Box','on');

        if f == 1
            set(ax,'YTickLabel',pair_labels);
            ax.YAxis.FontSize = 10;
        end

        title(ax, feat_names{f},'FontSize',11,'FontWeight','bold','Interpreter','none');
        ax.TitleFontSizeMultiplier = 1;  % prevent MATLAB from scaling title
    end

    % Fixed 0.30-inch band at the bottom for the x-axis label
    ann_h_q = 0.30 / fig_h;
    annotation(hfig_q,'textbox',[0, 0, 1, ann_h_q], ...
               'String','Quantile Difference (Q90)', ...
               'HorizontalAlignment','center','VerticalAlignment','middle', ...
               'FontSize',12,'FontWeight','bold','EdgeColor','none');

    set(hfig_q,'Visible','on'); drawnow;

    % ------------------------------------------------------------------
    %  Save figures
    % ------------------------------------------------------------------
    grp_str  = strjoin(sel_labels,'-');
    base_ks  = sprintf('KS_Statistics_%s_%s',   grp_str, timestamp_log);
    base_q   = sprintf('Quantile_Q90_%s_%s',     grp_str, timestamp_log);
    save_fig(HT_FLOG, hfig_ks, out_dir, base_ks);
    save_fig(HT_FLOG, hfig_q,  out_dir, base_q);

    % ------------------------------------------------------------------
    %  Export .xlsx  (two sheets: KS_D and Q90_diff)
    % ------------------------------------------------------------------
    xlsx_path = fullfile(out_dir, ...
        sprintf('KS_Quantile_Statistics_%s_%s.xlsx', grp_str, timestamp_log));
    export_xlsx(HT_FLOG, xlsx_path, pair_labels, feat_names, KS_D, Q90);

    % ------------------------------------------------------------------
    %  Close log
    % ------------------------------------------------------------------
    if HT_FLOG ~= -1
        fclose(HT_FLOG);
        fprintf('Log saved: %s\n', log_path);
    end

end  % end main function


%% ================================================================
%  LOCAL: ask_selections
%  Shows the full menu in one screen, collects group + feature
%  selections before returning — nothing appears one by one.
% ================================================================
function [sel_groups, sel_feats] = ask_selections(HT_FLOG, ht)

    all_labels = ht.labels;
    n_avail    = length(all_labels);
    feat_names = ht.pairwise.feat_names;
    n_feats    = length(feat_names);

    ht_fprintf(HT_FLOG, '\n============================================================\n');
    ht_fprintf(HT_FLOG, '  plot_KS_quantile_statistics — Selection\n');
    ht_fprintf(HT_FLOG, '============================================================\n');

    ht_fprintf(HT_FLOG, '\n  GROUPS\n');
    ht_fprintf(HT_FLOG, '  +-----+--------+-------------+\n');
    ht_fprintf(HT_FLOG, '  | idx | label  | n colonies  |\n');
    ht_fprintf(HT_FLOG, '  +-----+--------+-------------+\n');
    for k = 1:n_avail
        ht_fprintf(HT_FLOG, '  | %3d | %-6s | %11d |\n', ...
                   k, all_labels{k}, ht.groups(k).n_colonies);
    end
    ht_fprintf(HT_FLOG, '  +-----+--------+-------------+\n');

    ht_fprintf(HT_FLOG, '\n  FEATURES  (pre-computed in ht.pairwise)\n');
    ht_fprintf(HT_FLOG, '  +-----+---------------------+\n');
    for k = 1:n_feats
        ht_fprintf(HT_FLOG, '  | %3d | %-19s |\n', k, feat_names{k});
    end
    ht_fprintf(HT_FLOG, '  +-----+---------------------+\n');
    ht_fprintf(HT_FLOG, '\n============================================================\n');

    % --- Group selection ---
    while true
        fprintf('  Groups  — enter labels e.g. [ 3h 24h ]  or  ALL: ');
        raw_g = strtrim(input('','s'));
        if strcmpi(raw_g,'all')
            sel_groups = 1:n_avail;
            ht_fprintf(HT_FLOG, '  -> All %d groups selected.\n', n_avail);
            break;
        end
        tokens = regexp(raw_g,'[,\s]+','split');
        tokens = tokens(~cellfun(@isempty,tokens));
        sel_groups = zeros(1,length(tokens));
        ok = true;
        for t = 1:length(tokens)
            match = find(strcmpi(all_labels,tokens{t}));
            if isempty(match)
                fprintf('  WARNING  "%s" not recognised. Try again.\n',tokens{t});
                ok = false; break;
            end
            sel_groups(t) = match(1);
        end
        if ~ok, continue; end
        if length(sel_groups) < 2
            fprintf('  WARNING  Select at least 2 groups for pairwise comparison.\n');
            continue;
        end
        [~,ui]     = unique(sel_groups,'stable');
        sel_groups = sel_groups(ui);
        ht_fprintf(HT_FLOG, '  -> Groups: %s\n', strjoin(all_labels(sel_groups),', '));
        break;
    end

    % --- Feature selection ---
    while true
        fprintf('  Features — enter numbers e.g. [ 1 2 3 ]  or  ALL: ');
        raw_f = strtrim(input('','s'));
        if strcmpi(raw_f,'all')
            sel_feats = 1:n_feats;
            ht_fprintf(HT_FLOG, '  -> All %d features selected.\n', n_feats);
            return;
        end
        tokens = regexp(raw_f,'[,\s]+','split');
        nums   = str2double(tokens);
        if any(isnan(nums)) || any(nums<1) || any(nums>n_feats) || any(nums~=floor(nums))
            fprintf('  WARNING  Enter integers 1-%d or ALL.\n', n_feats);
            continue;
        end
        [~,ui]    = unique(nums,'stable');
        sel_feats = nums(ui);
        ht_fprintf(HT_FLOG, '  -> Features: %s\n', strjoin(feat_names(sel_feats),', '));
        return;
    end
end


%% ================================================================
%  LOCAL: nice_axis
%  Returns (abs_ceil, tick_step) for a symmetric axis half-width v.
% ================================================================
function [abs_ceil, tick_step] = nice_axis(v)
    if v <= 0, abs_ceil = 1; tick_step = 0.5; return; end
    raw_step   = v / 3;
    mag        = 10^floor(log10(raw_step));
    candidates = [1 2 2.5 5 10] * mag;
    for step = candidates
        abs_ceil = ceil(v / step) * step;
        if round(abs_ceil / step) <= 5
            tick_step = step;
            return;
        end
    end
    tick_step = candidates(end);
    abs_ceil  = ceil(v / tick_step) * tick_step;
end


%% ================================================================
%  LOCAL: export_xlsx
%  Writes KS_D and Q90 tables to two sheets in an .xlsx file.
% ================================================================
function export_xlsx(HT_FLOG, xlsx_path, pair_labels, feat_names, KS_D, Q90)
    n_pairs = length(pair_labels);
    n_feats = length(feat_names);
    header  = [{'Pair'}, feat_names];

    ks_table = cell(n_pairs+1, n_feats+1);
    q_table  = cell(n_pairs+1, n_feats+1);
    ks_table(1,:) = header;
    q_table(1,:)  = header;

    for r = 1:n_pairs
        ks_table{r+1,1} = pair_labels{r};
        q_table{r+1,1}  = pair_labels{r};
        for c = 1:n_feats
            ks_table{r+1,c+1} = KS_D(r,c);
            q_table{r+1,c+1}  = Q90(r,c);
        end
    end

    try
        writecell(ks_table, xlsx_path, 'Sheet','KS_D');
        writecell(q_table,  xlsx_path, 'Sheet','Q90_diff');
        ht_fprintf(HT_FLOG, 'Excel saved: %s\n', xlsx_path);
    catch ME
        ht_fprintf(HT_FLOG, 'WARNING: Could not save xlsx: %s\n', ME.message);
    end
end


%% ================================================================
%  LOCAL: save_fig
% ================================================================
function save_fig(HT_FLOG, hfig, out_dir, base_name)
    fig_path = fullfile(out_dir, base_name);
    savefig(hfig,  [fig_path '.fig']);
    exportgraphics(hfig, [fig_path '.pdf'], 'ContentType','vector');
    exportgraphics(hfig, [fig_path '.png'], 'Resolution',300);
    ht_fprintf(HT_FLOG, '  Saved: %s  (.fig / .pdf / .png)\n', base_name);
end


%% ================================================================
%  LOCAL: ht_fprintf
% ================================================================
function ht_fprintf(HT_FLOG, varargin)
    fprintf(varargin{:});
    if ~isempty(HT_FLOG) && HT_FLOG ~= -1
        fprintf(HT_FLOG, varargin{:});
    end
end
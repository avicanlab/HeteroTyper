%% HeteroTyper Pipeline for Bright Plates
% Plots doubling time histograms per selected sample group.
% Requires preprocess_pipeline_data(data) to have been run first.
%
%% USAGE:
%   plot_doublingTime(ht);

function plot_doublingTime(ht)
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
                        sprintf('plot_doublingTime_%s.txt', timestamp_log));
    HT_FLOG = fopen(log_path, 'w');
    if HT_FLOG == -1
        warning('Could not open log file: %s', log_path);
        HT_FLOG = -1;
    end
    fprintf('Log file: %s\n', log_path);


    % ------------------------------------------------------------------
    %  User prompt — which groups to plot
    % ------------------------------------------------------------------
    sel_idx  = ask_group_selection(HT_FLOG, ht);
    n_groups = length(sel_idx);

    ht_fprintf(HT_FLOG, '\n--- Doubling Time parameters ---\n');
    ht_fprintf(HT_FLOG, '  Groups selected : %s\n', strjoin(labels(sel_idx), ', '));
    ht_fprintf(HT_FLOG, '  Output directory: %s\n', p.out_dir);
    ht_fprintf(HT_FLOG, '--------------------------------\n\n');

    % ------------------------------------------------------------------
    %  Global x-axis range across selected groups (consistent axes)
    % ------------------------------------------------------------------
    all_vals = [];
    for k = 1:n_groups
        g     = sel_idx(k);
        valid = ht.groups(g).doublingTime(~isnan(ht.groups(g).doublingTime));
        all_vals = [all_vals; valid]; %#ok<AGROW>
    end

    if isempty(all_vals)
        ht_fprintf(HT_FLOG, 'No valid doubling time data found for selected groups.\n');
        return;
    end

    bin_width = 0.5;
    x_min     = floor(min(all_vals));
    % x_max is the true data maximum (not artificially capped) so the
    % x-axis always spans the full range of observed doubling times.
    x_max     = ceil(max(all_vals) / bin_width) * bin_width + bin_width;

    % ------------------------------------------------------------------
    %  Figure layout — manual position grid
    % ------------------------------------------------------------------
    panel_w  = 5.0;   panel_h  = 2.0;
    gap_y    = 0.60;
    margin_l = 0.90;  margin_r = 0.30;
    margin_b = 0.75;  margin_t = 0.25;

    fig_w = margin_l + panel_w + margin_r;
    fig_h = margin_b + n_groups*panel_h + (n_groups-1)*gap_y + margin_t;

    hfig = figure('Name','Doubling Time', ...
                  'Color','w', ...
                  'Units','inches', ...
                  'Position',[1 1 fig_w fig_h], ...
                  'PaperUnits','inches', ...
                  'PaperSize',[fig_w fig_h], ...
                  'PaperPosition',[0 0 fig_w fig_h]);

    combined_doublingT = [];

    % ------------------------------------------------------------------
    %  Plot loop
    % ------------------------------------------------------------------
    for k = 1:n_groups
        g      = sel_idx(k);
        grp    = ht.groups(g);
        validT = grp.doublingTime(~isnan(grp.doublingTime));

        % Axes position (row 1 = top)
        left   = margin_l / fig_w;
        bottom = (margin_b + (n_groups-k)*(panel_h+gap_y)) / fig_h;
        ax = axes('Position',[left, bottom, panel_w/fig_w, panel_h/fig_h]);

        if ~isempty(validT)
            histogram(ax, validT, ...
                      'BinWidth',   bin_width, ...
                      'BinLimits',  [x_min x_max], ...
                      'FaceColor',  colors(g,:), ...
                      'EdgeColor',  'none', ...
                      'FaceAlpha',  1);
            hold(ax, 'on');

            medVal = median(validT);
            yl     = ylim(ax);
            line(ax, [medVal medVal], [0 yl(2)], 'Color','r', 'LineWidth',1.5);

            % Smart corner: split x-range at its midpoint and compare
            % how much data mass falls on each side. Label goes to the
            % side with LESS mass (i.e. more empty space).
            x_mid      = (x_min + x_max) / 2;
            mass_left  = sum(validT <  x_mid);
            mass_right = sum(validT >= x_mid);

            if mass_left <= mass_right
                % emptier on left -> label top-left
                x_txt   = x_min + 0.02*(x_max-x_min);
                h_align = 'left';
            else
                % emptier on right -> label top-right
                x_txt   = x_max - 0.02*(x_max-x_min);
                h_align = 'right';
            end

            text(ax, x_txt, yl(2)*0.95, ...
                 sprintf('%s\nn = %d  |  median = %.2f h', ...
                         labels{g}, length(validT), medVal), ...
                 'FontSize',10, 'FontWeight','bold', 'Color','k', ...
                 'HorizontalAlignment',h_align, 'VerticalAlignment','top');
        else
            text(ax, 0.5, 0.5, sprintf('%s — no valid data', labels{g}), ...
                 'Units','normalized', 'HorizontalAlignment','center', ...
                 'FontSize',11, 'Color',[0.5 0.5 0.5]);
            axis(ax, 'off');
        end

        xlim(ax, [x_min x_max]);
        ylabel(ax, 'Frequency', 'FontSize',11, 'FontWeight','bold');

        % x-label only on bottom panel
        if k == n_groups
            xlabel(ax, 'Doubling Time (h)', 'FontSize',11, 'FontWeight','bold');
        end

        set(ax, 'FontSize',10, 'TickDir','out', 'Box','off');

        combined_doublingT = [combined_doublingT; validT]; %#ok<AGROW>
    end

    % ------------------------------------------------------------------
    %  Export to workspace
    % ------------------------------------------------------------------
    assignin('base', 'combined_doublingT', combined_doublingT);
    ht_fprintf(HT_FLOG, 'combined_doublingT written to workspace  (n = %d values).\n', ...
            length(combined_doublingT));

    % ------------------------------------------------------------------
    %  Save (.fig, .pdf, .png)
    % ------------------------------------------------------------------
    out_dir   = p.out_dir;
    timestamp = datestr(now, 'yyyymmdd_HHMMSS');
    grp_str   = strjoin(labels(sel_idx), '-');
    base_name = sprintf('DoublingTime_%s_%s', grp_str, timestamp);
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
    ht_fprintf(HT_FLOG, '  +-----+--------+------------------+\n');
    ht_fprintf(HT_FLOG, '  | idx | label  | n doubling times |\n');
    ht_fprintf(HT_FLOG, '  +-----+--------+------------------+\n');
    for k = 1:n_avail
        dt    = ht.groups(k).doublingTime;
        n_val = sum(~isnan(dt));
        ht_fprintf(HT_FLOG, '  | %3d | %-6s | %16d |\n', k, all_labels{k}, n_val);
    end
    ht_fprintf(HT_FLOG, '  +-----+--------+------------------+\n');
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
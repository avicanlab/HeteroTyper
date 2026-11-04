%% HeteroTyper Pipeline - Manual vs automatic colony counts over time
%
%%  PURPOSE
%  -------
%  Compares the automatic (pipeline) colony counts with manual counts
%  per plate over the imaging time course.
%
%   1) Asks for a colony size threshold (px). A colony is counted from
%      the first timepoint its smoothed size timecourse reaches >= the
%      threshold (same first-crossing rule as extract_cell_count_per_time
%      STEP 1). Counts are cumulative, i.e. the TOTAL number of colonies
%      on the plate at each timepoint.
%   2) Builds the automatic total count per timepoint for every plate from
%      data.processed{i}.colonies.new.timecourse_size_smoothed
%      (flag_colony_ok mask applied when available, as in the QC caching
%      of preprocess_pipeline_data). Time = time_info.elapsed_time_h + incTime.
%   3) Asks for the Excel file + sheet with the manual counts per time
%      course (long format: one row per plate x timepoint, e.g. columns
%      Position / Time / Total Counts) and which columns to use.
%   4) Combines automatic and manual counts by plate name (Pos....) and
%      time into one table (long + wide format), saved as .xlsx.
%   5) Line-dot plot: ONE main figure with one subplot per plate,
%      x = time (h, shared range), y = total colony count,
%      colour = Automatic / Manual (filled dots, no marker edge).
%
%  Plate names are taken from data.metadata.original.Position (e.g. Pos0073).
%  The output folder is always asked for (folder dialog, or typed path).
%
%%  USAGE
%  -----
%    plot_timecourse_manual_auto_counts(data);       % asks incTime
%    plot_timecourse_manual_auto_counts(data, ht);   % reuses ht.params.incTime
%                                               % (ht.params.out_dir = dialog start folder)
%
%  Outputs (in the chosen output folder):
%    plot_timecourse_manual_auto_counts_<thr>px_<timestamp>.xlsx         sheets Combined_long, Combined_wide
%    plot_timecourse_manual_auto_counts_<thr>px_<timestamp>.png/.pdf/.fig all plates as subplots
%    plot_timecourse_manual_auto_counts_log_<timestamp>.txt

function plot_timecourse_manual_auto_counts(data, ht)

    if nargin < 2
        ht = [];
    end
    if ~isstruct(data) || ~isfield(data, 'processed') || ~isfield(data, 'metadata')
        error('Input must be the pipeline data struct (with .processed and .metadata).');
    end

    AUTO_COLOR   = [0.165 0.471 0.839];   % #2a78d6  blue
    MANUAL_COLOR = [0.922 0.408 0.204];   % #eb6834  orange

    fprintf('\n============================================================\n');
    fprintf('  plot_timecourse_manual_auto_counts\n');
    fprintf('  Manual vs automatic total colony counts per plate over time\n');
    fprintf('============================================================\n');

    % ------------------------------------------------------------------
    %  1) Size threshold + time / output parameters
    % ------------------------------------------------------------------
    default_thr = 100;
    if isstruct(ht) && isfield(ht, 'params') && isfield(ht.params, 'size_threshold')
        default_thr = ht.params.size_threshold;
    end
    threshold = ask_positive_number( ...
        sprintf('\n  1) Count colonies of at least how many px? [default = %g]: ', default_thr), ...
        default_thr, false);
    fprintf('  -> Size threshold: >= %g px\n', threshold);

    % Incubation time is added to the AUTOMATIC time axis only (elapsed
    % imaging time from the first image). Manual times are used exactly
    % as given in the Excel sheet.
    if isstruct(ht) && isfield(ht, 'params') && isfield(ht.params, 'incTime')
        incTime = ht.params.incTime;
        fprintf('  -> RT incubation time from ht.params: %.4g h (automatic counts only)\n', incTime);
    else
        fprintf('\n  RT incubation time is added to the AUTOMATIC counts'' time only\n');
        fprintf('  (elapsed imaging time). Manual times are used as given in Excel.\n');
        incTime = ask_positive_number( ...
            '  RT incubation time (h) [default = 0]: ', 0, true);
        fprintf('  -> RT incubation time: %.4g h (automatic counts only)\n', incTime);
    end

    default_out = pwd;
    if isstruct(ht) && isfield(ht, 'params') && isfield(ht.params, 'out_dir') && isfolder(ht.params.out_dir)
        default_out = ht.params.out_dir;
    end
    out_dir = ask_out_dir(default_out);

    % ------------------------------------------------------------------
    %  Log file — mirrors everything printed to the terminal
    % ------------------------------------------------------------------
    timestamp = datestr(now, 'yyyymmdd_HHMMSS');
    log_path  = fullfile(out_dir, sprintf('plot_timecourse_manual_auto_counts_log_%s.txt', timestamp));
    flog = fopen(log_path, 'w');
    if flog == -1
        warning('Could not open log file: %s', log_path);
    end
    function tlog(varargin)
        fprintf(varargin{:});
        if flog ~= -1
            fprintf(flog, varargin{:});
        end
    end

    tlog('\n--- Parameters ---\n');
    tlog('  Size threshold       : >= %g px (first crossing of smoothed size)\n', threshold);
    tlog('  RT incubation time   : %.4g h (added to automatic times only; manual times as given)\n', incTime);
    tlog('  Output directory     : %s\n', out_dir);

    % ------------------------------------------------------------------
    %  2) Automatic total counts per plate per timepoint
    % ------------------------------------------------------------------
    plate_names = get_plate_names(data);
    auto = compute_auto_counts(data, plate_names, threshold, incTime);
    tlog('  Automatic counts     : %d plates with a size timecourse\n', numel(unique(auto.Plate)));

    % ------------------------------------------------------------------
    %  3) Manual counts from Excel
    % ------------------------------------------------------------------
    [manual, xls_path, sheet] = ask_manual_counts();
    tlog('  Manual counts file   : %s\n', xls_path);
    tlog('  Manual counts sheet  : %s\n', sheet);
    tlog('  Manual counts        : %d plates with at least one count\n', numel(unique(manual.Plate)));

    % ------------------------------------------------------------------
    %  Plates to plot
    % ------------------------------------------------------------------
    man_plates  = unique(manual.Plate, 'stable');
    auto_plates = unique(auto.Plate, 'stable');
    fprintf('\n  Which plates to combine / plot?\n');
    fprintf('    [1] Only plates with manual counts (default)\n');
    fprintf('    [2] All plates (manual and/or automatic)\n');
    raw = strtrim(input('  Choice [1/2]: ', 's'));
    if strcmp(raw, '2')
        plates = union(man_plates, auto_plates, 'stable');
    else
        plates = man_plates;
    end
    plates = sort_plates(plates);

    no_auto = setdiff(man_plates, auto_plates);
    if ~isempty(no_auto)
        tlog('  WARNING  No automatic counts for manual plate(s): %s\n', strjoin(no_auto, ', '));
    end
    tlog('  Plates combined      : %d (%s)\n', numel(plates), strjoin(plates, ', '));
    tlog('-------------------------\n\n');

    % ------------------------------------------------------------------
    %  4) Combine by plate + time
    % ------------------------------------------------------------------
    auto   = auto(ismember(auto.Plate, plates), :);
    manual = manual(ismember(manual.Plate, plates), :);

    auto.Source   = repmat({'Automatic'}, height(auto), 1);
    manual.Source = repmat({'Manual'},    height(manual), 1);
    long_tbl = [auto(:, {'Plate','Time_h','Count','Source'}); ...
                manual(:, {'Plate','Time_h','Count','Source'})];
    long_tbl = sortrows(long_tbl, {'Plate','Source','Time_h'});

    wide_tbl = combine_wide(auto, manual, plates);

    xlsx_path = fullfile(out_dir, sprintf('plot_timecourse_manual_auto_counts_%gpx_%s.xlsx', threshold, timestamp));
    writetable(long_tbl, xlsx_path, 'Sheet', 'Combined_long');
    writetable(wide_tbl, xlsx_path, 'Sheet', 'Combined_wide');
    tlog('  Combined table saved : %s\n\n', xlsx_path);

    % Per-plate summary: counts at the last manual timepoint
    tlog('  %-10s  %10s  %10s  %10s  %10s\n', 'Plate', 'Time (h)', 'Auto', 'Manual', 'Auto-Man');
    tlog('  %s\n', repmat('-', 1, 58));
    for k = 1:numel(plates)
        rows = wide_tbl(strcmp(wide_tbl.Plate, plates{k}) & ~isnan(wide_tbl.Manual_count), :);
        if isempty(rows)
            continue;
        end
        [~, j] = max(rows.Time_h);
        tlog('  %-10s  %10.2f  %10g  %10g  %10g\n', plates{k}, rows.Time_h(j), ...
            rows.Auto_count(j), rows.Manual_count(j), rows.Auto_count(j) - rows.Manual_count(j));
    end
    tlog('\n');

    % ------------------------------------------------------------------
    %  5) Line-dot plot: one main figure, one subplot (tile) per plate
    % ------------------------------------------------------------------
    n_pl  = numel(plates);
    n_col = ceil(sqrt(n_pl));
    n_row = ceil(n_pl / n_col);
    fig_w = 3.2 * n_col;
    fig_h = 2.6 * n_row + 0.9;

    % Shared x-range so plates are directly comparable; y is per plate
    all_t = [auto.Time_h; manual.Time_h];
    x_lim = [floor(min(all_t)), ceil(max(all_t))];
    if x_lim(1) == x_lim(2), x_lim = x_lim + [-1 1]; end

    hfig = figure('Color', 'w', 'Name', 'Manual vs automatic counts', 'Units', 'inches', ...
                  'Position', [1 1 fig_w fig_h], 'PaperUnits', 'inches', ...
                  'PaperSize', [fig_w fig_h], 'PaperPosition', [0 0 fig_w fig_h]);
    tl = tiledlayout(hfig, n_row, n_col, 'TileSpacing', 'compact', 'Padding', 'compact');

    h_auto = gobjects(0); h_man = gobjects(0);
    for k = 1:n_pl
        pl = plates{k};
        a  = sortrows(auto(strcmp(auto.Plate, pl), :), 'Time_h');
        m  = sortrows(manual(strcmp(manual.Plate, pl), :), 'Time_h');

        ax = nexttile(tl); hold(ax, 'on');
        if ~isempty(a)
            h_auto = plot(ax, a.Time_h, a.Count, '-o', 'Color', AUTO_COLOR, 'LineWidth', 1.2, ...
                'MarkerSize', 3.5, 'MarkerFaceColor', AUTO_COLOR, 'MarkerEdgeColor', 'none');
        end
        if ~isempty(m)
            h_man = plot(ax, m.Time_h, m.Count, '-o', 'Color', MANUAL_COLOR, 'LineWidth', 1.2, ...
                'MarkerSize', 3.5, 'MarkerFaceColor', MANUAL_COLOR, 'MarkerEdgeColor', 'none');
        end
        hold(ax, 'off');
        title(ax, pl, 'FontSize', 10, 'FontWeight', 'normal');
        set(ax, 'Box', 'off', 'TickDir', 'out', 'FontSize', 9);
        grid(ax, 'on'); ax.GridAlpha = 0.15;
        xlim(ax, x_lim);
        ylim(ax, [0 Inf]);
    end

    xlabel(tl, 'Time (h)');
    ylabel(tl, 'Total colony count');
    title(tl, sprintf('Manual vs automatic colony counts — colonies \\geq %g px', threshold), ...
          'Interpreter', 'tex', 'FontWeight', 'bold');

    h = [h_auto, h_man];
    lg = [repmat({'Automatic'}, 1, numel(h_auto)), repmat({'Manual'}, 1, numel(h_man))];
    if ~isempty(h)
        lgd = legend(h, lg, 'Orientation', 'horizontal', 'Box', 'off');
        lgd.Layout.Tile = 'north';
    end

    base = fullfile(out_dir, sprintf('plot_timecourse_manual_auto_counts_%gpx_%s', threshold, timestamp));
    savefig(hfig, [base '.fig']);
    print(hfig, [base '.png'], '-dpng', '-r300');
    print(hfig, [base '.pdf'], '-dpdf');
    tlog('  Figure saved (%d plates, %dx%d subplots): %s.png\n', n_pl, n_row, n_col, base);

    tlog('\nDone.\n');
    if flog ~= -1
        fclose(flog);
        fprintf('Log saved: %s\n', log_path);
    end
end


%% ================================================================
%  LOCAL FUNCTION: get_plate_names
%  Plate name per data.processed index from metadata Position column
%  (numeric positions are formatted as PosNNNN).
% ================================================================
function names = get_plate_names(data)
    n = numel(data.processed);
    names = cell(n, 1);
    meta = table();
    if isfield(data.metadata, 'original')
        meta = data.metadata.original;
    end
    has_pos = ~isempty(meta) && ismember('Position', meta.Properties.VariableNames);
    for i = 1:n
        if has_pos && i <= height(meta)
            v = meta.Position(i);
            if iscell(v), v = v{1}; end
            if isnumeric(v)
                names{i} = sprintf('Pos%04d', v);
            else
                names{i} = normalize_plate(char(string(v)));
            end
        else
            names{i} = sprintf('Pos%04d', i);
        end
    end
end


%% ================================================================
%  LOCAL FUNCTION: compute_auto_counts
%  Cumulative number of colonies whose smoothed size has reached
%  >= threshold at each timepoint, per plate. Long-format table.
% ================================================================
function T = compute_auto_counts(data, plate_names, threshold, incTime)
    plate_c = {}; time_c = []; count_c = [];
    for i = 1:numel(data.processed)
        pr = data.processed{i};
        if ~isstruct(pr) || ~isfield(pr, 'colonies') || ~isfield(pr.colonies, 'new')
            continue;
        end
        c_new = pr.colonies.new;
        if ~isfield(c_new, 'timecourse_size_smoothed') || isempty(c_new.timecourse_size_smoothed) || ...
           ~isfield(c_new, 'time_info') || ~isfield(c_new.time_info, 'elapsed_time_h')
            continue;
        end
        tc = c_new.timecourse_size_smoothed;          % [T x N]
        t  = c_new.time_info.elapsed_time_h(:) + incTime;

        % Valid-colony mask, same as the QC caching in preprocess_pipeline_data
        if isfield(pr.colonies, 'region_props') && isfield(pr.colonies.region_props, 'flag_colony_ok')
            flag = logical(pr.colonies.region_props.flag_colony_ok);
            if numel(flag) == size(tc, 2)
                tc = tc(:, flag);
            end
        end

        nT = min(size(tc, 1), numel(t));
        tc = tc(1:nT, :); t = t(1:nT);

        % First crossing per colony -> cumulative total per timepoint
        reached = cummax(double(tc >= threshold), 1);
        counts  = sum(reached, 2);

        plate_c = [plate_c; repmat(plate_names(i), nT, 1)]; %#ok<AGROW>
        time_c  = [time_c;  t]; %#ok<AGROW>
        count_c = [count_c; counts]; %#ok<AGROW>
    end
    if isempty(plate_c)
        error('No plate in data.processed has colonies.new.timecourse_size_smoothed + time_info.');
    end
    T = table(plate_c, time_c, count_c, 'VariableNames', {'Plate','Time_h','Count'});
end


%% ================================================================
%  LOCAL FUNCTION: ask_manual_counts
%  Excel file -> sheet -> plate / time / count columns.
% ================================================================
function [T, xls_path, sheet] = ask_manual_counts()
    fprintf('\n  3) Select the Excel file with the manual counts per time course.\n');
    try
        [f, d] = uigetfile({'*.xlsx;*.xls;*.xlsm', 'Excel files'}, 'Select manual counts Excel file');
    catch
        f = 0;   % no display available -> type the path instead
    end
    if isequal(f, 0)
        fprintf('  (dialog cancelled - type the path instead)\n');
        while true
            xls_path = strtrim(input('  Path to Excel file: ', 's'));
            xls_path = strrep(xls_path, '"', '');
            if isfile(xls_path), break; end
            fprintf('  WARNING  File not found. Try again.\n');
        end
    else
        xls_path = fullfile(d, f);
    end
    fprintf('  -> File: %s\n', xls_path);

    sheets = cellstr(sheetnames(xls_path));
    fprintf('\n  Sheets:\n');
    for k = 1:numel(sheets)
        fprintf('    [%d] %s\n', k, sheets{k});
    end
    k = ask_index('  Sheet with the manual counts (number): ', numel(sheets), []);
    sheet = sheets{k};
    fprintf('  -> Sheet: %s\n', sheet);

    opts = detectImportOptions(xls_path, 'Sheet', sheet, 'VariableNamingRule', 'preserve');
    raw  = readtable(xls_path, opts);
    cols = raw.Properties.VariableNames;

    fprintf('\n  Columns:\n');
    for k = 1:numel(cols)
        fprintf('    [%d] %s\n', k, cols{k});
    end
    d_plate = guess_col(cols, {'^position$', '^pos', 'plate'});
    d_time  = guess_col(cols, {'^time', 'time'});
    d_count = guess_col(cols, {'^total\s*count', 'total', 'count'});
    c_plate = ask_index(col_prompt('Plate name (Pos....) column', d_plate, cols), numel(cols), d_plate);
    c_time  = ask_index(col_prompt('Time (h) column', d_time, cols), numel(cols), d_time);
    c_count = ask_index(col_prompt('Total colony count column', d_count, cols), numel(cols), d_count);
    fprintf('  -> Plate: %s | Time: %s | Count: %s\n', cols{c_plate}, cols{c_time}, cols{c_count});

    plate = raw.(cols{c_plate});
    if isnumeric(plate)
        plate = arrayfun(@(v) sprintf('Pos%04d', v), plate, 'UniformOutput', false);
    else
        plate = cellfun(@(v) normalize_plate(char(string(v))), cellstr(string(plate)), 'UniformOutput', false);
    end
    time_h = to_numeric(raw.(cols{c_time}));
    count  = to_numeric(raw.(cols{c_count}));

    keep = ~cellfun(@isempty, plate) & ~isnan(time_h) & ~isnan(count);
    if ~any(keep)
        error('No rows with a plate name, time and count found in sheet "%s".', sheet);
    end
    T = table(plate(keep), time_h(keep), count(keep), 'VariableNames', {'Plate','Time_h','Count'});
end


%% ================================================================
%  LOCAL FUNCTION: combine_wide
%  One row per plate x timepoint with Auto_count and Manual_count.
%  Manual times are matched to the nearest automatic timepoint of the
%  same plate within half an imaging interval; unmatched times keep
%  their own row with NaN for the missing source.
% ================================================================
function W = combine_wide(auto, manual, plates)
    P = {}; Tm = []; A = []; M = [];
    for k = 1:numel(plates)
        a = sortrows(auto(strcmp(auto.Plate, plates{k}), :), 'Time_h');
        m = manual(strcmp(manual.Plate, plates{k}), :);
        ta = a.Time_h; ca = a.Count;
        cm = nan(size(ta));
        tol = 0.25;
        if numel(ta) > 1
            tol = 0.5 * median(diff(ta));
        end
        extra_t = []; extra_m = [];
        for r = 1:height(m)
            if isempty(ta)
                d = Inf;
            else
                [d, j] = min(abs(ta - m.Time_h(r)));
            end
            if d <= tol + 1e-9
                cm(j) = m.Count(r);
            else
                extra_t(end+1, 1) = m.Time_h(r); %#ok<AGROW>
                extra_m(end+1, 1) = m.Count(r); %#ok<AGROW>
            end
        end
        t_all = [ta; extra_t];
        P  = [P;  repmat(plates(k), numel(t_all), 1)]; %#ok<AGROW>
        Tm = [Tm; t_all]; %#ok<AGROW>
        A  = [A;  ca; nan(numel(extra_t), 1)]; %#ok<AGROW>
        M  = [M;  cm; extra_m]; %#ok<AGROW>
    end
    W = table(P, Tm, A, M, 'VariableNames', {'Plate','Time_h','Auto_count','Manual_count'});
    W.Auto_minus_Manual = W.Auto_count - W.Manual_count;
    W = sortrows(W, {'Plate','Time_h'});
end


%% ================================================================
%  Small helpers
% ================================================================
function s = normalize_plate(s)
    % 'Pos0073_48h_...' / 'pos73' -> 'Pos0073'; anything else kept as is
    s = strtrim(s);
    tok = regexp(s, '^[Pp]os(\d+)', 'tokens', 'once');
    if ~isempty(tok)
        s = sprintf('Pos%04d', str2double(tok{1}));
    end
end

function plates = sort_plates(plates)
    num = cellfun(@(s) str2double(regexp(s, '\d+', 'match', 'once')), plates);
    num(isnan(num)) = Inf;
    [~, o] = sort(num);
    plates = plates(o);
end

function v = to_numeric(x)
    if isnumeric(x)
        v = double(x);
    else
        v = str2double(string(x));
    end
    v = v(:);
end

function idx = guess_col(cols, patterns)
    idx = [];
    for p = 1:numel(patterns)
        hit = find(~cellfun(@isempty, regexpi(cols, patterns{p}, 'once')), 1);
        if ~isempty(hit)
            idx = hit; return;
        end
    end
end

function s = col_prompt(what, def, cols)
    if isempty(def)
        s = sprintf('  %s (number): ', what);
    else
        s = sprintf('  %s (number) [default = %d, %s]: ', what, def, cols{def});
    end
end

function idx = ask_index(prompt, n, def)
    while true
        raw = strtrim(input(prompt, 's'));
        if isempty(raw) && ~isempty(def)
            idx = def; return;
        end
        v = str2double(raw);
        if ~isnan(v) && v >= 1 && v <= n && floor(v) == v
            idx = v; return;
        end
        fprintf('  WARNING  Enter a whole number between 1 and %d.\n', n);
    end
end

function val = ask_positive_number(prompt, def, allow_zero)
    while true
        raw = strtrim(input(prompt, 's'));
        if isempty(raw)
            val = def; return;
        end
        val = str2double(raw);
        if ~isnan(val) && (val > 0 || (allow_zero && val == 0))
            return;
        end
        fprintf('  WARNING  Enter a positive number, or press Enter for the default.\n');
    end
end

function out_dir = ask_out_dir(def)
    fprintf('\n  OUTPUT FOLDER for the combined table, figure and log.\n');
    fprintf('  Select it in the dialog (cancel to type the path instead).\n');
    try
        sel = uigetdir(def, 'Select output folder');
    catch
        sel = 0;   % no display available -> type the path instead
    end
    if ~isequal(sel, 0)
        out_dir = sel;
        fprintf('  -> Output folder: %s\n', out_dir);
        return;
    end

    while true
        raw = strtrim(input(sprintf('  Output folder [Enter = %s]: ', def), 's'));
        raw = strrep(raw, '"', '');
        if isempty(raw)
            out_dir = def; break;
        end
        if isfolder(raw)
            out_dir = raw; break;
        end
        create_it = strtrim(input(sprintf('  Folder not found: %s\n  Create it now? (y/n): ', raw), 's'));
        if strcmpi(create_it, 'y')
            mkdir(raw);
            out_dir = raw; break;
        end
    end
    fprintf('  -> Output folder: %s\n', out_dir);
end

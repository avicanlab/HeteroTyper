%% HeteroTyper Pipeline - QC Growth Data Visualisation
%
%%  PURPOSE
%  -------
%  Four per-plate QC figures for every plate that has growth_quant == true.
%  All per-plate data (time vectors, size timecourses, lag times, early
%  doubling times, valid-colony flags) are read directly from ht.qc —
%  no access to the raw data struct is needed.
%
%%  USAGE
%  -----
%    plot_QC_growth_data(ht)

function plot_QC_growth_data(ht)

    % =========================================================
    %  Validate input
    % =========================================================
    if ~isstruct(ht) || ~isfield(ht, 'qc') || ~isfield(ht.qc, 'ix_growth')
        error(['plot_QC_growth_data: input must be the ht struct from preprocess_pipeline_data.\n' ...
               'Run:  preprocess_pipeline_data(data)  first.\n' ...
               'ht.qc.ix_growth not found — re-run the updated preprocess_pipeline_data.']);
    end

    % =========================================================
    %  Unpack from ht — no recomputation
    % =========================================================
    p          = ht.params;
    qc         = ht.qc;
    out_dir    = p.out_dir;

    ix_growth  = sort(qc.ix_growth(:));
    n_plates   = length(ix_growth);

    % Timestamp for filenames and log
    timestamp = datestr(now, 'yyyymmdd_HHMMSS');

    % =========================================================
    %  Open log file — mirrors everything printed to terminal
    % =========================================================
    log_path = fullfile(out_dir, sprintf('plot_QC_growth_data_%s.txt', timestamp));
    flog = fopen(log_path, 'w');
    if flog == -1
        warning('plot_QC_growth_data: Could not open log file: %s', log_path);
        flog = -1;
    end

    function tlog(varargin)
        fprintf(varargin{:});
        if flog ~= -1
            fprintf(flog, varargin{:});
        end
    end

    tlog('========================================\n');
    tlog('plot_QC_growth_data\n');
    tlog('Started: %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
    tlog('Output directory: %s\n', out_dir);
    tlog('========================================\n\n');
    tlog('Plates with growth data: %d\n\n', n_plates);

    if n_plates == 0
        tlog('WARNING: No plates with growth data found — nothing to plot.\n');
        if flog ~= -1, fclose(flog); end
        return;
    end

    % =========================================================
    %  Bin edges — from ht (no recomputation)
    % =========================================================
    xval  = 0 : 100 : ht.xbins.size_upper;            % colony size
    xval2 = p.incTime : 0.5 : p.max_lag;              % lag time
    xval3 = 0 : 0.5 : p.max_lag;                      % early doubling time

    % =========================================================
    %  Auto grid layout — most square grid fitting all plates
    % =========================================================
    n_cols_grid = ceil(sqrt(n_plates));
    n_rows_grid = ceil(n_plates / n_cols_grid);
    tlog('Grid: %d rows x %d cols (%d plates)\n\n', n_rows_grid, n_cols_grid, n_plates);

    % =========================================================
    %  Build plate -> group colour lookup
    % =========================================================
    plate_group = zeros(qc.nr_plates, 1);
    for g = 1:length(ht.groups)
        for pi = ht.groups(g).plate_indices(:)'
            if pi >= 1 && pi <= qc.nr_plates
                plate_group(pi) = g;
            end
        end
    end

    % =========================================================
    %  Screen size for figure positioning
    % =========================================================
    screen_px = get(0, 'ScreenSize');
    screen_w  = screen_px(3);
    screen_h  = screen_px(4);
    fig_scale = 0.85;

    % =========================================================
    %  PLOT 1 — Growth curves (size timecourse per plate)
    % =========================================================
    tlog('--- Plot 1: Growth curves ---\n');

    fig1 = figure('Name', sprintf('QC growth curves (%d plates)', n_plates), ...
                  'Units', 'pixels', ...
                  'Position', [screen_w*0.04, screen_h*0.04, ...
                               screen_w*fig_scale, screen_h*fig_scale], ...
                  'Color', 'w');

    for k = 1:n_plates
        i      = ix_growth(k);
        t_vec  = qc.plate_time{i};
        tc_mat = qc.plate_size_tc{i};   % [T x N_valid] — already flag-masked

        g_idx    = plate_group(i);
        fill_col = group_color(g_idx, ht.colors);

        ax = subplot(n_rows_grid, n_cols_grid, k);
        hold(ax, 'on');

        n_valid = 0;
        if ~isempty(t_vec) && ~isempty(tc_mat)
            n_valid = size(tc_mat, 2);
            plot(ax, t_vec, tc_mat, '-', 'Color', [fill_col, 0.45], 'LineWidth', 0.5);
            xlim(ax, [0, p.max_lag]);
            ylim(ax, [0, inf]);
        end

        % Colony count label — top-left, normalised coordinates
        text(ax, 0.02, 0.98, num2str(n_valid), ...
             'Units', 'normalized', ...
             'FontSize', 6, 'FontWeight', 'bold', ...
             'Color', fill_col * 0.6, ...
             'HorizontalAlignment', 'left', ...
             'VerticalAlignment', 'top');

        plate_title(ax, qc, i);
        set(ax, 'FontSize', 5);
        xlabel(ax, 'Time (h)', 'FontSize', 4);
        ylabel(ax, 'Size (px)', 'FontSize', 4);

        tlog('  Plate %2d: n_valid=%d\n', i, n_valid);
    end
    blank_unused(n_plates, n_rows_grid, n_cols_grid);

    save_figure(fig1, out_dir, sprintf('QC_growth_curves_%s', timestamp));
    tlog('  Saved: QC_growth_curves_%s\n\n', timestamp);

    % =========================================================
    %  PLOT 2 — Colony size distributions
    % =========================================================
    tlog('--- Plot 2: Colony size distributions ---\n');

    % Compute global x-max across all plates for shared axis
    global_size_max = 0;
    for k = 1:n_plates
        tc_k = qc.plate_size_tc{ix_growth(k)};
        if ~isempty(tc_k)
            global_size_max = max(global_size_max, max(tc_k(end, :), [], 'omitnan'));
        end
    end
    global_size_max = max(global_size_max, xval(end));
    tlog('  global_size_max = %.0f px\n', global_size_max);

    fig2 = figure('Name', sprintf('QC size distributions (%d plates)', n_plates), ...
                  'Units', 'pixels', ...
                  'Position', [screen_w*0.04, screen_h*0.04, ...
                               screen_w*fig_scale, screen_h*fig_scale], ...
                  'Color', 'w');

    for k = 1:n_plates
        i        = ix_growth(k);
        tc_k     = qc.plate_size_tc{i};
        g_idx    = plate_group(i);
        fill_col = group_color(g_idx, ht.colors);
        line_col = fill_col * 0.6;

        ax = subplot(n_rows_grid, n_cols_grid, k);
        hold(ax, 'on');

        n_valid = 0;
        if ~isempty(tc_k)
            final_sizes = tc_k(end, :)';
            final_sizes = final_sizes(isfinite(final_sizes) & final_sizes > 0);
            n_valid = length(final_sizes);

            if n_valid > 0
                [N_counts, edges] = histcounts(final_sizes, xval);
                plot_histogram(ax, edges, N_counts, fill_col, line_col);

                med_val = median(final_sizes);
                y_top   = max([N_counts, 1]);
                line(ax, [med_val med_val], [0, y_top], ...
                     'Color', 'r', 'LineWidth', 1.5);
                xlim(ax, [0, global_size_max]);
                ylim(ax, [0, max(1, ceil(y_top * 1.1))]);

                tlog('  Plate %2d: n=%d  size=[%.0f, %.0f]  median=%.0f\n', ...
                     i, n_valid, min(final_sizes), max(final_sizes), med_val);
            end
        end

        text(ax, 0.02, 0.98, num2str(n_valid), ...
             'Units', 'normalized', 'FontSize', 6, 'FontWeight', 'bold', ...
             'Color', line_col, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');

        plate_title(ax, qc, i);
        set(ax, 'FontSize', 5);
        xlabel(ax, 'Size (px)', 'FontSize', 4);
        ylabel(ax, 'Count', 'FontSize', 4);
    end
    blank_unused(n_plates, n_rows_grid, n_cols_grid);

    save_figure(fig2, out_dir, sprintf('QC_size_dist_%s', timestamp));
    tlog('  Saved: QC_size_dist_%s\n\n', timestamp);

    % =========================================================
    %  PLOT 3 — Lag time distributions
    % =========================================================
    tlog('--- Plot 3: Lag time distributions ---\n');

    fig3 = figure('Name', sprintf('QC lag time distributions (%d plates)', n_plates), ...
                  'Units', 'pixels', ...
                  'Position', [screen_w*0.04, screen_h*0.04, ...
                               screen_w*fig_scale, screen_h*fig_scale], ...
                  'Color', 'w');

    for k = 1:n_plates
        i        = ix_growth(k);
        lag_k    = qc.plate_lag{i};
        g_idx    = plate_group(i);
        fill_col = group_color(g_idx, ht.colors);
        line_col = fill_col * 0.6;

        ax = subplot(n_rows_grid, n_cols_grid, k);
        hold(ax, 'on');

        n_valid = 0;
        if ~isempty(lag_k)
            lag_finite = lag_k(isfinite(lag_k));
            n_valid    = length(lag_finite);

            if n_valid > 0
                [N_counts, edges] = histcounts(lag_finite, xval2);
                plot_histogram(ax, edges, N_counts, fill_col, line_col);

                med_val = median(lag_finite);
                y_top   = max([N_counts, 1]);
                line(ax, [med_val med_val], [0, y_top], ...
                     'Color', 'r', 'LineWidth', 1.5);
                xlim(ax, [p.incTime, p.max_lag]);
                ylim(ax, [0, max(1, ceil(y_top * 1.1))]);

                tlog('  Plate %2d: n=%d  lag=[%.1f, %.1f]  median=%.1f\n', ...
                     i, n_valid, min(lag_finite), max(lag_finite), med_val);
            end
        end

        text(ax, 0.02, 0.98, num2str(n_valid), ...
             'Units', 'normalized', 'FontSize', 6, 'FontWeight', 'bold', ...
             'Color', line_col, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');

        plate_title(ax, qc, i);
        set(ax, 'FontSize', 5);
        xlabel(ax, 'Lag time (h)', 'FontSize', 4);
        ylabel(ax, 'Count', 'FontSize', 4);
    end
    blank_unused(n_plates, n_rows_grid, n_cols_grid);

    save_figure(fig3, out_dir, sprintf('QC_lag_time_dist_%s', timestamp));
    tlog('  Saved: QC_lag_time_dist_%s\n\n', timestamp);

    % =========================================================
    %  PLOT 4 — Early doubling time distributions
    % =========================================================
    tlog('--- Plot 4: Early doubling time distributions ---\n');

    fig4 = figure('Name', sprintf('QC early doubling time distributions (%d plates)', n_plates), ...
                  'Units', 'pixels', ...
                  'Position', [screen_w*0.04, screen_h*0.04, ...
                               screen_w*fig_scale, screen_h*fig_scale], ...
                  'Color', 'w');

    for k = 1:n_plates
        i        = ix_growth(k);
        edt_k    = qc.plate_early_dt{i};
        g_idx    = plate_group(i);
        fill_col = group_color(g_idx, ht.colors);
        line_col = fill_col * 0.6;

        ax = subplot(n_rows_grid, n_cols_grid, k);
        hold(ax, 'on');

        n_valid = 0;
        if ~isempty(edt_k)
            edt_finite = edt_k(isfinite(edt_k) & edt_k > 0);
            n_valid    = length(edt_finite);

            if n_valid > 0
                [N_counts, edges] = histcounts(edt_finite, xval3);
                plot_histogram(ax, edges, N_counts, fill_col, line_col);

                med_val = median(edt_finite);
                y_top   = max([N_counts, 1]);
                line(ax, [med_val med_val], [0, y_top], ...
                     'Color', 'r', 'LineWidth', 1.5);
                xlim(ax, [0, p.max_lag]);
                ylim(ax, [0, max(1, ceil(y_top * 1.1))]);

                tlog('  Plate %2d: n=%d  early_DT=[%.1f, %.1f]  median=%.1f\n', ...
                     i, n_valid, min(edt_finite), max(edt_finite), med_val);
            end
        end

        text(ax, 0.02, 0.98, num2str(n_valid), ...
             'Units', 'normalized', 'FontSize', 6, 'FontWeight', 'bold', ...
             'Color', line_col, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');

        plate_title(ax, qc, i);
        set(ax, 'FontSize', 5);
        xlabel(ax, 'Early DT (h)', 'FontSize', 4);
        ylabel(ax, 'Count', 'FontSize', 4);
    end
    blank_unused(n_plates, n_rows_grid, n_cols_grid);

    save_figure(fig4, out_dir, sprintf('QC_early_DT_dist_%s', timestamp));
    tlog('  Saved: QC_early_DT_dist_%s\n\n', timestamp);

    % =========================================================
    %  Close log
    % =========================================================
    tlog('========================================\n');
    tlog('Completed: %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
    tlog('========================================\n');

    if flog ~= -1
        fclose(flog);
        fprintf('Log saved: %s\n', log_path);
    end

end  % end main function


%% ================================================================
%  LOCAL FUNCTION: group_color
%  Returns the group RGB colour, or a grey-blue fallback.
% ================================================================
function c = group_color(g_idx, colors)
    if g_idx >= 1 && g_idx <= size(colors, 1)
        c = colors(g_idx, :);
    else
        c = [0.4 0.6 0.8];
    end
end


%% ================================================================
%  LOCAL FUNCTION: plate_title
%  Sets the subplot title from ht.qc.fn or plate index.
% ================================================================
function plate_title(ax, qc, plate_idx)
    if ~isempty(qc.fn) && plate_idx <= length(qc.fn)
        t = strrep(qc.fn(plate_idx).name, '_', '-');
    else
        t = sprintf('Plate %d', plate_idx);
    end
    title(ax, t, 'FontWeight', 'Normal', 'Interpreter', 'none', 'FontSize', 5);
end


%% ================================================================
%  LOCAL FUNCTION: plot_histogram
%  Draws a filled-patch histogram + step outline, matching the
%  style used in plot_QC_full_dataset.
% ================================================================
function plot_histogram(ax, edges, counts, fill_col, line_col)
    n_bins = length(counts);
    % Filled patches (no edge)
    for b = 1:n_bins
        if counts(b) > 0
            patch(ax, ...
                  [edges(b), edges(b+1), edges(b+1), edges(b)], ...
                  [0, 0, counts(b), counts(b)], ...
                  fill_col, 'EdgeColor', 'none');
        end
    end
    % Step outline across all bins
    sx = zeros(1, 2*n_bins + 2);
    sy = zeros(1, 2*n_bins + 2);
    sx(1) = edges(1);  sy(1) = 0;
    for b = 1:n_bins
        sx(2*b)   = edges(b);     sy(2*b)   = counts(b);
        sx(2*b+1) = edges(b+1);   sy(2*b+1) = counts(b);
    end
    sx(end) = edges(end);  sy(end) = 0;
    plot(ax, sx, sy, '-', 'Color', line_col, 'LineWidth', 0.8);
end


%% ================================================================
%  LOCAL FUNCTION: blank_unused
%  Turns off empty subplot slots in a partial last row.
% ================================================================
function blank_unused(n_plates, nr, nc)
    for k = n_plates+1 : nr*nc
        ax_empty = subplot(nr, nc, k);
        axis(ax_empty, 'off');
    end
end


%% ================================================================
%  LOCAL FUNCTION: save_figure
%  Saves a figure as .fig, .pdf (vector), and .png (300 dpi).
% ================================================================
function save_figure(hfig, out_dir, base_name)
    fig_path = fullfile(out_dir, base_name);
    try
        savefig(hfig, [fig_path '.fig']);
    catch
        warning('save_figure: could not save .fig: %s', fig_path);
    end
    try
        exportgraphics(hfig, [fig_path '.pdf'], 'ContentType', 'vector');
    catch
        warning('save_figure: could not save .pdf: %s', fig_path);
    end
    try
        exportgraphics(hfig, [fig_path '.png'], 'Resolution', 300);
    catch
        warning('save_figure: could not save .png: %s', fig_path);
    end
end
%% HeteroTyper Pipeline - QC Visualisation for Full Dataset
%
%%  PURPOSE
%  -------
%  Quality-control plots for the full plate dataset.
%  All filtering, counts, area data, and manual counts are read from
%  ht (produced by preprocess_pipeline_data) — no recomputation needed.
%
%  data is only needed for Plot 3 (image montage).
%
%%  USAGE
%  -----
%    plot_QC_full_dataset(ht)            % Plots 1, 1A, 2 only
%    plot_QC_full_dataset(ht, data)      % + image montage (Plot 3)

function plot_QC_full_dataset(ht, data)

    % =========================================================
    %  Validate input
    % =========================================================
    if ~isstruct(ht) || ~isfield(ht, 'qc')
        error(['plot_QC_full_dataset: input must be the ht struct from preprocess_pipeline_data.\n' ...
               'Run:  preprocess_pipeline_data(data)  first.']);
    end

    have_data = nargin >= 2 && isstruct(data);

    % =========================================================
    %  Unpack from ht — no recomputation
    % =========================================================
    p   = ht.params;
    qc  = ht.qc;

    nr_plates          = qc.nr_plates;
    colony_count_all   = qc.colony_count_all;
    colony_count_clean = qc.colony_count_clean;
    colony_count_final = qc.colony_count_final;
    passes             = qc.passes;
    rel_plates         = find(passes);
    n_rel              = length(rel_plates);
    min_col            = p.min_col;
    max_col            = p.max_col;
    out_dir            = p.out_dir;

    % Timestamp for filenames and log
    timestamp = datestr(now, 'yyyymmdd_HHMMSS');

    % =========================================================
    %  Open log file — mirrors everything printed to terminal
    % =========================================================
    log_path = fullfile(out_dir, sprintf('plot_QC_full_dataset_%s.txt', timestamp));
    flog = fopen(log_path, 'w');
    if flog == -1
        warning('plot_QC_full_dataset: Could not open log file: %s', log_path);
        flog = -1;
    end

    function tlog(varargin)
        fprintf(varargin{:});
        if flog ~= -1
            fprintf(flog, varargin{:});
        end
    end

    tlog('========================================\n');
    tlog('plot_QC_full_dataset\n');
    tlog('Started: %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
    tlog('Output directory: %s\n', out_dir);
    tlog('========================================\n\n');

    tlog('QC filter: min_col=%d  max_col=%d  size_threshold=%d px\n', ...
         min_col, max_col, p.size_threshold);
    tlog('Plates passing QC: %d / %d\n\n', n_rel, nr_plates);

    % Log per-plate summary
    tlog('--- Per-plate colony counts ---\n');
    tlog('  %4s  %8s  %8s  %8s  %6s\n', 'Plate', 'All', 'Clean', 'Final', 'Pass');
    for i = 1:nr_plates
        tlog('  %4d  %8d  %8d  %8d  %6s\n', i, ...
             colony_count_all(i), colony_count_clean(i), colony_count_final(i), ...
             mat2str(passes(i)));
    end
    tlog('\n');

    % =========================================================
    %  Get screen size for figure sizing
    % =========================================================
    screen_px = get(0, 'ScreenSize');   % [left bottom width height] in px
    screen_w  = screen_px(3);
    screen_h  = screen_px(4);
    fig_scale = 0.85;                   % slightly smaller than screen

    % =========================================================
    %  PLOT 1 — Colony counts: all detected vs clean vs final
    % =========================================================
    tlog('--- Plot 1: Colony counts ---\n');

    fig1 = figure('Name', 'Colony count (all plates)', ...
                  'Units', 'pixels', ...
                  'Position', [screen_w*0.05, screen_h*0.05, ...
                               screen_w*fig_scale*0.7, screen_h*fig_scale*0.5]);

    plot_all   = colony_count_all;   plot_all(plot_all == 0)     = NaN;
    plot_clean = colony_count_clean; plot_clean(plot_clean == 0) = NaN;
    plot_final = colony_count_final; plot_final(plot_final == 0) = NaN;

    semilogy(plot_all,   'ok', 'DisplayName', 'All detected');  hold on;
    semilogy(plot_clean, 'or', 'DisplayName', 'Clean (flag\_colony\_ok)');
    semilogy(plot_final, 'sb', 'MarkerFaceColor', [0.2 0.5 0.9], ...
             'DisplayName', 'Final (size-filtered)');

    for i = 1:nr_plates
        if ~isnan(plot_all(i)) && ~isnan(plot_final(i))
            line([i i], [plot_final(i), plot_all(i)], 'Color', [0.7 0.7 0.7]);
        end
    end

    line([0.5, nr_plates+0.5], [min_col, min_col], 'Color', 'k', 'LineStyle', '--');
    line([0.5, nr_plates+0.5], [max_col, max_col], 'Color', 'k', 'LineStyle', '--');

    legend('Location', 'best');
    set(gca, 'FontSize', 10, 'XTick', 0:5:nr_plates);
    xlim([0, nr_plates+0.5]);
    xlabel('Plate index');
    ylabel('Colony count (log scale)');
    title(sprintf('Colony counts  —  thresholds: %d / %d', min_col, max_col));

    save_figure(fig1, out_dir, sprintf('QC_colony_counts_%s', timestamp));
    tlog('  Saved: QC_colony_counts_%s\n\n', timestamp);

    % =========================================================
    %  PLOT 1A — Automated (final, size-filtered) vs manual counts
    % =========================================================
    manual_counts = qc.manual_counts;
    count_col     = qc.manual_count_col;

    if ~isempty(count_col) && any(manual_counts > 0 & ~isnan(manual_counts))
        ix_original = find(manual_counts > 0 & ~isnan(manual_counts));
        n_ix        = length(ix_original);

        tlog('--- Plot 1A: Automated vs manual counts ---\n');
        tlog('  Column: %s\n', count_col);
        tlog('  Plates with manual counts: %d\n', n_ix);
        tlog('  %4s  %10s  %10s\n', 'Plate', 'Automated', 'Manual');
        for k = 1:n_ix
            tlog('  %4d  %10d  %10.0f\n', ix_original(k), ...
                 colony_count_final(ix_original(k)), manual_counts(ix_original(k)));
        end
        tlog('\n');

        fig1a = figure('Name', 'Automated vs manual colony count', ...
                       'Units', 'pixels', ...
                       'Position', [screen_w*0.05, screen_h*0.05, ...
                                    screen_w*fig_scale*0.75, screen_h*fig_scale*0.55]);

        % --- Left panel ---
        subplot(1, 2, 1);

        auto_plot   = colony_count_final(ix_original);
        manual_plot = manual_counts(ix_original);
        auto_plot(auto_plot == 0)     = NaN;
        manual_plot(manual_plot == 0) = NaN;

        semilogy(1:n_ix, auto_plot,   'or', 'DisplayName', 'Automated (size-filtered)');
        hold on;
        semilogy(1:n_ix, manual_plot, 'ob', 'DisplayName', 'Manual');
        legend('Location', 'best');

        % x-axis: start from 0, labels at multiples of 5, plates start at 1
        xtick_positions = 0:5:n_ix;
        if isempty(xtick_positions) || xtick_positions(end) < n_ix
            xtick_positions(end+1) = n_ix;
        end
        xtick_labels = cell(size(xtick_positions));
        for ti = 1:length(xtick_positions)
            if xtick_positions(ti) == 0
                xtick_labels{ti} = '0';
            else
                xtick_labels{ti} = num2str(ix_original(xtick_positions(ti)));
            end
        end
        set(gca, 'FontSize', 10, 'XTick', xtick_positions, 'XTickLabel', xtick_labels);
        xlim([0, n_ix + 0.5]);
        xlabel('Plate position');
        ylabel('Colony count');
        title('Count comparison per plate');
        xtickangle(-45);

        % --- Right panel: scatter (linear scale, sample-label colours) ---
        subplot(1, 2, 2);
        auto_vals   = colony_count_final(ix_original);
        manual_vals = manual_counts(ix_original);
        valid       = (auto_vals > 0) & (manual_vals > 0);

        % Build plate->group lookup (same logic as Plot 2)
        plate_group_1a = zeros(nr_plates, 1);
        for g = 1:length(ht.groups)
            for pi = ht.groups(g).plate_indices(:)'
                if pi >= 1 && pi <= nr_plates
                    plate_group_1a(pi) = g;
                end
            end
        end

        if any(valid)
            mv = manual_vals(valid);
            av = auto_vals(valid);
            px = ix_original(valid);   % plate indices for colour lookup

            % --- Statistics (compute before plotting regression line) ---
            % Pearson
            [r_pearson, p_pearson] = corr(mv, av, 'Type', 'Pearson');
            % Spearman
            [r_spearman, p_spearman] = corr(mv, av, 'Type', 'Spearman');
            % R^2 from Pearson r
            R2 = r_pearson^2;
            % Linear regression slope & intercept (least-squares)
            coeff     = polyfit(mv, av, 1);
            slope     = coeff(1);
            intercept = coeff(2);

            % Regression line spanning exact data range only (min to max)
            x_fit   = linspace(min(mv), max(mv), 200);
            y_fit   = polyval(coeff, x_fit);

            % Draw confidence band (95%) via residual std
            y_pred  = polyval(coeff, mv);
            resid   = av - y_pred;
            s_e     = std(resid);
            n_pts   = length(mv);
            x_bar   = mean(mv);
            ss_x    = sum((mv - x_bar).^2);
            se_band = s_e * sqrt(1/n_pts + (x_fit - x_bar).^2 / ss_x);
            t_crit  = tinv(0.975, n_pts - 2);
            y_upper = y_fit + t_crit * se_band;
            y_lower = y_fit - t_crit * se_band;

            hold on;
            % Shaded confidence band
            fill([x_fit, fliplr(x_fit)], [y_upper, fliplr(y_lower)], ...
                 [0.7 0.7 0.7], 'EdgeColor', 'none', 'FaceAlpha', 0.4);
            % Regression line
            plot(x_fit, y_fit, '-', 'Color', [0.4 0.4 0.4], 'LineWidth', 1.5);

            % Scatter dots on top
            for di = 1:length(px)
                g_idx = plate_group_1a(px(di));
                if g_idx >= 1 && g_idx <= size(ht.colors, 1)
                    dot_col = ht.colors(g_idx, :);
                else
                    dot_col = [0.3 0.3 0.3];
                end
                plot(mv(di), av(di), 'o', ...
                     'MarkerFaceColor', dot_col, 'MarkerEdgeColor', dot_col * 0.6, ...
                     'MarkerSize', 6);
            end

            % Axes start from 0, upper limit with 5% padding above data max
            x_lim_max = max(mv) * 1.05;
            y_lim_max = max(av) * 1.05;
            xlim([0, x_lim_max]);
            ylim([0, y_lim_max]);

            % Build equation string
            if intercept >= 0
                eq_str = sprintf('y = %.3fx + %.1f', slope, intercept);
            else
                eq_str = sprintf('y = %.3fx - %.1f', slope, abs(intercept));
            end

            % --- Print full stats to terminal / log ---
            tlog('\n  --- Automated vs Manual count statistics (n=%d) ---\n', sum(valid));
            tlog('  Pearson  r = %.4f,  p = %s\n', r_pearson, fmt_p(p_pearson));
            tlog('  Spearman r = %.4f,  p = %s\n', r_spearman, fmt_p(p_spearman));
            tlog('  R^2 (Pearson) = %.4f\n', R2);
            tlog('  Slope equation: %s\n\n', eq_str);

            % --- Annotation on plot: Spearman stats + equation ---
            ann_str = sprintf('\\rho = %.3f\np_{Sp} = %s\nR^2 = %.3f\n%s', ...
                              r_spearman, fmt_p(p_spearman), R2, eq_str);
            text(x_lim_max * 0.03, ...
                 y_lim_max * 0.97, ...
                 ann_str, 'FontSize', 8, 'VerticalAlignment', 'top', ...
                 'BackgroundColor', 'w', 'Margin', 3);
        end
        set(gca, 'FontSize', 10);
        xlabel(sprintf('Manual count  [%s]', strrep(count_col, '_', '\_')));
        ylabel('Automated count (size-filtered)');
        title('Agreement');

        save_figure(fig1a, out_dir, sprintf('QC_count_comparison_%s', timestamp));
        tlog('  Saved: QC_count_comparison_%s\n\n', timestamp);

    else
        tlog('  NOTE: No manual counts in ht.qc — skipping Plot 1A.\n\n');
    end

    % =========================================================
    %  PLOT 2 — Colony-size histograms (auto-sized grid)
    %  - Fixed bin width = 100 px, shared across all subplots
    %  - x-axis limit = global max area across all passing plates
    %  - Fill colour = sample-group colour (from ht.colors)
    %  - No bar edge; top of each bin drawn as a step line (same
    %    colour darkened by 40%) so individual bins are readable
    % =========================================================
    tlog('--- Plot 2: Colony size distributions ---\n');

    if n_rel == 0
        tlog('  WARNING: No plates pass QC filter — skipping.\n\n');
        warning('plot_QC_full_dataset: No plates pass QC filter — skipping size histograms.');
    else
        % ----------------------------------------------------------
        % Build plate -> group index lookup
        % ----------------------------------------------------------
        plate_group = zeros(nr_plates, 1);   % 0 = not assigned
        for g = 1:length(ht.groups)
            for pi = ht.groups(g).plate_indices(:)'
                if pi >= 1 && pi <= nr_plates
                    plate_group(pi) = g;
                end
            end
        end

        % ----------------------------------------------------------
        % Shared bin edges: bin_size = 100, global x-max across all
        % passing plates (rounded up to nearest bin boundary)
        % ----------------------------------------------------------
        bin_size   = 100;
        global_max = 0;
        for k = 1:n_rel
            col_area_k = qc.plate_area{rel_plates(k)};
            if ~isempty(col_area_k)
                global_max = max(global_max, max(col_area_k));
            end
        end
        global_max = ceil(global_max / bin_size) * bin_size;
        if global_max == 0, global_max = bin_size; end
        edges = 0:bin_size:global_max;

        % Nice x-axis ticks starting from 0
        x_tick_step = round_to_nice(global_max / 5);
        xtick_vals  = 0:x_tick_step:global_max;

        tlog('  bin_size=%d px  global_x_max=%d px  x_tick_step=%d px\n', ...
             bin_size, global_max, x_tick_step);

        % ----------------------------------------------------------
        % Grid layout
        % ----------------------------------------------------------
        n_cols_grid = ceil(sqrt(n_rel));
        n_rows_grid = ceil(n_rel / n_cols_grid);
        tlog('  Grid: %d rows x %d cols  (%d plates)\n', n_rows_grid, n_cols_grid, n_rel);

        fig2 = figure('Name', sprintf('Colony size distributions  (%d plates)', n_rel), ...
                      'Units', 'pixels', ...
                      'Position', [screen_w*0.04, screen_h*0.04, ...
                                   screen_w*fig_scale, screen_h*fig_scale]);

        for k = 1:n_rel
            i        = rel_plates(k);
            col_area = qc.plate_area{i};

            % Resolve group colour for this plate
            g_idx = plate_group(i);
            if g_idx >= 1 && g_idx <= size(ht.colors, 1)
                fill_col = ht.colors(g_idx, :);
            else
                fill_col = [0.4 0.6 0.8];   % fallback grey-blue
            end
            % Darken the fill colour by 40% for the step-line on top
            line_col = fill_col * 0.6;

            subplot(n_rows_grid, n_cols_grid, k);
            hold on;

            if ~isempty(col_area) && any(col_area > 0)
                n_colonies = length(col_area);
                counts_h   = histcounts(col_area, edges);
                y_max      = max(counts_h);
                n_bins     = length(counts_h);

                % --- Filled patch per bin (no EdgeColor) ---
                for b = 1:n_bins
                    if counts_h(b) > 0
                        x_left  = edges(b);
                        x_right = edges(b+1);
                        h_bar   = counts_h(b);
                        patch([x_left, x_right, x_right, x_left], ...
                              [0, 0, h_bar, h_bar], ...
                              fill_col, 'EdgeColor', 'none');
                    end
                end

                % --- Step line across the tops of all bins ---
                % Build staircase: left edge of bin 1 -> right edge of last bin
                stair_x = zeros(1, 2*n_bins + 2);
                stair_y = zeros(1, 2*n_bins + 2);
                stair_x(1) = edges(1);
                stair_y(1) = 0;
                for b = 1:n_bins
                    stair_x(2*b)   = edges(b);
                    stair_x(2*b+1) = edges(b+1);
                    stair_y(2*b)   = counts_h(b);
                    stair_y(2*b+1) = counts_h(b);
                end
                stair_x(end) = edges(end);
                stair_y(end) = 0;
                plot(stair_x, stair_y, '-', 'Color', line_col, 'LineWidth', 0.8);

                set(gca, 'XTick', xtick_vals);
                axis([0, global_max, 0, max(1, ceil(1.1 * y_max))]);

                % n = ... label upper-left
                text(global_max * 0.04, max(1, 0.90 * max(1, y_max)), ...
                     sprintf('n = %d', n_colonies), ...
                     'FontSize', 5, 'VerticalAlignment', 'top', 'Color', line_col);

                tlog('  Plate %2d (group %d): n=%d  area=[%.0f, %.0f]\n', ...
                     i, g_idx, n_colonies, min(col_area), max(col_area));
            else
                n_colonies = colony_count_final(i);
                text(0.5, 0.5, sprintf('n = %d\n(no area data)', n_colonies), ...
                     'Units', 'normalized', 'HorizontalAlignment', 'center', 'FontSize', 5);
                set(gca, 'XTick', xtick_vals, 'XLim', [0, global_max]);
                tlog('  Plate %2d: n=%d  (no area data)\n', i, n_colonies);
            end

            % Title from ht.qc.fn or plate index
            if ~isempty(qc.fn) && i <= length(qc.fn)
                t = strrep(qc.fn(i).name, '_', '-');
            else
                t = sprintf('Plate %d', i);
            end
            title(t, 'FontWeight', 'Normal', 'Interpreter', 'none');
            set(gca, 'FontSize', 5);
            xlabel('Area (px)', 'FontSize', 4);
            ylabel('Count',     'FontSize', 4);
        end

        save_figure(fig2, out_dir, sprintf('QC_size_distributions_%s', timestamp));
        tlog('  Saved: QC_size_distributions_%s\n\n', timestamp);
    end

    % =========================================================
    %  PLOT 3 — Image montage  (requires data argument)
    % =========================================================
    if ~have_data
        tlog('  NOTE: data not supplied — skipping image montage (Plot 3).\n');
        tlog('        Call as:  plot_QC_full_dataset(ht, data)\n\n');
    else
        tlog('--- Plot 3: Image montage ---\n');

        % ----------------------------------------------------------
        % Build all tile images first, then display in one figure.
        % Grid size is computed automatically from nr_plates.
        % ----------------------------------------------------------
        n_montage_cols = ceil(sqrt(nr_plates));
        n_montage_rows = ceil(nr_plates / n_montage_cols);

        tlog('  Montage grid: %d rows x %d cols (%d plates)\n', ...
             n_montage_rows, n_montage_cols, nr_plates);

        % Tile pixel size: imcrop rect [x0, y0, width, height] gives
        % a (height+1) x (width+1) image, so 1001x1001 px here.
        tile_w = 1000;   % imcrop width  argument
        tile_h = 1000;   % imcrop height argument

        % Figure: each tile rendered at ~150 px on screen for readability
        px_per_tile = 150;
        fig3_w = n_montage_cols * px_per_tile;
        fig3_h = n_montage_rows * px_per_tile;

        fig3 = figure('Name', 'Colony image montage', ...
                      'Units', 'pixels', ...
                      'Position', [screen_w*0.04, screen_h*0.04, fig3_w, fig3_h], ...
                      'Color', 'w');

        for i = 1:nr_plates
            mask_clean = data.processed{i}.colonies.debug.segmented;
            img_gray   = imadjust(rgb2gray(data.processed{i}.img_final), [0.025, 0.15], []);
            img_blend  = imfuse(mask_clean, img_gray, 'blend');
            tile       = imcrop(img_blend, [1500, 1500, tile_w, tile_h]);

            ax = subplot(n_montage_rows, n_montage_cols, i);
            imshow(tile, 'Parent', ax);
            hold(ax, 'on');

            % Colony count label: top-left corner of the tile in axes
            % normalised units (0=left/top, 1=right/bottom).
            % Use axes normalised coordinates so placement is
            % independent of tile pixel dimensions.
            n_count = colony_count_final(i);
            text(ax, 0.02, 0.98, num2str(n_count), ...
                 'Units',              'normalized', ...
                 'Color',              'k', ...
                 'FontSize',           7, ...
                 'FontWeight',         'bold', ...
                 'HorizontalAlignment','left', ...
                 'VerticalAlignment',  'top', ...
                 'BackgroundColor',    'w', ...
                 'Margin',             1);
        end

        % Remove any unused subplot slots (partial last row)
        for i = nr_plates+1 : n_montage_rows*n_montage_cols
            ax_empty = subplot(n_montage_rows, n_montage_cols, i);
            axis(ax_empty, 'off');
        end

        montage_name = sprintf('QC_montage_%s', timestamp);
        save_figure(fig3, out_dir, montage_name);
        tlog('  Saved: %s\n\n', montage_name);
    end

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


%% ================================================================
%  LOCAL FUNCTION: round_to_nice
%  Round x up to a "nice" tick step (1, 2, 5, 10, 20, 50, ...).
% ================================================================
function s = round_to_nice(x)
    if x <= 0, s = 1; return; end
    mag   = 10^floor(log10(x));
    frac  = x / mag;
    if     frac <= 1,  s = 1   * mag;
    elseif frac <= 2,  s = 2   * mag;
    elseif frac <= 5,  s = 5   * mag;
    else,              s = 10  * mag;
    end
end


%% ================================================================
%  LOCAL FUNCTION: fmt_p
%  Text display for a p-value. Pearson/Spearman p (via corr()) can
%  legitimately underflow to exact 0.0 in double precision for large
%  n / strong correlation — the true value is just too small to
%  represent, not actually zero — shown as "<1e-300" instead.
% ================================================================
function s = fmt_p(p, fmt)
    if nargin < 2, fmt = '%.4g'; end
    if isnan(p)
        s = 'NaN';
    elseif p == 0
        s = '<1e-300';
    else
        s = sprintf(fmt, p);
    end
end
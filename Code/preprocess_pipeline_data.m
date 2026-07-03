%% HeteroTyper Pipeline
%
%%  PURPOSE
%  -------
%  Performs every calculation that is currently repeated across all
%  plotting/export scripts in one single pass and deposits the results
%  into the MATLAB base workspace so any downstream script can use them
%  directly without re-computing.
%
%%  USAGE
%  -----
%    preprocess_pipeline_data(data);          % uses data.params for col limits
%    preprocess_pipeline_data(data, cfg);     % pass optional configuration struct
%
%  After running, the following variables are available in the workspace:
%
%    ht                  - master struct (all outputs, see below)
%    ht.params           - all user/pipeline parameters
%    ht.ix               - {4x1} cell: plate indices per sample group
%    ht.labels           - discovered group labels, naturally sorted (e.g. {'2h','7h','24h'})
%    ht.colors           - [4x3] RGB color matrix
%    ht.groups(g)        - per-group processed data struct (see below)
%    ht.global           - global maxima across all groups (for axis scaling)
%    ht.xbins            - x-axis bin edge vectors for all 10 features
%    ht.pairwise         - pairwise KS D / Q90 statistics + normality-aware
%                          hypothesis test (Welch t-test or Wilcoxon rank-sum,
%                          see .p/.test_name/.normal_A/.normal_B) for every
%                          sample-group pair and all 10 features (STEP 6).
%                          .p_fdr is the Benjamini-Hochberg FDR-corrected
%                          p-value, corrected per feature across all group
%                          pairs tested for that feature.
%    ht.biorep_group_tests(r) - same hypothesis test as ht.pairwise, but run
%                          separately per biological replicate (bio_rep_col),
%                          comparing sample groups using ONLY that replicate's
%                          colonies (STEP 6i) — see field list at that step.
%
%  Per-group struct  ht.groups(g)  contains:
%    .label              - group label string as found in metadata (e.g. '24h')
%    .plate_indices      - vector of plate indices in this group
%    .n_colonies         - total valid colonies accumulated
%    .lag_time           - [Nx1] lag time + incTime (h)
%    .size               - [Nx1] final colony size (px)
%    .area               - [Nx1] area (px)
%    .intensity          - [Nx1] total intensity
%    .mean_intensity     - [Nx1] mean gray intensity
%    .int_per_size       - [Nx1] intensity / size
%    .perimeter          - [Nx1] perimeter (px)
%    .circularity        - [Nx1] circularity
%    .eccentricity       - [Nx1] eccentricity
%    .solidity           - [Nx1] solidity
%    .centroid           - [Nx2] [x y] centroid coordinates
%    .size_timecourse    - [T x N] smoothed size timecourse matrix
%    .time_h             - [Tx1] elapsed time vector (h, incl. incTime)
%    .doublingTime       - [Nx1] colony doubling times (h)  (NaN if not computable)
%    .hist               - struct of normalised histogram counts for all 10 features
%    .gini               - struct of Gini indices for all 10 features
%    .stats              - struct of median / IQR / n per feature
%
%  Additionally, these workspace variables are created for backward
%  compatibility with scripts that reference them by name:
%    combined_lag_time_<label>  (one per discovered group)
%
%  A log file is automatically generated in the output directory:
%    preprocess_log_YYYYMMDD_HHMMSS.txt

function preprocess_pipeline_data(data, cfg)

    % =========================================================
    %  STEP 0 - Default config + optional override
    % =========================================================
    if nargin < 2
        cfg = struct();
    end

    % --- Colony count limits (prefer data.params, then cfg, then defaults) ---
    if isfield(data, 'params') && isfield(data.params, 'min_colony_nr')
        default_min_col = data.params.min_colony_nr;
        default_max_col = data.params.max_colony_nr;
    else
        default_min_col = 5;
        default_max_col = 700;
    end

    p.min_col        = get_cfg(cfg, 'min_col',        default_min_col);
    p.max_col        = get_cfg(cfg, 'max_col',        default_max_col);
    p.size_threshold = get_cfg(cfg, 'size_threshold', 100);
    p.ecc_threshold  = get_cfg(cfg, 'ecc_threshold',  0.70);
    p.lag_step       = get_cfg(cfg, 'lag_step',       0.5);

    nr_plates = length(data.processed);

    % =========================================================
    %  STEP 1 - Interactive parameter form (terminal)
    % =========================================================
    fprintf('\n============================================================\n');
    fprintf('HeteroTyper Pipeline - Parameter Configuration\n');
    fprintf('============================================================\n');
    
    user = get_user_params(data, nr_plates);

    p.incTime   = user.incTime;
    p.img_int   = user.img_int;
    p.max_lag   = user.max_lag;
    p.out_dir   = user.out_dir;
    p.group_col = user.group_col;
    p.bio_rep_col  = user.bio_rep_col;
    p.tech_rep_col = user.tech_rep_col;

    % =========================================================
    %  STEP 1B - Initialize logging system (AFTER all user input)
    % =========================================================
    log_filename = sprintf('preprocess_log_%s.txt', datetime('now','Format','yyyyMMdd_HHmmss'));
    log_filepath = fullfile(p.out_dir, log_filename);
    
    fid = fopen(log_filepath, 'w');
    if fid == -1
        error('Could not create log file at: %s', log_filepath);
    end
    
    % Write header to log
    log_message(fid, '========================================');
    log_message(fid, 'HeteroTyper Pipeline - Preprocessing Log');
    log_message(fid, '========================================');
    log_message(fid, sprintf('Started: %s', datetime('now','Format','yyyy-MM-dd HH:mm:ss')));
    log_message(fid, '');

    fprintf('\n--- Parameters confirmed ---\n');
    log_message(fid, '--- Parameters confirmed ---');
    
    fprintf('  RT incubation time  : %.4g h\n', p.incTime);
    log_message(fid, sprintf('  RT incubation time  : %.4g h', p.incTime));
    
    fprintf('  Imaging interval    : %.4g h\n', p.img_int);
    log_message(fid, sprintf('  Imaging interval    : %.4g h', p.img_int));
    
    fprintf('  Max lag time        : %.4g h\n', p.max_lag);
    log_message(fid, sprintf('  Max lag time        : %.4g h', p.max_lag));
    
    fprintf('  Min colonies/plate  : %d\n',     p.min_col);
    log_message(fid, sprintf('  Min colonies/plate  : %d', p.min_col));
    
    fprintf('  Max colonies/plate  : %d\n',     p.max_col);
    log_message(fid, sprintf('  Max colonies/plate  : %d', p.max_col));
    
    fprintf('  Size threshold      : %d px\n',  p.size_threshold);
    log_message(fid, sprintf('  Size threshold      : %d px', p.size_threshold));
    
    fprintf('  Ecc threshold       : %.2f\n',   p.ecc_threshold);
    log_message(fid, sprintf('  Ecc threshold       : %.2f', p.ecc_threshold));
    
    fprintf('  Grouping column     : %s\n',     p.group_col);
    log_message(fid, sprintf('  Grouping column     : %s', p.group_col));
    
    bio_rep_display  = '(not set)';
    tech_rep_display = '(not set)';
    if isfield(p, 'bio_rep_col')  && ~isempty(p.bio_rep_col),  bio_rep_display  = p.bio_rep_col;  end
    if isfield(p, 'tech_rep_col') && ~isempty(p.tech_rep_col), tech_rep_display = p.tech_rep_col; end
    fprintf('  Biological rep col  : %s\n',     bio_rep_display);
    log_message(fid, sprintf('  Biological rep col  : %s', bio_rep_display));
    fprintf('  Technical rep col   : %s\n',     tech_rep_display);
    log_message(fid, sprintf('  Technical rep col   : %s', tech_rep_display));
    
    fprintf('  Output directory    : %s\n',     p.out_dir);
    log_message(fid, sprintf('  Output directory    : %s', p.out_dir));
    
    fprintf('----------------------------\n\n');
    log_message(fid, '----------------------------');
    log_message(fid, '');

    % =========================================================
    %  STEP 2 - Filter plates: growth flag + final colony count
    % =========================================================
    growth_available = false(nr_plates, 1);
    fprintf('\n--- STEP 2: QC Filtering (based on FINAL segmented colonies) ---\n');
    log_message(fid, '--- STEP 2: QC Filtering (based on FINAL segmented colonies) ---');
    
    for i = 1:nr_plates
        growth_flag = data.processed{i}.growth_quant;
        
        % Count FINAL valid colonies (after all processing/segmentation)
        % Look at the actual segmented data structure
        if isfield(data.processed{i}.colonies, 'new') && ...
           isfield(data.processed{i}.colonies.new, 'lag_time')
            % Count non-NaN lag times = final valid colonies
            n_col = sum(~isnan(data.processed{i}.colonies.new.lag_time));
        elseif isfield(data.processed{i}.colonies, 'region_props') && ...
               isfield(data.processed{i}.colonies.region_props, 'size')
            % Fallback: count region properties
            n_col = length(data.processed{i}.colonies.region_props.size);
        else
            n_col = 0;
        end
        
        growth_available(i) = growth_flag;
        
        % Check filtering criteria based on FINAL colony count
        col_count_pass = (n_col >= p.min_col) && (n_col <= p.max_col);
        passes_qc = growth_flag && col_count_pass;
        
        if ~passes_qc
            reason = '';
            if ~growth_flag
                reason = 'growth_quant=False';
            end
            if ~col_count_pass
                if ~isempty(reason)
                    reason = [reason ' AND '];
                end
                if n_col < p.min_col
                    reason = sprintf('%scolony_count=%d < min=%d', reason, n_col, p.min_col);
                else
                    reason = sprintf('%scolony_count=%d > max=%d', reason, n_col, p.max_col);
                end
            end
            detail_msg = sprintf('  Plate %2d final_colonies=%3d: REJECTED - %s', i, n_col, reason);
            fprintf('%s\n', detail_msg);
            log_message(fid, detail_msg);
        else
            detail_msg = sprintf('  Plate %2d final_colonies=%3d: PASSED', i, n_col);
            fprintf('%s\n', detail_msg);
            log_message(fid, detail_msg);
        end
    end
    
    ix_growth = find(growth_available);
    log_message(fid, '');
    
    qc_msg = sprintf('Plates passing QC filter: %d / %d', length(ix_growth), nr_plates);
    fprintf('\n%s\n', qc_msg);
    log_message(fid, qc_msg);

    % =========================================================
    %  STEP 3 - Split into sample groups (detected dynamically)
    % =========================================================
    all_times = data.metadata.original.(p.group_col);

    % Guard: metadata table may be shorter than data.processed
    n_times = length(all_times);
    out_of_range = ix_growth(ix_growth > n_times);
    if ~isempty(out_of_range)
        warn_msg = sprintf(['  WARNING  %d plate(s) in ix_growth exceed the metadata ' ...
                 'table length (%d): indices [%s].\n' ...
                 '           These plates will be excluded from group assignment.'], ...
                length(out_of_range), n_times, num2str(out_of_range(:)'));
        fprintf('%s\n', warn_msg);
        log_message(fid, warn_msg);
        ix_growth = ix_growth(ix_growth <= n_times);
    end

    % Discover labels dynamically from the data
    present_raw = unique(all_times(ix_growth));
    present_raw = present_raw(~cellfun(@isempty, present_raw));
    % Strip surrounding apostrophes that some importers embed in strings
    % (e.g. "'3h'" stored as a 4-char string -> strip to "3h")
    present_raw = cellfun(@strip_quotes, present_raw, 'UniformOutput', false);

    % Regime-aware stable sort
    num_prefix = regexp(present_raw, '^\d+(\.\d+)?', 'match', 'once');
    has_num    = ~cellfun(@isempty, num_prefix);

    % Numeric-prefix group: sort by numeric value
    idx_num  = find(has_num);
    vals_num = cellfun(@str2double, num_prefix(idx_num));
    [~, o]   = sort(vals_num);
    idx_num  = idx_num(o);

    % Pure-text group: sort alphabetically (case-insensitive)
    idx_txt  = find(~has_num);
    [~, o]   = sort(lower(present_raw(idx_txt)));
    idx_txt  = idx_txt(o);

    labels = present_raw([idx_num; idx_txt]);

    sort_msg = sprintf('  Label sort: %d numeric-prefix, %d text-only  -> [%s]', ...
            length(idx_num), length(idx_txt), strjoin(labels, ', '));
    fprintf('%s\n', sort_msg);
    log_message(fid, sort_msg);

    n_groups = length(labels);
    colors   = generate_group_colors(n_groups);

    % Strip apostrophes from all_times for consistent matching
    all_times_clean = cellfun(@strip_quotes, all_times, 'UniformOutput', false);

    ix = cell(n_groups, 1);
    for g = 1:n_groups
        ix{g} = intersect(find(strcmp(all_times_clean, labels{g})), ix_growth);
        group_msg = sprintf('  Group %-10s : %d plates  color [%.2f %.2f %.2f]', ...
                labels{g}, length(ix{g}), colors(g,1), colors(g,2), colors(g,3));
        fprintf('%s\n', group_msg);
        log_message(fid, group_msg);
    end
    log_message(fid, '');

    % =========================================================
    %  STEP 4 - Accumulate all per-group data in one pass
    % =========================================================
    accum_msg = 'Accumulating colony data across all groups...';
    fprintf('%s\n', accum_msg);
    log_message(fid, accum_msg);

    % Global maxima accumulators (for x-axis scaling shared across groups)
    gmax.size     = 0;
    gmax.area     = 0;
    gmax.int      = 0;
    gmax.mean_int = 0;
    gmax.int_size = 0;
    gmax.peri     = 0;
    gmax.img_count = 0;

    groups(n_groups) = struct();

    % Per-plate size cache — populated inside the plate loop below.
    % per_plate_size{i} holds the size-threshold-passing colony sizes
    % (rp_clean.Area, fail_size masked) for plate i.  Used by
    % plot_morphology_colonies_Gini without needing raw data.
    per_plate_size = cell(nr_plates, 1);

    % Loop over each group and accumulate colony data
    % Use colonies.new which contains ONLY the final filtered colonies from detection (Step 4)
    % These are the exact colonies shown in the colored segmented images - NO additional filtering
    for g = 1:n_groups
        group_label = labels{g};
        groups(g).label = group_label;
        groups(g).plate_indices = ix{g};
        groups(g).n_colonies = 0;
        
        % Pre-allocate arrays for this group
        max_colonies_est = 10000;
        size_acc        = nan(max_colonies_est, 1);
        area_acc        = nan(max_colonies_est, 1);
        int_acc         = nan(max_colonies_est, 1);
        mean_int_acc    = nan(max_colonies_est, 1);
        int_size_acc    = nan(max_colonies_est, 1);
        peri_acc        = nan(max_colonies_est, 1);
        circ_acc        = nan(max_colonies_est, 1);
        ecc_acc         = nan(max_colonies_est, 1);
        solid_acc       = nan(max_colonies_est, 1);
        centroid_acc    = nan(max_colonies_est, 2);
        lag_time_acc    = nan(max_colonies_est, 1);
        doubling_acc    = nan(max_colonies_est, 1);
        % Pass/fail mask for the final-timepoint size threshold
        % (p.size_threshold, default 100 px). Entries that fail are
        % NaN-masked on every per-colony feature accumulator below so
        % they are excluded from all histograms, medians, IQR, Gini, and
        % CDF curves while still being counted in n_colonies.
        pass_size_acc    = true(max_colonies_est, 1);
        n_size_rejected  = 0;

        % Timecourse accumulators — grown lazily once the time-grid length
        % is known from the first plate encountered in this group.
        % tc_acc  : [T x max_colonies_est] — NaN-initialised, filled per plate.
        % tc_T    : scalar, number of imaging timepoints (determined on first plate).
        % tc_time : [T x 1] elapsed-time vector (h) = incTime + (0:T-1)*img_int.
        tc_acc  = [];   % allocated on first plate with valid timecourse
        tc_T    = 0;
        tc_time = [];

        colony_idx = 0;

        % Per-plate colony counts, aligned with ix{g} — lets downstream
        % steps (e.g. per-bio-rep pooling) slice out exactly the colonies
        % belonging to a single plate from the concatenated groups(g).*
        % vectors below. Plates skipped by any of the `continue`s stay 0.
        plate_ncol = zeros(length(ix{g}), 1);

        % Loop over all plates in this group
        for plate_idx = 1:length(ix{g})
            i = ix{g}(plate_idx);

            if ~isfield(data.processed{i}, 'colonies')
                continue;
            end

            colonies = data.processed{i}.colonies;

            % Use colonies.new - these are the FINAL filtered colonies (Step 4 of detection)
            if ~isfield(colonies, 'new')
                continue;
            end

            c_new = colonies.new;
            n_colonies_here = length(c_new.lag_time);

            if n_colonies_here == 0
                continue;
            end
            plate_ncol(plate_idx) = n_colonies_here;
            
            % For geometry: map back to region_props_clean if available
            rp_clean = [];
            if isfield(colonies, 'region_props_clean')
                rp_clean = colonies.region_props_clean;
            elseif isfield(colonies, 'region_props')
                % Use flag_colony_ok to find which ones to use
                rp_all = colonies.region_props;
                if isfield(colonies, 'flag_colony_ok') && length(colonies.flag_colony_ok) == height(rp_all)
                    clean_idx = find(colonies.flag_colony_ok);
                    if length(clean_idx) == n_colonies_here
                        rp_clean = rp_all(clean_idx, :);
                    end
                end
            end
            
            % Extract geometry from region_props_clean
            if ~isempty(rp_clean) && height(rp_clean) == n_colonies_here
                size_acc(colony_idx + (1:n_colonies_here))     = rp_clean.Area;
                area_acc(colony_idx + (1:n_colonies_here))     = rp_clean.Area;
                ecc_acc(colony_idx + (1:n_colonies_here))      = rp_clean.Eccentricity;
                peri_acc(colony_idx + (1:n_colonies_here))     = rp_clean.Perimeter;
                circ_acc(colony_idx + (1:n_colonies_here))     = rp_clean.Circularity;
                solid_acc(colony_idx + (1:n_colonies_here))    = rp_clean.Solidity;
                centroid_acc(colony_idx + (1:n_colonies_here), :) = rp_clean.Centroid;
            end
            
            % Store the final filtered colony data from colonies.new.
            %
            % Two corrections applied here:
            %
            % 1) Add RT incubation time so that lag_time represents total
            %    elapsed time from inoculation, not just imaging-phase time.
            %    This matches the .lag_time docstring contract (line 31)
            %    and the histogram bin edges below
            %    (p.incTime : p.lag_step : p.max_lag).
            %
            % 2) Right-censoring: colonies that never grew within the
            %    imaging window are stored upstream with lag_time pinned
            %    at (or fractionally above) max_lag. These are NOT real
            %    measured lag times -- they only mean "the colony did not
            %    start growing before imaging ended". Including them
            %    massively biases the median and produces a fake spike at
            %    the right edge of the histogram. We convert them to NaN
            %    so they are excluded from the histogram, median, IQR,
            %    and Gini index, while still being counted in n_colonies.
            %    The number excluded per group is logged below.
            lt_raw = c_new.lag_time(:) + p.incTime;
            cens_tol = 0.5 * p.lag_step;
            is_censored = ~isnan(lt_raw) & (lt_raw >= p.max_lag - cens_tol);
            lt_raw(is_censored) = NaN;
            lag_time_acc(colony_idx + (1:n_colonies_here)) = lt_raw;
            
            % Extract intensity metrics - must be careful about indexing!
            % colonies.new only has data for colonies that passed detection filters
            
            % Method 1: Use final_col_size which should be available
            % This is the final intensity from the timecourse (last timepoint)
            if isfield(c_new, 'final_col_size') && ~isempty(c_new.final_col_size)
                % final_col_size might actually be size, not intensity
                % Check if intensity data exists
            end
            
            % Method 2: Extract from timecourse_intensity_smoothed
            % The last row contains the final intensity for each colony
            if isfield(colonies, 'timecourse_intensity_smoothed') && ...
               ~isempty(colonies.timecourse_intensity_smoothed)
                tc_int = colonies.timecourse_intensity_smoothed;
                if size(tc_int, 2) == n_colonies_here
                    % Last row = final intensity for each colony
                    final_intensity = tc_int(end, :)';
                    int_acc(colony_idx + (1:n_colonies_here)) = final_intensity;
                end
            end
            
            % Method 3: Mean intensity from region_props_clean
            % Match it with the clean colonies that passed filters
            if ~isempty(rp_clean) && height(rp_clean) == n_colonies_here
                % Try different field names for mean intensity
                if isfield(rp_clean, 'MeanIntensity')
                    mean_int_acc(colony_idx + (1:n_colonies_here)) = rp_clean.MeanIntensity;
                elseif isfield(rp_clean, 'mean_intensity')
                    mean_int_acc(colony_idx + (1:n_colonies_here)) = rp_clean.mean_intensity;
                else
                    % Compute mean intensity from total intensity / area
                    int_vals = int_acc(colony_idx + (1:n_colonies_here));
                    size_vals = size_acc(colony_idx + (1:n_colonies_here));
                    valid_mask = ~isnan(int_vals) & ~isnan(size_vals) & size_vals > 0;
                    if any(valid_mask)
                        mean_int_temp = int_vals;
                        mean_int_temp(valid_mask) = int_vals(valid_mask) ./ size_vals(valid_mask);
                        mean_int_acc(colony_idx + (1:n_colonies_here)) = mean_int_temp;
                    end
                end
            end
            
            % Method 4: Compute intensity/size ratio if both available
            size_vals = size_acc(colony_idx + (1:n_colonies_here));
            int_vals = int_acc(colony_idx + (1:n_colonies_here));
            valid_mask = ~isnan(size_vals) & ~isnan(int_vals) & size_vals > 0;
            if any(valid_mask)
                int_per_size_temp = int_vals;
                int_per_size_temp(valid_mask) = int_vals(valid_mask) ./ size_vals(valid_mask);
                int_size_acc(colony_idx + (1:n_colonies_here)) = int_per_size_temp;
            end

            % ----------------------------------------------------------
            % FINAL-TIMEPOINT SIZE FILTER
            % Any colony whose final segmented size is below
            % p.size_threshold (default 100 px) is treated as a non-pass
            % detection (specks, debris, fragments). Such slots get
            % NaN-masked across every per-colony feature so they are
            % excluded from histograms, medians, IQR, Gini, CDF curves,
            % and intensity-derived metrics. They are still counted in
            % n_colonies; their tally is exposed as groups(g).n_size_rejected
            % and logged in the per-group accumulation summary.
            slot_range = colony_idx + (1:n_colonies_here);
            size_slot  = size_acc(slot_range);
            fail_size  = ~isnan(size_slot) & (size_slot < p.size_threshold);
            if any(fail_size)
                fail_idx = slot_range(fail_size);
                size_acc(fail_idx)        = NaN;
                area_acc(fail_idx)        = NaN;
                int_acc(fail_idx)         = NaN;
                mean_int_acc(fail_idx)    = NaN;
                int_size_acc(fail_idx)    = NaN;
                peri_acc(fail_idx)        = NaN;
                circ_acc(fail_idx)        = NaN;
                ecc_acc(fail_idx)         = NaN;
                solid_acc(fail_idx)       = NaN;
                centroid_acc(fail_idx, :) = NaN;
                lag_time_acc(fail_idx)    = NaN;
                doubling_acc(fail_idx)    = NaN;
                pass_size_acc(fail_idx)   = false;
                n_size_rejected = n_size_rejected + sum(fail_size);
            end

            % Exclude ALL morphological features for lag-censored colonies.
            % A censored colony never completed its lag phase within the
            % imaging window.  Its final size, shape, and intensity reflect
            % partial (ongoing) growth, not a settled phenotype, so including
            % it in any feature distribution or Gini/stats calculation would
            % mix two fundamentally different biological states.
            % lag_time is already NaN for these colonies (set above).
            % We NaN every other feature here so all downstream analyses
            % (histograms, medians, IQR, Gini, correlations, bio-rep stats)
            % automatically exclude them via the standard ~isnan() filter.
            if any(is_censored)
                cens_idx = slot_range(is_censored);
                size_acc(cens_idx)        = NaN;
                area_acc(cens_idx)        = NaN;
                int_acc(cens_idx)         = NaN;
                mean_int_acc(cens_idx)    = NaN;
                int_size_acc(cens_idx)    = NaN;
                peri_acc(cens_idx)        = NaN;
                circ_acc(cens_idx)        = NaN;
                ecc_acc(cens_idx)         = NaN;
                solid_acc(cens_idx)       = NaN;
                centroid_acc(cens_idx, :) = NaN;
                doubling_acc(cens_idx)    = NaN;
                if tc_T > 0
                    tc_acc(:, cens_idx)   = NaN;
                end
            end

            % Cache per-plate final colony sizes (size-threshold filtered).
            % Captured here — after fail_size masking — so only colonies
            % that passed the size threshold are stored.
            sz_plate = size_acc(slot_range);
            per_plate_size{i} = sz_plate(~isnan(sz_plate));

            % ----------------------------------------------------------
            % TIMECOURSE ACCUMULATION
            % Collect timecourse_size_smoothed from colonies.new.
            % The field is [T x n_colonies_here] (rows = timepoints,
            % cols = colonies), matching the structure expected by
            % plot_combined_samples_growth_curves.
            % ----------------------------------------------------------
            if isfield(c_new, 'timecourse_size_smoothed') && ...
               ~isempty(c_new.timecourse_size_smoothed)

                tc_plate = c_new.timecourse_size_smoothed;  % [T_plate x n_colonies_here]

                % On first plate with valid timecourse: fix the time grid
                % and allocate the accumulator matrix.
                if tc_T == 0
                    tc_T    = size(tc_plate, 1);
                    tc_time = p.incTime + (0 : tc_T-1)' * p.img_int;   % [T x 1]
                    tc_acc  = nan(tc_T, max_colonies_est);
                end

                T_plate = size(tc_plate, 1);

                if T_plate == tc_T
                    % Same length as the established grid — copy directly.
                    tc_acc(:, colony_idx + (1:n_colonies_here)) = tc_plate;
                elseif T_plate < tc_T
                    % Shorter plate: pad with NaN at the end.
                    tc_acc(1:T_plate,  colony_idx + (1:n_colonies_here)) = tc_plate;
                else
                    % Longer plate: truncate to the established grid length.
                    tc_acc(:, colony_idx + (1:n_colonies_here)) = tc_plate(1:tc_T, :);
                end

                % NaN-mask the same colonies that failed the size threshold.
                if any(fail_size)
                    tc_acc(:, slot_range(fail_size)) = NaN;
                end
            end

            colony_idx = colony_idx + n_colonies_here;

            % Update global maxima
            gmax.size     = max(gmax.size,     nanmax(size_acc(1:colony_idx)));
            gmax.area     = max(gmax.area,     nanmax(area_acc(1:colony_idx)));
            gmax.int      = max(gmax.int,      nanmax(int_acc(1:colony_idx)));
            gmax.mean_int = max(gmax.mean_int, nanmax(mean_int_acc(1:colony_idx)));
            gmax.int_size = max(gmax.int_size, nanmax(int_size_acc(1:colony_idx)));
            gmax.peri     = max(gmax.peri,     nanmax(peri_acc(1:colony_idx)));
        end
        
        % Trim arrays and store in group struct
        groups(g).n_colonies   = colony_idx;
        groups(g).size         = size_acc(1:colony_idx);
        groups(g).area         = area_acc(1:colony_idx);
        groups(g).intensity    = int_acc(1:colony_idx);
        groups(g).mean_intensity = mean_int_acc(1:colony_idx);
        groups(g).int_per_size = int_size_acc(1:colony_idx);
        groups(g).perimeter    = peri_acc(1:colony_idx);
        groups(g).circularity  = circ_acc(1:colony_idx);
        groups(g).eccentricity = ecc_acc(1:colony_idx);
        groups(g).solidity     = solid_acc(1:colony_idx);
        groups(g).centroid     = centroid_acc(1:colony_idx, :);
        groups(g).lag_time     = lag_time_acc(1:colony_idx);
        groups(g).doubling_time = doubling_acc(1:colony_idx);

        % Per-plate segment boundaries into the concatenated vectors above.
        % plate_indices(k) contributed rows plate_seg_start(k):plate_seg_end(k)
        % (seg_end < seg_start for a plate with 0 colonies, i.e. an empty range).
        groups(g).plate_ncol      = plate_ncol;
        groups(g).plate_seg_end   = cumsum(plate_ncol);
        groups(g).plate_seg_start = groups(g).plate_seg_end - plate_ncol + 1;

        % ------------------------------------------------------------------
        % TIMECOURSE: trim accumulator to actual colony count and store.
        %   .size_timecourse  — [T x N] smoothed colony-size matrix
        %                        (NaN for colonies rejected by size threshold)
        %   .time_h           — [T x 1] elapsed time in hours (incl. incTime)
        % If no plate in this group had a valid timecourse field the fields
        % are stored as empty so downstream scripts can detect the gap.
        % ------------------------------------------------------------------
        if tc_T > 0 && colony_idx > 0
            groups(g).size_timecourse = tc_acc(1:tc_T, 1:colony_idx);
            groups(g).time_h          = tc_time;
        else
            groups(g).size_timecourse = [];
            groups(g).time_h          = [];
        end

        % ------------------------------------------------------------------
        % DOUBLING TIME: fit log-linear model to each colony's size timecourse
        %   S(t) = S0 * exp(r*t)  =>  ln S = ln S0 + r*t
        %   doubling time = ln(2) / r  (NaN if r <= 0 or fit not possible)
        % Only the post-lag portion is used: frames where t > lag_time.
        % Falls back to all-NaN if no timecourse is available.
        % ------------------------------------------------------------------
        dt_vec = nan(colony_idx, 1);
        if tc_T > 0 && colony_idx > 0 && ~isempty(tc_time)
            tc_mat  = tc_acc(1:tc_T, 1:colony_idx);  % [T x N]
            lt_vec  = lag_time_acc(1:colony_idx);     % [N x 1] lag times (h)
            for ci = 1:colony_idx
                tc_col = tc_mat(:, ci);               % [T x 1] size timecourse
                lt_ci  = lt_vec(ci);
                if isnan(lt_ci), continue; end        % censored / size-failed
                % Use frames strictly after lag time and with positive size
                post_mask = (tc_time > lt_ci) & (tc_col > 0) & ~isnan(tc_col);
                if sum(post_mask) < 2, continue; end  % need >=2 points to fit
                t_fit = tc_time(post_mask);
                s_fit = log(tc_col(post_mask));
                % Linear regression: s_fit = a + r*t_fit
                X = [ones(length(t_fit),1), t_fit];
                coeffs = X \ s_fit;
                r = coeffs(2);                        % growth rate (h^-1)
                if r > 0
                    dt_vec(ci) = log(2) / r;
                end
            end
        end
        groups(g).doublingTime  = dt_vec;   % camelCase — matches plot_doublingTime
        groups(g).doubling_time = dt_vec;   % snake_case alias kept for compatibility

        % Final-timepoint size filter bookkeeping (see filter block above).
        % pass_size       : logical vector aligned with all per-colony fields;
        %                   true = passed size threshold, false = NaN-masked
        %                   as a non-pass detection.
        % n_size_rejected : count of NaN-masked colonies in this group.
        groups(g).pass_size       = pass_size_acc(1:colony_idx);
        groups(g).n_size_rejected = n_size_rejected;

        % Right-censoring bookkeeping (see censoring block above).
        % n_lag_censored : colonies whose lag_time was pinned at max_lag
        %                  and converted to NaN above. They are still
        %                  counted in n_colonies, but excluded from all
        %                  lag-time stats, histograms, and CDFs.
        % NOTE: count only true right-censored colonies; entries that
        % became NaN because they failed the size threshold are tallied
        % separately as n_size_rejected.
        % n_lag_grew     : colonies with a real (non-censored) lag time.
        is_censored_lag = isnan(groups(g).lag_time) & groups(g).pass_size;
        groups(g).n_lag_censored = sum(is_censored_lag);
        groups(g).n_lag_grew     = sum(~isnan(groups(g).lag_time));
        % Compute Gini indices and stats for all 10 features.
        % Histograms are computed after xbins is finalised (see below).
        groups(g).hist  = struct();
        groups(g).gini  = struct();
        groups(g).stats = struct();

        feature_names = {'lag_time', 'size', 'area', 'intensity', 'mean_intensity', ...
                        'int_per_size', 'perimeter', 'circularity', 'eccentricity', 'solidity'};

        for feat_idx = 1:length(feature_names)
            feat     = feature_names{feat_idx};
            data_vec = groups(g).(feat);

            if isempty(data_vec)
                groups(g).gini.(feat)  = 0;
                groups(g).hist.(feat)  = [];
                continue;
            end

            groups(g).gini.(feat) = compute_gini(data_vec);

            data_clean = data_vec(~isnan(data_vec));
            if ~isempty(data_clean)
                groups(g).stats.(feat) = struct( ...
                    'median', median(data_clean), ...
                    'iqr', iqr(data_clean), ...
                    'n', length(data_clean));
            else
                groups(g).hist.(feat) = [];
            end
        end
        
        % Print distribution statistics (median, IQR, Gini) for this group
        stat_hdr = sprintf('  %-16s  %10s  %10s  %8s  %8s', ...
            'Feature', 'Median', 'IQR', 'Gini', 'n');
        fprintf('%s\n', stat_hdr);
        log_message(fid, stat_hdr);
        stat_sep = sprintf('  %s', repmat('-', 1, 58));
        fprintf('%s\n', stat_sep);
        log_message(fid, stat_sep);
        feat_labels_print = {'Lag Time','Size','Area','Intensity','Mean Intensity', ...
                             'Int/Size','Perimeter','Circularity','Eccentricity','Solidity'};
        for feat_idx = 1:length(feature_names)
            fn = feature_names{feat_idx};
            if isfield(groups(g).stats, fn)
                s     = groups(g).stats.(fn);
                gv    = groups(g).gini.(fn);
                sline = sprintf('  %-16s  %10.4g  %10.4g  %8.4f  %8d', ...
                    feat_labels_print{feat_idx}, s.median, s.iqr, gv, s.n);
                fprintf('%s\n', sline);
                log_message(fid, sline);
            end
        end
        fprintf('\n');
        log_message(fid, '');

        n_cens = groups(g).n_lag_censored;
        n_grew = groups(g).n_lag_grew;
        n_rej  = groups(g).n_size_rejected;
        if colony_idx > 0
            cens_frac = 100 * n_cens / colony_idx;
            rej_frac  = 100 * n_rej  / colony_idx;
        else
            cens_frac = 0;
            rej_frac  = 0;
        end
        accum_detail = sprintf( ...
            ['  Group %s: %d colonies  (grew=%d, censored=%d [%.1f%%], ' ...
             'size<%d px rejected=%d [%.1f%%])'], ...
            group_label, colony_idx, n_grew, n_cens, cens_frac, ...
            p.size_threshold, n_rej, rej_frac);
        fprintf('%s\n', accum_detail);
        log_message(fid, accum_detail);
        if cens_frac > 25
            warn_msg = sprintf([ ...
                '  WARNING  Group %s: %.1f%% of colonies are right-censored.\n' ...
                '           Reported lag-time median/IQR reflect only growers and\n' ...
                '           understate the true central tendency of the population.'], ...
                group_label, cens_frac);
            fprintf('%s\n', warn_msg);
            log_message(fid, warn_msg);
        end
    end
    
    log_message(fid, '');

    % Print summary of STEP 4 results
    log_message(fid, '');
    log_message(fid, '--- STEP 4 SUMMARY ---');
    for g = 1:n_groups
        summary_msg = sprintf('Group %s: %d colonies', groups(g).label, groups(g).n_colonies);
        fprintf('%s\n', summary_msg);
        log_message(fid, summary_msg);
    end
    log_message(fid, '');
    
    % =========================================================
    %  FINALIZATION - Write summary and close log
    % =========================================================
    log_message(fid, '========================================');
    log_message(fid, sprintf('Completed: %s', datetime('now','Format','yyyy-MM-dd HH:mm:ss')));
    log_message(fid, '========================================');
    
    % Close the log file
    fclose(fid);
    
    % Confirm log location (use fprintf directly to avoid Windows path escape issues)
    fprintf('\nLog file saved to: %s\n', log_filepath);
    
    % =========================================================
    %  BUILD PER-PLATE QC ARRAYS  (used by plot_QC_full_dataset)
    %
    %  For every plate we compute three things:
    %    colony_count_all   - total detected colonies (length of flag_colony_ok)
    %    colony_count_clean - colonies passing flag_colony_ok
    %    colony_count_final - colonies surviving ALL filters applied in STEP 4:
    %                         colonies.new entries whose rp_clean.Area >= size_threshold
    %                         This is the count that matches ht.groups data exactly.
    %    plate_area         - cell array of per-plate Area vectors (rp_clean, size-filtered)
    %                         used by plot_QC_full_dataset Plot 2 without needing raw data
    % =========================================================
    qc_colony_count_all   = zeros(nr_plates, 1);
    qc_colony_count_clean = zeros(nr_plates, 1);
    qc_colony_count_final = zeros(nr_plates, 1);
    qc_passes             = false(nr_plates, 1);
    qc_plate_area         = cell(nr_plates, 1);   % per-plate Area (size-filtered)

    for i = 1:nr_plates
        if ~isfield(data.processed{i}, 'colonies')
            continue;
        end
        colonies_i = data.processed{i}.colonies;

        % --- raw / clean counts from flag_colony_ok ---
        if isfield(colonies_i, 'region_props') && ...
           isfield(colonies_i.region_props, 'flag_colony_ok')
            flags = colonies_i.region_props.flag_colony_ok;
            qc_colony_count_all(i)   = length(flags);
            qc_colony_count_clean(i) = sum(flags ~= 0);
        end

        % --- final count and area: mirror STEP 4 filtering ---
        if ~isfield(colonies_i, 'new') || ...
           ~isfield(colonies_i.new, 'lag_time')
            continue;
        end
        n_new = length(colonies_i.new.lag_time);
        if n_new == 0, continue; end

        % Resolve rp_clean (same logic as STEP 4)
        rp_c = [];
        if isfield(colonies_i, 'region_props_clean')
            rp_c = colonies_i.region_props_clean;
        elseif isfield(colonies_i, 'region_props')
            rp_all_i = colonies_i.region_props;
            if isfield(colonies_i, 'flag_colony_ok') && ...
               length(colonies_i.flag_colony_ok) == height(rp_all_i)
                cidx = find(colonies_i.flag_colony_ok);
                if length(cidx) == n_new
                    rp_c = rp_all_i(cidx, :);
                end
            end
        end

        if isempty(rp_c) || height(rp_c) ~= n_new
            % No geometry available — use n_new as best estimate
            qc_colony_count_final(i) = n_new;
            continue;
        end

        % Apply size threshold (same as STEP 4)
        area_i    = rp_c.Area;
        pass_mask = area_i >= p.size_threshold;

        qc_colony_count_final(i) = sum(pass_mask);
        qc_plate_area{i}         = area_i(pass_mask);   % size-filtered areas

        % Plate passes QC if growth_quant is true AND final count in [min_col, max_col]
        n_final_i = qc_colony_count_final(i);
        qc_passes(i) = data.processed{i}.growth_quant && ...
                       (n_final_i >= p.min_col) && (n_final_i <= p.max_col);
    end

    % Store plate filenames
    if isfield(data.metadata, 'fn')
        qc_fn = data.metadata.fn;
    else
        qc_fn = [];
    end

    % Ask user for manual count column
    qc_manual_col    = ask_manual_count_column(data);
    qc_manual_counts = get_manual_counts(data, qc_manual_col, nr_plates);

    % =========================================================
    %  BUILD PER-PLATE QC GROWTH ARRAYS  (used by plot_QC_growth_data)
    %
    %  For every plate in ix_growth (growth_quant==true) we store:
    %    qc_plate_time{i}     - [T x 1] elapsed time vector (h, incl. incTime)
    %    qc_plate_size_tc{i}  - [T x N] smoothed colony size timecourse
    %    qc_plate_valid{i}    - [N x 1] logical valid-colony mask (flag_colony_ok)
    %    qc_plate_lag{i}      - [N x 1] lag times (h, incl. incTime)
    %    qc_plate_early_dt{i} - [N x 1] early doubling times (h)
    %  All cell arrays are nr_plates x 1; non-growth plates hold empty cells.
    % =========================================================
    fprintf('\nSTEP QC-growth: Caching per-plate growth data for plot_QC_growth_data...\n');

    qc_plate_time     = cell(nr_plates, 1);
    qc_plate_size_tc  = cell(nr_plates, 1);
    qc_plate_valid    = cell(nr_plates, 1);
    qc_plate_lag      = cell(nr_plates, 1);
    qc_plate_early_dt = cell(nr_plates, 1);

    for i = ix_growth(:)'
        if ~isfield(data.processed{i}, 'colonies') || ...
           ~isfield(data.processed{i}.colonies, 'new')
            continue;
        end
        c_new_i = data.processed{i}.colonies.new;

        % --- Time vector ---
        if isfield(c_new_i, 'time_info') && isfield(c_new_i.time_info, 'elapsed_time_h')
            qc_plate_time{i} = c_new_i.time_info.elapsed_time_h(:) + p.incTime;
        end

        % --- Valid-colony mask (flag_colony_ok from region_props) ---
        flag_i = [];
        if isfield(data.processed{i}.colonies, 'region_props') && ...
           isfield(data.processed{i}.colonies.region_props, 'flag_colony_ok')
            flag_i = logical(data.processed{i}.colonies.region_props.flag_colony_ok);
        end
        qc_plate_valid{i} = flag_i;

        % --- Smoothed size timecourse [T x N_valid] ---
        if isfield(c_new_i, 'timecourse_size_smoothed') && ...
           ~isempty(c_new_i.timecourse_size_smoothed)
            tc_full = c_new_i.timecourse_size_smoothed;   % [T x N_all]
            if ~isempty(flag_i) && size(tc_full, 2) == length(flag_i)
                qc_plate_size_tc{i} = tc_full(:, flag_i);
            else
                qc_plate_size_tc{i} = tc_full;
            end
        end

        % --- Lag times (include incTime, apply flag mask) ---
        if isfield(c_new_i, 'lag_time')
            lt_i = c_new_i.lag_time(:) + p.incTime;
            if ~isempty(flag_i) && length(lt_i) == length(flag_i)
                qc_plate_lag{i} = lt_i(flag_i);
            else
                qc_plate_lag{i} = lt_i;
            end
        end

        % --- Early doubling time (apply flag mask) ---
        if isfield(c_new_i, 'early_doublingtime')
            edt_i = c_new_i.early_doublingtime(:);
            if ~isempty(flag_i) && length(edt_i) == length(flag_i)
                qc_plate_early_dt{i} = edt_i(flag_i);
            else
                qc_plate_early_dt{i} = edt_i;
            end
        end
    end
    fprintf('  Done.\n');

    % =========================================================
    %  BUILD XBINS FOR HISTOGRAM BINNING
    % =========================================================
    % Create bin edge vectors for all 10 features
    % These are used by plot_combined_samples for histogram display
    
    xbins = struct();
    
    % Lag time bins
    xbins.lag_step = p.lag_step;
    xbins.lag = p.incTime:p.lag_step:p.max_lag;
    
    % Size bins - round to nice number
    size_step = round_to_sig(gmax.size / 100, 1);  % ~100 bins
    xbins.size_step = size_step;
    xbins.size = 0:size_step:gmax.size*1.1;
    xbins.size_upper = max(xbins.size);
    
    % Area bins
    area_step = round_to_sig(gmax.area / 100, 1);
    xbins.area_step = area_step;
    xbins.area = 0:area_step:gmax.area*1.1;
    xbins.area_upper = max(xbins.area);
    
    % Perimeter bins
    peri_step = round_to_sig(gmax.peri / 100, 1);
    xbins.peri_step = peri_step;
    xbins.perimeter = 0:peri_step:gmax.peri*1.1;
    xbins.peri_upper = max(xbins.perimeter);
    
    % Circularity, Eccentricity, Solidity: already bounded [0-1]
    xbins.circularity = 0:0.01:1;
    xbins.eccentricity = 0:0.01:1;
    xbins.solidity = 0:0.01:1;
    

    % =========================================================
    %  COMPUTE HISTOGRAMS  (after xbins is finalised)
    %  All bin edges are now fixed to their global-maximum values,
    %  so every group's hist.* vector is the same length as the
    %  corresponding xp_* midpoint vector in plot_combined_samples.
    % =========================================================
    % feature_names_h / bin_map are reused as-is by STEP 6i below (per-bio-rep
    % histograms), so they are computed once here rather than per group.
    feature_names_h = {'lag_time', 'size', 'area', 'intensity', 'mean_intensity', ...
                        'int_per_size', 'perimeter', 'circularity', 'eccentricity', 'solidity'};
    bin_map = struct( ...
        'lag_time',       xbins.lag, ...
        'size',           xbins.size, ...
        'area',           xbins.area, ...
        'intensity',      (0 : round_to_sig(gmax.int      / 100, 1) : ceil(gmax.int      / round_to_sig(gmax.int      / 100, 1)) * round_to_sig(gmax.int      / 100, 1) + round_to_sig(gmax.int      / 100, 1)), ...
        'mean_intensity', (0 : round_to_sig(gmax.mean_int / 100, 2) : ceil(gmax.mean_int / round_to_sig(gmax.mean_int / 100, 2)) * round_to_sig(gmax.mean_int / 100, 2) + round_to_sig(gmax.mean_int / 100, 2)), ...
        'int_per_size',   (0 : round_to_sig(gmax.int_size / 100, 1) : ceil(gmax.int_size / round_to_sig(gmax.int_size / 100, 1)) * round_to_sig(gmax.int_size / 100, 1) + round_to_sig(gmax.int_size / 100, 1)), ...
        'perimeter',      xbins.perimeter, ...
        'circularity',    xbins.circularity, ...
        'eccentricity',   xbins.eccentricity, ...
        'solidity',       xbins.solidity);

    for g = 1:n_groups
        for fi = 1:length(feature_names_h)
            feat      = feature_names_h{fi};
            data_vec  = groups(g).(feat);
            bin_edges = bin_map.(feat);
            if isempty(data_vec) || isempty(bin_edges)
                groups(g).hist.(feat) = zeros(1, length(bin_edges) - 1);
                continue;
            end
            data_clean = data_vec(~isnan(data_vec));
            counts = histcounts(data_clean, bin_edges);
            total  = sum(counts);
            if total > 0
                groups(g).hist.(feat) = counts / total;
            else
                groups(g).hist.(feat) = counts;  % all-zero
            end
        end
    end
    % =========================================================
    %  BUILD ht.pairwise  — KS D / Q90 and a normality-aware two-sample
    %  hypothesis test (Welch t-test or Wilcoxon rank-sum) for all group
    %  pairs across all 10 features, computed once here so that
    %  plot_KS_quantile_statistics / plot_combined_samples can read them
    %  without raw data.
    %
    %  Hypothesis test selection per pair per feature:
    %    1. Test both groups' raw (unbinned) colony values for normality
    %       with the Lilliefors test (lillietest).
    %    2. If both are consistent with a normal distribution, use a
    %       two-sample Welch t-test (ttest2, unequal variances).
    %    3. Otherwise use the Wilcoxon rank-sum test (ranksum, i.e. the
    %       Mann-Whitney U test) — appropriate for the typically skewed,
    %       right-censored colony-size/lag-time distributions here.
    %  Raw per-colony values are used (not the binned probability-density
    %  histogram), since binning would discard within-bin variation the
    %  rank-based test relies on.
    % =========================================================
    pw_feat_names = {'Lag time','Colony size','Area','Intensity', ...
                     'Mean intensity','Intensity / size','Perimeter', ...
                     'Circularity','Eccentricity','Solidity'};
    pw_feat_fields = {'lag_time','size','area','intensity', ...
                      'mean_intensity','int_per_size','perimeter', ...
                      'circularity','eccentricity','solidity'};
    n_pw_feats = length(pw_feat_fields);

    % Enumerate all unique ordered pairs (i < j)
    n_pw_pairs  = n_groups * (n_groups - 1) / 2;
    pw_pair_idx    = zeros(n_pw_pairs, 2);
    pw_pair_labels = cell(n_pw_pairs, 1);
    pw_KS_D        = nan(n_pw_pairs, n_pw_feats);
    pw_Q90         = nan(n_pw_pairs, n_pw_feats);
    pw_p           = nan(n_pw_pairs, n_pw_feats);
    pw_stat        = nan(n_pw_pairs, n_pw_feats);
    pw_test_name   = cell(n_pw_pairs, n_pw_feats);
    pw_normal_A    = false(n_pw_pairs, n_pw_feats);
    pw_normal_B    = false(n_pw_pairs, n_pw_feats);
    pw_n_A         = nan(n_pw_pairs, n_pw_feats);
    pw_n_B         = nan(n_pw_pairs, n_pw_feats);
    pw_median_A    = nan(n_pw_pairs, n_pw_feats);
    pw_median_B    = nan(n_pw_pairs, n_pw_feats);

    pi_pw = 0;
    for ga = 1:n_groups
        for gb = ga+1:n_groups
            pi_pw = pi_pw + 1;
            pw_pair_idx(pi_pw,:)    = [ga, gb];
            pw_pair_labels{pi_pw}   = sprintf('%s vs %s', labels{ga}, labels{gb});

            for fi = 1:n_pw_feats
                fld = pw_feat_fields{fi};
                va  = groups(ga).(fld);
                vb  = groups(gb).(fld);
                va  = va(~isnan(va));
                vb  = vb(~isnan(vb));
                if isempty(va) || isempty(vb), continue; end

                % KS statistic D — computed as the max absolute difference
                % between the two empirical CDFs (same definition as kstest2).
                all_x  = sort([va(:); vb(:)]);
                cdf_a  = arrayfun(@(x) mean(va <= x), all_x);
                cdf_b  = arrayfun(@(x) mean(vb <= x), all_x);
                pw_KS_D(pi_pw, fi) = max(abs(cdf_a - cdf_b));

                % Q90 difference: 90th percentile of b minus a
                pw_Q90(pi_pw, fi)  = prctile(vb, 90) - prctile(va, 90);

                % Normality-aware two-sample hypothesis test
                res = two_sample_test(va, vb);
                pw_p(pi_pw, fi)         = res.p;
                pw_stat(pi_pw, fi)      = res.stat;
                pw_test_name{pi_pw, fi} = res.test_name;
                pw_normal_A(pi_pw, fi)  = res.normal_a;
                pw_normal_B(pi_pw, fi)  = res.normal_b;
                pw_n_A(pi_pw, fi)       = res.n_a;
                pw_n_B(pi_pw, fi)       = res.n_b;
                pw_median_A(pi_pw, fi)  = res.median_a;
                pw_median_B(pi_pw, fi)  = res.median_b;
            end
        end
    end

    % Benjamini-Hochberg FDR correction, applied per feature (i.e. within
    % the family of all group-pair comparisons run for that one feature —
    % the same convention as a post-hoc test after an omnibus comparison).
    pw_p_fdr = nan(size(pw_p));
    for fi = 1:n_pw_feats
        pw_p_fdr(:, fi) = bh_fdr(pw_p(:, fi));
    end

    pairwise.feat_names  = pw_feat_names;
    pairwise.pair_idx    = pw_pair_idx;
    pairwise.pair_labels = pw_pair_labels;
    pairwise.KS_D        = pw_KS_D;
    pairwise.Q90         = pw_Q90;
    pairwise.p           = pw_p;            % chosen test's uncorrected p-value
    pairwise.p_fdr       = pw_p_fdr;        % Benjamini-Hochberg FDR-adjusted p, per feature
    pairwise.stat        = pw_stat;         % chosen test's statistic (t or z/ranksum)
    pairwise.test_name   = pw_test_name;    % 't-test (Welch)' | 'Wilcoxon rank-sum'
    pairwise.normal_A    = pw_normal_A;     % Lilliefors: group A consistent with normal
    pairwise.normal_B    = pw_normal_B;
    pairwise.n_A         = pw_n_A;
    pairwise.n_B         = pw_n_B;
    pairwise.median_A    = pw_median_A;
    pairwise.median_B    = pw_median_B;
    pairwise.median_diff = abs(pw_median_A - pw_median_B);


    % =========================================================
    %  STEP 6d - Gompertz growth rate fitting and trade-off statistics
    %
    %  Fits S(t) = Amax / exp( exp( k*(tm - t) ) )  to each colony's
    %  smoothed size timecourse.  Amin = 0.
    %
    %  Maximum growth rate:  mu_max = (Amax / e) * k   [px/h]
    %
    %  Results stored per group in groups(g).growth:
    %    .mu_max          - per-colony maximum growth rate (NaN if fit failed)
    %    .Amax / .k / .tm - fitted Gompertz parameters
    %    .fit_ok          - logical: true where fit converged
    %    .local_density   - 2D bin-count density (40x40 bins, mu_max vs lag_time)
    %    .spearman_rho    - Spearman rho (mu_max vs lag_time)
    %    .spearman_p      - two-tailed p-value for Spearman rho
    %    .quantile_slope  - median quantile regression slope (lag ~ mu_max, tau=0.5)
    %    .quantile_intcpt - intercept of the median quantile regression
    %    .wt_spearman     - density-weighted Spearman rho
    %    .wt_qr_slope     - density-weighted quantile regression slope
    % =========================================================
    fprintf('\nSTEP 6d: Fitting Gompertz growth models...\n');

    gompertz_opts = optimoptions('lsqcurvefit', ...
        'Display',            'off', ...
        'MaxIterations',      500,   ...
        'FunctionTolerance',  1e-8);

    gompertz_fun = @(b, t) b(1) ./ exp(exp(b(2) .* (b(3) - t)));

    for g = 1:n_groups
        tc_mat  = groups(g).size_timecourse;   % [T x N]
        time_h  = groups(g).time_h;            % [T x 1]
        lag     = groups(g).lag_time;          % [N x 1]

        if isempty(tc_mat) || isempty(time_h)
            groups(g).growth = empty_growth_struct(length(lag));
            fprintf('  %s: no timecourse — skipped.\n', groups(g).label);
            continue;
        end

        [~, n_c] = size(tc_mat);
        mu_max_vec  = NaN(n_c, 1);
        Amax_vec    = NaN(n_c, 1);
        k_vec       = NaN(n_c, 1);
        tm_vec      = NaN(n_c, 1);
        fit_ok_vec  = false(n_c, 1);

        MIN_PTS = 6;

        for ci = 1:n_c
            s = tc_mat(:, ci);
            valid = s > 0 & isfinite(s);
            if sum(valid) < MIN_PTS, continue; end

            t_fit = time_h(valid);
            a_fit = s(valid);

            Amax0 = max(a_fit);
            k0    = 0.3;
            da    = diff(a_fit);
            [~, imax] = max(da);
            tm0   = t_fit(min(imax+1, length(t_fit)));

            lb = [0,   0, min(time_h)];
            ub = [Inf, 5, max(time_h)];

            try
                b = lsqcurvefit(gompertz_fun, [Amax0, k0, tm0], ...
                                t_fit, a_fit, lb, ub, gompertz_opts);
                Amax_vec(ci)   = b(1);
                k_vec(ci)      = b(2);
                tm_vec(ci)     = b(3);
                mu_max_vec(ci) = (b(1) / exp(1)) * b(2);
                fit_ok_vec(ci) = true;
            catch
                % fit failed — leave NaN
            end
        end

        ok      = fit_ok_vec & isfinite(mu_max_vec) & isfinite(lag);
        mu_ok   = mu_max_vec(ok);
        lag_ok  = lag(ok);
        n_ok    = sum(ok);

        fprintf('  %s: %d/%d fits converged.\n', groups(g).label, n_ok, n_c);

        % --- 2D local density (40×40 bin counts) — vectorised --------------
        % Uses histcounts to bin all colonies at once (O(N) not O(N^2)).
        local_density = NaN(n_c, 1);
        if n_ok >= 2
            n_bins    = 40;
            mu_edges  = linspace(min(mu_ok),  max(mu_ok)  + eps, n_bins+1);
            lag_edges = linspace(min(lag_ok), max(lag_ok) + eps, n_bins+1);

            % Assign each colony to a bin index (clamp to [1, n_bins])
            bm = max(1, min(n_bins, discretize(mu_ok,  mu_edges)));
            bl = max(1, min(n_bins, discretize(lag_ok, lag_edges)));

            % Build 2D count matrix and look up count for each colony
            cnt_mat = accumarray([bl, bm], 1, [n_bins, n_bins]);
            lin_idx = sub2ind([n_bins, n_bins], bl, bm);
            local_density_ok = max(1, cnt_mat(lin_idx));

            local_density(ok) = local_density_ok;
        end

        % --- Spearman correlation (mu_max vs lag_time) -------------------
        sp_rho = NaN;  sp_p = NaN;
        if n_ok >= 3
            [sp_rho, sp_p] = corr(mu_ok, lag_ok, 'Type', 'Spearman');
        end

        % --- Median quantile regression: lag ~ mu_max, tau = 0.5 --------
        qr_slope = NaN;  qr_intcpt = NaN;
        if n_ok >= 4
            [qr_slope, qr_intcpt] = quantile_regression(lag_ok, mu_ok, 0.5, []);
        end

        % --- Density-weighted Spearman and QR ----------------------------
        wt_sp      = NaN;
        wt_qr_slope = NaN;
        if n_ok >= 3 && any(local_density_ok > 0)
            w      = local_density_ok(:);
            w_norm = w / sum(w);
            % Weighted rank correlation approximation via rank transform
            [~, rmu]  = sort(mu_ok);   rank_mu  = zeros(n_ok,1); rank_mu(rmu)  = (1:n_ok)';
            [~, rlag] = sort(lag_ok);  rank_lag = zeros(n_ok,1); rank_lag(rlag) = (1:n_ok)';
            wt_sp = sum(w_norm .* rank_mu .* rank_lag) / ...
                    sqrt(sum(w_norm .* rank_mu.^2) * sum(w_norm .* rank_lag.^2) + eps);
        end
        if n_ok >= 4 && any(local_density_ok > 0)
            w_vec = local_density_ok(:);
            [wt_qr_slope, ~] = quantile_regression(lag_ok, mu_ok, 0.5, w_vec);
        end

        groups(g).growth.mu_max         = mu_max_vec;
        groups(g).growth.Amax           = Amax_vec;
        groups(g).growth.k              = k_vec;
        groups(g).growth.tm             = tm_vec;
        groups(g).growth.fit_ok         = fit_ok_vec;
        groups(g).growth.local_density  = local_density;
        groups(g).growth.spearman_rho   = sp_rho;
        groups(g).growth.spearman_p     = sp_p;
        groups(g).growth.quantile_slope = qr_slope;
        groups(g).growth.quantile_intcpt= qr_intcpt;
        groups(g).growth.wt_spearman    = wt_sp;
        groups(g).growth.wt_qr_slope    = wt_qr_slope;
    end


    % =========================================================
    %  STEP 6e - Lag time vs morphology: Pearson + Spearman correlations
    %            and linear regression per feature per group.
    %
    %  For each group and each morphological feature, computes:
    %    .r          - Pearson r (lag_time vs feature)
    %    .r_sq       - R^2 = r^2
    %    .p_pearson  - two-tailed p-value for Pearson r
    %    .rho        - Spearman rho
    %    .p_spearman - two-tailed p-value for Spearman rho
    %    .slope      - linear regression slope (feature ~ lag_time)
    %    .intercept  - linear regression intercept
    %    .eq_str     - equation string  "y = m*x + b"
    %
    %  Stored in groups(g).morph_corr.(feature_field)
    % =========================================================
    fprintf('\nSTEP 6e: Computing lag time vs morphology correlations...\n');

    morph_fields = {'size','area','intensity','mean_intensity', ...
                    'int_per_size','perimeter','circularity', ...
                    'eccentricity','solidity'};

    for g = 1:n_groups
        lag = groups(g).lag_time;   % [N x 1], NaN for censored colonies
        groups(g).morph_corr = struct();

        for fi = 1:length(morph_fields)
            fld  = morph_fields{fi};
            feat = groups(g).(fld);

            % Default: all NaN
            mc = struct('r',NaN,'r_sq',NaN,'p_pearson',NaN, ...
                        'rho',NaN,'p_spearman',NaN, ...
                        'slope',NaN,'intercept',NaN,'eq_str','');

            if isempty(feat) || isempty(lag)
                groups(g).morph_corr.(fld) = mc;
                continue;
            end

            ok = isfinite(lag) & isfinite(feat);
            x  = lag(ok);
            y  = feat(ok);
            n_ok = sum(ok);

            if n_ok < 4
                groups(g).morph_corr.(fld) = mc;
                continue;
            end

            % Pearson r and p-value
            [mc.r, mc.p_pearson] = corr(x, y, 'Type','Pearson');
            mc.r_sq = mc.r^2;

            % Spearman rho and p-value
            [mc.rho, mc.p_spearman] = corr(x, y, 'Type','Spearman');

            % Linear regression: y = slope*x + intercept
            X = [ones(n_ok,1), x];
            b = X \ y;
            mc.intercept = b(1);
            mc.slope     = b(2);

            % Equation string (sign-aware)
            if mc.intercept >= 0
                mc.eq_str = sprintf('y = %.3g x + %.3g', mc.slope, mc.intercept);
            else
                mc.eq_str = sprintf('y = %.3g x - %.3g', mc.slope, abs(mc.intercept));
            end

            groups(g).morph_corr.(fld) = mc;
        end

        fprintf('  %s: morph correlations done (n=%d valid colonies)\n', ...
            groups(g).label, sum(isfinite(lag)));
    end

    % =========================================================
    %  STEP 6f - Biological Replicate Correlation & Statistical Analysis
    %
    %  For each sample group, identifies every unique biological replicate
    %  (using the bio_rep_col chosen at parameter step 6), then pools ALL
    %  colonies from ALL plates belonging to the same bio-rep label — so
    %  technical replicates are merged, not compared.  Statistics are then
    %  computed between the resulting bio-rep colony pools.
    %
    %  Tests computed per replicate-pair per feature:
    %
    %    Pearson r           - Linear correlation between the 20-quantile profiles
    %                          of two replicates. Parametric complement to Spearman;
    %                          sensitive to linear shifts in the distribution shape.
    %                          r = 1 means perfect linear agreement.
    %
    %    Pearson p-value     - Two-tailed t-test p-value for Pearson r (df = 18).
    %
    %    Spearman rho        - Rank correlation of the same 20-quantile profiles.
    %                          Non-parametric; robust to skewed colony distributions.
    %                          rho = 1 means perfect rank agreement.
    %
    %    Spearman p-value    - Two-tailed p-value for Spearman rho (df = 18).
    %
    %    KS D statistic      - Kolmogorov-Smirnov max CDF distance between
    %                          the two replicates (0 = identical distribution,
    %                          1 = completely non-overlapping). Captures
    %                          distributional differences beyond central tendency.
    %
    %    KS p-value          - Two-sample KS p-value (H0: same distribution).
    %                          p < 0.05 indicates significant distributional
    %                          difference between replicates.
    %
    %    Jensen-Shannon Div  - Symmetric, bounded [0, 1] divergence between
    %                          the normalised histograms of two replicates.
    %                          JSD = 0 means identical distributions; JSD = 1
    %                          means maximally different.
    %
    %    |Delta Gini|        - Absolute difference in Gini indices between
    %                          the two replicates. A small |ΔGini| means both
    %                          replicates have similar within-replicate colony
    %                          heterogeneity levels.
    %
    %    Median difference   - Absolute difference in medians between replicates.
    %                          Simple, interpretable effect size.
    %
    %    CV of medians (%)   - Coefficient of variation of the per-replicate
    %                          medians within a group. Measures how consistent
    %                          the central tendency is across all replicates.
    %                          Computed at the group level (not pair level).
    %
    %  Per-replicate (not per-pair) quantities — stored as [n_reps x 10]:
    %
    %    Gini index          - Classic Gini coefficient on the raw colony values
    %                          per replicate. Measures within-replicate phenotypic
    %                          heterogeneity at the single-colony level.
    %                          0 = all colonies identical, 1 = maximum inequality.
    %
    %  Results stored in:
    %    ht.biorep_corr(g)   - per-group struct
    %      .group_label      - label string
    %      .plate_labels     - {n_reps x 1} unique biological replicate label strings
    %      .rep_n_plates     - [n_reps x 1] number of plates pooled per bio-rep
    %      .n_reps           - number of unique biological replicates in this group
    %      .features         - {1 x 10} feature name strings
    %      .pair_labels      - {n_pairs x 1} 'RepA vs RepB' strings
    %      .pearson_r        - [n_pairs x 10] Pearson r on quantile profiles
    %      .pearson_p        - [n_pairs x 10] Pearson p-value (two-tailed, df=18)
    %      .spearman_rho     - [n_pairs x 10] Spearman rho on quantile profiles
    %      .spearman_p       - [n_pairs x 10] Spearman p-value (two-tailed, df=18)
    %      .ks_D             - [n_pairs x 10] KS D statistic
    %      .ks_p             - [n_pairs x 10] KS p-value
    %      .jsd              - [n_pairs x 10] Jensen-Shannon divergence
    %      .delta_gini       - [n_pairs x 10] |Gini_A - Gini_B| on raw values
    %      .median_diff      - [n_pairs x 10] |median(A) - median(B)|
    %      .gini_values      - [n_reps  x 10] classic Gini on raw colony values
    %      .cv_medians       - [1 x 10] CV% of medians across all replicates
    % =========================================================
    fprintf('\nSTEP 6f: Biological replicate correlation analysis...\n');
    fprintf('  Strategy: plates are pooled by biological replicate label.\n');
    fprintf('  All colonies from all technical replicate plates sharing the\n');
    fprintf('  same bio-rep label are concatenated into one pool before any\n');
    fprintf('  statistics are computed.  Technical replicates are NOT compared\n');
    fprintf('  against each other.\n\n');

    biorep_feat_fields = {'lag_time','size','area','intensity', ...
                          'mean_intensity','int_per_size','perimeter', ...
                          'circularity','eccentricity','solidity'};
    biorep_feat_labels = {'Lag Time','Colony Size','Area','Intensity', ...
                          'Mean Intensity','Int/Size','Perimeter', ...
                          'Circularity','Eccentricity','Solidity'};
    n_brf = length(biorep_feat_fields);

    biorep_corr(n_groups) = struct();

    % Verify bio_rep_col is usable — abort with clear message if not.
    % A silent plate-level fallback would produce meaningless statistics
    % (e.g. 276 pairs for 24 plates) so we stop here instead.
    has_bio_col = isfield(p, 'bio_rep_col') && ~isempty(p.bio_rep_col) && ...
                  ~isempty(data.metadata.original) && ...
                  ismember(p.bio_rep_col, data.metadata.original.Properties.VariableNames);

    if ~has_bio_col
        if ~isfield(p, 'bio_rep_col') || isempty(p.bio_rep_col)
            fprintf('\n  STEP 6f SKIPPED: no biological replicate column was selected\n');
            fprintf('  (you pressed 0 at parameter step 6).\n');
            fprintf('  To run bio-rep statistics, re-run preprocess_pipeline_data\n');
            fprintf('  and select the metadata column that identifies each biological\n');
            fprintf('  replicate (e.g. "Set", "Strain", "BioRep").\n\n');
        else
            fprintf('\n  STEP 6f SKIPPED: bio_rep_col "%s" not found in metadata table.\n', p.bio_rep_col);
            fprintf('  Available columns: %s\n', ...
                strjoin(data.metadata.original.Properties.VariableNames, ', '));
            fprintf('  Re-run preprocess_pipeline_data and select the correct column.\n\n');
        end
        % Store an empty marker so plot_combined_samples prints a helpful message
        for g_skip = 1:n_groups
            biorep_corr(g_skip).group_label  = labels{g_skip};
            biorep_corr(g_skip).n_reps       = 0;
            biorep_corr(g_skip).features     = {};
            biorep_corr(g_skip).plate_labels = {};
            biorep_corr(g_skip).pair_labels  = {};
            biorep_corr(g_skip).skipped      = true;
            biorep_corr(g_skip).skip_reason  = 'bio_rep_col not set or not found in metadata';
        end
        fprintf('STEP 6f skipped.\n\n');
    else

    for g = 1:n_groups
        plate_ix = ix{g};   % all QC-passing plate indices in this sample group
        n_plates = length(plate_ix);

        biorep_corr(g).group_label = labels{g};
        biorep_corr(g).features    = biorep_feat_labels;

        % ------------------------------------------------------------------
        %  DISCOVER UNIQUE BIOLOGICAL REPLICATES within this sample group.
        %
        %  If a bio_rep_col was supplied by the user (e.g. "Strain", "Set",
        %  "BioRep"), we read it from the metadata table to assign each plate
        %  to a bio-rep label.  All plates sharing the same bio-rep label —
        %  regardless of how many technical replicates they represent — are
        %  pooled together.
        %
        %  If no bio_rep_col was supplied (user pressed 0 at step 6), we fall
        %  back to treating each plate as its own bio-rep.  In that case the
        %  user will see a warning and the comparison is plate-level, not
        %  bio-rep-level.
        % ------------------------------------------------------------------
        plate_bio_labels = cell(n_plates, 1);   % bio-rep label for each plate

        meta_tbl = data.metadata.original;
        for k = 1:n_plates
            pi_k = plate_ix(k);
            if pi_k <= height(meta_tbl)
                plate_bio_labels{k} = val2str(meta_tbl.(p.bio_rep_col)(pi_k));
            else
                plate_bio_labels{k} = sprintf('Plate%d', pi_k);
            end
        end

        % Diagnostic: show first 4 plate->bio-rep assignments so the
        % user can immediately verify the correct column is being used
        if g == 1
            fprintf('  [DIAG] Bio-rep column "%s", first 4 plate assignments:\n', p.bio_rep_col);
            for k_diag = 1:min(4, n_plates)
                fprintf('         Plate %2d  ->  "%s"\n', plate_ix(k_diag), plate_bio_labels{k_diag});
            end
        end

        % Unique bio-rep labels, preserving first-occurrence order
        [unique_bio_labels, first_occ] = unique(plate_bio_labels, 'stable');
        n_bioreps = length(unique_bio_labels);

        biorep_corr(g).n_reps      = n_bioreps;
        biorep_corr(g).plate_labels = unique_bio_labels;   % one entry per bio-rep

        fprintf('  Group %s: %d plates -> %d biological replicate(s): %s\n', ...
            labels{g}, n_plates, n_bioreps, strjoin(unique_bio_labels, ' | '));

        % ------------------------------------------------------------------
        %  POOL COLONIES PER BIOLOGICAL REPLICATE.
        %
        %  rep_data{r}{fi} = [N_pooled x 1] cleaned vector of all colony
        %  values for bio-rep r, feature fi, concatenated across every plate
        %  that belongs to that bio-rep (technical replicates included).
        %  rep_n_plates(r) = how many plates were pooled into bio-rep r.
        % ------------------------------------------------------------------
        rep_data    = cell(n_bioreps, 1);
        rep_n_plates = zeros(n_bioreps, 1);

        for r = 1:n_bioreps
            rep_data{r} = cell(1, n_brf);
            for fi = 1:n_brf
                rep_data{r}{fi} = [];   % will grow by concatenation below
            end
        end

        for k = 1:n_plates
            pi_k    = plate_ix(k);
            bio_lbl = plate_bio_labels{k};

            % Map this plate to its bio-rep index
            r_idx = find(strcmp(unique_bio_labels, bio_lbl), 1);
            if isempty(r_idx), continue; end
            rep_n_plates(r_idx) = rep_n_plates(r_idx) + 1;

            % ---- Extract colony data from this plate (same logic as STEP 4) ----
            if ~isfield(data.processed{pi_k}, 'colonies') || ...
               ~isfield(data.processed{pi_k}.colonies, 'new')
                continue;
            end

            c_new_k = data.processed{pi_k}.colonies.new;
            n_col_k = length(c_new_k.lag_time);
            if n_col_k == 0, continue; end

            col_k = data.processed{pi_k}.colonies;

            % Resolve rp_clean (same logic as STEP 4)
            rp_k = [];
            if isfield(col_k, 'region_props_clean')
                rp_k = col_k.region_props_clean;
            elseif isfield(col_k, 'region_props')
                rp_all_k = col_k.region_props;
                if isfield(col_k, 'flag_colony_ok') && ...
                   length(col_k.flag_colony_ok) == height(rp_all_k)
                    cidx_k = find(col_k.flag_colony_ok);
                    if length(cidx_k) == n_col_k
                        rp_k = rp_all_k(cidx_k, :);
                    end
                end
            end

            % Build feature vectors for this plate
            lt_k = c_new_k.lag_time(:) + p.incTime;
            cens_tol_k = 0.5 * p.lag_step;
            lt_k(lt_k >= p.max_lag - cens_tol_k) = NaN;

            sz_k = nan(n_col_k,1);  ar_k = nan(n_col_k,1);
            ec_k = nan(n_col_k,1);  pe_k = nan(n_col_k,1);
            ci_k = nan(n_col_k,1);  so_k = nan(n_col_k,1);
            in_k = nan(n_col_k,1);  mi_k = nan(n_col_k,1);
            ip_k = nan(n_col_k,1);

            if ~isempty(rp_k) && height(rp_k) == n_col_k
                sz_k = rp_k.Area;         ar_k = rp_k.Area;
                ec_k = rp_k.Eccentricity; pe_k = rp_k.Perimeter;
                ci_k = rp_k.Circularity;  so_k = rp_k.Solidity;
            end

            if isfield(col_k, 'timecourse_intensity_smoothed') && ...
               ~isempty(col_k.timecourse_intensity_smoothed) && ...
               size(col_k.timecourse_intensity_smoothed, 2) == n_col_k
                in_k = col_k.timecourse_intensity_smoothed(end, :)';
            end

            if ~isempty(rp_k) && height(rp_k) == n_col_k
                if isfield(rp_k, 'MeanIntensity')
                    mi_k = rp_k.MeanIntensity;
                elseif isfield(rp_k, 'mean_intensity')
                    mi_k = rp_k.mean_intensity;
                end
            end
            nan_mi_k = isnan(mi_k);
            if any(nan_mi_k) && any(~isnan(in_k)) && any(~isnan(sz_k))
                valid_mi = ~isnan(in_k) & ~isnan(sz_k) & sz_k > 0 & nan_mi_k;
                mi_k(valid_mi) = in_k(valid_mi) ./ sz_k(valid_mi);
            end

            valid_ip_k = ~isnan(in_k) & ~isnan(sz_k) & sz_k > 0;
            ip_k(valid_ip_k) = in_k(valid_ip_k) ./ sz_k(valid_ip_k);

            % Size threshold (same as STEP 4)
            fail_k = ~isnan(sz_k) & (sz_k < p.size_threshold);
            lt_k(fail_k) = NaN;  sz_k(fail_k) = NaN;  ar_k(fail_k) = NaN;
            in_k(fail_k) = NaN;  mi_k(fail_k) = NaN;  ip_k(fail_k) = NaN;
            pe_k(fail_k) = NaN;  ci_k(fail_k) = NaN;  ec_k(fail_k) = NaN;
            so_k(fail_k) = NaN;

            % Lag-censoring exclusion (same rule as STEP 4):
            % NaN all morphology for colonies with undetectable lag time.
            cens_k = isnan(lt_k);
            sz_k(cens_k) = NaN;  ar_k(cens_k) = NaN;
            in_k(cens_k) = NaN;  mi_k(cens_k) = NaN;  ip_k(cens_k) = NaN;
            pe_k(cens_k) = NaN;  ci_k(cens_k) = NaN;
            ec_k(cens_k) = NaN;  so_k(cens_k) = NaN;

            all_vecs_k = {lt_k, sz_k, ar_k, in_k, mi_k, ip_k, pe_k, ci_k, ec_k, so_k};

            % Concatenate into the bio-rep pool for this group
            for fi = 1:n_brf
                v = all_vecs_k{fi};
                rep_data{r_idx}{fi} = [rep_data{r_idx}{fi}; v(~isnan(v))];
            end
        end % plates loop

        % Log pooling summary
        for r = 1:n_bioreps
            n_col_pool = length(rep_data{r}{1});   % use feature 1 as proxy count
            fprintf('    Bio-rep %-20s : %d plate(s) pooled, ~%d colonies\n', ...
                unique_bio_labels{r}, rep_n_plates(r), n_col_pool);
        end
        fprintf('\n');

        biorep_corr(g).plate_labels  = unique_bio_labels;
        biorep_corr(g).rep_n_plates  = rep_n_plates;   % [n_bioreps x 1] plate counts

        % Enumerate all unique bio-rep pairs (i < j)
        n_pairs = n_bioreps * (n_bioreps - 1) / 2;
        pair_labels  = cell(n_pairs, 1);
        pearson_r_mat  = nan(n_pairs, n_brf);
        pearson_p_mat  = nan(n_pairs, n_brf);
        sp_rho_mat     = nan(n_pairs, n_brf);
        sp_p_mat       = nan(n_pairs, n_brf);
        ks_D_mat       = nan(n_pairs, n_brf);
        ks_p_mat       = nan(n_pairs, n_brf);
        jsd_mat        = nan(n_pairs, n_brf);
        delta_gini_mat = nan(n_pairs, n_brf);
        meddiff_mat    = nan(n_pairs, n_brf);

        pi_br = 0;
        for ra = 1:n_bioreps
            for rb = ra+1:n_bioreps
                pi_br = pi_br + 1;
                pair_labels{pi_br} = sprintf('%s vs %s', unique_bio_labels{ra}, unique_bio_labels{rb});

                for fi = 1:n_brf
                    va = rep_data{ra}{fi};
                    vb = rep_data{rb}{fi};

                    if length(va) < 4 || length(vb) < 4
                        continue;
                    end

                    % --- Quantile profile (shared basis for Pearson & Spearman) ---
                    % Both replicates may have different N, so we compare their
                    % distributional shape via 20 evenly spaced percentiles
                    % (5th through 95th). This N-invariant fingerprint is used
                    % for both Pearson r and Spearman rho.
                    q_probs = linspace(5, 95, 20);   % 20 percentile points
                    qa = prctile(va, q_probs);        % [1 x 20] quantile profile A
                    qb = prctile(vb, q_probs);        % [1 x 20] quantile profile B

                    if std(qa) > 0 && std(qb) > 0
                        % --- Pearson r (linear correlation of quantile profiles) ---
                        [pearson_r_mat(pi_br, fi), pearson_p_mat(pi_br, fi)] = ...
                            corr(qa(:), qb(:), 'Type', 'Pearson');

                        % --- Spearman rho (rank correlation of quantile profiles) ---
                        [sp_rho_mat(pi_br, fi), sp_p_mat(pi_br, fi)] = ...
                            corr(qa(:), qb(:), 'Type', 'Spearman');
                    end

                    % --- KS test -----------------------------------------
                    all_x   = sort([va(:); vb(:)]);
                    cdf_a   = arrayfun(@(x) mean(va <= x), all_x);
                    cdf_b   = arrayfun(@(x) mean(vb <= x), all_x);
                    ks_D_mat(pi_br, fi) = max(abs(cdf_a - cdf_b));

                    % Approximate KS p-value (Kolmogorov distribution series)
                    n_a = length(va);  n_b = length(vb);
                    n_eff = sqrt(n_a * n_b / (n_a + n_b));
                    lambda = (n_eff + 0.12 + 0.11/n_eff) * ks_D_mat(pi_br, fi);
                    j_vec = (1:100)';
                    ks_p_mat(pi_br, fi) = 2 * sum((-1).^(j_vec-1) .* exp(-2 * lambda^2 * j_vec.^2));
                    ks_p_mat(pi_br, fi) = max(0, min(1, ks_p_mat(pi_br, fi)));

                    % --- Jensen-Shannon Divergence -----------------------
                    % Computed on shared 50-bin edges spanning the union
                    % range of both replicates.
                    lo = min(min(va), min(vb));
                    hi = max(max(va), max(vb));
                    if hi > lo
                        edges_jsd = linspace(lo, hi, 51);
                        ca = histcounts(va, edges_jsd);
                        cb = histcounts(vb, edges_jsd);
                        pa = ca / sum(ca);
                        pb = cb / sum(cb);
                        pm = 0.5 * (pa + pb);
                        eps_jsd = 1e-12;
                        kl_a = sum(pa(pa > 0) .* log2(pa(pa > 0) ./ (pm(pa > 0) + eps_jsd)));
                        kl_b = sum(pb(pb > 0) .* log2(pb(pb > 0) ./ (pm(pb > 0) + eps_jsd)));
                        jsd_mat(pi_br, fi) = max(0, 0.5 * (kl_a + kl_b));

                        % |ΔGini| on raw colony values (not PDF bins)
                        delta_gini_mat(pi_br, fi) = abs( ...
                            compute_gini(va(:)) - compute_gini(vb(:)));
                    end

                    % --- Absolute median difference ----------------------
                    meddiff_mat(pi_br, fi) = abs(median(va) - median(vb));
                end % features loop
            end % rb loop
        end % ra loop

        % --- Per-replicate Gini indices (classic Gini on raw colony values) ---
        %
        % gini_values(r, fi) — Gini coefficient on the raw colony-level values.
        % Measures within-replicate phenotypic inequality: 0 = all colonies
        % identical, 1 = maximum inequality (one colony dominates).
        gini_values_mat = nan(n_bioreps, n_brf);

        for r = 1:n_bioreps
            for fi = 1:n_brf
                v = rep_data{r}{fi};
                if length(v) < 2, continue; end
                gini_values_mat(r, fi) = compute_gini(v(:));
            end
        end

        % --- CV of medians across all replicates (group-level summary) ---
        cv_med = nan(1, n_brf);
        if n_bioreps >= 2
            for fi = 1:n_brf
                rep_medians = nan(n_bioreps, 1);
                for r = 1:n_bioreps
                    v = rep_data{r}{fi};
                    if ~isempty(v)
                        rep_medians(r) = median(v);
                    end
                end
                valid_m = rep_medians(~isnan(rep_medians));
                if length(valid_m) >= 2 && mean(valid_m) ~= 0
                    cv_med(fi) = 100 * std(valid_m) / abs(mean(valid_m));
                elseif length(valid_m) >= 2
                    cv_med(fi) = 0;
                end
            end
        end

        biorep_corr(g).pair_labels    = pair_labels;
        biorep_corr(g).pearson_r      = pearson_r_mat;
        biorep_corr(g).pearson_p      = pearson_p_mat;
        biorep_corr(g).spearman_rho   = sp_rho_mat;
        biorep_corr(g).spearman_p     = sp_p_mat;
        biorep_corr(g).ks_D           = ks_D_mat;
        biorep_corr(g).ks_p           = ks_p_mat;
        biorep_corr(g).jsd            = jsd_mat;
        biorep_corr(g).delta_gini     = delta_gini_mat;
        biorep_corr(g).median_diff    = meddiff_mat;
        biorep_corr(g).gini_values    = gini_values_mat;
        biorep_corr(g).cv_medians     = cv_med;

        % Console summary
        fprintf('  Group %s: %d bio-replicate(s), %d pair(s)\n', labels{g}, n_bioreps, n_pairs);

        % Per-replicate Gini summary
        fprintf('    Per-bio-rep Gini (classic, raw colony values):\n');
        for r = 1:n_bioreps
            gini_strs = '';
            for fi = 1:n_brf
                gv = gini_values_mat(r, fi);
                if ~isnan(gv)
                    gini_strs = [gini_strs, sprintf('%s: %.3f  ', ...
                        biorep_feat_labels{fi}, gv)]; %#ok<AGROW>
                end
            end
            fprintf('      [%s]  %s\n', unique_bio_labels{r}, gini_strs);
        end

        if n_pairs > 0 && n_bioreps >= 2
            for pi_pr = 1:n_pairs
                fprintf('    %s\n', pair_labels{pi_pr});
                for fi = 1:n_brf
                    pr_v  = pearson_r_mat(pi_pr, fi);
                    rho_v = sp_rho_mat(pi_pr, fi);
                    ks_v  = ks_D_mat(pi_pr, fi);
                    jsd_v = jsd_mat(pi_pr, fi);
                    dg_v  = delta_gini_mat(pi_pr, fi);
                    if ~isnan(pr_v)
                        fprintf('      %-20s  Pearson_r=%+.3f  Spearman_rho=%+.3f  KS_D=%.3f  JSD=%.3f  |dGini|=%.3f\n', ...
                            biorep_feat_labels{fi}, pr_v, rho_v, ks_v, jsd_v, dg_v);
                    end
                end
            end
        else
            fprintf('    (only 1 biological replicate — no pairwise comparison possible)\n');
        end
    end % groups loop

    fprintf('STEP 6f complete.\n\n');
    end  % if has_bio_col ... else ... end

    % =========================================================
    %  STEP 6i - Per-Biological-Replicate, Cross-Group Hypothesis Tests
    %
    %  STEP 6f above compares REPLICATES against each other WITHIN one
    %  sample group (e.g. Rep1 vs Rep2, both "7h"). This step does the
    %  opposite: for EACH biological replicate found in the metadata
    %  (bio_rep_col), it compares SAMPLE GROUPS against each other using
    %  ONLY that replicate's own colonies — e.g. for "Rep1": 3h vs 7h,
    %  3h vs 24h, 7h vs 24h, ...; then the same set of comparisons again
    %  for "Rep2", independently. This directly answers the reviewer
    %  question "are colonies from the 7h sample really larger than
    %  <other time>?" on a per-replicate basis rather than pooling every
    %  replicate together.
    %
    %  Uses the SAME normality-aware test as STEP 6/ht.pairwise above
    %  (Welch t-test if both groups' raw values pass a Lilliefors
    %  normality test, else Wilcoxon rank-sum) on RAW per-colony values
    %  (not the binned probability-density histogram).
    %
    %  Colonies are pooled per (replicate, sample group) using the
    %  plate_seg_start/plate_seg_end boundaries recorded during STEP 4 —
    %  this correctly excludes colonies NaN-masked by the final-size
    %  filter / lag-censoring, unlike naively slicing by raw plate
    %  colony count.
    %
    %  Results stored in ht.biorep_group_tests(r):
    %    .rep_label       - biological replicate label string
    %    .groups_present  - {1 x n_gp} sample-group labels present in this replicate
    %    .features        - {1 x 10} feature display names (matches pw_feat_names)
    %    .hist            - [n_gp x 10] cell of normalised histogram counts
    %                       (same bin edges as groups(g).hist.*), for plotting
    %    .n               - [n_gp x 10] valid (non-NaN) colony count per group/feature
    %    .median          - [n_gp x 10] median per group/feature
    %    .gini            - [n_gp x 10] Gini index per group/feature
    %    .pair_labels     - {n_pairs x 1} 'GroupA vs GroupB' strings
    %    .p / .p_fdr / .stat / .test_name / .normal_A / .normal_B / .n_A / .n_B /
    %    .median_A / .median_B / .median_diff  — same meaning as ht.pairwise,
    %                       but computed on this replicate's colonies only.
    %                       .p_fdr is Benjamini-Hochberg corrected per feature
    %                       across only this replicate's own pairs (STEP 6i
    %                       replicates are NOT pooled into one correction family).
    %    .skipped         - true if bio_rep_col was not set/found (see ht.biorep_corr)
    % =========================================================
    fprintf('\nSTEP 6i: Per-biological-replicate, cross-group hypothesis tests...\n');

    if ~has_bio_col
        biorep_group_tests = struct( ...
            'rep_label', {}, 'groups_present', {}, 'features', {}, ...
            'hist', {}, 'n', {}, 'median', {}, 'gini', {}, ...
            'pair_labels', {}, 'p', {}, 'p_fdr', {}, 'stat', {}, 'test_name', {}, ...
            'normal_A', {}, 'normal_B', {}, 'n_A', {}, 'n_B', {}, ...
            'median_A', {}, 'median_B', {}, 'median_diff', {}, ...
            'skipped', {}, 'skip_reason', {});
        biorep_group_tests(1).skipped     = true;
        biorep_group_tests(1).skip_reason = 'bio_rep_col not set or not found in metadata';
        fprintf('  STEP 6i SKIPPED: %s\n\n', biorep_group_tests(1).skip_reason);
    else
        meta_tbl = data.metadata.original;

        % Per-group, per-plate bio-rep label (aligned with groups(g).plate_indices)
        group_plate_biolbl = cell(n_groups, 1);
        for g = 1:n_groups
            n_pl_g = length(groups(g).plate_indices);
            lbls_g = cell(n_pl_g, 1);
            for k = 1:n_pl_g
                pi_k = groups(g).plate_indices(k);
                if pi_k <= height(meta_tbl)
                    lbls_g{k} = val2str(meta_tbl.(p.bio_rep_col)(pi_k));
                else
                    lbls_g{k} = sprintf('Plate%d', pi_k);
                end
            end
            group_plate_biolbl{g} = lbls_g;
        end

        all_bio_labels_global = unique(vertcat(group_plate_biolbl{:}), 'stable');
        n_bioreps_global = length(all_bio_labels_global);

        fprintf('  Biological replicates found (%d): %s\n', ...
            n_bioreps_global, strjoin(all_bio_labels_global, ', '));

        biorep_group_tests(n_bioreps_global) = struct();

        for r = 1:n_bioreps_global
            rep_lbl = all_bio_labels_global{r};

            % Pool raw per-colony vectors for this replicate, kept separate
            % per sample group (unlike STEP 6f, which pools across groups).
            pooled          = cell(n_groups, n_pw_feats);
            group_has_data  = false(n_groups, 1);

            for g = 1:n_groups
                match_k = find(strcmp(group_plate_biolbl{g}, rep_lbl));
                if isempty(match_k), continue; end
                group_has_data(g) = true;

                for fi = 1:n_pw_feats
                    acc = [];
                    for kk = 1:length(match_k)
                        k  = match_k(kk);
                        s0 = groups(g).plate_seg_start(k);
                        e0 = groups(g).plate_seg_end(k);
                        if e0 < s0, continue; end   % plate contributed 0 colonies
                        seg = groups(g).(pw_feat_fields{fi})(s0:e0);
                        acc = [acc; seg(~isnan(seg))]; %#ok<AGROW>
                    end
                    pooled{g, fi} = acc;
                end
            end

            groups_present_idx = find(group_has_data);
            n_gp = length(groups_present_idx);

            biorep_group_tests(r).rep_label      = rep_lbl;
            biorep_group_tests(r).groups_present = labels(groups_present_idx);
            biorep_group_tests(r).features       = pw_feat_names;
            biorep_group_tests(r).skipped        = false;

            % Per-group/feature descriptives + histogram (for plotting)
            hist_c   = cell(n_gp, n_pw_feats);
            n_mat    = nan(n_gp, n_pw_feats);
            med_mat  = nan(n_gp, n_pw_feats);
            gini_mat = nan(n_gp, n_pw_feats);
            for a = 1:n_gp
                g = groups_present_idx(a);
                for fi = 1:n_pw_feats
                    v = pooled{g, fi};
                    n_mat(a, fi) = length(v);
                    if isempty(v), continue; end
                    med_mat(a, fi)  = median(v);
                    gini_mat(a, fi) = compute_gini(v);
                    bin_edges = bin_map.(feature_names_h{fi});
                    counts = histcounts(v, bin_edges);
                    total  = sum(counts);
                    if total > 0
                        hist_c{a, fi} = counts / total;
                    else
                        hist_c{a, fi} = counts;
                    end
                end
            end
            biorep_group_tests(r).hist   = hist_c;
            biorep_group_tests(r).n      = n_mat;
            biorep_group_tests(r).median = med_mat;
            biorep_group_tests(r).gini   = gini_mat;

            % Pairwise hypothesis tests among groups present in this replicate
            n_pairs_r = n_gp * (n_gp - 1) / 2;
            pair_labels_r = cell(n_pairs_r, 1);
            p_mat    = nan(n_pairs_r, n_pw_feats);
            stat_mat = nan(n_pairs_r, n_pw_feats);
            tn_mat   = cell(n_pairs_r, n_pw_feats);
            nrmA_mat = false(n_pairs_r, n_pw_feats);
            nrmB_mat = false(n_pairs_r, n_pw_feats);
            nA_mat   = nan(n_pairs_r, n_pw_feats);
            nB_mat   = nan(n_pairs_r, n_pw_feats);
            medA_mat = nan(n_pairs_r, n_pw_feats);
            medB_mat = nan(n_pairs_r, n_pw_feats);

            pi_r = 0;
            for a = 1:n_gp
                for b = a+1:n_gp
                    pi_r = pi_r + 1;
                    ga = groups_present_idx(a);
                    gb = groups_present_idx(b);
                    pair_labels_r{pi_r} = sprintf('%s vs %s', labels{ga}, labels{gb});

                    for fi = 1:n_pw_feats
                        res = two_sample_test(pooled{ga, fi}, pooled{gb, fi});
                        p_mat(pi_r, fi)    = res.p;
                        stat_mat(pi_r, fi) = res.stat;
                        tn_mat{pi_r, fi}   = res.test_name;
                        nrmA_mat(pi_r, fi) = res.normal_a;
                        nrmB_mat(pi_r, fi) = res.normal_b;
                        nA_mat(pi_r, fi)   = res.n_a;
                        nB_mat(pi_r, fi)   = res.n_b;
                        medA_mat(pi_r, fi) = res.median_a;
                        medB_mat(pi_r, fi) = res.median_b;
                    end
                end
            end

            % BH-FDR correction, per feature, across this replicate's own
            % set of pairs only (replicates are corrected independently —
            % consistent with "Rep1 A vs B, A vs C" not being pooled with
            % "Rep2 A vs B, A vs C").
            p_fdr_mat = nan(size(p_mat));
            for fi = 1:n_pw_feats
                p_fdr_mat(:, fi) = bh_fdr(p_mat(:, fi));
            end

            biorep_group_tests(r).pair_labels = pair_labels_r;
            biorep_group_tests(r).p           = p_mat;
            biorep_group_tests(r).p_fdr       = p_fdr_mat;
            biorep_group_tests(r).stat        = stat_mat;
            biorep_group_tests(r).test_name   = tn_mat;
            biorep_group_tests(r).normal_A    = nrmA_mat;
            biorep_group_tests(r).normal_B    = nrmB_mat;
            biorep_group_tests(r).n_A         = nA_mat;
            biorep_group_tests(r).n_B         = nB_mat;
            biorep_group_tests(r).median_A    = medA_mat;
            biorep_group_tests(r).median_B    = medB_mat;
            biorep_group_tests(r).median_diff = abs(medA_mat - medB_mat);

            fprintf('  Replicate %-20s : %d/%d sample group(s) present (%s), %d pair(s)\n', ...
                rep_lbl, n_gp, n_groups, strjoin(labels(groups_present_idx), ', '), n_pairs_r);
        end
    end
    fprintf('STEP 6i complete.\n\n');

    % =========================================================
    %  BUILD AND OUTPUT HT STRUCTURE TO WORKSPACE
    % =========================================================
    % Assemble the final ht structure with all results
    % This is what plot_combined_samples(ht) needs
    ht = struct();
    ht.params  = p;
    ht.labels  = labels;
    ht.colors  = colors;
    ht.groups  = groups;
    ht.ix      = ix;
    ht.global  = gmax;
    ht.xbins   = xbins;  % ← Bin edges for plotting histograms
    ht.pairwise     = pairwise;      % pairwise KS D / Q90 / hypothesis test for all group pairs and features
    ht.biorep_corr  = biorep_corr;  % biological replicate correlation / stats (STEP 6f)
    ht.biorep_group_tests = biorep_group_tests;  % per-replicate, cross-group hypothesis tests (STEP 6i)

    % Metadata table — stored so downstream plot scripts (e.g.
    % plot_lagTime_scatter_violin) can access all metadata columns
    % without needing the raw  data  struct in the workspace.
    if isfield(data, 'metadata') && isfield(data.metadata, 'original')
        ht.metadata = data.metadata.original;   % full MATLAB table
    else
        ht.metadata = table();
    end

    % Per-plate biological + technical replicate labels
    % Stored as ht.plate_labels{i} = 'BioRep_TechRep' (or just 'BioRep')
    % so plot_lagTime_scatter_violin can label the x-axis without raw data.
    ht.plate_labels = build_plate_labels(ht.metadata, p, nr_plates);

    % Per-plate QC data — used by plot_QC_full_dataset (no raw data access needed)
    ht.qc.nr_plates            = nr_plates;
    ht.qc.colony_count_all     = qc_colony_count_all;    % length(flag_colony_ok)
    ht.qc.colony_count_clean   = qc_colony_count_clean;  % sum(flag_colony_ok)
    ht.qc.colony_count_final   = qc_colony_count_final;  % after size_threshold filter (= plot_combined_samples count)
    ht.qc.plate_area           = qc_plate_area;          % cell: per-plate Area (size-filtered)
    ht.qc.passes               = qc_passes;              % logical: passes both gates
    ht.qc.fn                   = qc_fn;                  % plate filenames
    ht.qc.manual_count_col     = qc_manual_col;          % chosen metadata column name
    ht.qc.manual_counts        = qc_manual_counts;       % numeric [nr_plates x 1]
    % Per-plate growth data — used by plot_QC_growth_data (no raw data access needed)
    ht.qc.ix_growth            = ix_growth;              % plate indices with growth_quant==true
    ht.qc.plate_time           = qc_plate_time;          % cell: per-plate time vector (h)
    ht.qc.plate_size_tc        = qc_plate_size_tc;       % cell: per-plate size timecourse [T x N]
    ht.qc.plate_valid          = qc_plate_valid;         % cell: per-plate valid-colony logical mask
    ht.qc.plate_lag            = qc_plate_lag;           % cell: per-plate lag times (h)
    ht.qc.plate_early_dt       = qc_plate_early_dt;      % cell: per-plate early doubling times (h)

    % Per-plate final colony size — used by plot_morphology_colonies_Gini.
    % per_plate_size{i} = size-threshold-passing rp_clean.Area values for plate i.
    ht.per_plate_size = per_plate_size;

    % Output to base workspace so other scripts can use it
    assignin('base', 'ht', ht);
    
    fprintf('\n========================================\n');
    fprintf('Preprocessing complete!\n');
    fprintf('ht structure assigned to workspace.\n');
    fprintf('  ht.qc  — per-plate QC counts, area arrays, manual counts, growth data\n');
    fprintf('Run: plot_combined_samples(ht);\n');
    fprintf('Run: plot_lagTime_colonies_Gini(ht, data);        %% lag time violin + Gini\n');
    fprintf('Run: plot_morphology_colonies_Gini(ht, data);     %% lag time + colony size violin + Gini\n');
    fprintf('Run: plot_lagTime_vs_growthRate(ht);  %% Gompertz trade-off plot\n');
    fprintf('Run: plot_QC_full_dataset(ht);        %% Plots 1, 1A, 2\n');
    fprintf('Run: plot_QC_full_dataset(ht, data);  %% + image montage\n');
    fprintf('Run: plot_QC_growth_data(ht);          %% growth QC plots\n');
    fprintf('========================================\n');

end


%% ================================================================
%  LOCAL FUNCTION: log_message
%  Writes a message to the log file.
% ================================================================
function [slope, intercept] = quantile_regression(y, x, tau, w)
    % Iteratively reweighted quantile regression at quantile tau.
    n = length(x);
    if nargin < 4 || isempty(w), w = ones(n,1); end
    w = w(:) / sum(w);
    b = [median(y); 0];
    X = [ones(n,1), x(:)];
    y = y(:);
    for iter = 1:200
        r  = y - X*b;
        qw = tau * (r >= 0) + (tau-1) * (r < 0);
        rw = max(abs(r), 1e-8);
        ww = w .* abs(qw) ./ rw;
        W  = diag(ww);
        b_new = (X'*W*X + 1e-10*eye(2)) \ (X'*W*y);
        if norm(b_new - b) < 1e-8, break; end
        b = b_new;
    end
    intercept = b(1);
    slope     = b(2);
end


function gr = empty_growth_struct(n)
    gr.mu_max          = NaN(n,1);
    gr.Amax            = NaN(n,1);
    gr.k               = NaN(n,1);
    gr.tm              = NaN(n,1);
    gr.fit_ok          = false(n,1);
    gr.local_density   = NaN(n,1);
    gr.spearman_rho    = NaN;
    gr.spearman_p      = NaN;
    gr.quantile_slope  = NaN;
    gr.quantile_intcpt = NaN;
    gr.wt_spearman     = NaN;
    gr.wt_qr_slope     = NaN;
end


function log_message(fid, msg)
    fprintf(fid, '%s\n', msg);
end


%% ================================================================
%  LOCAL FUNCTION: two_sample_test
%  Normality-aware two-sample hypothesis test on raw (unbinned)
%  per-colony values.
%    - Both samples consistent with a normal distribution (Lilliefors
%      test) -> Welch's two-sample t-test (unequal variances).
%    - Otherwise                                        -> Wilcoxon
%      rank-sum test (= Mann-Whitney U).
%  Returns NaN/'' fields if either sample has fewer than 3 valid
%  (non-NaN) observations.
% ================================================================
function res = two_sample_test(va, vb)
    va = va(~isnan(va));
    vb = vb(~isnan(vb));

    res.n_a         = length(va);
    res.n_b         = length(vb);
    res.median_a    = NaN;
    res.median_b    = NaN;
    res.p           = NaN;
    res.stat        = NaN;
    res.test_name   = '';
    res.normal_a    = false;
    res.normal_b    = false;

    if res.n_a < 3 || res.n_b < 3
        return;
    end

    res.median_a = median(va);
    res.median_b = median(vb);
    res.normal_a = local_is_normal(va);
    res.normal_b = local_is_normal(vb);

    if res.normal_a && res.normal_b
        [~, p, ~, stats] = ttest2(va, vb, 'Vartype', 'unequal');
        res.p         = p;
        res.stat      = stats.tstat;
        res.test_name = 't-test (Welch)';
    else
        [p, ~, stats] = ranksum(va, vb);
        res.p = p;
        if isfield(stats, 'zval')
            res.stat = stats.zval;
        else
            res.stat = stats.ranksum;
        end
        res.test_name = 'Wilcoxon rank-sum';
    end
end


%% ================================================================
%  LOCAL FUNCTION: local_is_normal
%  Lilliefors test for normality (built into the Statistics and
%  Machine Learning Toolbox — no external/file-exchange dependency).
%  h==1 rejects normality at alpha=0.05; treated conservatively as
%  "not normal" (-> non-parametric test) on error or insufficient data.
% ================================================================
function tf = local_is_normal(x)
    tf = false;
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
%  Benjamini-Hochberg false discovery rate correction for one family
%  of p-values (NaN entries are ignored and returned as NaN — they
%  correspond to comparisons that were never run, e.g. too little
%  data, and must not count towards the family size m).
% ================================================================
function q = bh_fdr(pvals)
    q = nan(size(pvals));
    valid = ~isnan(pvals);
    m = sum(valid);
    if m == 0
        return;
    end

    idx           = find(valid);
    [ps, order]   = sort(pvals(idx));
    ranks         = (1:m)';
    q_sorted      = ps(:) .* m ./ ranks;

    % Enforce monotonicity: q_(i) = min(q_(i), q_(i+1), ..., q_(m))
    for i = m-1:-1:1
        q_sorted(i) = min(q_sorted(i), q_sorted(i+1));
    end
    q_sorted = min(q_sorted, 1);

    q(idx(order)) = q_sorted;
end


%% ================================================================
%  LOCAL FUNCTION: compute_gini
%  Gini index (0 = perfectly equal, 1 = maximally unequal).
% ================================================================
function g = compute_gini(x)
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
%  LOCAL FUNCTION: generate_group_colors
%  Returns an [n x 3] RGB matrix with one perceptually distinct
%  colour per sample group.
% ================================================================
function C = generate_group_colors(n)
    base = [ ...
        0.37  0.21  0.65;   %  1  purple
        0.12  0.69  0.70;   %  2  teal
        0.87  0.71  0.00;   %  3  amber
        0.90  0.52  0.10;   %  4  orange
        0.20  0.63  0.17;   %  5  green
        0.84  0.15  0.16;   %  6  red
        0.12  0.47  0.71;   %  7  blue
        0.58  0.40  0.74;   %  8  lavender
        0.00  0.62  0.45;   %  9  sea green
        0.94  0.39  0.63;   % 10  rose / hot pink
        0.50  0.33  0.17;   % 11  brown
        0.00  0.45  0.70;   % 12  steel blue
        0.80  0.47  0.74;   % 13  orchid
        0.34  0.71  0.91;   % 14  sky blue
        0.50  0.50  0.00;   % 15  olive
        0.96  0.65  0.14;   % 16  golden yellow
        0.40  0.76  0.65;   % 17  sage / mint
        0.70  0.19  0.38;   % 18  crimson
    ];

    if n <= size(base, 1)
        C = base(1:n, :);
    else
        hsv_colors = hsv(n);
        hsv_colors(:,2) = 0.75;
        hsv_colors(:,3) = min(0.90, max(0.55, hsv_colors(:,3)));
        C = hsv2rgb(hsv_colors);
    end
end


%% ================================================================
%  LOCAL FUNCTION: ask_manual_count_column
%  Asks the user which metadata column holds manual colony counts.
%  Called at the end of preprocessing so plot_QC_full_dataset
%  needs no direct access to data.metadata.
% ================================================================
function col_name = ask_manual_count_column(data)

    if ~isfield(data.metadata, 'original')
        warning('preprocess_pipeline_data: data.metadata.original not found — skipping manual count column.');
        col_name = '';
        return;
    end

    meta_tbl  = data.metadata.original;
    col_names = meta_tbl.Properties.VariableNames;

    fprintf('\n--- 6) Metadata column for manual colony counts (QC Plot 1A) ---\n');
    fprintf('   This column is used to compare automated vs. manual counts.\n');
    fprintf('   Enter 0 to skip this plot.\n\n');
    fprintf('   Available columns:\n');

    for k = 1:length(col_names)
        try
            vals = unique(meta_tbl.(col_names{k}));
            if iscell(vals)
                preview = strjoin(vals(1:min(3,end)), ', ');
            elseif isnumeric(vals)
                preview = strtrim(num2str(vals(1:min(3,end))'));
            else
                preview = '(preview unavailable)';
            end
        catch
            preview = '(preview unavailable)';
        end
        fprintf('   [%2d] %-25s  e.g.  %s\n', k, col_names{k}, preview);
    end

    fprintf('\n');
    while true
        raw = strtrim(input('   Enter column number (or 0 to skip): ', 's'));
        idx = str2double(raw);
        if ~isnan(idx) && idx == 0
            col_name = '';
            fprintf('   -> Skipping manual count comparison.\n');
            return;
        end
        if ~isnan(idx) && idx >= 1 && idx <= length(col_names) && floor(idx) == idx
            col_name = col_names{idx};
            fprintf('   -> Manual count column: %s\n\n', col_name);
            return;
        end
        fprintf('   WARNING  Please enter a whole number between 0 and %d.\n', length(col_names));
    end
end


%% ================================================================
%  LOCAL FUNCTION: get_manual_counts
%  Extracts manual count vector from the chosen metadata column.
%  Returns a numeric [nr_plates x 1] vector (NaN where unavailable).
% ================================================================
function counts = get_manual_counts(data, col_name, nr_plates)
    counts = nan(nr_plates, 1);

    if isempty(col_name) || ~isfield(data.metadata, 'original')
        return;
    end

    meta_tbl = data.metadata.original;
    if ~ismember(col_name, meta_tbl.Properties.VariableNames)
        warning('preprocess_pipeline_data: column "%s" not found in metadata.', col_name);
        return;
    end

    raw = meta_tbl.(col_name);
    n   = min(nr_plates, length(raw));

    if isnumeric(raw)
        counts(1:n) = raw(1:n);
    elseif iscell(raw)
        counts(1:n) = str2double(raw(1:n));
    end
end


%% ================================================================
%  LOCAL FUNCTION: round_to_sig
%  Round x to n significant figures (for bin step sizing).
% ================================================================
function y = round_to_sig(x, n)
    if x == 0, y = 1; return; end
    d = ceil(log10(abs(x)));
    p = n - d;
    y = round(x * 10^p) / 10^p;
    if y <= 0, y = 10^(-p); end
end


%% ================================================================
%  LOCAL FUNCTION: get_cfg
%  Read field from cfg struct if present, else return default.
% ================================================================
function val = get_cfg(cfg, field, default)
    if isfield(cfg, field)
        val = cfg.(field);
    else
        val = default;
    end
end


%% ================================================================
%  LOCAL FUNCTION: get_user_params
%  Collects interactive user input for pipeline parameters.
%  INCLUDES METADATA COLUMN SELECTION
% ================================================================
function params = get_user_params(data, nr_plates)
    
    % 1. Incubation time (RT)
    fprintf('\n--- 1) Room Temperature (RT) incubation time ---\n');
    while true
        raw = input('   Enter incubation time in hours [e.g. 20 or 20.5]: ', 's');
        val = str2double(strtrim(raw));
        if ~isnan(val) && val >= 0
            incTime = val;
            fprintf('   -> RT incubation time: %.4g h\n', incTime);
            break;
        end
        fprintf('   WARNING  Invalid input. Please enter a positive number.\n');
    end

    % 2. Imaging interval
    fprintf('\n--- 2) Imaging interval ---\n');
    while true
        raw_unit = input('   Time unit - enter  1  for minutes  or  2  for hours: ', 's');
        unit_choice = str2double(strtrim(raw_unit));
        if ismember(unit_choice, [1 2]), break; end
        fprintf('   WARNING  Please enter 1 (minutes) or 2 (hours).\n');
    end
    while true
        raw_period = input('   Interval value [e.g. 30, 0.5]: ', 's');
        period_val = str2double(strtrim(raw_period));
        if ~isnan(period_val) && period_val > 0, break; end
        fprintf('   WARNING  Invalid input. Please enter a positive number.\n');
    end
    if unit_choice == 1
        img_int = period_val / 60;
        fprintf('   -> Imaging interval: %.4g min = %.4g h\n', period_val, img_int);
    else
        img_int = period_val;
        fprintf('   -> Imaging interval: %.4g h\n', img_int);
    end

    % 3. Lag time
    fprintf('\n--- 3) Lag Time ---\n');
    max_img_count = 0;
    for k = 1:nr_plates
        if isfield(data.processed{k}.colonies.new, 'timecourse_size_smoothed')
            n = size(data.processed{k}.colonies.new.timecourse_size_smoothed, 1);
            max_img_count = max(max_img_count, n);
        end
    end
    max_lag_calc = incTime + (max_img_count - 1) * img_int;

    fprintf('   Max image count found: %d images\n', max_img_count);
    h_part = floor(max_lag_calc);
    m_part = round((max_lag_calc - h_part) * 60);
    if m_part == 0
        fprintf('   Formula: %.4g + (%d-1) x %.4g = %.4g h\n', incTime, max_img_count, img_int, max_lag_calc);
    else
        fprintf('   Formula: %.4g + (%d-1) x %.4g = %.4g h  (%d h %d min)\n', ...
                incTime, max_img_count, img_int, max_lag_calc, h_part, m_part);
    end
    fprintf('   Press ENTER to use calculated value (%.4g h), or enter override: ', max_lag_calc);
    raw_lag = strtrim(input('', 's'));
    if isempty(raw_lag)
        max_lag = max_lag_calc;
        fprintf('   -> Using calculated max lag time: %.4g h\n', max_lag);
    else
        custom = str2double(raw_lag);
        if ~isnan(custom) && custom > 0
            max_lag = custom;
            fprintf('   -> Using custom max lag time: %.4g h\n', max_lag);
        else
            fprintf('   WARNING  Invalid — falling back to calculated value.\n');
            max_lag = max_lag_calc;
        end
    end

    % 4. Output directory
    fprintf('\n--- 4) Output directory ---\n');
    fprintf('   Enter full path to existing folder,\n');
    fprintf('   or press ENTER for current directory (%s).\n', pwd);
    while true
        raw_dir = strtrim(input('   Output directory: ', 's'));
        if isempty(raw_dir)
            out_dir = pwd;
            fprintf('   -> Saving to: %s\n', out_dir);
            break;
        end
        if isunix && strncmp(raw_dir, '~', 1)
            raw_dir = [getenv('HOME') raw_dir(2:end)];
        end
        if isfolder(raw_dir)
            out_dir = raw_dir;
            fprintf('   -> Output directory: %s\n', out_dir);
            break;
        else
            fprintf('   WARNING  Directory not found: %s\n', raw_dir);
            create_it = strtrim(input('   Create it now? (y/n): ', 's'));
            if strcmpi(create_it, 'y')
                mkdir(raw_dir);
                out_dir = raw_dir;
                fprintf('   -> Directory created: %s\n', out_dir);
                break;
            else
                fprintf('   Please enter a valid path or press ENTER for current directory.\n');
            end
        end
    end

    % 5. METADATA COLUMN SELECTION (with preview)
    fprintf('\n--- 5) Metadata column for sample grouping ---\n');
    
    if isfield(data.metadata, 'original')
        meta_tbl  = data.metadata.original;
        col_names = meta_tbl.Properties.VariableNames;
    else
        error('ERROR: No metadata.original table found in data structure.');
    end
    
    fprintf('   Available columns:\n');
    for k = 1:length(col_names)
        try
            vals = unique(meta_tbl.(col_names{k}));
            if iscell(vals)
                preview = strjoin(vals(1:min(3, end)), ', ');
            else
                preview = strtrim(num2str(vals(1:min(3, end))'));
            end
        catch
            preview = '(preview unavailable)';
        end
        fprintf('   [%d] %-20s  e.g. %s\n', k, col_names{k}, preview);
    end
    
    while true
        raw_col = strtrim(input('   Enter column number: ', 's'));
        col_idx = str2double(raw_col);
        if ~isnan(col_idx) && col_idx >= 1 && col_idx <= length(col_names) && floor(col_idx) == col_idx
            group_col = col_names{col_idx};
            fprintf('   -> Using column: %s\n', group_col);
            break;
        end
        fprintf('   WARNING  Please enter a whole number between 1 and %d.\n', length(col_names));
    end

    % 6. BIOLOGICAL REPLICATE COLUMN
    fprintf('\n--- 6) Biological replicate column ---\n');
    fprintf('   Identifies the biological sample / strain / condition\n');
    fprintf('   (e.g. Set, Strain, BioRep, SampleName).\n');
    fprintf('   This will be used as the X-axis label in violin plots.\n');
    fprintf('   Enter 0 to use plate position numbers only.\n\n');
    fprintf('   Available columns:\n');
    for k = 1:length(col_names)
        try
            vals = unique(meta_tbl.(col_names{k}));
            if iscell(vals)
                preview = strjoin(vals(1:min(4, end)), ', ');
            else
                preview = strtrim(num2str(vals(1:min(4, end))'));
            end
        catch
            preview = '(preview unavailable)';
        end
        fprintf('   [%d] %-20s  e.g. %s\n', k, col_names{k}, preview);
    end

    while true
        raw_bio = strtrim(input('   Enter column number (0 to skip): ', 's'));
        bio_idx = str2double(raw_bio);
        if ~isnan(bio_idx) && bio_idx == 0
            bio_rep_col = '';
            fprintf('   -> Skipped — will use plate numbers.\n');
            break;
        end
        if ~isnan(bio_idx) && bio_idx >= 1 && bio_idx <= length(col_names) && floor(bio_idx) == bio_idx
            bio_rep_col = col_names{bio_idx};
            fprintf('   -> Biological replicate column: %s\n', bio_rep_col);
            break;
        end
        fprintf('   WARNING  Please enter a whole number between 0 and %d.\n', length(col_names));
    end

    % 7. TECHNICAL REPLICATE COLUMN
    fprintf('\n--- 7) Technical replicate column ---\n');
    fprintf('   Identifies replicate plates of the same biological sample\n');
    fprintf('   (e.g. Replicate, TechRep, PlateRep, R1/R2/R3).\n');
    fprintf('   This is appended to the bio-rep label: BioRep_TechRep.\n');
    fprintf('   Enter 0 to skip (bio-rep label only).\n\n');

    while true
        raw_tech = strtrim(input('   Enter column number (0 to skip): ', 's'));
        tech_idx = str2double(raw_tech);
        if ~isnan(tech_idx) && tech_idx == 0
            tech_rep_col = '';
            fprintf('   -> Skipped.\n');
            break;
        end
        if ~isnan(tech_idx) && tech_idx >= 1 && tech_idx <= length(col_names) && floor(tech_idx) == tech_idx
            tech_rep_col = col_names{tech_idx};
            fprintf('   -> Technical replicate column: %s\n', tech_rep_col);
            break;
        end
        fprintf('   WARNING  Please enter a whole number between 0 and %d.\n', length(col_names));
    end

    fprintf('\n============================================================\n\n');

    params.incTime      = incTime;
    params.img_int      = img_int;
    params.max_lag      = max_lag;
    params.out_dir      = out_dir;
    params.group_col    = group_col;
    params.bio_rep_col  = bio_rep_col;
    params.tech_rep_col = tech_rep_col;

end

%% ================================================================
%  LOCAL FUNCTION: build_plate_labels
%  Builds a {nr_plates x 1} cell of x-axis label strings using the
%  bio-rep and tech-rep columns chosen during parameter setup.
%  Label format:
%    bio+tech:  "<BioRep>_<TechRep>"   e.g.  "Set1_R2"
%    bio only:  "<BioRep>"             e.g.  "Set1"
%    neither:   "Plate<i>"             e.g.  "Plate7"
% ================================================================
function labels_out = build_plate_labels(meta_tbl, p, nr_plates)

    labels_out = cell(nr_plates, 1);

    has_bio  = isfield(p, 'bio_rep_col')  && ~isempty(p.bio_rep_col)  && ...
               ~isempty(meta_tbl) && ismember(p.bio_rep_col,  meta_tbl.Properties.VariableNames);
    has_tech = isfield(p, 'tech_rep_col') && ~isempty(p.tech_rep_col) && ...
               ~isempty(meta_tbl) && ismember(p.tech_rep_col, meta_tbl.Properties.VariableNames);

    n_meta = height(meta_tbl);

    for i = 1:nr_plates
        if i > n_meta || (~has_bio && ~has_tech)
            labels_out{i} = sprintf('Plate%d', i);
            continue;
        end

        if has_bio
            bio_lbl = val2str(meta_tbl.(p.bio_rep_col)(i));
        else
            bio_lbl = sprintf('Plate%d', i);
        end

        if has_tech
            tech_lbl = val2str(meta_tbl.(p.tech_rep_col)(i));
            labels_out{i} = [bio_lbl '_' tech_lbl];
        else
            labels_out{i} = bio_lbl;
        end
    end

    % Confirm to console
    fprintf('\n  Bio/tech rep labels — first 6 plates:\n');
    for i = 1:min(6, nr_plates)
        fprintf('    Plate %2d  ->  %s\n', i, labels_out{i});
    end
    if nr_plates > 6
        fprintf('    ... (%d more)\n', nr_plates - 6);
    end
    fprintf('\n');
end


%% ================================================================
%  LOCAL FUNCTION: strip_quotes
%  Remove surrounding single-quote apostrophes from a string if
%  present (e.g. "'Set1'" -> "Set1", "'3h'" -> "3h").
%  Some metadata importers store cell-string values with the
%  apostrophes baked in as literal characters.
% ================================================================
function s = strip_quotes(s)
    if ischar(s) && length(s) >= 2 && s(1) == '''' && s(end) == ''''
        s = s(2:end-1);
    end
end


%% ================================================================
%  LOCAL FUNCTION: val2str
%  Convert a single metadata table value to a plain string.
%  Also strips surrounding apostrophes that some importers leave
%  baked into the string value (e.g. "'Set1'" -> "Set1").
% ================================================================
function s = val2str(v)
    if iscell(v),             s = strtrim(char(v{1}));
    elseif isnumeric(v),      s = strtrim(num2str(v));
    elseif ischar(v),         s = strtrim(v);
    elseif isstring(v),       s = strtrim(char(v));
    else,                     s = '?';
    end
    % Strip surrounding single-quote apostrophes that some CSV/Excel
    % importers leave embedded in the string value, e.g. "'Set1'" -> "Set1"
    if length(s) >= 2 && s(1) == '''' && s(end) == ''''
        s = s(2:end-1);
    end
end
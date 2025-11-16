%% 11.11.2025 - HeteroTyper Pipeline for Dark Plates 

function plot_QC_full_dataset(data)

    nr_plates = length(data.processed);
    
    % Define min-max colony numbers
    min_col = 10;
    max_col = 350;
    
    %% Plot 1: number of detected colonies and number of "clean" colonies
    colony_count = zeros(nr_plates,1);
    colony_count_clean = zeros(nr_plates,1);
    
    for i = 1:nr_plates
        colony_count(i) = length(data.processed{i}.colonies.region_props.flag_colony_ok);
        colony_count_clean(i) = length(find(data.processed{i}.colonies.region_props.flag_colony_ok));
    
    end
    
    figure('Name','colony count');
    semilogy(colony_count,'ok');
    hold on;
    semilogy(colony_count_clean,'or');
    for i = 1:nr_plates
        line([i i],[colony_count(i) colony_count_clean(i)],'Color','r');
    end
    line([0.5 nr_plates+0.5],[50 50],'Color','k','LineStyle','--');
    line([0.5 nr_plates+0.5],[500 500],'Color','k','LineStyle','--');
    
    set(gca,'FontSize',10);
    axis([0.5 nr_plates+0.5 1 inf]);
    
    %% Plot 1A: comparison with manual counts
    ix_original = find(data.metadata.original.Count > 0);
    figure('Name','colony count');
    subplot(1,2,1),...
    semilogy(colony_count_clean(ix_original),'or');
    hold on;
    semilogy(data.metadata.original.Count(ix_original),'ob');
    set(gca,'FontSize',10,'XTick',[1:length(ix_original)],'XTickLabel',ix_original);
    axis([0.5 length(ix_original)+0.5 1 inf]);
    xlabel('plate positions');
    ylabel('colony count');
    xtickangle(-45);
    
    subplot(1,2,2),...
    loglog(data.metadata.original.Count(ix_original),colony_count_clean(ix_original),'ok');
    hold on;
    line([1 1000],[1 1000],'Color','k');
    set(gca,'FontSize',10);
    xlabel('manual count');
    ylabel('automated count');
    
    
    %% Plot 2: distribution of colony sizes
    figure('Name','colony size');
    runIx = 0;
    rel_plates = [];
    
    for i = 1:nr_plates
        flag_ok  = find(data.processed{i}.colonies.region_props.flag_colony_ok);
        col_area = data.processed{i}.colonies.region_props.Area;
        if (~isempty(flag_ok))
            colony_size_distr = col_area(flag_ok);
            colony_size_median(i) = median(colony_size_distr);
        else
            colony_size_distr = [];
            colony_size_median(i) = 0;
        end
        if(length(colony_size_distr) > min_col)&&(length(colony_size_distr) < max_col)
            rel_plates = [rel_plates;i];
            runIx = runIx + 1;
    
            edges   = 0:100:max(col_area);
            counts  = histcounts(colony_size_distr, edges);
            y_max   = max(counts);
    
            subplot(3,8,runIx),...      % subplot(nrow,ncol,..) ->  nrow: number of sample groups | ncol: number of plates in sample groups (take max across all)
            histogram(colony_size_distr, edges);
            hold on;
            t = data.metadata.fn(i).name;
            t = strrep(t,'_','-');
            title(t,'FontWeight','Normal');
            set(gca,'FontSize',6);
            axis([0 max(col_area) 0 max(1, ceil(1.1*y_max))]);
            text(0, max(1, 0.85*y_max), num2str(length(colony_size_distr)));
        end
    end
    
    
    
    %% Plot 3: visual inspection of colonies

    img_montage = {};
    img_montage_all = {};
    runIx = 0;
    runIx_all = 0;
    
    text_x = repmat([25:590:2800],1,5);
    text_y = [repmat(text_x(1),1,6),repmat(text_x(2),1,6),repmat(text_x(3),1,6),repmat(text_x(4),1,6),repmat(text_x(5),1,6)];
    
    x_shift = 20;
    text_x_shifted = text_x + x_shift;
    
    for i = 1:nr_plates
        runIx = runIx + 1;
    
        % get mask and overlay
        mask_clean = data.processed{i}.colonies.debug.segmented;
        img_tmp2   = imadjust(rgb2gray(data.processed{i}.img_final), [0.025 0.15], []);
        img_tmp    = imfuse(mask_clean, img_tmp2, 'blend');
    
        % crop region of interest (adapt if needed)
        img_montage{runIx} = imcrop(img_tmp, [1500 1500 1000 1000]);
    
        % show montage when we collected 30 images or reached the last plate
        if runIx == 30 || i == nr_plates
            runIx_all = runIx_all + 1;
            figure('Name', sprintf('colony image: %d', runIx_all));
            montage(img_montage, 'Size', [5 6], 'BorderSize', 5, 'BackgroundColor', 'w');
            hold on;
    
            % starting plate index for this batch
            group_start = i - runIx + 1;
            group_idx   = group_start:i;
    
            % annotate each tile with the clean colony count for its plate
            for j = 1:runIx
                rectangle('Position', [text_x_shifted(j)-40, text_y(j)-15, 80, 40], ...
                          'FaceColor', 'w', 'EdgeColor', 'none');
                text(text_x_shifted(j), text_y(j), ...
                     num2str(colony_count_clean(group_idx(j))), ...
                     'Color', 'k', 'FontSize', 12, ...
                     'HorizontalAlignment', 'center', ...
                     'VerticalAlignment', 'middle');
            end
    
            % reset for next batch
            runIx       = 0;
            img_montage = {};
        end
    end
end
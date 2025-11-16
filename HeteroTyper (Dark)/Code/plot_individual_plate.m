%% 11.11.2025 - HeteroTyper Pipeline for Dark Plates

function plot_individual_plate(data)

    max_lag = 55;

    % Define Room Temperature Incubation Time
    incTime = 24;
    
    % identify plate position
    plate_id = inputdlg('Plate position (1-104)?');
    plate_id = str2num(plate_id{1});

    if plate_id < 1 || plate_id > length(data.processed)
        error('Invalid plate index.');
    end
    
    % color-mode
    color_mode = inputdlg('color by (0) segmentation, (1) lag-time, (2), col size?');
    color_mode = str2num(color_mode{1});
    
    % load mask
    mask_clean = data.processed{plate_id}.colonies.mask_clean;
    % load final img
    img_tmp2 = imadjust(rgb2gray(data.processed{plate_id}.img_final),[0.025 0.15],[]);
    % colony flag
    flag_col = find(data.processed{plate_id}.colonies.region_props.flag_colony_ok);
    
    max_val = max(data.processed{plate_id}.colonies.new.timecourse_size_smoothed(end,:));

    fig_w = 1200;
    fig_h = 400;
    
    figure('Name',strcat('Plate:',num2str(plate_id)),'Position',[100 100 fig_w fig_h]);
    subplot(1,3,1:2),...
        
    if(color_mode == 0)
        mask_bin = label2rgb(mask_clean,'jet','k','shuffle');
        mask_bin_r = mask_bin(:,:,1);
        mask_bin_g = mask_bin(:,:,2);
        mask_bin_b = mask_bin(:,:,3);
        
        img_tmp = imfuse(mask_bin,img_tmp2,'blend');
        
        imshow(img_tmp);
        hold on;
    elseif(color_mode == 1)
        mask_bin = zeros(size(mask_clean));
        for i = 1:length(flag_col)
            lag_time = data.processed{plate_id}.colonies.new.lag_time(flag_col(i)) + incTime;
            mask_bin(find(mask_clean == flag_col(i))) = lag_time;
            
            feat_ratio(flag_col(i),1) = lag_time./max_lag;
            
        end
        imshow(mask_bin);
        c = colormap(parula(max_lag-1));
        c = [0 0 0;c];
        colormap(c);
        caxis([0 max_lag]);
        colorbar('eastoutside');
        colorrange_nr = size(c,1);
        colorrange = [0:1/(colorrange_nr-1):1];
        feat_ratio(find(feat_ratio>1)) = 1;
        
    elseif(color_mode == 2)
        mask_bin = zeros(size(mask_clean));
        for i = 1:length(flag_col)
            final_col_size = data.processed{plate_id}.colonies.new.timecourse_size_smoothed(end,flag_col(i));
            mask_bin(find(mask_clean == flag_col(i))) = final_col_size;
            feat_ratio(flag_col(i),1) = final_col_size./max_val;
        end
        c = colormap(parula(15));
        c = [0 0 0;c];
        imshow(mask_bin);
        colormap(c);
        caxis([0 max_val]);
        colorbar('eastoutside');
        
        colorrange_nr = size(c,1);
        colorrange = [0:1/(colorrange_nr-1):1];
        feat_ratio(find(feat_ratio>1)) = 1;
    end
    
    
    
    
    for i = 1:length(flag_col)
        col_center = data.processed{plate_id}.colonies.region_props.Centroid(flag_col(i),:);
        if(color_mode == 0)
            text(col_center(1),col_center(2),num2str(flag_col(i)),'Color','r','FontSize',8);
        else
            %
            ixt = find(colorrange>=feat_ratio(flag_col(i),1));
            
            curve_col(flag_col(i),:) = c(ixt(1),:);
            
            text(col_center(1),col_center(2),num2str(flag_col(i)),'Color',[0.8 0.8 0.8],'FontSize',8);
        end
    end
    
    if(data.processed{plate_id}.growth_quant == 1)
        time = data.processed{plate_id}.time.elapsed_time_h + incTime;
        
        
        for j = 1:length(flag_col)
            y_val = data.processed{plate_id}.colonies.new.timecourse_size_smoothed(:,flag_col(j));
            
            
            if(color_mode == 0)
                r_props_pxls = data.processed{plate_id}.colonies.region_props.PixelIdxList(flag_col(j));
                colony_color(1) = mask_bin_r(r_props_pxls{1}(1));
                colony_color(2) = mask_bin_g(r_props_pxls{1}(1));
                colony_color(3) = mask_bin_b(r_props_pxls{1}(1));
            
                subplot(1,3,3),...
                plot(time,y_val,'-k','Color',colony_color);
                hold on;
                
            else
                colony_color = curve_col(flag_col(j),:);
                subplot(1,3,3),...
                plot(time,y_val,'-k','Color',colony_color);
                hold on;
            end
        end
        subplot(1,3,3),...
        axis([0 max_lag 0 max_val]);
        xlabel('time [h]');
        ylabel('colony size [px]');
        
    end
end
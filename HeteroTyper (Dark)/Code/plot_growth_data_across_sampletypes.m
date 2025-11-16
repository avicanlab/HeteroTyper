%% 11.11.2025 - HeteroTyper Pipeline for Dark Plates 

function data = plot_growth_data_across_sampletypes(data)

    nr_plates = length(data.processed);
    
    % Define maximum lag time
    max_lag = 55;

    % Define Room Temperature Incubation Time
    incTime = 24;
    
    % Define min-max colony numbers
    min_col = 10;
    max_col = 350;

    n_groups = 3;
    max_plate = 8;

    %% extract plates where we did quantify growth parameters 
    for i = 1:nr_plates 
        growth_available(i,1) = data.processed{i}.growth_quant; 
        
        % --- Skip plates with too many or too few colonies ---
        % manually adapt if needed: 
        n_col = length(find(data.processed{i}.colonies.region_props.flag_colony_ok));
        if n_col > max_col || n_col < min_col
            growth_available(i,1) = 0;
        end
    end 
    
    ix_growth = find(growth_available); 
    

    % Following parameters (max_val, x_val, etc.) must be defined after identifying the maximum values for each parameter. 
    % For these parameters to be defined, code must be run once.
    max_val         = 10800; 
    xval            = [0:50:max_val]; % colony size bins 
    xval2           = [incTime:0.5:max_lag]; % lag time bins 
    xval3           = [0:1:max_lag]; % early DT bins 
    
    %% split things by organ
    ix_MLN = find(strcmpi(data.metadata.original.Organ,'MLN')); 
    ix_Spleen = find(strcmpi(data.metadata.original.Organ,'Spleen')); 
    ix_Liver = find(strcmpi(data.metadata.original.Organ,'Liver'));  
    labels_type = {'MLN','Spleen','Liver'}; 
    
    ix1 = ix_MLN; 
    ix1 = intersect(ix1,ix_growth); 
    ix2 = ix_Spleen; 
    ix2 = intersect(ix2,ix_growth); 
    ix3 = ix_Liver; 
    ix3 = intersect(ix3,ix_growth); 
    
    ix{1} = ix1; 
    ix{2} = ix2; 
    ix{3} = ix3; 
    
    colors = [0.37 0.21 0.65;0.12 0.69 0.70;0.87 0.71 0];

    %% Define filter parameters
    size_threshold = 200; 
    ecc_threshold = 0.70;

    figure('Name','final col size vs col intensity');
    for i = 1:length(ix)
        ix_tmp = ix{i};
        for j = 1:length(ix_tmp)
            col_size = data.processed{ix_tmp(j)}.colonies.new.timecourse_size_smoothed(end,:);
            col_ecc  = data.processed{ix_tmp(j)}.colonies.region_props.Eccentricity;
            col_int  = data.processed{ix_tmp(j)}.colonies.new.timecourse_intensity_smoothed(end,:);
            col_flag = data.processed{ix_tmp(j)}.colonies.region_props.flag_colony_ok;
            
            valid_idx = (col_flag == 1 & col_size > size_threshold & col_ecc < ecc_threshold);

            subplot(1,length(ix),i),...
            plot(col_size(valid_idx),col_int(valid_idx),'.k');
            hold on;
        end
        axis([0 inf 0 inf]);
    end
    
    figure('Name','split distributions, size');
    for i = 1:length(ix)
        ix_tmp = ix{i};
        for j = 1:length(ix_tmp)
            col_size = data.processed{ix_tmp(j)}.colonies.new.timecourse_size_smoothed(end,:);
            col_ecc  = data.processed{ix_tmp(j)}.colonies.region_props.Eccentricity;
            col_flag = data.processed{ix_tmp(j)}.colonies.region_props.flag_colony_ok;
            
            valid_idx = (col_flag == 1 & col_size > size_threshold & col_ecc < ecc_threshold);
    
            subplot(n_groups,max_plate,j+(i-1)*max_plate),...
            histogram(col_size(valid_idx),xval,'FaceColor',colors(i,:),'EdgeColor','none','FaceAlpha',1);
            axis([0 inf 0 inf]);
    
            %t = strcat('Pos:',num2str(ix_tmp(j)),'-',data.metadata.original.Time(ix_tmp(j)),'-(',num2str(length(find(col_flag))),')');
            t = data.metadata.original.Time(ix_tmp(j));
            t = strrep(t,'_',' ');
            % t = strrep(t,'_','-');
            title(t,'FontWeight','normal');
            set(gca,'FontSize',8);
        end
    end
    
    figure('Name','split distributions, lag time');
    for i = 1:length(ix)
        ix_tmp = ix{i};
        for j = 1:length(ix_tmp)
            lag_time = data.processed{ix_tmp(j)}.colonies.new.lag_time + incTime;
            col_size = data.processed{ix_tmp(j)}.colonies.new.timecourse_size_smoothed(end,:);
            col_ecc  = data.processed{ix_tmp(j)}.colonies.region_props.Eccentricity;
            col_flag = data.processed{ix_tmp(j)}.colonies.region_props.flag_colony_ok;
            
            valid_idx = (col_flag == 1 & col_size > size_threshold & col_ecc < ecc_threshold);
    
            subplot(n_groups,max_plate,j+(i-1)*max_plate),...
            histogram(lag_time(valid_idx),xval2,'FaceColor',colors(i,:),'EdgeColor','none','FaceAlpha',1);
            axis([0 inf 0 inf]);
    
            %t = strcat('Pos:',num2str(ix_tmp(j)),'-',data.metadata.original.Time(ix_tmp(j)),'-(',num2str(length(find(col_flag))),')');
            %text(0,0.8.*max(N),num2str(length(find(col_flag))));
            t = data.metadata.original.Time(ix_tmp(j));
            t = strrep(t,'_',' ');
            title(t,'FontWeight','normal');
            set(gca,'FontSize',8);
        end
    end
    
    figure('Name','split distributions, early DT');
    for i = 1:length(ix)
        ix_tmp = ix{i};
        for j = 1:length(ix_tmp)
            earlyDT = data.processed{ix_tmp(j)}.colonies.new.early_doublingtime;
            col_size = data.processed{ix_tmp(j)}.colonies.new.timecourse_size_smoothed(end,:);
            col_ecc  = data.processed{ix_tmp(j)}.colonies.region_props.Eccentricity;
            col_flag = data.processed{ix_tmp(j)}.colonies.region_props.flag_colony_ok;
            
            valid_idx = (col_flag == 1 & col_size > size_threshold & col_ecc < ecc_threshold);
    
            subplot(n_groups,max_plate,j+(i-1)*max_plate),...
            histogram(earlyDT(valid_idx),xval3,'FaceColor',colors(i,:),'EdgeColor','none','FaceAlpha',1);
            axis([0 inf 0 inf]);
    
            %t = strcat('Pos:',num2str(ix_tmp(j)),'-',data.metadata.original.Time(ix_tmp(j)),'-(',num2str(length(find(col_flag))),')');
            t = data.metadata.original.Time(ix_tmp(j));
            t = strrep(t,'_',' ');
            title(t,'FontWeight','normal');
            set(gca,'FontSize',8);
        end
    end
    
    %%Scatter plots with 2D distributions
    figure('Name','split distributions, col size vs lag time with density');
    for i = 1:length(ix)
        ix_tmp = ix{i};
        for j = 1:length(ix_tmp)
            lag_time = data.processed{ix_tmp(j)}.colonies.new.lag_time + incTime;
            col_size = data.processed{ix_tmp(j)}.colonies.new.timecourse_size_smoothed(end,:);
            col_ecc  = data.processed{ix_tmp(j)}.colonies.region_props.Eccentricity;
            col_flag = data.processed{ix_tmp(j)}.colonies.region_props.flag_colony_ok;
            
            valid_idx = (col_flag == 1 & col_size > size_threshold & col_ecc < ecc_threshold);
    
            subplot(n_groups,max_plate,j+(i-1)*max_plate),...
            dscatter(col_size(valid_idx)',lag_time(valid_idx));
            m1 = median(col_size(valid_idx)');
            m2 = median(lag_time(valid_idx));
    
            hold on;
            line([m1 m1],[0 m2],'Color','r');
            line([0 m1],[m2 m2],'Color','r');
    
            histogram(lag_time(find(col_flag)),xval2,'FaceColor',colors(i,:),'EdgeColor','none','FaceAlpha',1);
            axis([0 inf 0 inf]);
    
            %t = strcat('Pos:',num2str(ix_tmp(j)),'-',data.metadata.original.Time(ix_tmp(j)),'-(',num2str(length(find(col_flag))),')');
            t = data.metadata.original.Time(ix_tmp(j));
            t = strrep(t,'_',' ');
            title(t,'FontWeight','normal');
            set(gca,'FontSize',8);
        end
    end
    
    
    figure('Name','split distributions, col size vs early DT with density');
    for i = 1:length(ix)
        ix_tmp = ix{i};
        for j = 1:length(ix_tmp)
            earlyDT = data.processed{ix_tmp(j)}.colonies.new.early_doublingtime; 
            col_size = data.processed{ix_tmp(j)}.colonies.new.timecourse_size_smoothed(end,:);
            col_ecc  = data.processed{ix_tmp(j)}.colonies.region_props.Eccentricity;
            col_flag = data.processed{ix_tmp(j)}.colonies.region_props.flag_colony_ok;
            
            valid_idx = (col_flag == 1 & col_size > size_threshold & col_ecc < ecc_threshold);
    
            subplot(n_groups,max_plate,j+(i-1)*max_plate),...
            dscatter(col_size(valid_idx)',earlyDT(valid_idx));
            m1 = median(col_size(valid_idx)');
            m2 = median(earlyDT(valid_idx));
    
            hold on;
            line([m1 m1],[0 m2],'Color','r');
            line([0 m1],[m2 m2],'Color','r');
    
            %histogram(lag_time(find(col_flag)),xval3,'FaceColor',colors(i,:),'EdgeColor','none','FaceAlpha',1);
            axis([0 inf 0 inf]);
    
            t = data.metadata.original.Time(ix_tmp(j));
            t = strrep(t,'_',' ');
            title(t,'FontWeight','normal');
            set(gca,'FontSize',8);
        end
    end
    
    figure('Name','split distributions, lag time vs early DT with density');
    for i = 1:length(ix)
        ix_tmp = ix{i};
        for j = 1:length(ix_tmp)
            lag_time = data.processed{ix_tmp(j)}.colonies.new.lag_time + incTime;
            earlyDT = data.processed{ix_tmp(j)}.colonies.new.early_doublingtime;
            col_size = data.processed{ix_tmp(j)}.colonies.new.timecourse_size_smoothed(end,:);
            col_ecc  = data.processed{ix_tmp(j)}.colonies.region_props.Eccentricity;
            col_flag = data.processed{ix_tmp(j)}.colonies.region_props.flag_colony_ok;
            
            valid_idx = (col_flag == 1 & col_size > size_threshold & col_ecc < ecc_threshold);
    
            subplot(n_groups,max_plate,j+(i-1)*max_plate),...
            dscatter(lag_time(valid_idx),earlyDT(valid_idx));
            m1 = median(lag_time(valid_idx)');
            m2 = median(earlyDT(valid_idx));
    
            hold on;
            line([m1 m1],[0 m2],'Color','r');
            line([0 m1],[m2 m2],'Color','r');
    
            %histogram(lag_time(find(col_flag)),xval2,'FaceColor',colors(i,:),'EdgeColor','none','FaceAlpha',1);
            axis([0 inf 0 inf]);
    
            %t = strcat('Pos:',num2str(ix_tmp(j)),'-',data.metadata.original.Time(ix_tmp(j)),'-(',num2str(length(find(col_flag))),')');
            t = data.metadata.original.Time(ix_tmp(j));
            t = strrep(t,'_',' ');
            title(t,'FontWeight','normal');
            set(gca,'FontSize',8);
        end
    end
    
    
    figure('Name','split distributions, col size vs lag time');
    for i = 1:length(ix)
        ix_tmp = ix{i};
        for j = 1:length(ix_tmp)
            lag_time = data.processed{ix_tmp(j)}.colonies.new.lag_time + incTime;
            col_size = data.processed{ix_tmp(j)}.colonies.new.timecourse_size_smoothed(end,:);
            col_ecc  = data.processed{ix_tmp(j)}.colonies.region_props.Eccentricity;
            col_flag = data.processed{ix_tmp(j)}.colonies.region_props.flag_colony_ok;
            
            valid_idx = (col_flag == 1 & col_size > size_threshold & col_ecc < ecc_threshold);
    
            subplot(n_groups,max_plate,j+(i-1)*max_plate),...
            plot(col_size(valid_idx)',lag_time(valid_idx),'.k');
            m1 = median(col_size(valid_idx)');
            m2 = median(lag_time(valid_idx));
    
            hold on;
            line([m1 m1],[0 m2],'Color','r');
            line([0 m1],[m2 m2],'Color','r');
    
            %histogram(lag_time(find(col_flag)),xval2,'FaceColor',colors(i,:),'EdgeColor','none','FaceAlpha',1);
            axis([0 inf 0 inf]);
    
            %t = strcat('Pos:',num2str(ix_tmp(j)),'-',data.metadata.original.Time(ix_tmp(j)),'-(',num2str(length(find(col_flag))),')');
            t = data.metadata.original.Time(ix_tmp(j));
            t = strrep(t,'_',' ');
            title(t,'FontWeight','normal');
            set(gca,'FontSize',8);
        end
    end
    
    figure('Name','split distributions, col size vs early DT');
    for i = 1:length(ix)
        ix_tmp = ix{i};
        for j = 1:length(ix_tmp)
            earlyDT = data.processed{ix_tmp(j)}.colonies.new.early_doublingtime;
            col_size = data.processed{ix_tmp(j)}.colonies.new.timecourse_size_smoothed(end,:);
            col_ecc  = data.processed{ix_tmp(j)}.colonies.region_props.Eccentricity;
            col_flag = data.processed{ix_tmp(j)}.colonies.region_props.flag_colony_ok;
            
            valid_idx = (col_flag == 1 & col_size > size_threshold & col_ecc < ecc_threshold);
    
            subplot(n_groups,max_plate,j+(i-1)*max_plate),...
            plot(col_size(valid_idx)',earlyDT(valid_idx),'.k');
            m1 = median(col_size(valid_idx)');
            m2 = median(earlyDT(valid_idx));
    
            hold on;
            line([m1 m1],[0 m2],'Color','r');
            line([0 m1],[m2 m2],'Color','r');
    
            %histogram(lag_time(find(col_flag)),xval2,'FaceColor',colors(i,:),'EdgeColor','none','FaceAlpha',1);
            axis([0 inf 0 inf]);
    
            %t = strcat('Pos:',num2str(ix_tmp(j)),'-',data.metadata.original.Time(ix_tmp(j)),'-(',num2str(length(find(col_flag))),')');
            t = data.metadata.original.Time(ix_tmp(j));
            t = strrep(t,'_',' ');
            title(t,'FontWeight','normal');
            set(gca,'FontSize',8);
        end
    end
    
    figure('Name','split distributions, lag time vs early DT');
    for i = 1:length(ix)
        ix_tmp = ix{i};
        for j = 1:length(ix_tmp)
            lag_time = data.processed{ix_tmp(j)}.colonies.new.lag_time + incTime;
            earlyDT = data.processed{ix_tmp(j)}.colonies.new.early_doublingtime; 
            col_size = data.processed{ix_tmp(j)}.colonies.new.timecourse_size_smoothed(end,:);
            col_ecc  = data.processed{ix_tmp(j)}.colonies.region_props.Eccentricity;
            col_flag = data.processed{ix_tmp(j)}.colonies.region_props.flag_colony_ok;
            
            valid_idx = (col_flag == 1 & col_size > size_threshold & col_ecc < ecc_threshold);
    
            subplot(n_groups,max_plate,j+(i-1)*max_plate),...
            plot(lag_time(valid_idx),earlyDT(valid_idx),'.k');
            m1 = median(lag_time(valid_idx)');
            m2 = median(earlyDT(valid_idx));
    
            hold on;
            line([m1 m1],[0 m2],'Color','r');
            line([0 m1],[m2 m2],'Color','r');
    
            % histogram(lag_time(find(col_flag)),xval2,'FaceColor',colors(i,:),'EdgeColor','none','FaceAlpha',1);
            axis([0 inf 0 inf]);
    
            %t = strcat('Pos:',num2str(ix_tmp(j)),'-',data.metadata.original.Time(ix_tmp(j)),'-(',num2str(length(find(col_flag))),')');
            t = data.metadata.original.Time(ix_tmp(j));
            t = strrep(t,'_',' ');
            title(t,'FontWeight','normal');
            set(gca,'FontSize',8);
        end
    end
    
    
    for i = 1:length(ix)
        ix_tmp = ix{i};
    
        median_size_tmp = [];
        median_lag_tmp = [];
        median_earlyDT_tmp = [];
        col_nr_tmp = [];
        for j = 1:length(ix_tmp)
            lag_time = data.processed{ix_tmp(j)}.colonies.new.lag_time + incTime;
            col_size = data.processed{ix_tmp(j)}.colonies.new.timecourse_size_smoothed(end,:);
            early_DT = data.processed{ix_tmp(j)}.colonies.new.early_doublingtime;   
            col_flag = data.processed{ix_tmp(j)}.colonies.region_props.flag_colony_ok;
            
            valid_idx = (col_flag == 1 & col_size > size_threshold & col_ecc < ecc_threshold);
    
            median_size_tmp(j) = median(col_size(valid_idx));
            median_lag_tmp(j) = median(lag_time(valid_idx));
            median_earlyDT_tmp(j) = median(early_DT(valid_idx));  
            col_nr_tmp(j) = length(col_flag);
        end
        median_size{i} = median_size_tmp;
        median_lagtime{i} = median_lag_tmp;
        median_earlyDT{i} = median_earlyDT_tmp;
        col_nr{i} = col_nr_tmp;
    end
    
    figure('Name','size and lag time distr comparison');
    for i = 1:length(ix)
        ix_tmp = ix{i};
        for j = 1:length(ix_tmp)
            lag_time = data.processed{ix_tmp(j)}.colonies.new.lag_time + incTime;
            col_size = data.processed{ix_tmp(j)}.colonies.new.timecourse_size_smoothed(end,:);
            col_ecc  = data.processed{ix_tmp(j)}.colonies.region_props.Eccentricity;
            early_DT = data.processed{ix_tmp(j)}.colonies.new.early_doublingtime;  
            col_flag = data.processed{ix_tmp(j)}.colonies.region_props.flag_colony_ok;
            
            valid_idx = (col_flag == 1 & col_size > size_threshold & col_ecc < ecc_threshold);
    
            [N ,plot_edges,~] = histcounts_clean(col_size(valid_idx),xval);
            [N2 ,plot_edges2,~] = histcounts_clean(lag_time(valid_idx),xval2);
            [N3 ,plot_edges3,~] = histcounts_clean(early_DT(valid_idx),xval3);
    
    
            % plot size dist
            subplot(4,3,1+(i-1)*3),...
            %bar(plot_edges,N,'EdgeColor',colors(i,:),'FaceColor','none');
            plot(plot_edges,N,'-k','Color',colors(i,:));
            %histogram(col_size(find(col_flag)),xval,'DisplayStyle','bar','FaceColor',colors(i,:));
            hold on;
            axis([0 inf 0 inf]);
            line([median(col_size(valid_idx)) median(col_size(valid_idx))],[0 max(N)],'Color','b','LineStyle','--');
    
            subplot(4,3,2+(i-1)*3),...
            %bar(plot_edges2,N2,'EdgeColor','none','FaceColor',colors(i,:));
            plot(plot_edges2,N2,'-k','Color',colors(i,:));
            hold on;
            axis([0 inf 0 inf]);
            line([median(lag_time(valid_idx)) median(lag_time(valid_idx))],[0 max(N2)],'Color','b','LineStyle','--');
    
            subplot(4,3,3+(i-1)*3),...
            plot(plot_edges3,N3,'-k','Color',colors(i,:));
            %bar(plot_edges3,N3,'EdgeColor','none','FaceColor',colors(i,:));
            hold on;
            axis([0 inf 0 inf]);
            line([median(early_DT(valid_idx)) median(early_DT(valid_idx))],[0 max(N3)],'Color','b','LineStyle','--');
        end
    end
    
    
    
    %% plot overarching summary
    figure('Name','extracted params-median');
    for i = 1:length(median_size)
        median_size_tmp = median_size{i};
        median_lag_tmp = median_lagtime{i};
        median_earlyDT_tmp = median_earlyDT{i};
    
        col_nr_tmp = col_nr{i};
    
        x_wiggle = (0.5-rand(size(median_size_tmp)))./5;
        for j = 1:length(median_size_tmp)
    
            subplot(1,4,1),...
            plot(i+x_wiggle,median_size_tmp,'ok','Color',colors(i,:));
            hold on;
    
            subplot(1,4,2),...
            plot(i+x_wiggle,median_lag_tmp,'ok','Color',colors(i,:));
            hold on;
    
            subplot(1,4,3),...
            plot(i+x_wiggle,median_earlyDT_tmp,'ok','Color',colors(i,:));
            hold on;
    
            subplot(1,4,4),...
            plot(i+x_wiggle,col_nr_tmp,'ok','Color',colors(i,:));
            hold on;
        end
    end
    
    subplot(1,4,1),...
    set(gca,'XTick',[1:5],'XTickLabel',labels_type);
    xtickangle(-45);
    axis([0.5 n_groups+0.5 0 inf]);
    ylabel('median col size [px]');
    subplot(1,4,2),...
    set(gca,'XTick',[1:5],'XTickLabel',labels_type);
    xtickangle(-45);
    ylabel('median lag time [h]');
    axis([0.5 n_groups+0.5 0 inf]);
    subplot(1,4,3),...
    set(gca,'XTick',[1:5],'XTickLabel',labels_type);
    xtickangle(-45);
    ylabel('early DT [h]');
    axis([0.5 n_groups+0.5 0 inf]);
    subplot(1,4,4),...
    set(gca,'XTick',[1:5],'XTickLabel',labels_type);
    xtickangle(-45);
    ylabel('colony nr [-]');
    axis([0.5 n_groups+0.5 0 inf]);    
    
    % calculate p-values
    [~,p(1,1)] = ttest2(median_size{1},median_size{1});
    [~,p(1,2)] = ttest2(median_size{2},median_size{2});
    [~,p(1,3)] = ttest2(median_size{3},median_size{3});
    
    [~,p(2,1)] = ttest2(median_lagtime{1},median_lagtime{1});
    [~,p(2,2)] = ttest2(median_lagtime{2},median_lagtime{2});
    [~,p(2,3)] = ttest2(median_lagtime{3},median_lagtime{3});
    
    [~,p(3,1)] = ttest2(col_nr{1},col_nr{1});
    [~,p(3,2)] = ttest2(col_nr{2},col_nr{2});
    [~,p(3,3)] = ttest2(col_nr{3},col_nr{3});
    
    [~,p(4,1)] = ttest2(median_earlyDT{1},median_earlyDT{1});
    [~,p(4,2)] = ttest2(median_earlyDT{2},median_earlyDT{2});
    [~,p(4,3)] = ttest2(median_earlyDT{3},median_earlyDT{3});
    
    %% 
    data.across_sample_types.median_size = median_size;
    data.across_sample_types.median_lagtime = median_lagtime;
    data.across_sample_types.median_earlyDT = median_earlyDT;
    data.across_sample_types.col_nr = col_nr;
    data.across_sample_types.sample_types = labels_type;
    
    data.across_sample_types.p_values = p;
    data.across_sample_types.p_values_columns = {'1:MLN','2:Spleen','3:Liver'};
    data.across_sample_types.p_values_rows = {'1:median size','2:median lag time','3:colony number','4:median early DT'};

end
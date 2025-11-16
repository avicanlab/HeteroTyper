%% 11.11.2025 - HeteroTyper Pipeline for Bright Plates 

% plot summary of growth data

function plot_QC_growth_data(data)

    nr_plates = length(data.processed);
    
    % Define maximum lag time
    max_lag = 52;

    % Define Room Temperature Incubation Time
    incTime = 20;
    
    % Define min-max colony numbers
    min_col = 10;
    max_col = 650;

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
    
    % depending on the number of plates we will use different layouts
    if(length(ix_growth)<= 24)
        nr_rows = 4;
        nr_colums = 6;
    elseif(length(ix_growth)<= 30)
        nr_rows = 5;
        nr_colums = 6;
    elseif(length(ix_growth)<= 48)
        nr_rows = 6;
        nr_colums = 8;
    else
        nr_rows = 10;
        nr_colums = 11;
    end
    
    
    %% plot growth curves
    figure('Name','QC - all growth curves');
    max_val = 4100;
    
    for i = 1:length(ix_growth)
        subplot(nr_rows,nr_colums,i),...
            
        t1 = data.processed{ix_growth(i)}.colonies.new.time_info.elapsed_time_h + incTime;
        col_size = data.processed{ix_growth(i)}.colonies.new.timecourse_size_smoothed;
        col_flag = data.processed{ix_growth(i)}.colonies.region_props.flag_colony_ok;
       
        plot(t1,col_size(:,find(col_flag)),'-');
        hold on;
        axis([0 inf 0 inf]);

        text(0,0.5*max(col_size(end,:)),num2str(length(find(col_flag))));
        t = data.metadata.fn(ix_growth(i)).name;
        t = strrep(t,'_',' ');
        title(t,'FontWeight','normal');
        set(gca,'FontSize',7); % size was 6
    end
    
    %% plot distributions of final cell size, lag time, and early DT
    figure('Name','QC - size distr');
    xval = [0:100:max_val];
    
    for i = 1:length(ix_growth)
        subplot(nr_rows,nr_colums,i),...
            
        t = data.processed{ix_growth(i)}.colonies.new.time_info.elapsed_time_h + incTime;
        col_size = data.processed{ix_growth(i)}.colonies.new.timecourse_size_smoothed(end,:);
        col_flag = data.processed{ix_growth(i)}.colonies.region_props.flag_colony_ok;
        
        %[N,~] = histcounts(col_size(find(col_flag)),xval);
        [N ,plot_edges,~] = histcounts_clean(col_size(find(col_flag)),xval);
        
        bar(plot_edges,N,'FaceColor','k');
        axis([0 inf 0 inf]);
        hold on;
        line([median(col_size(find(col_flag))) median(col_size(find(col_flag)))],[0 max(N)],'Color','r','LineWidth',2);
        
        text(0,0.8.*max(N),num2str(length(find(col_flag))));
        t = data.metadata.fn(ix_growth(i)).name;
        t = strrep(t,'_',' ');
        title(t,'FontWeight','normal');
        set(gca,'FontSize',8);
    end
    
    figure('Name','QC - lag time distr');
    xval2 = [incTime:0.5:max_lag];
    
    for i = 1:length(ix_growth)
        subplot(nr_rows,nr_colums,i),...
        
        lag_time = data.processed{ix_growth(i)}.colonies.new.lag_time + incTime;
        col_flag = data.processed{ix_growth(i)}.colonies.region_props.flag_colony_ok;
        
        %[N,~] = histcounts(lag_time(find(col_flag)),xval2);
        [N ,plot_edges,~] = histcounts_clean(lag_time(find(col_flag)),xval2);
        bar(plot_edges,N,'FaceColor','k');
        axis([incTime inf 0 inf]);
        hold on;
        line([median(lag_time(find(col_flag))) median(lag_time(find(col_flag)))],[0 max(N)],'Color','r','LineWidth',2);
        
        text(0,0.8.*max(N),num2str(length(find(col_flag))));
        t = data.metadata.fn(ix_growth(i)).name;
        t = strrep(t,'_',' ');
        title(t,'FontWeight','normal');
        set(gca,'FontSize',8);
    end
    
    figure('Name','QC - early DT distr');
    xval3 = [0:0.5:max_lag];
    
    for i = 1:length(ix_growth)
        subplot(nr_rows,nr_colums,i),...
        
        early_DT = data.processed{ix_growth(i)}.colonies.new.early_doublingtime;
        col_flag = data.processed{ix_growth(i)}.colonies.region_props.flag_colony_ok;
        
        %[N,~] = histcounts(lag_time(find(col_flag)),xval2);
        [N ,plot_edges,~] = histcounts_clean(early_DT(find(col_flag)),xval3);
        bar(plot_edges,N,1,'FaceColor','k');
        axis([-inf inf 0 inf]);
        hold on;
        line([median(early_DT(find(col_flag))) median(early_DT(find(col_flag)))],[0 max(N)],'Color','r','LineWidth',2);
        
        text(0,0.8.*max(N),num2str(length(find(col_flag))));
        t = data.metadata.fn(ix_growth(i)).name;
        t = strrep(t,'_',' ');
        title(t,'FontWeight','normal');
        set(gca,'FontSize',8);
    end

end
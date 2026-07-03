%% 11.11.2025 - HeteroTyper Pipeline for Bright Plates

function data = extract_early_growth(data)

nr_plates = length(data.processed);
mue_mult = data.params.early_DT_range;

for i = 1:nr_plates
    if(data.processed{i}.growth_quant == 1)
        col_size = data.processed{i}.colonies.new.timecourse_size_smoothed;
        time = data.processed{i}.colonies.new.time_info.elapsed_time_h;
        
        lag_time_ix = [];
        lag_time_ix_2 = [];
        
        for j = 1:size(col_size,2)
            tmp = find(col_size(:,j)>data.params.lag_time_thresh,1);
            tmp2 = find(col_size(:,j)>data.params.lag_time_thresh*mue_mult,1);
            
            if(~isempty(tmp))
                lag_time_ix(j) = tmp;
            else
                lag_time_ix(j) = length(time);
            end
            if(~isempty(tmp2))
                lag_time_ix_2(j) = tmp2;
            else
                lag_time_ix_2(j) = length(time);
            end
        end
        delta_tmp = time(lag_time_ix_2) - time(lag_time_ix);
        
        data.processed{i}.colonies.new.early_doublingtime = delta_tmp;
    end
end

end
%% 11.11.2025 - HeteroTyper Pipeline for Bright Plates

% extract actual time information from filenames

function out = extract_time_info(fn)

tp_info = zeros(size(fn,1),6);
tp_datenum = zeros(size(fn,1),1);

ix = [18 14 12 9 7 5]; % hard-coded position of y,m,d,hour,min,sec
for i = 1:length(fn)
    tmp = fn{i};
    
    for j = 1:length(ix)
        
        if(j == 1)
            t1 = tmp((end-ix(j)):(end-ix(j)+3));
        else
            t1 = tmp((end-ix(j)):(end-ix(j)+1));
        end
        tp_info(i,j) = str2num(t1);
    end
    tp_datenum(i,1) = datenum(tp_info(i,:));
    
end

%% 
out.date_info = tp_info;
out.date_info_ids = {'y','m','d','hour','min','sec'};
out.elapsed_time_h = (tp_datenum - tp_datenum(1))*24;

end
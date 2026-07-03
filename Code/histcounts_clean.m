%% 11.11.2025 - HeteroTyper Pipeline for Bright Plates

function [N ,plot_edges,label] = histcounts_clean(Y,bin_edges)


for i = 1:length(bin_edges)
    
    if(i == 0) % 
        ix = find(Y < bin_edges(1));
        plot_edges(i,1) = bin_edges(i) - (bin_edges(i+1)-bin_edges(i))./2;
        N(i,1) = length(ix);
        label(ix,1) = i;
    elseif(i>0)&&(i<length(bin_edges))
        ix1 = find(Y >= bin_edges(i));
        ix2 = find(Y < bin_edges(i+1));
        ix = intersect(ix1,ix2);
        
        plot_edges(i+1,1) = bin_edges(i) + (bin_edges(i+1)-bin_edges(i))./2;
        N(i+1,1) = length(ix);
        label(ix,1) = i+1;
        
    elseif(i == length(bin_edges))
        ix = find(Y >= bin_edges(i));
        N(i+1,1) = length(ix);
        label(ix,1) = i+1;
        plot_edges(i+1,1) = bin_edges(i) + (bin_edges(i)-bin_edges(i-1))./2;
    end
end


% ix = find(Y < bin_edges(1));
% plot_edges_tmp = bin_edges(1) - (bin_edges(2)-bin_edges(1))./2;
% N(i,1) = length(ix);
% label(ix,1) = i;
end
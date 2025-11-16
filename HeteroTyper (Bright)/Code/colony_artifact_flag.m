%% 11.11.2025 - HeteroTyper Pipeline for Bright Plates

function [colony_flag,mask_clean] = colony_artifact_flag(img,colonies,radius,centers,params,plot_flag)

    %% criterion 2: touches border region
    % Accept either 'mask' (old name) or 'mask_filled' (new name)
    if isfield(colonies, 'mask')
        m = colonies.mask;
    elseif isfield(colonies, 'mask_filled')
        m = colonies.mask_filled;
    else
        error('colony_artifact_flag: colonies structure does not contain mask or mask_filled field.');
    end
    
    
    % Pre-allocate colony_flag with the number of region_props entries (if available)
    nRegions = 0;
    if isfield(colonies,'region_props') && ~isempty(colonies.region_props)
        nRegions = height(colonies.region_props);
    end
    colony_flag = zeros(nRegions,1);
    
    
    [rows cols] = meshgrid(1:size(m,2),1:size(m,1));
    circle_pxl = (rows - centers(1)).^2 + (cols - centers(2)).^2 <= round(radius).^2;  
    %exclude 
    m2 = m;
    m2(m2>0) = 1;
    m2(~circle_pxl) = 1;
    m2 = imclearborder(m2); 
    m(m2==0) = 0;
    ix_include_2 = double(unique(m));
    ix_include_2 = ix_include_2(ix_include_2 > 0);
    
    %% merge two lists
    % ix_include_merged = intersect(ix_include_1,ix_include_2);
    colony_flag(ix_include_2) = 1;   % colony_flag(ix_include_merged) = 1;
    
    %% set output to results
    m_cleaned = m;
    ix_colony_not_ok = find(colony_flag == 0);
    
    for i = 1:length(ix_colony_not_ok)
        m_cleaned(m_cleaned==ix_colony_not_ok(i)) = 0;
    end
    mask_clean = m_cleaned;
    
    % visualize results
    if(plot_flag == 1)
        m_cleaned_plot = mask_clean;
        m_cleaned_plot(m_cleaned_plot>0) = 1;
        m_original_plot = m;
        m_original_plot(m_original_plot>0) = 1;
    
        figure('Name','colony cleanup: colonies in purple flagged as outliers');
        imshowpair(m_cleaned_plot,m_original_plot,'falsecolor');
    end
end

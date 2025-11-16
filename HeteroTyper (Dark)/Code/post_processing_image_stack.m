%% 11.11.2025 - HeteroTyper Pipeline for Dark Plates
    
function out = post_processing_image_stack(inp, params, plot_flag)

%% 1. try to find center of plate
start_img = inp.img{1};

center_pos = find_plate_center(start_img);
max_radius = size(start_img,1)./2;
rel_radius = max_radius * params.border_range;

% check if the found center_pos is off compared to the median --> if yes,
% set to median
if(abs(center_pos(1) - params.center_median(1))./max_radius >= 0.1)
    center_pos = params.center_median;
end

% if(plot_flag == 1)
%     % plot original image
%     figure('Name',strcat(inp.position_name,'-timepoint 1 with identified center'));
%     imshow(start_img);
%     hold on;
%     viscircles(center_pos,max_radius,'Color','b');
%     viscircles(center_pos,rel_radius,'Color','r');
% end


out.center_pos = center_pos;
out.radius = rel_radius;

%% 2. image postprocessing
fprintf(strcat('processing...',inp.position_name,'\n'));

for i = 1:size(inp.img,1)
    img_tmp = inp.img{i};
    img_tmp = imsubtract(img_tmp,start_img);
    
    img_tmp2 = img_post_processing(img_tmp,center_pos,rel_radius);
    
    out.img{i,1} = img_tmp2;
end

% if(plot_flag == 1)
%     figure('Name',strcat(inp.position_name,'-final (purple) vs first (green) timepoint after processing'));
%     imshowpair(out.img{1},out.img{end},'falsecolor');
%     hold on;
%     viscircles(center_pos,rel_radius,'Color','r');
% end
fprintf(strcat('processing...',inp.position_name,':done\n'));

end

%% This function attempts to find plate center
function center_pos = find_plate_center(img)

img_gray = rgb2gray(img);
img_binarized = ~imbinarize(img_gray);
im_filled = imfill(img_binarized,'holes'); 
im_filled = bwareaopen(im_filled,100000); % remove all small objects
% find all objects, identify the biggest one (presumably the plate)
rp = regionprops(im_filled);
[~,ix] = max([rp.Area]);

center_pos = rp(ix).Centroid;

end

%% This function performs the actual post-processing
function img_out = img_post_processing(img_in,centers,radius)

img_out = img_in;
% Step 1: crop border regions
[rows cols] = meshgrid(1:size(img_in,2),1:size(img_in,1));
circle_pxl = (rows - centers(1)).^2 + (cols - centers(2)).^2 <= round(radius).^2;  

for i = 1:3
    tmp = img_out(:,:,i);
    
    tmp(~circle_pxl) = 0;
    img_out(:,:,i) = tmp;
end


end
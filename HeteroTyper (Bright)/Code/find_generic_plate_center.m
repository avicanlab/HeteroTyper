%% 11.11.2025 - HeteroTyper Pipeline for Bright Plates

function data = find_generic_plate_center(data)

fn = data.metadata.fn;
nr_plates = 4;  % Must be changed to the number of samples
%% load first image of every first image (for first 35 images)
for i = 1:nr_plates
    %i
    
    t1 = strcat(fn(i).folder,'\',fn(i).name); % get pointer
    
    fn_tmp = dir(t1);
    fn_tmp = fn_tmp(4,:);
    
    t2 = strcat(fn_tmp.folder,'\',fn_tmp.name);
    img{i,1} = imread(t2);
    
    center_pos(i,:) = find_center_individual(img{i,1});
end

center_median = median(center_pos,1);
radius = size(img{1},1)./2;

data.params.center_median = center_median;
data.params.radius = radius;

%% plot results
figure('Name','find plate centers');
for i = 1:nr_plates
    subplot(4,6,i),...  % CHANGE
    imshow(img{i});
    hold on;    
    viscircles(center_median,radius*0.9,'Color','b');
    
    if(abs(center_pos(i,1)-center_median(1))./radius > 0.1)
        viscircles(center_pos(i,:),radius*0.9,'Color','r');
    else
        viscircles(center_pos(i,:),radius*0.9,'Color','g');
    end
end

end


%% This function attempts to find plate center
function center_pos = find_center_individual(img)

img_gray = rgb2gray(img);
img_binarized = ~imbinarize(img_gray);
im_filled = imfill(img_binarized,'holes'); 
im_filled = bwareaopen(im_filled,100000); % remove all small objects
% find all objects, identify the biggest one (presumably the plate)
rp = regionprops(im_filled);
[~,ix] = max([rp.Area]);

center_pos = rp(ix).Centroid;

end
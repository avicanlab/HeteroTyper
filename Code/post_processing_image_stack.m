%% 11.11.2025 - HeteroTyper Pipeline for Bright Plates

% post-processing of image stack
        
function out = post_processing_image_stack(inp, params, plot_flag)

%% 1. Use per-plate center and radius from params (set by find_generic_plate_center)
%
% Previously this function ran its own find_plate_center() internally and
% used image_height/2 as the radius.  Both are now replaced by the
% per-plate median values computed across the image stack in
% find_generic_plate_center, which are injected as plate_center_current
% and plate_r_current before this function is called.
%
% Fallback chain (should not normally be needed):
%   1. plate_center_current / plate_r_current  (per-plate median — preferred)
%   2. center_median(1,:) / radius(1)          (first row, last resort)
%   3. image centre / image_height/2           (blind fallback)

start_img = inp.img{1};
[img_h, img_w, ~] = size(start_img);

if isfield(params, 'plate_center_current') && ~any(isnan(params.plate_center_current)) ...
        && isfield(params, 'plate_r_current') && ~isnan(params.plate_r_current) && params.plate_r_current > 0
    center_pos = params.plate_center_current;          % [cx, cy]
    max_radius = params.plate_r_current;
else
    % Last-resort fallback — should not occur if main script is correct
    warning('post_processing_image_stack: plate_center_current not set, using image centre.');
    center_pos = [img_w/2, img_h/2];
    max_radius = img_h / 2;
end

rel_radius = max_radius * params.border_range;

out.center_pos = center_pos;
out.radius     = rel_radius;

%% 2. Image postprocessing
fprintf(strcat('processing...', inp.position_name, '\n'));

for i = 1:size(inp.img, 1)
    img_tmp  = inp.img{i};
    img_tmp  = imsubtract(img_tmp, start_img);
    img_tmp2 = img_post_processing(img_tmp, center_pos, rel_radius);
    out.img{i,1} = img_tmp2;
end

fprintf(strcat('processing...', inp.position_name, ':done\n'));

end


%% Crop image to plate disk — pixels outside the circle set to zero
function img_out = img_post_processing(img_in, centers, radius)

img_out = img_in;
[rows, cols] = meshgrid(1:size(img_in,2), 1:size(img_in,1));
circle_pxl = (rows - centers(1)).^2 + (cols - centers(2)).^2 <= round(radius).^2;

for i = 1:3
    tmp = img_out(:,:,i);
    tmp(~circle_pxl) = 0;
    img_out(:,:,i) = tmp;
end

end
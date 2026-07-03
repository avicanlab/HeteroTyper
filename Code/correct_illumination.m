%% 11.11.2025 - HeteroTyper Pipeline for Bright Plates

% correct for drifts in illumination

function out = correct_illumination(inp,params)

%% get plate center:
centers = inp.processed.center_pos;
radius = inp.processed.radius;

[rows cols] = meshgrid(1:size(inp.processed.img{1},2),1:size(inp.processed.img{1},1));
circle_pxl = (rows - centers(1)).^2 + (cols - centers(2)).^2 <= round(radius).^2;  
circle_pxl(find(inp.processed.colonies.mask_clean)) = 0;

for i = 1:length(inp.raw.img)
    x1 = rgb2gray(inp.raw.img{i});
    x2 = x1(circle_pxl);
    out.median_signal(i,1) = mean(double(x2));
    out.mean_signal(i,1) = mean(double(x2));
    out.percentile_signal(i,1) = prctile(double(x2),75);
end

figure;
subplot(3,1,1),...
    plot(inp.processed.time.elapsed_time_h,out.median_signal-out.median_signal(1),'-ok');
axis([0 25 -20 20]);
subplot(3,1,2),...
    plot(inp.processed.time.elapsed_time_h,out.mean_signal-out.mean_signal(1),'-ok');
axis([0 25 -20 20]);

subplot(3,1,3),...
    plot(inp.processed.time.elapsed_time_h,out.percentile_signal-out.percentile_signal(1),'-ok');
axis([0 25 -20 20]);

end
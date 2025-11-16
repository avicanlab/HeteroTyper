%% 11.11.2025 - HeteroTyper Pipeline for Bright Plates

% load all images within input folder

function out = load_image_stack(inp,position_name)

fn = dir(inp);
fn = fn(3:end,:);

fprintf(strcat('loading...',position_name,'\n'));

for i = 1:length(fn)
    out.filename{i,1} = fn(i).name;
    t = strcat(inp,'\',fn(i).name);
    out.img{i,1} = imread(t);
end
out.position_name = position_name;

end
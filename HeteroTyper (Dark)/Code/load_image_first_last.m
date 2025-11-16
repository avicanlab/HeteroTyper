%% 11.11.2025 - HeteroTyper Pipeline for Dark Plates

% load first and last image from each image stack

function out = load_image_first_last(inp,position_name)

fn = dir(inp);
fn = fn(3:end,:);

fprintf(strcat('loading...',position_name,'\n'));

id_first_last = [1 length(fn)];

for i = 1:length(id_first_last)
    out.filename{i,1} = fn(id_first_last(i)).name;
    t = strcat(inp,'\',fn(id_first_last(i)).name);
    out.img{i,1} = imread(t);
end
out.position_name = position_name;

end
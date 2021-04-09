function add_path()
path0 = mfilename('fullpath');
path_num = strfind(path0, '\');
path_str = path0(1 : path_num(end));

new_path = ['\Proposed_schemes', '\Basic', '\Hy_BD', '\Two_stage'];
new_path_num = strfind([new_path '\'], '\');
for path_i = 1 : length(new_path_num) - 1
    path(path, [path_str new_path(new_path_num(path_i) : new_path_num(path_i + 1) - 1)]);
end

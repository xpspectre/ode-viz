function init_ode_viz
% Add needed paths to use this library
ode_viz_path = fileparts(mfilename('fullpath'));
addpath(ode_viz_path);
addpath(fullfile(ode_viz_path, 'jsonlab-1.5'));
end


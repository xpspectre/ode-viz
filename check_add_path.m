function added = check_add_path(directory)
% Check of folder is on path and adds it if it isn't already
%
% Inputs:
%   directory [ string ]
%       Directory to check and add to path
%
% Outputs:
%   added [ true | false ]
%       Whether the directory was added (true if directory wasn't already
%       on the path)

pathCell = regexp(path, pathsep, 'split');
if ispc  % Windows is not case-sensitive
  onPath = any(strcmpi(directory, pathCell));
else
  onPath = any(strcmp(directory, pathCell));
end

if onPath
    added = false;
    return
else
    addpath(get_full_path(directory));
    added = true;
end

end


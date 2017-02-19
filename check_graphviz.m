function isPresent = check_graphviz()
% Detect whether Graphviz executbles are present
%
% Outputs:
%   isPresent [ true | false ]
%       Returns true if OK

exes = {'dot', 'neato'};
for i = 1:length(exes)
    exe = exes{i};
    [status, ~] = system(['where ' exe]);
    if status == 1
        isPresent = false;
        return
    end
end
isPresent = true;
end


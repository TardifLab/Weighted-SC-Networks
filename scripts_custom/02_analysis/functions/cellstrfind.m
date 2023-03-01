function outi = cellstrfind(targ,str)
% Searches cell arrays for a given string occurrence. Returns indices of
% cells in which target string found.
%
% Input:
%   targ    : Cell array to search for string
%   str     : target string
%
% Output:
%   outi    : indices of string
%
% 2021 Mark C Nelson, McConnell Brain Imaging Centre, MNI, McGill
%--------------------------------------------------------------------------
opt1='ignorecase';
opt2='UniformOutput';

if ischar(str)                                                              % search input is string vector
    tmp                 = cellfun(@(c) regexp(c,str,opt1),targ,opt2,0);
    outi                = find(cellfun(@(c) ~isempty(c),tmp));
elseif iscell(str)                                                          % search input may be a cell array of strings
    outi                = cell(size(str));
    for ii = 1:numel(str)
        if ischar(str{ii})
            tmp         = cellfun(@(c) regexp(c,str{ii},opt1),targ,opt2,0);
            outi{ii}    = find(cellfun(@(c) ~isempty(c),tmp));
        end
    end
else
    error('Unknown input data type: CELLSTRFIND accepts strings or cell arrays of strings')
end
        
%--------------------------------------------------------------------------
end
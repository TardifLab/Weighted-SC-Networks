function [offset] = dirOffset(sourceDir, skipMe)
%  This function accepts any directory & determines an offset value at which
%  any irrelevant items in directory are skipped. 
%
% NOTE: Irrelevant files must be manually indicated by skipMe & should be
% found at the top of the directory.
%
% Input:
%   sourceDir   : Target directory
%   skipMe      : Things to skip in directory
%
% Output:
%   Offset      : Start point of relevant objects in dir
%
% 2021 Mark C Nelson, McConnell Brain Imaging Centre, MNI, McGill
%--------------------------------------------------------------------------

% Get target dir structure
dirTrg          = dir(sourceDir);

% Iterate through names of dir components to find 1st relevant file (could
% be modified to index 
for ii = 1 : length(dirTrg)
    if any(strcmp( dirTrg(ii).name, skipMe ))
    else
        offset = ii;
        break
    end
end
%--------------------------------------------------------------------------
end
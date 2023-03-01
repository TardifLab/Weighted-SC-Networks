function [D,labels] = conn_selectData(D,labels,select,dim)
%  Selects the indicated row(s) or column(s) from the overall 
%  W x P data array. 
%
% Input:
%   D       : IxJ cell array of stacked subject-level connectomes
%   labels  : Cell array of data labels (must match one of the dimensions
%             of the data array D)
%   select  : cell array of strings indicating data to select
%   dim     : dim of data array for selection (1=rows,2=cols)
%
% Output:
%   D       : Reduced data array
%   labels  : updated labels (i.e. select)
%
% 2021 Mark C Nelson, McConnell Brain Imaging Centre, MNI, McGill
%--------------------------------------------------------------------------
%% Check inputs
if ischar(select)
    select                  = {select};
end
if isempty(select)
    error('No label for selection given!')
else
    indices                 = cellstrfind(labels,select);
end
if all(cellfun(@(c) isempty(c),indices))
    if all(cellfun(@(c) isempty(c),select))
        disp('* *  ------------- All labels retained ------------- * *')
        disp('          (no input given to conn_selectData)           ')
        return
    else
        error('Unable to select your data: Check your SELECT input to conn_selectData')
    end
end

%% Deal with ambiguous labels

% Manual list of commonly mixed up labels. Add your problem labels to this
ambiglbls                   = {'schaefer-100' 'schaefer-1000' 'ad' 'adc' 'COMMIT' 'icvf-commit' 'aparc' 'aparc-a2009s' 'icvf-commit' 'icvf-noddi'};

% Check for labels in list and use more rigorous search criteria to differentiate
if any(cellfun(@(c) numel(c),indices)>1) && any(ismember(select,ambiglbls) | ismember(strrep(select,' ','-'),ambiglbls))                       
    ambigi                  = find(cellfun(@(c) numel(c), indices) > 1);    % indices of labels returning multiple matches
    for aa = 1 : numel(ambigi)
        it                  = ambigi(aa);
        lblt                = select(it);
        if ismember(lblt,ambiglbls) || ismember(strrep(lblt,' ','-'),ambiglbls)
            indices{it}     = find(cellfun(@(c) strcmpi(lblt,c),labels));   % more rigid-search criteria
            if length(indices{it}) > 1
                disp('* * * WARNING: Unable to resolve ambiguity!')
                keyboard                                                    
            end
        end
    end
end

%% Select the data corresponding to these labels
disp(['Selecting: ' select])
SXi                         = indices;

% Convert SXi to array and deal with possible multiple indices per-label (sometimes multiple indices is desirable)
if any(cellfun(@(c) numel(c), SXi) > 1)
    tmpsxi                  = [];
    for ii = 1 : numel(SXi)
        tmpsxi              = [tmpsxi, SXi{ii}];
    end
    SXi                     = tmpsxi;
else
    SXi                     = cell2mat(SXi);
end
SXi                         = unique(SXi);                                  % deal with possibility of duplicates & sort

% Select data
if ~any(isnan(SXi))
    RMi                     = setdiff(1:numel(labels), SXi);
    labels(RMi)             = [];
    switch dim
        case 1
            D(RMi,:)        = [];
        case 2
            D(:,RMi)        = [];
        otherwise
            error('Requires valid input for dimension')
    end
else
    disp('ERROR: one or more labels not found. Check variable: select')
end
%--------------------------------------------------------------------------
end
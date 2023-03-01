function D = conn_normalize(D,meth)
% Normalizes all datasets in the cell array D according to METH
%
% Input:
%   D       : IxJ cell array, each cell containing an NxNxS stack of
%             2D connectivity matrices
%   meth    : 1x2 cell array with contents:
%               1. Char vector indicating scaling method. Options include:
%               'normalize' rescale all weight magnitudes by largest
%                           (sign preserved)
%               'range' data scaled to indicated range
%
%               2. Necessary inputs (depends on method)
%
%           e.g. meth = {'range', [0 1]}
%
% Output:
%   D       : Cell array reflecting desired modification   
%
% 2021 Mark C Nelson, McConnell Brain Imaging Centre, MNI, McGill
%--------------------------------------------------------------------------
[I,J]               = size(D);
if ~iscell(meth); meth = {meth}; end
methname            = meth{1};

%% --------------------------- Normalize data -------------------------- %%
switch lower(methname)
    
    case 'normalize'                                                        % scale all weight magnitudes by largest, preserve sign
        
        for ii = 1 : I
            for jj = 1 : J
                dt              = D{ii,jj};                                 % stack of subject-level connectomes
%                 dt             = groupavg(dt,3);                          % Option to use Group max (may lead to some values > 1)
                xmax            = max(abs(dt(dt~=0)));                      % Max of pooled subject data
                Sn              = size(dt,3);                               % number of subjects
                for kk = 1 : Sn
                    dts         = dt(:,:,kk);                               % single subjects data
                    dts         = dts./xmax;                                % rescale
                    dt(:,:,kk)  = dts;
                end
                D{ii,jj}        = dt;
            end
        end
    
    case 'range'                                                            % min-max scaling to indicated range

        [a,b]                   = deal(meth{2}(1), meth{2}(2));
        for ii = 1 : I
            for jj = 1 : J
                dt              = D{ii,jj};                                 % stack of subject-level connectomes
                Sn              = size(dt,3);                               % number of subjects
                xmin            = min(dt(dt~=0));                           % pooled subject min (want absolute minimum)
                dtg             = groupavg(dt,3);
                xmax            = max(dtg(dtg~=0));                         % group average max (to avoid spuriously high values)
                dt(dt==0)       = nan;                                      % to allow preservation of 0s
                for kk = 1 : Sn
                    dts         = dt(:,:,kk);                               % single subjects data
                    it          = find(triu(dts)~=0);                       % get upper tri indices
                    dtn         = nan(length(dts));                         % empty matrix
                    
                    % scale to range (vectorized & subject-level much faster!)
                    dtn(it)     = arrayfun(@(x)  a+(((x-xmin)*(b-a))/(xmax-xmin)), dts(it));
                    dt(:,:,kk)  = dtn;
                end
                dt(isnan(dt))   = 0;                                        % retain original 0s
                D{ii,jj}        = dt;
            end
        end
        
    otherwise
        error('ERROR: unknown normalization option!')
end
%--------------------------------------------------------------------------
end
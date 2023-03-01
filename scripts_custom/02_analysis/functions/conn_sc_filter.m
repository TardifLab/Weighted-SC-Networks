function [D] = conn_sc_filter(D,W8s,filtMeth,varargin)
%  Nulls out edges in connectomes based on input criteria. Works on
%  structural connectomes only unless otherwise indicated (FCrm).
%
%                       + Filtering methods +
%   sift2       : edges with SIFT2<1 removed 
%                 * interpreted as < 1 streamline
%                 * DO NOT USE IF NORMALIZING BY NODE VOLUME!
%                 * Smith R. et al., NeuroImage. 2015
%   commit      : edges with COMMIT weight <1e-12 removed from all SC 
%                 * machine precision 0
%                 * Daducci A. et al., IEEE. 2015
%   nz          : edges = 0 in any subject removed from all subjects
%   uniform     : only edges present in > a given threshold of subjects 
%                 are retained 
%
% Input:
%                           + Required +
%   D           : IxJ cell array, each cell containing an NxNxS stack of
%                 connectivity matrices
%   W8s         : Cell array of strings describing edge weights
%   filtMeth    : character vector indicating filter method (switch) 
%                 e.g. {'sift2' 'commit' 'nz' 'uniform'}
%
%                           + Optional +
%   Uthr        : uniform threshold; numeric value [1 99]
%   filtALL     : []/1 option to apply filter to all structural connectomes
%                 (Cant be used if multiple filtering methods given)
%                 (1=default)
%   FCrm        : []/1 option to filter FC data too (1=filter, default=[])
%
% Output:
%   D           : filtered data
%
%
% Note: nz & uniform filtering methods are relative to dataset i.e.
%       densities of filtered networks may not match across weights. 
%       This could be acheived in a 2-stage filtering process e.g.
%           1. Cs = conn_sc_filter(Cs,Ws,'commit')
%           2. Cs = conn_sc_filter(Cs,Ws,'nz')
%
%
% 2021 Mark C Nelson, McConnell Brain Imaging Centre, MNI, McGill
%--------------------------------------------------------------------------

%% Optional inputs
nin=max(nargin,1) - 3; 
defaults                = {[] 1 'cor' []};
[Uthr,filtALL,FCrm]     = INhandler(varargin,nin,defaults);

%% Setup
[I,J]                       = size(D);
stradd                      = '';

if isempty(FCrm)
    FCrm                    = 0;
elseif FCrm==1
    stradd                  = '& FUNCTIONAL ';
end
if ~isempty(filtALL)
    if ~isnumeric(filtALL) || filtALL~=0 && filtALL~=1
        error('** FILTER FAIL: filtALL should be: [], 0, or 1')             % Acceptable values for filtALL: {0 1 []}
    end
end

% filtMeth should be a character vector
if ~ischar(filtMeth)
    error('** FILTER FAIL: Check your input to filtMeth')
end

% W8s should be a 1D cell array of strings
if ~iscell(W8s)
    W8s                     = {W8s};                       
end
if any(size(W8s) == 1)                                                      % should be 1D
    % Check contents of cells
    if any(cellfun(@(c) ~ischar(c),W8s))
        error('** FILTER FAIL: W8s input should contain character labels!')
    end
else
    error('** FILTER FAIL: W8s input should be 1D!')
end

iw8                         = cellstrfind(W8s,filtMeth);
if isempty(iw8); iw8        = 1; end                                        % allows Uthr & nz filtering to proceed, no effect             

%% Filter Option 1

% Only filter indicated dataset(s) i.e. do not apply to all
if isempty(filtALL)
    
    disp(['Filtering: ' W8s{iw8(kk)} ' only as ' filtMeth])
    keyboard                                                                % This section is old, not sure if it still works
%     for jj = 1 : J
%         dselt                               = D{iw8(kk),jj};
%         switch lower(filtMeth)
%             case 'sift2'
%                 disp('Setting all edges <1 to 0 in SIFT2 weighted connectomes only')
%                 dselt(dselt<1)              = 0;                            % edges with streamline count < 1 --> go to 0
%             case 'commit'
%                 disp('No change to COMMIT weighted connectomes!')
%             otherwise
%                 error('** FILTER FAIL: unable to use your filtNAME input!')
%         end
%         D{iw8(kk),jj}                      = dselt;
%     end

    
%% Filter Option 2 (standard)
% Filter all connectomes

elseif ~isempty(filtALL) && filtALL==1

    switch lower(filtMeth)
        case 'sift2'
            disp(['Applying ' upper(filtMeth) ' filter to all STRUCTURAL ' stradd 'connectomes '])
            disp( '* IMPORTANT * DO NOT USE IF NORMALIZING BY NODE VOLUME!!! ')
        case 'commit'
            disp(['Applying ' upper(filtMeth) ' filter to all STRUCTURAL ' stradd 'connectomes '])
        case 'nz'
            disp(['Filtering: Removing all edges which = 0 for any subject from all STRUCTURAL ' stradd 'connectomes'])
        case 'uniform'
            if isempty(Uthr) || ~isnumeric(Uthr) || Uthr < 1 || Uthr > 99
                error('Check your input to Uthr')
            end
            disp('*---------------------- Applying uniform threshold -----------------------* ')
            disp([' Will remove edges which are not present in at least ' num2str(Uthr) '% of subjects'])
        otherwise
            error('** FILTER FAIL: unable to use your filtMeth input!')
    end

    % Filter data
    for jj = 1 : J
        dselt                           = D{iw8,jj};                        % dataset to be used as filter reference
        for ii = 1 : I
            if isempty(cellstrfind(W8s(ii), 'fc')) || FCrm==1               % Only apply filter to FC if desired
                dt                      = D{ii,jj};                         % dataset to be filtered
                % Apply filter
                switch lower(filtMeth)
                    case 'sift2'
                        dt(dselt<1)         = 0;
                    case 'commit'
                         dt(dselt<1e-12)    = 0;
                    case 'nz'
                        mask            = sum(double(dt==0),3)==0;          % 1=edges which have no 0 entries across all subjects
                        for ss = 1:size(dt,3)
                            dt(:,:,ss)  =  dt(:,:,ss).*mask;                % apply mask to all subjects
                        end
                    case 'uniform'
                        C               = Uthr*size(dt,3)/100;              % number of subjects corresponding to threshold
                        mask            = sum(double(dt~=0),3)>C;           % counts non-zero edges over 3rd dim & sets to 1 in mask if > C
                        for ss = 1:size(dt,3)
                            dt(:,:,ss)  =  dt(:,:,ss).*mask;                % apply mask to all subjects
                        end
                end
                D{ii,jj}                = dt;
            end 
        end 
    end 
else
    error('** FILTER FAIL: filtALL input should be [], 0 or 1')
end
%--------------------------------------------------------------------------
end
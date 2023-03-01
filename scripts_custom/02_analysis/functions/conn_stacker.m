function [Cs, Ps, Ws, Ss] = conn_stacker(sourceDir, varargin)
%  Iterates over a directory of subject level connectomes and stacks all
%  files into a cell array.
%
%  Expects file tree to be organized as follows: 
% 
%   connectomes/filter-type/parcellation/edge-weight/subject-connectomes.txt
%
%  Connectomes are stacked into a WEIGHT x PARCELLATION cell array (WxP)
%  where each cell contains an NxNxS matrix (N=Nodes, S=Subjects).
%
% NOTES:
%  (1) Dir structure must be connectomes/filter/parcellation/weight/files
%  (2) Subject labels on files must include "HC" tag  (e.g. HC001) and must
%      be consistent across weights and parcellations. If subjects are not
%      the same, should call this function separately for each subset with
%      consistent subjects.
%  (3) Depending on directory structure, offsets may need to be modified
%      i.e. position of first parcellation, weight & connectome may vary by
%      directory and system.
%
% Input:
%                           + Required +
%   sourceDir   : Target directory (see above for necessary structure)
%
%                           + Optional +
%   outDir      : Save here
%
% Output:
%   Cs      : WEIGHT x PARC cell array. Each cell containing a stack of
%             subject-level connectomes of size NxNxS (S=subs)
%   Ps      : 1xP cell array of parcellation labels
%   Ws      : 1xW cell array edge weight labels
%   Ss      : 1xS cell array of subject labels
%
% 2021 Mark C Nelson, McConnell Brain Imaging Centre, MNI, McGill
%--------------------------------------------------------------------------

%% Optional inputs
nin=max(nargin,1)-1; 
defaults    = {[] };
[outDir]    = INhandler(varargin,nin,defaults);

%% Settings

% get 1st pass of directory dimensions
dirPrc  = dir(sourceDir); 
dirW8   = dir([sourceDir '/' dirPrc(end).name]);
dirSub  = dir([sourceDir '/' dirPrc(end).name '/' dirW8(end).name]);


% Get 1st pass at offsets  (Alternatively: find(strcmp({dirTrg.name}, 'glasser-360')==1))
dSkip   = {'.', '..','.DS_Store'};                                          % set names of dir entries of non-interest
L1os    = dirOffset(sourceDir, dSkip);
L2os    = dirOffset([sourceDir '/' dirPrc(L1os).name], dSkip);
L3os    = dirOffset([sourceDir '/' dirPrc(L1os).name '/' dirW8(L2os).name], dSkip);

% Loop vars
Pn      = length(dirPrc) - L1os + 1;
Wn      = length(dirW8)  - L2os + 1;
Sn      = length(dirSub) - L3os + 1;

Ps      = cell( 1 , Pn );                                                   % Parcellations
Ws      = cell( 1 , Wn );                                                   % Weights
Ss      = cell( 1 , Sn );                                                   % Subjects
Cs      = cell( Wn, Pn );                                                   % Connectomes

%% iterate through directory to extract individual data
for ii = L1os : length(dirPrc)

    Pi                  = ii - L1os + 1;
    PARC                = dirPrc(ii).name;
    Ps{Pi}              = PARC;                                             % extract parcellation names

    disp('----------------------------------------------------------------')
    disp(['Parcellation scheme: ' PARC])


    % Get this specific dir & double check offset
    dirW8               = dir([sourceDir '/' PARC]);
    L2os                = dirOffset([sourceDir '/' PARC], dSkip);
    Wi                  = 0;

    for jj = L2os : length(dirW8)                                           % assumes weights same for all parcs

        Wi              = Wi + 1;
        W8              = dirW8(jj).name;
        Ws{Wi}          = W8;                                               % extract weight name

        fprintf(['\nConnection weight: ' W8 '\n'])

        %Get this specific dir & offset
        dirSub          = dir([sourceDir '/' PARC '/' W8]);                 % active dir = subj connectomes for this weight & parcell
        L3os            = dirOffset([sourceDir '/' PARC '/' W8], dSkip);
        Si              = 0;

        % Ensure SUBS array is large enough: number of subjects may not be uniform 
        tmpsn           = length(dirSub) - L3os + 1;
        while tmpsn > length(Ss)
            Sn          = Sn + 1;
            Ss{Sn}      = [];
        end

        for kk = L3os : length(dirSub)                                      % loop over individual weighted connectome files (.txt)

            Si          = Si + 1;
            FILE        = dirSub(kk).name;

            % Get sub ID tag
            tmpcell     = strsplit(FILE,'_');
            ind         = cellstrfind(tmpcell,'HC');
            if ~isempty(ind) && length(ind)==1
                SUB     = tmpcell{ind};
            else
                disp(['WARNING: unable to parse subject tag: ' FILE])
            end


            % Control HC tags
            if isempty(Ss{Si});
                Ss{Si}              = SUB;                                  % Extract subject id from file name on 1st pass
            elseif ~strcmp(Ss{Si}, SUB);                                    % else confirm match with file
                disp(['WARNING: Filename & Subject tag DO NOT match for ' SUB ' and ' Ss{Si}])
                keyboard
                Dtmp                = nan(size(Cs{Wi,Pi},1));
                Cs{Wi,Pi}(:,:,Si)   = Dtmp; clear Dtmp                      % Insert slice of nans to hold position of missing sub
                Si                  = Si + 1;
            end

            disp( ['ALL GOOD: Extracting : ' FILE] )
            try
                Dtmp    = dlmread([sourceDir '/' PARC '/' W8 '/' FILE]);    % ACTIVE FILE
                Cs{Wi,Pi}(:,:,Si) = Dtmp; clear Dtmp
            catch
                keyboard
            end
        end
    end
end                                                                         % outer loop

% %% Save data
if ~isempty(outDir)
    fprintf('\nSaving data...\n')
    save([outDir '/1_parcelations.mat'], 'Ps')
    save([outDir '/2_weights.mat'],      'Ws')
    save([outDir '/3_subjects.mat'],     'Ss')
    save([outDir '/0_connectomes.mat'],  'Cs')
end
%--------------------------------------------------------------------------
end


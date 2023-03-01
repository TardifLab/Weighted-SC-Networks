function [D] = conn_cut(D,W8s,seg,pinfo)
%  Cuts out matrix entries corresponding to excluded ROIs. 
%  These may include: the L & R medial wall, subcortex, cerebellum. 
%
% Input:
%   D       : IxJ cell array, each cell containing an NxNxS stack of
%             2D connectivity matrices
%   W8s     : Cell array of labels describing edge weights in D
%             (used to identify SC & FC data)
%   seg     : desired output type: {
%             'full' : cortex + subcortex + cerebellum
%             'cor'  : cortex only
%             'sub'  : cortex + subcortex only}
%   pinfo   : 1xJ structure with parcellation info (see conn_getParcInfo.m)
%
% Output:
%   D       : Trimmed dataset
% 
% 2021 Mark C Nelson, McConnell Brain Imaging Centre, MNI, McGill
%--------------------------------------------------------------------------

% Setup
[I,J]       = size(D);
mwfc        = 49;                                                           % medial wall in all FC datasets DOUBLE CHECK! May depend on how they were registered...)
mwsc        = {pinfo.medwallSC};                                            % indices of SC medial wall
exi         = {pinfo.exclude};                                              % indices of additional excluded ROIS from SC & FC


% Remove medial wall
for ii = 1 : I
    if isempty(cellstrfind(W8s(ii), 'fc'))
        
        % SC data: 2 medial wall parcels
        for jj = 1 : J
            
            % Add indices for additional exclusions
            if ~isempty(exi{jj})
                rmi = [mwsc{jj}; exi{jj}(:)];
            else
                rmi = mwsc{jj};
            end
            rmi     = rmi + 48;                                             % 48 nodes accounts for subcortex & cerebellum
            
            dtmp    = D{ii,jj};                                             % Data
            dbin    = double(dtmp~=0);
            mwit    = find(sum(sum(dbin,2),3)==0);                          % Medial wall parcels should be 0
            
            if numel(mwit)==2 && all(mwit == rmi(1:2))                      % if perfect match for medial wall indices
                dtmp(rmi,:,:)=[]; dtmp(:,rmi,:)=[];                         % Remove these rows/columns (might include additional exclusions!)
                D{ii,jj}=dtmp;                    
            elseif numel(mwit)>2 && any(mwit == rmi(1)) && any(mwit == rmi(2)) % many rows sum to 0, including the mw indices taken from pinfo
                dtmp(rmi,:,:)=[]; dtmp(:,rmi,:)=[];                         % Remove these rows/columns
                D{ii,jj}=dtmp;
            else
                disp('* * * WARNING: Double check these medial wall indices!!')
                keyboard
            end
        end
    else
        % FC data: 1 medial wall parcel
        for jj = 1 : J
            if ~isempty(exi{jj})
                EXIt        = exi{jj} + 48;                                 % + 48 to account for subcortex & cerebellum
                EXIt(2,:)   = EXIt(2,:) - 1;                                % -1 from the 2nd row to account for lack of R medial wall in FC (assumes pairs are ipsilateral ROIs)
                rmi         = [mwfc; EXIt(:)];
            else
                rmi         = mwfc;
            end
            dtmp            = D{ii,jj};
            dtmp(rmi,:,:)   = []; 
            dtmp(:,rmi,:)   = [];                                           % Remove left medial wall from FC data
            D{ii,jj}        = dtmp;
        end
    end
end

% Perform desired segmentation
switch seg
    case 'full'
        disp('Eliminating L & R medial wall only --> CORTICAL-SUBCORTICAL-CEREBELLAR connectomes')
        
    case 'cor'
        % eliminate subcortical structures as well (should be the same in all datasets...)
        disp('Eliminating L & R medial wall, subcortex & cerebellum --> CORTICAL connectomes')
        for ii = 1 : I
            for jj = 1 : J
                dtmp=D{ii,jj};
                dtmp(1:48,:,:)=[]; dtmp(:,1:48,:)=[];
                D{ii,jj}=dtmp;
            end
        end
             
    case 'sub'
        % eliminate cerebellum & keep subcortex
        disp('Eliminating L & R medial wall & cerebellum --> CORTICAL-SUBCORTICAL connectomes')
        for ii = 1 : I
            for jj = 1 : J
                dtmp=D{ii,jj};
                dtmp(15:48,:,:)=[]; dtmp(:,15:48,:)=[];
                D{ii,jj}=dtmp;
            end
        end
    otherwise
        error('Please include a string indicating your desired output type...')
end

% Double check result (imperfect solution, won't always have coordinates)
tmp     = cellfun(@(c) length(c),D);                                        % retained node count in D by parcellation & weight
tmp2    = arrayfun(@(c) length(c.coor),pinfo);                              % node coordinates by parcellation
if any(range(tmp) ~= 0) || any(tmp(1,:) ~= tmp2)
    disp('* * * WARNING: Check dimensions of D')
    keyboard
end

%--------------------------------------------------------------------------
end
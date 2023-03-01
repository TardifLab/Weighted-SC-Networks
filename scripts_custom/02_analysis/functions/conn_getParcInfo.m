function [pinfo] = conn_getParcInfo(parcs,ppath,seg)
%  loads look-up-table & conte69 labels by parcellation to get relevant info
%
%   NOTES: 
%       1. Primarily tested on aparc, glasser and Schaefer parcellations. 
%          May need to be expanded to meet your needs.
%       2. aparc (desikan-killiany) labels have problems with some plotting
%          software. Micapipe aparc labels may have an additional ROI not
%          contained in the original DK atlas 
%          (Insula? corpus callosum?)
%       3. COR & SUB segmentations have been tested extensively. FULL 
%          (including the cerebellum) may require expansion/adpatation.
%
% Input:
%   parcs   : 1xP cell array of strings describing parcellations
%              (should match filenames e.g. 'glasser-360')
%   ppath   : location of parcellation info files
%   seg     : string indicating if subcortex/cerbellum retained 
%             (see conn_cut.m)
%
% Output:
%   pinfo   : 1xP structure with fields:
%
%       pinfo.name          : parcellation name
%       pinfo.ncor          : cortical node count
%       pinfo.nsub          : subcortical node count
%       pinfo.labels        : character labels for rois/nodes
%       pinfo.conte69i      : indices to map nodes to conte69 surface
%       pinfo.coor          : x,y,z node coordinates
%       pinfo.cis           : indices for community assignments
%       pinfo.cirois        : nodes relabeled by community assignments
%       pinfo.clabels       : text labels for communities
%       pinfo.eucliddist    : geodesic distance between all parcel centroids
%       pinfo.medwallSC     : indices for medial wall in SC data
%       pinfo.exclude       : indices for additioanl excluded ROIs from SC & FC
%       pinfo.seg           : brain coverage e.g. cor, sub or full
%                             (see conn_cut.m)
%       pinfo.InodeHemiL    : Indices for LEFT hemipshere nodes (nx1 double)
%       pinfo.InodeHemiR    : Indices for RIGHT hemipshere nodes (nx1 double)
%       masks               : Edge masks (1x3 structure)
%                             - hemisphere: Within/between hemisphere
%                             - community : Within/between community
%                             - gradient  : Within/between unimodal/transmodal
%
% 2021 Mark C Nelson, McConnell Brain Imaging Centre, MNI, McGill
%--------------------------------------------------------------------------

% setup
pn      = numel(parcs);
lutpath = [ppath '/lut'];
opt     = 'UniformOutput';

pinfo   = struct('name',       parcs,   ...
                 'ncor',       [],      ...
                 'nsub',       [],      ...
                 'labels',     [],      ...
                 'conte69i',   [],      ...
                 'coor',       [],      ...
                 'cis',        [],      ...
                 'cirois',     [],      ...
                 'clabels',    [],      ...
                 'eucliddist', [],      ...
                 'medwallSC',  [],      ...
                 'exclude',    [],      ...
                 'seg',        seg,     ...
                 'InodeHemiL', [],      ...
                 'InodeHemiR', [],      ...
                 'masks', []);

%% Community names: parcellation specific (Manually input, may need to be expanded/adapted)
% Yeah, I did all this manually... You're welcome :)
C                       = []; 

% Subcortex only
C.sub.Subcortex         = {'Thalamus' 'Caudate' 'Putamen' 'Pallidum' 'Hippocampus' 'Amygdala' 'Accumbens'};

% Subcortex & Cerebellum
C.full.Subcortex        = {'Thalamus' 'Caudate' 'Putamen' 'Pallidum' 'Hippocampus' 'Amygdala' 'Accumbens'};
C.full.Cerebellum       = {};

% Glasser 360
C.glasser.Visual        = {'V1' 'V2' 'V3' 'V4' 'V3A' 'V3B' 'V6' 'V6A' 'V7' 'IPS1' 'V8' 'VVC' 'PIT' 'FFC' 'VMV1' 'VMV2' 'VMV3' 'V3CD' 'LO1' 'LO2' 'LO3' 'V4t' 'FST' 'MT' 'MST' 'PH'};
C.glasser.SomMot        = {'4' '3a' '3b' '1' '2' '24dd' '24dv' '6mp' '6ma' 'SCEF' '5m' '5L' '5mv' '55b' '6d' '6a' 'FEF' '6v' '6r' 'PEF' '43' 'FOP1' 'OP1' 'OP2-3' 'OP4' 'PFcm'};
C.glasser.EarlyAud      = {'A1' 'LBelt' 'MBelt' 'PBelt' 'RI'};
C.glasser.AssocAud      = {'A4' 'A5' 'STSdp' 'STSda' 'STSvp' 'STSva' 'STGa' 'TA2'};
C.glasser.InsulFOC      = {'52' 'PI' 'Ig' 'PoI1' 'PoI2' 'FOP2' 'FOP3' 'MI' 'AVI' 'AAIC' 'Pir' 'FOP4' 'FOP5'};
C.glasser.medTemp       = {'H' 'PreS' 'EC' 'PeEc' 'PHA1' 'PHA2' 'PHA3'};
C.glasser.latTemp       = {'PHT' 'TE1p' 'TE1m' 'TE1a' 'TE2p' 'TE2a' 'TGv' 'TGd' 'TF'};
C.glasser.TPOJ          = {'TPOJ1' 'TPOJ2' 'TPOJ3' 'STV' 'PSL'};
C.glasser.supParie      = {'LIPv' 'LIPd' 'VIP' 'AIP' 'MIP' '7PC' '7AL' '7Am' '7PL' '7Pm'};
C.glasser.infParie      = {'PGp' 'PGs' 'PGi' 'PFm' 'PF' 'PFt' 'PFop' 'IP0' 'IP1' 'IP2'};
C.glasser.postCing      = {'DVT' 'ProS' 'POS1' 'POS2' 'RSC' 'v23ab' 'd23ab' '31pv' '31pd' '31a' '23d' '23c' 'PCV' '7m'};
C.glasser.antCingMPFC   = {'33pr' 'p24pr' 'a24pr' 'p24' 'a24' 'p32pr' 'a32pr' 'd32' 'p32' 's32' '8BM' '9m' '10v' '10r' '25'};
C.glasser.OrbPolFC      = {'47s' '47m' 'a47r' '11l' '13l' 'a10p' 'p10p' '10pp' '10d' 'OFC' 'pOFC'};
C.glasser.infFC         = {'44' '45' 'IFJp' 'IFJa' 'IFSp' 'IFSa' '47l' 'p47r'};
C.glasser.dlPFC         = {'8C' '8Av' 'i6-8' 's6-8' 'SFL' '8BL' '9p' '9a' '8Ad' 'p9-46v' 'a9-46v' '46' '9-46d'};


% All schaefer parcellations
C.schaefer.Visual           = {'Vis'};
C.schaefer.SomMot           = {'SomMot'};
C.schaefer.DorsAtt          = {'DorsAttn'};
C.schaefer.SalVentAtt       = {'SalVentAttn'};
C.schaefer.Limbic           = {'Limbic'};
C.schaefer.Control          = {'Cont'};
C.schaefer.Default          = {'Default'};


% Destrieux
C.aparca2009s.Frontal   = {'S_precentral-sup-part' 'S_precentral-inf-part' 'S_central' 'G_precentral' 'S_front_inf' 'S_orbital-H_Shaped' 'G_rectus' ... 
                           'G_front_inf-Triangul' 'G_front_inf-Opercular' 'G_front_inf-Orbital' 'G_orbital' 'G_and_S_subcentral' 'G_front_middle' ... 
                           'S_front_sup' 'S_front_middle' 'G_front_sup' 'S_suborbital' 'S_orbital_med-olfact' 'G_and_S_frontomargin' ... 
                           'G_and_S_transv_frontopol' };

C.aparca2009s.Temporal  = {'S_temporal_inf' 'S_temporal_sup' 'G_temp_sup-Lateral' 'G_temp_sup-G_T_transv' 'G_temp_sup-Plan_tempo' 'Pole_temporal' ... 
                           'G_temp_sup-Plan_polar' 'S_temporal_transverse' 'G_temporal_middle' 'S_collat_transv_ant' 'S_collat_transv_post' };

C.aparca2009s.Parietal  = {'S_parieto_occipital' 'S_postcentral' 'S_intrapariet_and_P_trans' 'G_and_S_paracentral' 'S_subparietal' 'G_precuneus' ... 
                           'G_postcentral' 'S_interm_prim-Jensen' 'G_pariet_inf-Supramar' 'G_pariet_inf-Angular' 'G_parietal_sup' };

C.aparca2009s.Occipit   = {'G_and_S_occipital_inf' 'S_occipital_ant' 'S_oc_sup_and_transversal' 'G_occipital_middle' 'S_oc_middle_and_Lunatus' ... 
                           'S_calcarine' 'G_cuneus' 'Pole_occipital' 'G_occipital_sup'};
                           
C.aparca2009s.TempOcc   = {'S_oc-temp_lat' 'S_oc-temp_med_and_Lingual' 'G_oc-temp_lat-fusifor' 'G_oc-temp_med-Lingual' 'G_oc-temp_med-Parahip' ... 
                           'G_temporal_inf'};

C.aparca2009s.Insul     = {'S_circular_insula_sup' 'S_circular_insula_inf' 'S_circular_insula_ant' 'S_orbital_lateral' 'G_Ins_lg_and_S_cent_ins' ... 
                           'Lat_Fis-ant-Horizont' 'Lat_Fis-ant-Vertical' 'Lat_Fis-post' 'G_insular_short' };
                           
C.aparca2009s.Limbic    = {'S_cingul-Marginalis' 'G_subcallosal' 'S_pericallosal' 'G_and_S_cingul-Ant' 'G_and_S_cingul-Mid-Ant' ... 
                           'G_and_S_cingul-Mid-Post' 'G_cingul-Post-dorsal' 'G_cingul-Post-ventral'};
                       

% Desikan-Killiany
C.aparc.Temporal        = {'superiortemporal' 'middletemporal' 'inferiortemporal' 'transversetemporal' 'bankssts' 'entorhinal' 'parahippocampal' 'temporalpole' 'fusiform'};
C.aparc.Frontal         = {'superiorfrontal' 'rostralmiddlefrontal' 'caudalmiddlefrontal' 'frontalpole' 'parsopercularis' 'parstriangularis' 'parsorbitalis' 'lateralorbitofrontal' 'medialorbitofrontal' 'precentral' 'paracentral'};
C.aparc.Insula          = { 'insula' };
C.aparc.Parietal        = {'supramarginal' 'superiorparietal' 'inferiorparietal' 'postcentral' 'precuneus'};
C.aparc.Occipital       = {'lingual' 'pericalcarine' '-cuneus' 'lateraloccipital'};
C.aparc.Cingulate       = {'rostralanteriorcingulate' 'caudalanteriorcingulate' 'posteriorcingulate' 'isthmuscingulate'};


%% ----------------------- Load & extract data ----------------------- %%
warning off                                                                 % supress readtable's warnings

%% Subcortex/cerebellum
sub_coor                    = [];
sub_labels                  = [];
sub_Ci                      = [];
sub_Cl                      = [];
sub_rois                    = [];
sub_n                       = 0;

if ~strncmpi('cor',seg,3)
    disp('* NOTE: Subcortical info will be appended to all cortical parcellations.')
    disp('        If this is not desired, check your input to SEG')
    lutfile                 = [lutpath '/lut_subcortical-cerebellum_mics.csv'];
    luttable                = readtable(lutfile);                           % read csv into table
    lutheaders              = luttable.Properties.VariableNames;            % get headers
    lutcell                 = table2cell(luttable);                         % convert to cell array

    if strncmpi('sub',seg,3)
        lutcell(15:48,:)    = [];                                           % eliminate cerebelum if only subcortex desired
    end

    % Get coords & labels
    coori                   = cellstrfind(lutheaders,'coor');
    labeli                  = cellstrfind(lutheaders,'label');              % column indices for coords & labels in LUT
    sub_coor                = cell2mat(lutcell(:,coori));
    sub_labels              = lutcell(:,labeli);
    sub_n                   = length(sub_labels);                           % number of subcortical parcels retained
    
    % Get community assignments if they exist
    S                       = C.(seg);
    sub_Cl                  = fieldnames(S);                                % manual community labels
    sub_Ci                  = cell(1,numel(sub_Cl));                        % community assignment indices
    Lall                    = sub_labels;                                   % node labels (from .csv)
    for li = 1 : numel(sub_Cl)
        Lt                  = S.(sub_Cl{li});                               % subset of labels for this community
        % Try rigid search first
        tmp                 = cellfun(@(c) find(strcmp(Lall,c)),Lt,opt,0);  % node indices for this community
        % Loosen search criteria if necessary
        if all(cellfun(@(c) isempty(c),tmp))
            tmp             = cellfun(@(c) cellstrfind(Lall,c),Lt,opt,0);
        end
        try
            tmp             = cell2mat(tmp);
        catch
            disp('** Check your manual labels!')
            keyboard
        end
        sub_Ci{li}          = sort(tmp(:));
    end
    sub_rois                = 1:sub_n;                                      % Starting node order

    % Relabel nodes by community assignments
    for cc = 1 : numel(sub_Ci)
        mask                        = 1:sub_n;
        cit                         = sub_Ci{cc};
        mask(~ismember(mask,cit))   = 0;
        sub_rois(mask~=0)           = cc;
    end
end
clear S

%% Cortical parcellations
for pp = 1 : pn
    disp(' ')
    disp(['------------------------------ ' parcs{pp} ' ------------------------------'])
    % load conte69 csv
    contefile               = [ppath   '/' parcs{pp} '_conte69.csv'];
    if exist(contefile,'file')
        pinfo(pp).conte69i  = load(contefile);
    else
        disp(['* * * WARNING: Unable to load CONTE69 csv for ' parcs{pp}])
    end 
    
    % Load & process lut csv
    lutfile                     = [lutpath '/lut_' parcs{pp} '_mics.csv'];    
    if exist(lutfile,'file')
        % load
        luttable                = readtable(lutfile);                           % read csv into table
        lutheaders              = luttable.Properties.VariableNames;            % get headers
        lutcell                 = table2cell(luttable);                         % convert to cell array
        
        % EXCLUSIONS
        
        % Get indices for Medial wall
        labeli                  = cellstrfind(lutheaders,'label');
        rmi                     = cellstrfind(lutcell(:,labeli),'medial_wall');
        if isempty(rmi);  rmi   = cellstrfind(lutcell(:,labeli),'medial wall'); end
        pinfo(pp).medwallSC     = rmi;                                      % store for use in cutting connetomes
        
        
        % Manually exclude corpus callosum from aparc to get 68 node atlas
        if strncmpi(parcs{pp},'aparc',4)
            exli                = cellstrfind(lutcell(:,labeli),'corpuscallosum');
            pinfo(pp).exclude   = exli;
            rmi                 = [rmi; exli];
        end
        
        % Exclude these labels
        lutcell(rmi,:)              = [];
        
        
        % Extract coords & labels
        coori                   = cellstrfind(lutheaders,'coor');
        coort                   = cell2mat(lutcell(:,coori));
        labelt                  = lutcell(:,labeli);
        N                       = numel(labelt);
        pinfo(pp).ncor          = N;                                        % number of cortical parcels
        
        
        % If coords are missing, get them from conte69 surface
        if any(any(arrayfun(@(c) isnan(c),coort))) || isempty(coort)
            disp(['* NOTE: Some coords may be missing for: ' parcs{pp}])
            
            if ~isempty(pinfo(pp).conte69i)
                disp('        Getting coordinates from CONTE69 surface')
                [Sl,Sr]             = load_conte69();                       % Load surface from brainspace
                S                   = combine_surfaces(Sl,Sr);              % combine to single structure
                L                   = pinfo(pp).conte69i;                   % conte69 labels
                ct_c69              = nan(N,3);
                for nn = 1 : N
                    ct_c69(nn,1)    = mean(S.coord(1,L==nn));               % get avg coords from surface for each node
                    ct_c69(nn,2)    = mean(S.coord(2,L==nn));
                    ct_c69(nn,3)    = mean(S.coord(3,L==nn));
                end
                coort               = ct_c69;
                clear S Sl Sr L ct_c69
            else
                disp(['* * * WARNING: NO CONTE69 labels found for ' parcs{pp} ' * * * '])
            end
        end
        
        
        % get community assignments from LUT character labels
        shortname           = strsplit(parcs{pp},'-');
        if length(shortname) == 2
            if strncmpi('a2',shortname{2},2)
                shortname   = [shortname{1} shortname{2}];
            else
                shortname   = shortname{1};
            end
        else 
            shortname       = shortname{1};
        end
        
        if ~isempty(cellstrfind(fieldnames(C),shortname))
            S               = C.(shortname);                                % communities
        else
            disp(['* * * WARNING: Unable to find community labels for: ' parcs{pp} ' * * *'])
            keyboard
            S               = [];                                           % If community labels dont exist for this parc
        end
        if ~isempty(S)
            Cl              = fieldnames(S);                                % community labels (manual input)
            Ci              = cell(1,numel(Cl));
            Lall            = labelt;                                       % labels for cortical parcels (from .csv)
            for li = 1 : numel(Cl)
                Lt          = S.(Cl{li});                                   % subset of labels for this community
                % Try rigid search first
                tmp         = cellfun(@(c) find(strcmp(Lall,c)),Lt,opt,0);  % node indices for this community
                % Loosen search criteria if necessary
                if all(cellfun(@(c) isempty(c),tmp))
                    tmp     = cellfun(@(c) cellstrfind(Lall,c),Lt,opt,0);
                end
                try
                    tmp     = cell2mat(tmp);
                catch
                    disp('** Check your manual labels!')
                    keyboard
                end
                Ci{li}      = sort(tmp(:)) + sub_n;                         % make sure to shift by subcortical node count!
            end
            % Double check that community assignments & node counts match
            if sum(cellfun(@(c) numel(c),Ci)) ~= N
                disp(['WARNING: check labels & cis for: ' parcs{pp}])
                keyboard
            end
            rois                            = 1:N;
            
            % Relabel nodes by community assignments
            for cc = 1 : numel(Ci)
                mask                        = 1:N;
                cit                         = Ci{cc} - sub_n;               % make sure to account for the shift from subcortex
                mask(~ismember(mask,cit))   = 0;
                rois(mask~=0)               = cc;
            end
            if ~isempty(max(sub_rois))
                rois                        = rois + max(sub_rois);         % adjust cortical rois to account for subcortex
            end
            
            % Store info for cortex & subcortex together
            pinfo(pp).coor                  = [sub_coor; coort];            % Node centroid coordinates in conte69 space
            pinfo(pp).labels                = [sub_labels ; labelt];        % Node labels from LUT csv
            pinfo(pp).cirois                = [sub_rois, rois];             % cis projected onto rois
            pinfo(pp).cis                   = [sub_Ci, Ci];                 % community assignment indices
            pinfo(pp).clabels               = [sub_Cl; Cl];                 % community labels
            pinfo(pp).nsub                  = sub_n;                        % number of cortical parcels appended
        end
        
        
        % Compute Euclidean Distance between each node centroid (quick loop solution)
        Ntot                        = length(pinfo(pp).labels);
        if any(any(arrayfun(@(c) isnan(c),pinfo(pp).coor))) || isempty(pinfo(pp).coor)
            disp([' Coords still empty, skipping ED for ' parcs{pp}])
        else
            EDmat                   = zeros(Ntot);
            coort                   = pinfo(pp).coor;
            for ii = 1 : Ntot-1
                for jj = ii+1 : Ntot
                    EDmat(ii,jj)    = sqrt(sum((coort(ii,:)-coort(jj,:)).^2));  % 3D Euclidean distance
                end
            end
            pinfo(pp).eucliddist    = EDmat + EDmat';
        end
        
        
        % Get indices for nodes by hemisphere
        lbls    = pinfo(pp).labels;
        
        Ilh     = cellstrfind(lbls,'left');                                 % Subcortical
        Irh     = cellstrfind(lbls,'right');
        
        switch shortname                                                    % Cortex: syntax varies by atlas
            
            case 'schaefer'
                Ilh     = [Ilh; cellstrfind(lbls,'_LH_')];                            
                Irh     = [Irh; cellstrfind(lbls,'_RH_')];
            case 'aparc'
                Ilh     = [Ilh; cellstrfind(lbls,'ctx-lh-')];                            
                Irh     = [Irh; cellstrfind(lbls,'ctx-rh-')];
            otherwise
                disp('CODE UNFINISHED FOR THIS PARCELLATION!')
                keyboard
        end
        
        dblchk  = [isempty(intersect(Ilh,Irh)),...                          % L & R hemi indices should not overlap
                   length(Ilh)==length(Irh) & ~isempty(Ilh)];               % number of indices should be equal & nonzero across L & R hemi
        if all(dblchk)
            pinfo(pp).InodeHemiL = Ilh;
            pinfo(pp).InodeHemiR = Irh;
        else
            disp('WARNING: check sorting of nodes by hemisphere in conn_getParcInfo.m')
            keyboard
        end
        
        %% Create edge group masks for this parcellation (Only works for Schaefer?)
        if ~strncmpi(shortname,'scha',4)                                    % IF NOT using one of the schaefer parcellations
            disp('WARNING: this section might only work with the Schaefer parcellations!')
            disp(' See approx line 390 of conn_GetParcInfo.m')
            keyboard
        end
        
        maskut          = logical(triu(ones(Ntot)));                        % mask for upper triangle WITH main diagonal
        tmpmat          = ones(Ntot);
        tmpmat(maskut)  = 0;
        [tmpa,tmpb]     = find(tmpmat);                                     % row & column indices of all lower tri edges
        tmpind          = find(tmpmat);                                     % linear indices of all lower tri edges
        Ncl             = length(pinfo(pp).cis);
        Clbls           = pinfo(pp).clabels;
        
        % Hemisphere within/between
        tmpid           = double(arrayfun(@(x,y) [ismember(x,Ilh) & ismember(y,Ilh) | ... 
                                                  ismember(x,Irh) & ismember(y,Irh)], tmpa, tmpb)); % logical TRUE for WITHIN hemispheric edges
        maskh           = nan(Ntot);
        maskh(tmpind)   = tmpid;
        maskh(maskh==0) = 2;                                                % BETWEEN hemispheric edges set to 2
        
        
        % Networks within/between
        maskc           = 2*ones(Ntot);                                     % BETWEEN network edges will be left as 2s
        for cc = 1 : Ncl
            it          = pinfo(pp).cis{cc};
            maskc(it,it)= 1;                                                % Set WITHIN network edges to 1
        end
        maskc(maskut)   = nan;

        
        % Unimodal/Transmodal within/between
        if ~strncmpi(shortname,'scha',4)
            maskg=[];
        else
            Ic_exl          = cell2mat(cellstrfind(Clbls,{'Limbic'}));              % indices for excluded networks in community labels
            if ~strncmpi('cor',seg,3)
                Ic_exl      = [Ic_exl cell2mat(cellstrfind(Clbls,{'Subcortex'}))];  % indices for excluded networks in community labels
            end
            Ic_uni          = cell2mat(cellstrfind(Clbls,{'Visual' 'SomMot'}));     % unimodal communities
            Ic_trn          = setdiff(1:Ncl,[Ic_uni Ic_exl]);                       % transmodal communities
            maskg           = ones(Ntot).*2;                                        % 2s will be used for between networks

            it              = [];
            for cc = Ic_uni
                it          = [it; pinfo(pp).cis{cc}];
            end
            maskg(it,it)    = 1;                                                    % Within unimodal network edges

            it              = [];
            for cc = Ic_trn
                it          = [it; pinfo(pp).cis{cc}];
            end
            maskg(it,it)    = 3;                                                    % Within transmodal network edges

            it              = [];
            for cc = Ic_exl
                it          = [it; pinfo(pp).cis{cc}];
            end
            maskg(it,:)     = 0;                                                    % All excluded connections
            maskg(:,it)     = 0;

            maskg(maskut)   = nan;
            clear tmp*
        end
        
        % Save masks
        pinfo(pp).masks.hemisphere  = maskh;                                % 1 = within hemisphere, 2 = Between hemispheres
        pinfo(pp).masks.community   = maskc;                                % 1 = within community,  2 = Between community
        pinfo(pp).masks.gradient    = maskg;                                % 1 = unimodal, 2 = between, 3 = transmodal, 0 = excluded (subcortex & limbic)
        
        
    else
        disp(['NOTE: Unable to load LUT csv for ' parcs{pp} ' from ' lutpath])
    end
end
%--------------------------------------------------------------------------
end
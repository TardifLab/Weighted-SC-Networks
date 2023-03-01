function [Osub_r,varargout] = conn_controlW8(D,W8s,target)
% Controls for variance related to the dataset indicated by TARGET in all
% other datasets in D.
%
% Variance is removed using linear regression. 
%
% Subject & group level connectomes are returned with edge values 
% equivalent to the residuals after linear regression.
%
% Input
%                           + Required +
%   D       : IXJ cell array, each cell containing an NxNxS stack of
%             connectivity matrices (S=subs) 
%   W8s     : 1XI cell array of strings describing edge weights in D
%   target  : character vector indicating which dataset in D to control
%
%
% Outputs
%                           + Required +
%   Dsub_r  : I-1 X J cell array of residuals 
%             (Same shape as D less the target dataset)
%
%                           + Optional +
%   W8s_r   : W8s input modified to show regression of target
%   Ogrp_r  : Group level residuals of data in D
%
% 2022 Mark C Nelson, McConnell Brain Imaging Centre, MNI, McGill
%--------------------------------------------------------------------------
[I,J]                       = size(D);

%% Optional inputs
% nin                         = max(nargin,1) - 3;
% defaults                    = {};
% []   = INhandler(varargin,nin,defaults);

%% Setup
Itrg                        = cellstrfind(W8s,target);                      % index of control dataset
Iothr                       = setdiff(1:I,Itrg);
Osub_r                      = cell(I-1,J);
Ogrp_r                      = cell(I-1,J);
W8s_r                       = cell(1,I-1);

for jj = 1 : J
    trg_s                   = D{Itrg,jj};
    trg_g                   = groupavg(trg_s,3,'nz');
    Nnode                   = length(trg_s);    
    Nsub                    = size(trg_s,3);
    maskut                  = logical(triu(ones(Nnode)));                   % upper tri WITH main diagonal (assuming symmetry!)
    
    trg_g(maskut)           = nan;                                          % Isolate lower triangle
    trg_g(trg_g==0)         = nan;
    tmpind                  = find(~isnan(trg_g));                          % These indices will be used to extract from target & dataset of interest
    Ntmpind                 = length(tmpind);
        
    
    % Loop over all other datasets
    it=0;
    for ii = Iothr
        it=it+1;
        TRUE=0;
        
        D_s             = D{ii,jj};
        
        %% Group Level
        D_g             = groupavg(D_s,3,'nz');
        D_g(maskut)     = nan;                                              % Isolate lower triangle
        D_g(D_g==0)     = nan;                                              % Isolate non-zero values in lower tri
                
        vec_g           = D_g(tmpind);                                      % Vectorize non-zero edges in lower tri
        
        % Log transform if necessary
        if abs(skewness(vec_g)) > 3
            TRUE            = 1;
            vec_g           = log(vec_g);
        end
        
        % Extract target data
        trgvec_g            = trg_g(tmpind);

        % Regress target data variance
        offsetg             = ones(Ntmpind,1);                              % intercept
        P                   = [offsetg, trgvec_g];
        [~,~,vec_rsd_g]     = regress(vec_g,P);
        
        
        % Restore original distribution     (SMART?)
        if TRUE==1;
            vec_rsd_g       = exp(vec_rsd_g);
        end
        
        % Store group level data
        D_g(:)              = 0;                                            % Clear original group matrix
        D_g(tmpind)         = vec_rsd_g;                                    % Map residuals to connectome edges
        D_g                 = D_g + D_g';                                   % Restore upper tri
        Ogrp_r{it,jj}       = D_g;


        %% Subject-Level (Would it be better to run the regression on all subject data together???)
        for ss = 1 : Nsub
            
            % Extract target data
            trgsub              = trg_s(:,:,ss);                            % Extract target data
            trgsub(maskut)      = nan;                                      % Isolate lower triangle
            trgsub(trgsub==0)   = nan;                                      % Isolate non-zero values in lower tri
            tmpind              = find(~isnan(trgsub));
            Ntmpind             = length(tmpind);
            trgsubvec           = trgsub(tmpind);

            % Extract dataset of interest
            Dsub                = D_s(:,:,ss);
            Dsub(maskut)        = nan;                                      % Isolate lower triangle
            Dsub(Dsub==0)       = nan;                                      % Isolate non-zero values in lower tri
            Dsubvec             = Dsub(tmpind);                             % Vectorize non-zero edges in lower tri

            % Log transform if necessary
            if TRUE==1
                Dsubvec         = log(Dsubvec);
            end
            

            % Regress LoS
            offsets             = ones(Ntmpind,1);                          % intercept
            P                   = [offsets, trgsubvec];
            [~,~,Dsubvec_rsd]   = regress(Dsubvec,P);

            % Restore original distribution     (SMART?)
            if TRUE==1;
                Dsubvec_rsd     = exp(Dsubvec_rsd);
            end
            
            % Store results
            Dsub(:)             = 0;                                        % Clear original subject matrix
            Dsub(tmpind)        = Dsubvec_rsd;                              % Map residuals to connectome edges
            Dsub                = Dsub + Dsub';                             % Restore upper tri
            D_s(:,:,ss)         = Dsub;
        end
        Osub_r{it,jj}           = D_s;
        W8s_r{it}               = [W8s{ii} ' (rsd)'];
    end
end



%% Optional output
nout=nargout-1;
if nout > 0
    varargout = cell(1,nout);
    for oo = 1 : nout
        switch oo
            case 1
                varargout{oo}   = W8s_r;
            case 2
                varargout{oo}   = Ogrp_r;
        end
    end
end


%--------------------------------------------------------------------------
end
function plot_conn_histSmooth(D,W8s,varargin)
%  Plots smoothed histograms of input data using a kernal smoothing
%  function
%
% Input:
%                           + Required +
%   D       :  IxJ cell array, each cell containing an NxNxS stack of
%              connectivity matrices (S=subs)
%   W8s     : Cell array of strings describing edge weights
%
%                           + Optional +
%   PARCS   : Cell array of strings describing parcellations
%   str     : character vector for plot title
%   cmap    : Ix3 numerical vector indicating desired color for each weight
%   ppos    : indicates plot size & position
%   Zdo     : 0/1 zscore here (0=default)
%
% 2022 Mark C Nelson, McConnell Brain Imaging Centre, MNI, McGill
%--------------------------------------------------------------------------
%% Setup
[I,J]   = size(D);
fntsz   = 25;


%% Optional inputs
nin                                         = max(nargin,1) - 2;
if nin < 5 || isempty(varargin{5}); Zdo     = 0; end
if nin < 4 || isempty(varargin{4}); ppos    = [-2559 839 1276 367]; end
if nin < 3 || isempty(varargin{3}); cmap    = colormap(hsv(I)); end
if nin < 2 || isempty(varargin{2}); str     = ''; end
if nin < 1 || isempty(varargin{1}) 
    PARCS = cell(1,J);
    for jj = 1 : J; PARCS{jj} = ['unknown parcellation ' num2str(jj)]; end 
end

if nin > 0
    for kk = 1:nin
        if ~isempty(varargin{kk})
            switch kk
                case 1
                    PARCS   = varargin{kk};
                case 2
                    str     = varargin{kk};
                case 3
                    cmap    = varargin{kk};
                case 4
                    ppos    = varargin{kk};
                case 5
                    Zdo     = varargin{kk};
            end  
        end
    end
end


%% Plots
% Over parcellations
for jj = 1 : J
    strplt  = ['Subject level ' upper(str) ' distributions: ' PARCS{jj} ' (normalized)'];
    myfig(strplt,ppos);
    W8plt   = W8s;
  
    % Over weights
    for ii = 1 : I
        d           = D{ii,jj}; 
        d(isnan(d)) = 0;
        
        % switch to zscore data if necessary
        switch Zdo
            case 1; d = zscore(d(d~=0));
            case 0; d = d(d~=0);
            otherwise
                disp('ERROR: Check Zdo input...')
                return
        end
        

        % Compute density based on skewness
        if abs(skewness(d)) > 3
            [f,xi]      = ksdensity(log10(d));
            W8plt{ii}   = [W8s{ii} ' (log10)'];
        else
            [f,xi]      = ksdensity(d);
            f=f./max(f);
        end
        
        plot(xi,f,'Color',cmap(ii,:),'LineWidth',8); hold on

    end
    title(strplt); 
    ylabel('Probability Density Estimate');
    legend(W8plt,'FontSize', fntsz);
    set(gca,'FontSize', fntsz);
end
%--------------------------------------------------------------------------
end

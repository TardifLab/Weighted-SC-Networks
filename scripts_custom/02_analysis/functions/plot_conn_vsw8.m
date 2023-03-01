function [r] = plot_conn_vsw8(D,W8s,Wvs,varargin)
%  Plots group- or pooled-subject-level heat scatter plots of the value of
%  a given edge in one connectome vs another. Only non-zero values are
%  retained.
%
%  * WARNING * visualizing at subject level very slow!!! Recommend only
%  doing this with reduced quantities of data i.e. a 2x1 cell array 
%  or S<10. N<400 also recommended.
%
% Input:
%                           + Required +
%   D       : IxJ cell array, each cell containing an NxNxS stack of
%             connectivity matrices (S=subs)
%   W8s     : Cell array of strings describing edge weights
%   Wvs     : Weight to plot on x axis i.e. all weights vs this one
%
%                           + Optional +
%   PARCS   : Cell array of strings describing parcellations
%   str     : character vector for plot title
%   level   : string indicating level of visualization {'subject' 'group'}
%              
%   avg     : String indicating method for group-averaging (see groupavg.m)
%             default = 'nz'
%   plotON  : []/1 toggle plotting off/on (default=on)
%
% Output:
%   r       : correlations
%
% 2021 Mark C Nelson, McConnell Brain Imaging Centre, MNI, McGill
%--------------------------------------------------------------------------
[I,J]                       = size(D);
fntsz                       = 25;
Type                        = 'Spearman';
Copt                        = {'Type',Type,'Rows','complete'};              % Correlation call flags

switch lower(Type)
    case 'pearson'
        strtyp          = 'r';
    case 'spearman'
        strtyp          = '\rho';
end

%% Optional inputs
nin                         = max(nargin,1) - 3;
if nin < 5; plotON          = 1;                    end
if nin < 4; avg             = 'nz';                 end
if nin < 3 
    if all(all(cellfun(@(c) length(size(c)),D) == 3))                       % assume group average desired if data 3D
        level               = 'group';
    else
        level               = 'subject';
    end
end
if nin < 2; str             = 'data';               end
if nin < 1 
    PARCS                   = cell(1,J);
    for jj = 1:J
        PARCS{jj}           = ['PARC-' num2str(jj)]; 
    end       
end

if nin > 0
    for kk = 1:nin
        switch kk
            case 1
                PARCS       = varargin{kk};
            case 2
                str         = varargin{kk};
            case 3
                level       = varargin{kk};
            case 4
                avg         = varargin{kk};
            case 5
                plotON      = varargin{kk};
        end  
    end
end

%% Setup
maxsp           = 20;
if (I-1) > maxsp
    n           = 5;                                                        % dims for subplot
    m           = 4;
    ppos        = [1442 -539 2559 1336];                                    % plot size & location [1442 -539 2350 1336]
else
    n           = ceil(sqrt(I-1));
    m           = ceil((I-1)/n);
    ppos        = [1442 -457 2226 1254];
end
vsi             = cellstrfind(W8s, Wvs);
r               = nan(I-1,J);

if isempty(vsi); error('Unable to find an index for your input: "Wvs"'); end
if numel(vsi)>1; error(['* Ambiguous label! Unable to plot heat scatter vs ' Wvs]); end

%% Process
for jj = 1 : J
    strplt  = [upper(Wvs) ' vs ' str ': ' level '-level : ' PARCS{jj}];
    sp      = 0;
    in      = 0;

    if ~isempty(plotON) && plotON==1
        myfig(strplt,ppos)
    end

   % Extract & plot each weight vs the target weight
    for ii = setdiff(1:I, vsi)
        in                      = in + 1;                                   % index for storing correlations
        Dx                      = D{vsi,jj};                                % re-extract X for each weight (exclusions vary)
        Dy                      = D{ii,jj};

        switch level
            case 'subject'
                Dx(isnan(Dx))   = 0;
                Dy(isnan(Dy))   = 0;
                Dx              = Dx(Dy~=0);                                % vectorize by non-zero elements of Dy
                Dy              = Dy(Dy~=0);
                Dy              = Dy(Dx~=0);                                % further index non-zero elements of Dx
                Dx              = Dx(Dx~=0);
            case 'group'
                Dy              = groupavg(Dy,3,avg);                       % Group-average
                Dx              = groupavg(Dx,3,avg);
                V               = conn2vec({Dy Dx});
                Dy              = V{1};
                Dx              = V{2};
            otherwise
                disp('Issue with your "level" input... using group-level data by default')
                Dy              = groupavg(Dy,3,avg);                       % Group-average
                Dx              = groupavg(Dx,3,avg);
                V               = conn2vec({Dy Dx});
                Dy              = V{1};
                Dx              = V{2};
        end

       % Prepare data for plotting
        Cx                      = [0 max(Dx)];                              % Determine limits for this weight
        Cy                      = [min(Dy) max(Dy)];
        r(in,jj)                = corr(Dx, Dy, Copt{:});
        rstr                    = {sprintf(' = %.3f', r(in,jj))};

        % plot
        if ~isempty(plotON) && plotON==1
            sp                  = sp + 1;
            if sp > maxsp
                sp              = 1;
                myfig(strplt,ppos)
            end
            subplot(m,n,sp)
            heatscatter(Dy,Dx)
            title([strtyp rstr{:}])
            ylabel(W8s{vsi}); xlabel(W8s{ii}); ylim(Cx); xlim(Cy);
            set(gca,'FontSize', fntsz)
        end
   end
end
%--------------------------------------------------------------------------
end
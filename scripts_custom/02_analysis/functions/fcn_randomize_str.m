function [B,energymin, init_energy] = fcn_randomize_str(A,varargin)
%FCN_RANDOMIZE_STR  rewire to preserve degree and strength
%
%   [B,ENERGYMIN] = FCN_RANDOMIZE_STR(A,VARARGIN) rewires the matrix A,
%   preserving degree sequence exactly and strength sequence approximately.
%
%   Algorithm: Uses Maslov & Sneppen rewiring model to produce a binary
%   matrix with the same degree sequence as A. Next, the original edge
%   weights are randomly assigned to existing edges. The weights of the
%   edges are then permuted and optimized using a simulated annealing
%   routine, which tries to match the strength sequence of the randomized
%   network to that of A.
%
%
%       'maxswap' - number of proposed swaps for the edge rewiring
%       algorithm.
%
%       'nstage' - number of stages for annealing step.
%
%       'niter' - number of iterations per stage.
%
%       'temp' - initial temperature.
%
%       'frac' - fractional decrease in temperature per stage.
%
%       'energy' - the energy function to minimze. Choices are:
%           'euclidean' - Euclidean distance between strength sequence
%                         vectors of the original network and the
%                         randomized network.
%           'max' - the single largest value by which the strength
%                         sequences deviate.
%
%       'verbose' - print status to screen at the end of every stage.
%
%   Input:      A,          symmetric/weighted adjacency matrix
%
%   Outputs:    B,          randomized network
%               ENERGYMIN,  minimum energy obtained by annealing
%
%
%   Richard Betzel, Indiana University, 2014

narg = length(varargin);
for i = 1:2:narg
    str = varargin{i};
    switch str
        case 'maxswap'
            maxswap = varargin{i + 1};
        case 'nstage'
            nstage = varargin{i + 1};
        case 'niter'
            niter = varargin{i + 1};
        case 'temp'
            temp = varargin{i + 1};
        case 'frac'
            f = varargin{i + 1};
        case 'verbose'
            verbose = varargin{i + 1};
        case 'energy'
            energytype = varargin{i + 1};
    end
end

n = length(A);
if ~exist('maxswap','var')
    maxswap = 2^10;
end
if ~exist('nstage','var')
    nstage = 1;
end
if ~exist('niter','var')
    niter = 50000;
end
if ~exist('temp','var')
    temp = 0;
end
if ~exist('frac','var')
    f = 0.8;
end
if ~exist('verbose','var')
    verbose = false;
end
if ~exist('energytype','var')
    energytype = 'euclidean';
end

s = sum(A,2);

B = subfcn_randomize_graph(A,maxswap);

[u,v,wts] = find(triu(B,1));
ind = (v - 1)*n + u;
m = length(wts);
sb = sum(B,2);

switch energytype
    case 'euclidean'
        energy = sum((s - sb).^2);
    case 'max'
        energy = max(abs(s - sb));
end
init_energy = energy;
energymin = energy;
wtsmin = wts;

for istage = 1:nstage
    
    naccept = 0;
    for iter = 1:niter
        
        e1 = randi(m);
        e2 = randi(m);
        
        a = u(e1); b = v(e1);
        c = u(e2); d = v(e2);
        
        sb_prime = sb;
        sb_prime([a,b]) = sb_prime([a,b]) - wts(e1) + wts(e2);
        sb_prime([c,d]) = sb_prime([c,d]) + wts(e1) - wts(e2);
        
        switch energytype
            case 'euclidean'
                energy_prime = sum((sb_prime - s).^2);
            case 'max'
                energy_prime = max(abs(sb_prime - s));
        end
        if energy_prime < energy || rand < exp(-(energy_prime - energy)/temp);
            sb = sb_prime;
            wts([e1,e2]) = wts([e2,e1]);
            energy = energy_prime;
            if energy < energymin
                energymin = energy;
                wtsmin = wts;
            end
            naccept = naccept + 1;
        end
        
    end
    temp = temp*f;
    if verbose
        fprintf('stage %i, temp %.5f, best energy %.5f, frac of accepted moves %.3f\n',istage,temp,energymin,naccept/niter);
    end
end
B = zeros(n);
B(ind) = wtsmin;
B = B + B';

function A = subfcn_randomize_graph(A,maxswap)
%SUBFCN_RANDOMIZE_GRAPH    swap edges/preserve degree sequence
%
%   A = FCN_RANDOMIZE_GRAPH(A,MAXSWAP) takes adjacency matrix A and
%       performs MAXSWAP rewirings, returning the randomized matrix A.
%
%   Inputs:       A,      adjacency matrix
%           MAXSWAP,      number of rewirings
%
%   Outputs:      A,      randomized matrix
%
%   Richard Betzel, Indiana University, 2013
%
%   Notes:
%   1. Based on script written by Olaf Sporns (randmio_und.m).
%   2. Graph may become disconnected as a result of rewiring. Always
%      important to check.
%   3. A can be weighted, though the weighted degree sequence will not be
%      preserved.
%   4. A must be undirected.
%

[i,j] = find(triu(A,1));
m = length(i);
nswap = 0;
while nswap < maxswap
    while 1
        e1 = randi(m); e2 = randi(m);
        while e2 == e1
            e2 = randi(m);
        end
        a = i(e1); b = j(e1);
        c = i(e2); d = j(e2);
        if all(a~=[c,d]) && all(b~=[c,d])
            break
        end
    end
    if rand > 0.5
        i(e2) = d; j(e2) = c;
        c = i(e2); d = j(e2);
    end
    if ~(A(a,d) || A(c,b))
        A(a,d) = A(a,b); A(a,b) = 0;
        A(d,a) = A(b,a); A(b,a) = 0;
        A(c,b) = A(c,d); A(c,d) = 0;
        A(b,c) = A(d,c); A(d,c) = 0;
        j(e1) = d;
        j(e2) = b;
        nswap = nswap + 1;
    end
end
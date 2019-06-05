% Copyright (c) 2017, Evangelos Evangelou
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
%
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the distribution
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.


function [ maxeig, pfe, ncore, naftercore, acssize, nodecore, ...
           fitness, Af, sedge ] = ...
  getstatyr5( year, fdr, folder, alllb )
%GETSTATYR Summary of this function goes here
%   Detailed explanation goes here

%folder = '../Received/zoomed';
pv=dlmread([folder,'/percentile_',int2str(year),'_4.tsv']);
A=dlmread([folder,'/matrix_',int2str(year),'_4.tsv']);
labels=textread([folder,'/labels_',int2str(year),'_4.tsv'],'%s');
A(logical(eye(size(A)))) = 0;
pv(logical(eye(size(pv)))) = 1;
%alllb = textread([folder,'/all_labels.tsv'],'%s');
N=numel(alllb);
Af=zeros(N);
[ll,ii]=ismember(labels,alllb);
ii = ii(ii > 0); % Exclude those labels not found in alllb
A = A(ll,ll); pv = pv(ll,ll); % Exclude undesired labels
h = fdr_bh(pv, fdr);
Ap=A.*h;
% Af(ii,ii) = Ap;
Af(ii,ii) = Ap > 0; % Using binary matrix
sedge = zeros(N);
sedge(ii,ii) = h;

[blocks,growth,core,mingrowth,ismin,maxgrowth,ismax] = ...
  nodestrength(Af);

nodecore = zeros([N,1]);
for b = 1:numel(blocks)
  nodecore(blocks{b}) = core(b);
end

% ACS
maxeig = maxgrowth(1);
ncore = sum(cellfun(@numel,blocks(core == 2)));
naftercore = sum(cellfun(@numel,blocks(core == 1)));
acssize = ncore + naftercore;

% Fitness
fitness = zeros(N,1);
for b = 1:numel(blocks)
    fitness(blocks{b}) = growth(1,b);
end

pfe = pfecalc(Af);
end

function v = pfecalc(A)
A = sparse(A);
[v,~] = eigs(A',1);
if (max(v) < abs(min(v)))               % All negative
  v = -v;
end
v(v<0) = 0;
v = v./sum(v);
end


function [blocks,growth,core,...
   mingrowth,ismin,maxgrowth,ismax] = ...
   nodestrength(G)
%Given a graph G, finds the strongest nodes.

s = size(G,1);                          %  Number of nodes
% Get irreducible blocks (Tarjan's algorithm)
E = sparse(G); % MATLAB expects a sparse matrix
[n,idx] = graphconncomp(E); % Num of blocks and Block id
blocks = cell(n,1);
for i = 1:n
  blocks{i} = find(idx == i);
end
% iscycle = false(n,1); % Is the block a cycle?
% for i = 1:n
%   iscycle(i) = numel(blocks{i}) > 1 || E(blocks{i},blocks{i});
% end

% Construct irreducible blocks matrix
D = sparse(n,n);
for i = 1:n-1
  for j = i+1:n
    if(any(any(E(blocks{j},blocks{i})))) % Link from block j to block i
      D(j,i) = 1;
    elseif(any(any(E(blocks{i},blocks{j})))) % Link from block i to block j
      D(i,j) = 1;
    end
  end
end

core = zeros(n,1); % Is this block a core? 0 = not strong; 1 = periphery; 2 = core
mingrowth = [sum(sum(abs(G))); 0];
maxgrowth = [0; 0];
growth = repmat([0;-1],[1,n]); % The growth rate of each block
order = graphtopoorder(D); % Order downstream

for i = 1:n
  ii = order(i); % Current block
  iparents = find(D(:,ii) ~= 0); % Immediate parents of ii
  np = numel(iparents);
  for j = 1:np
    jj = iparents(j); % jth parent block
    growth(:,ii) = updtgrowth(growth(:,ii),growth(:,jj));
  end
  [growth(:,ii), core(ii)] = growthrate(E(blocks{ii},blocks{ii}),growth(:,ii));
  if np == 0 % Can find lowest growth among head nodes
    if growth(1,ii) < mingrowth(1) || ...
        ((growth(1,ii) - mingrowth(1) < 1e-6) && growth(2,ii) < mingrowth(2))
      mingrowth = growth(:,ii);
    end
  end
  if growth(1,ii) > maxgrowth(1) || ...
      ((maxgrowth(1) - growth(1,ii) < 1e-6) && growth(2,ii) > maxgrowth(2))
    maxgrowth = growth(:,ii);
  end
end

% Find the block with the min growth
lmin = (growth(1,:) - mingrowth(1) < 1e-6) & ...
  (growth(2,:) - mingrowth(2) < 1e-6);
weakest = cat(2,blocks{lmin});
ismin = false(s,1);
ismin(weakest) = true;

% Find the block with the max growth
if 0 == maxgrowth(1) && 0 == maxgrowth(2)
  lmax = [];
  strongest = [];
  core = zeros(n,1);
  ismax = false(s,1);
else
  lmax = (maxgrowth(1) - growth(1,:) < 1e-6) & ...
    (maxgrowth(2) - growth(2,:) < 1e-6);
  if all(lmax) && maxgrowth(1) > 1 + 1e-6 % lamdaPF > 1 so we need the eigenvector
%     leafs = find(sum(~D,2)); % Nodes without outgoing link
%     nleafs = numel(leafs);
%     heads = find(sum(~D,1)); % Nodes without incoming link (core)
%     nheads = numel(heads);
%     dleaf = n*ones(nleafs,1); % How far from the core is each leaf?
%     for i = 1:nleafs
%       for j = 1:nheads
%         dst = graphshortestpath(D,heads(j),leafs(i),'Method','Acyclic');
%         dleaf(i) = min(dleaf(i),dst);
%       end
%     end
%     imin = leafs(dleaf == min(dleaf));
%     weakest = cat(2,blocks{imin});
%     lmin = false(n,1);
%     lmin(imin) = true;
%     lmax = ~lmin;
%     strongest = cat(2,blocks{lmax});
    [ismin] = eigenvecmax ( E, maxgrowth(1) );
    weakest = find(ismin);
    ismax = ~ismin;
    strongest = find(ismax);
  else
    strongest = cat(2,blocks{lmax});
    core(~lmax) = 0;
    ismax = false(s,1);
    ismax(strongest) = true;
  end
end
end


function growth = updtgrowth(old,new)
%updtgrowth - Check if new growth is better than old and update
%
% Syntax: growth = updtgrowth(old,new)
%
% old and new are vectors of length 2. The first element is the eigenvalue exponent and the second is the polynomial power. The first element is tested first and if it is found to be higher, then the function returns new, otherwise if they are equal, the second element is tested, otherwise it returns old.
  if new(1) > old(1)
    growth = new;
  elseif new(1) < old(1)
    growth = old;
  elseif new(2) > old(2)
    growth = new;
  else
    growth = old;
  end
end

function [ growth, core ] = growthrate( graph, growth )
%growthrate For a given graph, compute its growth rate. Input growth is the
%current growth in the path
%   Detailed explanation goes here
lambda = maxeig(graph);
if (growth(1) == lambda) % This node is growing with the same exponent
  growth(2) = growth(2) + 1;
  core = 1 + (lambda > 0); % Core node
elseif (growth(1) < lambda) % This node is growing faster
  growth = [lambda, 0];
  core = 2; % Core node
else
  core = 1; % Not a core node
end
end

function lambda = maxeig( graph )
cgraph = graph ~= 0;
csum = sum(cgraph);
if (all(csum == 1)) % Number of links = number of nodes
  lambda = 1;
elseif (all(csum == 0)) % No links
  lambda = 0;
else
  ee = eig(full(graph));
  lambda = max(abs(ee));
  %lambda = abs(eigs(graph,1));
end
end

function [lmin, lmax] = eigenvecmax ( graph, lambda )
[ V, L ] = eig(full(graph)');
evals = diag(L);
evalsr = real(evals);
evalsi = imag(evals);
ii = lambda - evalsr < 1e-6;
lmult = sum(ii);
if lmult == 1
  X = abs(V(:,ii));
else
  VV = V(:,ii);
  VI = imag(VV);
  VR = real(VV);
  ok = all(VI < 1e-6) & (all(VR <= 0) | all(VR >= 0));
  X = max(abs(VR(:,ok)),[],2);
end
Xmin = min(X);
lmin = X - Xmin < 1e-6;
Xmax = max(X);
lmax = Xmax - X < 1e-6;
end

% We use the MUNKRES (aka Hungarian Algorithm) implementation by Yi Ciao,
% downloaded from 
% https://www.mathworks.com/matlabcentral/fileexchange/20652-hungarian-algorithm-for-linear-assignment-problems-v2-3
% See "munkres" folder for the code and the disclaimer. 
% NOTE: We also tried LAPJV. LAPJV is sometimes faster, but it is 
% not stable numerically. (for some instances, it takes considerably more 
% time than MUNKRES, or sometimes the experiment stops prematurely)

function P = ProjPermMatrix(P)
% projects onto set of permutation matrices
n = size(P,1);
order = munkres(max(P(:))-P); % equivalent to >order = munkres(-P);
% order = lapjv(max(P(:))-P); % equivalent to >order = lapjv(-P);
P = sparse(1:n,order,ones(n,1)); 
% for small problems it can be more efficient not to use sparse
end
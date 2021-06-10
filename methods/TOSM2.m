function [Z, info, Zrdd] = TOSM2(A,B,varargin)
%% 3 operator splitting (Split 2)
% Implements TOS (Split2) for our QAP experiments
%
% [YMS21] A. Yurtsever, V. Mangalick, and S. Sra,
% "Three Operator Splitting with a Nonconvex Loss Function"
% International Conference on Machine Learning, 2021
% 
% contact information: https://github.com/alpyurtsever

n = size(A,1);

tol = -inf;
maxit = 1e5;
X = ones(n,n)./n;
if ~isempty(varargin)
    for tt = 1:2:length(varargin)
        switch lower(varargin{tt})
            case 'tol'
                tol = varargin{tt+1};
            case 'maxit'
                maxit = varargin{tt+1};
            case 'x0'
                X = varargin{tt+1};
            otherwise
                warning(['Unknown option: ',varargin{tt}]);
        end
    end
end

if ~issymmetric(A) || ~issymmetric(B)
    grad = @(X) 0.5.*(A*X*B + A'*X*B');
else
    grad = @(X) A*X*B;
end
proj1 = @(X) max(X,0);
proj2 = @(X) ProjCSSum(X);

info.iter = [];
info.time = [];
info.gap = [];
info.feas = [];
info.obj = [];
info.rddobj = [];

if issparse(A), nA = svds(A,1); else, nA = norm(A); end
if issparse(B), nB = normest(B); %nB = svds(B,1);
else, nB = norm(B); end
L = nA*nB;
if L == 0, L = 1e-6; end
s = 1/L; 
l = 1; 
Zrdd = [];

tic;
Y = X;
for t = 1:maxit
    Z = proj2(Y);
    D = grad(Z);
    
    X = proj1(2*Z - Y - s*D);
    Y = Y + l*(X - Z);
    
    if (2^floor(log2(t)) == t) || (t == maxit)
        
        info.iter(end+1,1) = t;
        H = LAP(D);
        info.gap(end+1,1) = D(:)'*(Z(:)-H(:));
        info.feas(end+1,1) = norm(X - Z,'fro');
        info.obj(end+1,1) = iprod(A'*Z,Z*B);
        Zrdd = ProjPermMatrix(Z);
        info.rddobj(end+1,1) = iprod(A'*Zrdd,Zrdd*B);
        info.time(end+1,1) = toc;
        
        % check stopping criterion
        stopErr1 = abs(info.gap(end)) / max(1,abs(info.obj(end)));
        stopErr2 = info.feas(end)/sqrt(n);
        fprintf('Iteration: %d | Tol: %e, Err1: %e, Err2: %e \n', t, tol, stopErr1, stopErr2);
        if max(stopErr1,stopErr2) <= tol
            break;
        end
    end
    
end

end

function X = ProjCSSum(X)
n = size(X,1);
sX1 = sum(X,1);
sX2 = sum(X,2);
sXa = sum(sX1);
X = X + (1/n + sXa/n^2) - ((1/n)*sX2) - ((1/n)*sX1);
end

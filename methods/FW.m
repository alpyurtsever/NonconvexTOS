function [X,info,Xrdd] = FW(A,B,varargin)
%% Frank-Wolfe
% Implements the Frank-Wolfe algorithm for our QAP experiments
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
% grad = @(X) 0.5.*(A*X*B' + A'*X*B);
obj = @(X) iprod(A'*X,X*B); % = trace(X'*A*X*B); % These are same but iprod is more efficient

info.iter = [];
info.time = [];
info.obj = [];
info.rddobj = [];
info.gap = [];
info.failFlag = false;

tic;
for t = 1:maxit
    
    % find the gradient
    G = grad(X);
    % solve the linear assignment subproblem
    H = LAP(G);
    % update direction
    D = H-X;
            
    if (2^floor(log2(t)) == t) || t == maxit
        
        % compute error
        Xrdd = ProjPermMatrix(X);
        info.iter(end+1,1) = t;
        info.gap(end+1,1) = -iprod(G,D);
        info.obj(end+1,1) = obj(X);
        info.rddobj(end+1,1) = obj(Xrdd);
        info.time(end+1,1) = toc;
        
        % check stopping criterion
        stopErr = abs(info.gap(end)) / max(1,abs(info.obj(end)));
        fprintf('Iteration: %d | Tol: %e, Err: %e \n', t, tol, stopErr);
        if stopErr <= tol
            break;
        end
    end
    
    % find the step-size
    % Here, "obj = iprod(X+eta*D, A*(X+eta*D)*B)"
    % We can rewrite this as
    % eta^2 <D,ADB> + eta <D,AXB> + eta <X,ADB> + <X,AXB>
    % if "<D,ADB> > 0", this quadratic is increasing, to find the minimum,
    % we can take the gradient wrt eta and equate to 0.
    % ==> 2 eta <D,ADB> + <D,AXB> + <X,ADB> = 0;
    % otherwise, the quadratic is decreasing so the minima occurs either at
    % "eta = 0" or "eta = 1".
    % Note, however, if "eta = 0", this would imply X is a stationary point
    % and "gap = 0". Since we already checked gap, this cannot be the case,
    % so we choose "eta = 1".
    tA = iprod(D,A*D*B);
    if tA > 0
        tB = iprod(X,A*D*B) + iprod(D,A*X*B);
        eta = (-tB)/(2*tA);
        eta = max(min(eta,1),0);
    else
        eta = 1;
        % you can check eta = 0, but it will not happen unless we have a bug!
    end
    if eta == 0 
        % compute error
        Xrdd = ProjPermMatrix(X);
        info.iter(end+1,1) = t;
        info.gap(end+1,1) = -iprod(G,D);
        info.obj(end+1,1) = obj(X);
        info.rddobj(end+1,1) = obj(Xrdd);
        info.time(end+1,1) = toc;
        stopErr = abs(info.gap(end)) / max(1,abs(info.obj(end)));
        fprintf('Iteration: %d | Tol: %e, Err: %e \n', t, tol, stopErr);
        info.failFlag = true;
        break;
    end
    
    % update the decision variable
    X = (1-eta)*X + eta*H;
    
end

info.totaliter = t;

end


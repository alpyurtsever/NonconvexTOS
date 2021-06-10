function X = ApproxProjBirkhoff(Xin,maxit)
% Approximates the projection of Xin onto the Birkhoff Polytope.
Y = Xin;
s = 1; l = 1;
for iter = 1:maxit
    Z = Proj1(Y);
    X = Proj2(2*Z - Y - s*(Z - Xin));
    Y = Y + l*(X - Z);
end
end

function X = Proj1(X)
X = max(X,0);
end

function X = Proj2(X)
n = size(X,1);
e = ones(n,1);
n = size(X,1);
Z = (1/n + sum(X(:))/n^2)*eye(n) - (1/n)*X;
X = X + sum(Z,2)*e' - ((1/n)*e)*(X'*e)';
end
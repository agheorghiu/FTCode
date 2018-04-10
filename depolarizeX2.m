function out = depolarizeX2(p, n, k, rho)

X = sparse([0 1; 1 0]);
bigX = kron(speye(2^(k-1)), kron(X, speye(2^(n - k))));

out = (1-p) * rho + p * (bigX * rho * bigX);

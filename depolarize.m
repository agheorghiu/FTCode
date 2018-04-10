function out = depolarize(p, n, k, rho)

X = sparse([0 1; 1 0]);
Y = sparse([0 -sqrt(-1); sqrt(-1) 0]);
Z = sparse([1 0; 0 -1]);

bigX = kron(speye(2^(k-1)), kron(X, speye(2^(n - k))));
bigY = kron(speye(2^(k-1)), kron(Y, speye(2^(n - k))));
bigZ = kron(speye(2^(k-1)), kron(Z, speye(2^(n - k))));

out = (1-3*p/4) * rho + (p/4) * ((bigX * rho * bigX) + (bigY * rho * bigY) + (bigZ * rho * bigZ));

function [chi] = chi_Gaussian(x, Sigma)

Nmode = size(Sigma,1)/2;
% Matrix Omega
In = eye(Nmode);
Ome0 = [0,1;-1,0];
Ome = kron(In, Ome0);

chi = exp(-1/2 * x' * Ome * Sigma * Ome' * x);

end
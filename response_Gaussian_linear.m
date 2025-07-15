function [response] = response_Gaussian_linear(xv, pv, Sigma)

Nx = max(size(xv));
Np = max(size(pv));
response = zeros(Nx, Np);
% Matrix Omega
Omega = [0, 1;-1, 0];

% Number of modes
K = size(Sigma,1)/2;

for nx = 1:Nx
    for np = 1:Np
        x = xv(nx);
        p = pv(np);
        xtable = [-x,p,x,-p];
        chi_in = zeros(2*K, 1);

        for i = 2:K-1
            chi_in(2*i-1) = xtable(mod(i-2,4)+1);
        end

        chi_in(1) = -p;
        chi_in(2) = -x;
        chi_in(2*K-1:2*K) = [x,p]*Omega^K;

        response(nx, np) = chi_Gaussian(chi_in, Sigma);
    end
end

end
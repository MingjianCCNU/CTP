function [response] = response_Gaussian(xv, pv, Sigma)

Nx = max(size(xv));
Np = max(size(pv));
response = zeros(Nx, Np);

N = size(Sigma,1)/2;
% Matrix Omega

for nx = 1:Nx
    for np = 1:Np
        x = xv(nx);
        p = pv(np);
        chi_in = zeros(2*N, 1);

        for i = 2:N-1
            chi_in(2*i-1) = x/(N-2);
        end

        chi_in(1) = p;
        chi_in(2) = x;
        chi_in(2*N-1) = -p;
        chi_in(2*N) = x;

        % TMSV case: testing
        % chi_in(1) = p;
        % chi_in(2) = x;
        % chi_in(2*N-1) = x;
        % chi_in(2*N) = p;

        response(nx, np) = chi_Gaussian(chi_in, Sigma);
    end
end


end
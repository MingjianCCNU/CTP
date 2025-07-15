function [CF] = CF_Gaussian(xi, Sigma)
% input: xi_a, xi_b, xi_c, xi_d, ...

Nmode = size(Sigma,1)/2;

Chi_in = zeros(2*Nmode, 1);

for i = 1:Nmode
    Chi_in(2*i-1) = real(xi(i));
    Chi_in(2*i) = imag(xi(i));
end

CF =  chi_Gaussian(Chi_in, Sigma);

end
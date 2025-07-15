clc;clear
K = 4;

M = rand(2*K);
nphotons = 5; % number for photons subtracted for the two modes

psv = zeros(2*K, 1);
psv(1:2) = nphotons;
psv(end-1:end) = nphotons;


loopNumber = (nphotons+1)^8;
p = zeros(8, 1);

true_mat = false(2*K);
true_mat(1,1) = true;
true_mat(1,2) = true;
true_mat(2,2) = true;
true_mat(1,2*K-1) = true;
true_mat(2,2*K) = true;
true_mat(2*K-1,2*K) = true;
true_mat(2*K,2*K) = true;
true_mat(2*K-1,2*K-1) = true;

tic
Hsum = 0;
for k = 1:loopNumber
    % fit p into matrix
    n2 = zeros(2*K, 2*K);
    n2(true_mat) = p;

    Hsum = Hsum + Hsingle(n2, psv, rand(2*K, 1), M, K);

    % Increase the index vector:
    for ip = 1:8
      if p(ip) < nphotons
        p(ip) = p(ip) + 1;
        break;           % Stop "for ip" loop
      end
      p(ip) = 0;         % Reset this index
   end
end
toc

function [H] = Hsingle(n2, n, x, M, K)

H = 1;

q = zeros(2*K, 1);

for i = 1:2*K
    q(i) = n(i)-n2(i,i);
    for j = 1:2*K
        q(i) = q(i) - n2(min(i,j), max(i,j));
    end
end

for i = 1 : (2*K)
    if q(i) >= 0
        H = H*n(i)*x(i)^q(i)/factorial(q(i));
        for j = (i+1) : (2*K)
            H = H * (2*M(i,j))^n2(i,j)/factorial(n2(i,j));
        end
    else
        H = 0;
    end

end






end
function [H] = HMatIn(n2, n, x, M, K)
% x: 输入自变量矩阵

% 判断哪些ij值需要去计算H
q = zeros(2*K, 1);
for i = 1:2*K
    q(i) = n(i)-n2(i,i);
    for j = 1:2*K
        q(i) = q(i) - n2(min(i,j), max(i,j));
    end
end

Nplot = size(x,2);
H = ones(Nplot);
if sum(q<0) == 0 % 若输入的{n_ij}不导致{q_i}<0
    for i = 1 : (2*K)
        H = H.*factorial(n(i)).*reshape(x(i,:,:),[Nplot,Nplot]).^q(i)/factorial(q(i))...
            * M(i,i)^n2(i,i)/factorial(n2(i,i));
        for j = (i+1) : (2*K)
            H = H .* (2*M(i,j))^n2(i,j)/factorial(n2(i,j));
        end
    end
else
    H = 0;
end


end
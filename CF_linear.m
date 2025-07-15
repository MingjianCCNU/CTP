function [response] = CF_linear(x,p,Vt,Matin,K,psv,nphotons)
% 用于求归一化系数 - 废弃

In = eye(K);
X = kron(In, [0,1;1,0]);

xtable = [-x,p,x,-p];
Ftxtable = [-p,-x,p,x];
Ftptable = [x,-p,-x,p];
xin = zeros(2*K, 1);
xin(1:2) = [-p-1i*x, -p+1i*x];
for i = 2:K-1
    xin(2*i+1:2*i+2) = [xtable(mod(i-2,4)+1),...
        xtable(mod(i-2,4)+1)];
end
xin(end-1:end) = [...
    Ftxtable(mod(K+1,4)+1)-1i*Ftptable(mod(K+1,4)+1),...
    Ftxtable(mod(K+1,4)+1)+1i*Ftptable(mod(K+1,4)+1)];

% TMSV: testing case
% xin = [p+1i*x; p-1i*x; x+1i*p; x-1i*p];

xin = X*Vt*xin;

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


response = 0;
for k = 1:loopNumber
    % fit p into matrix
    n2 = zeros(2*K, 2*K);
    n2(true_mat) = p;

    response = response + Hsingle(n2, psv, xin, Matin, K);

    % Increase the index vector:
    for ip = 1:8
        if p(ip) < nphotons
            p(ip) = p(ip) + 1;
            break;           % Stop "for ip" loop
        end
        p(ip) = 0;         % Reset this index
    end
end

end
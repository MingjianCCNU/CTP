function [response] = CFMatIn(xin,Vt,Matin,K,psv,nphotons)
% xin: 相空间坐标向量输入

In = eye(K);
X = kron(In, [0,1;1,0]);

Nx = size(xin,2);
Np = size(xin,2);

xin = reshape(pagemtimes(X*Vt, reshape(xin, [2*K, 1, Nx, Np])), [2*K, Nx, Np]);


% 考虑哪些元素取决于对哪些模施加PS
Dmat = abs(Matin);
Dmattf = Dmat>1e-10;


% psmat = zeros(2*K);
% psmat(:, true&psv) = true;
% psmat(true&psv, :) = true;
% Dmattf = psmat&triu(Dmattf);

Dmattf = triu(Dmattf);
Dmattf = Dmattf&psv*psv';


loopNumber = (nphotons+1)^sum(Dmattf,'all');
p = zeros(sum(Dmattf,'all'), 1);

response = 0;
for k = 1:loopNumber
    % fit p into matrix
    n2 = zeros(2*K, 2*K);
    n2(Dmattf) = p;

    response = response + HMatIn(n2, psv, xin, Matin, K);

    % Increase the index vector:
    for ip = 1:sum(Dmattf,'all')
        if p(ip) < nphotons
            p(ip) = p(ip) + 1;
            break;           % Stop "for ip" loop
        end
        p(ip) = 0;         % Reset this index
    end
end

end
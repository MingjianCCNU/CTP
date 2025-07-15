clc;clear

K = 5;
nphotons = 1;
psv = zeros(2*K, 1);
psv(1:2) = nphotons;
psv(end-1:end) = nphotons;
n = psv;

load("C:\Users\PS\Documents\BaiduSyncdisk\0. 科研\0. 文章\2025\Matlab codes\Mtest.mat")
Matin = M;

% 考虑哪些元素取决于对哪些模施加PS
Dmat = abs(Matin);
Dmattf = Dmat>1e-10;
Dmattf = triu(Dmattf);

Dmattf = Dmattf&psv*psv';


loopNumber = (nphotons+1)^sum(Dmattf,'all');
p = zeros(sum(Dmattf,'all'), 1);

x = rand(10,1);


response = 0;
for k = 1:loopNumber
    % fit p into matrix
    n2 = zeros(2*K, 2*K);
    n2(Dmattf) = p;

    q = zeros(2*K, 1);
    for i = 1:2*K
        q(i) = n(i)-n2(i,i);
        for j = 1:2*K
            q(i) = q(i) - n2(min(i,j), max(i,j));
        end
    end

    H = 1;
    if sum(q<0) == 0 % 若输入的{n_ij}不导致{q_i}<0
        for i = 1 : (2*K)
            H = H*factorial(n(i))*x(i)^q(i)/factorial(q(i))...
                * M(i,i)^n2(i,i)/factorial(n2(i,i));
            for j = (i+1) : (2*K)
                H = H * (2*M(i,j))^n2(i,j)/factorial(n2(i,j));
            end
        end
    else
        H = 0;
    end

    response = response + H;

    % Increase the index vector:
    for ip = 1:sum(Dmattf,'all')
        if p(ip) < nphotons
            p(ip) = p(ip) + 1;
            break;           % Stop "for ip" loop
        end
        p(ip) = 0;         % Reset this index
    end
end






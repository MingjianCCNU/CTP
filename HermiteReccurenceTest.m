clc;clear all
% 艾米尔特函数递归性质测试

K = 5;
N = 2;
Matin = linspace(1,2,2*K);
Matin = Matin + Matin';
M = Matin;
xin = linspace(1,2,2*K);

% global M x;
x = zeros([2*K, N, N]);

for n1 = 1:N
    for n2 = 1:N
        x(:,n1,n2) = xin'+0.001*n1;
    end
end

psv = zeros(2*K, 1);
psv(1) = 1;
psv(2) = 0;

nphotons = 10;

% response_2110 = GetH(xin,Matin,K,[10 10 10 10]',nphotons);
% toc
% response_2100 = GetH(xin,Matin,K,[2 1 0 0]',nphotons);
% response_1000 = GetH(xin,Matin,K,[1 0 0 0]',nphotons);
% response_0100 = GetH(xin,Matin,K,[0 1 0 0]',nphotons);
% response_1100 = GetH(xin,Matin,K,[1 1 0 0]',nphotons);
% response_0000 = GetH(xin,Matin,K,[0 0 0 0]',nphotons);

nin = [10 10 10 10 0 0 0 0 0 0];
tic
for i = 1:N^2
    t1 = multiHermite(nin, xin+0.001, Matin+0.001);
    clear multiHermite
end
toc

% tic
% for i = 1:N^2
%     t2 = GetH(xin,Matin,K,nin',5);
% end
% toc

tic
t3 = multiHermite_vectorized(nin, x, M);
toc

t3(1)


% response_1000_rec = xin(1) * response_0000;
% response_0100_rec = xin(2) * response_0000;
% response_1100_rec = ...
%     x(2) * response_1000 +...
%     2 * M(2,1) * 1 * response_0000;
% response_2100_rec = ...
%     x(1) * response_1100 +...
%     2 * M(1,1) * 1 * response_0100 +...
%     2 * M(2,1) * 1 * response_1000
    
% response_t_2100 = getHermite(1.01, 1.1, [2,1,0,0], 1)


function [response] = GetH(xin,Matin,K,psv,nphotons)
% 以下代码段用于通过解析式计算一个H
% 考虑哪些元素取决于对哪些模施加PS
Dmat = abs(Matin);
Dmattf = Dmat>1e-10;

Dmattf = triu(Dmattf);
Dmattf = Dmattf&psv*psv';

% Dmattf = triu(Matin)&ones(2*K);

loopNumber = (nphotons+1)^sum(Dmattf,'all');
p = zeros(sum(Dmattf,'all'), 1);

response = 0;
for k = 1:loopNumber
    % fit p into matrix
    n2 = zeros(2*K, 2*K);
    n2(Dmattf) = p;

    response = response + Hsingle(n2, psv, xin, Matin, K);

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

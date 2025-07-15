clc;clear

Nr = 1;
K = 5;
nphotons = 3;

psv = zeros(2*K, 1);
psv(1:2) = 5;
psv(end-1:end) = 5;
rdBbin = -4;

xplim = 4;
xsum = 51;
x = linspace(-xplim,xplim,xsum);
p = x;
dx = diff(x);
[xm, pm] = meshgrid(x, p);

% CF for coherent states
CF_coh = exp(-0.5*(xm.^2+pm.^2));
% sum(CF_coh.*CF_coh',"all")*dx(1)^2; % fidelity test

% Matrix Omega and some constant
In = eye(K);
Ome0 = [0,1;-1,0];
Ome = kron(In, Ome0);

Z = [1,0;0,-1];
Zk = kron(In, Z);
J = kron(In, [1,1i;1,-1i])/2;
X = kron(In, [0,1;1,0]);

% Adjacent matrix for the linear cluster states
Uadj = zeros(K);
Uadj(1,2) = 1;
Uadj(K,K-1) = 1;
for i = 2:K-1
    Uadj(i,i-1) = 1;
    Uadj(i,i+1) = 1;
end

% 文献[1]中U的格式是按照x1,x2,...,p1,p2排列，需转换 [1]PHYSICAL REVIEW A 91, 032314 (2015)
U = (In + 1i*Uadj)/(sqrtm(Uadj*Uadj+In));

Ur = zeros(K);
for i = 1:K
    Ur(2*i-1, i) = 1;
    Ur(2*i, i+K) = 1;
end
% CM for the cluster state made via linear optics
S = [real(U), -imag(U);imag(U), real(U)];
Sr = Ur * S * Ur';

% Calsulate the x,p matrices


Nx = size(x,2);
Np = size(p,2);

xtable =   [-x,p,x,-p];
Ftxtable = [-p,-x,p,x];
Ftptable = [x,-p,-x,p];
[Mx, Mp] = meshgrid(x, p);
% 注意此处xin输入为复变量形式，所以会跟公式稍有不同
% 坐标变量x,p的矩阵化，与r值无关
xin = zeros(2*K, Nx, Np);

xin(1,:,:) = -Mp-1i*Mx;
xin(2,:,:) = -Mp+1i*Mx;
for i = 1:K-2
    switch mod(i,4)
        case 1
            xin(2*i+1,:,:) = -Mx;
            xin(2*i+2,:,:) = -Mx;
        case 2
            xin(2*i+1,:,:) = Mp;
            xin(2*i+2,:,:) = Mp;
        case 3
            xin(2*i+1,:,:) = Mx;
            xin(2*i+2,:,:) = Mx;
        case 0
            xin(2*i+1,:,:) = -Mp;
            xin(2*i+2,:,:) = -Mp;
    end
end

switch mod(K,4)
    case 1
        xin(2*K-1,:,:) = -Mp+1i*Mx;
        xin(2*K,:,:) = -Mp-1i*Mx;
    case 2
        xin(2*K-1,:,:) = -Mx-1i*Mp;
        xin(2*K,:,:) = -Mx+1i*Mp;
    case 3
        xin(2*K-1,:,:) = Mp-1i*Mx;
        xin(2*K,:,:) = Mp+1i*Mx;
    case 0
        xin(2*K-1,:,:) = Mx+1i*Mp;
        xin(2*K,:,:) = Mx-1i*Mp;
end

tic
Fps = zeros(Nr, 1);
F = Fps;
for pindex = 1:Nr
    r = rdBbin(pindex)/8.69;
    % Initial covariance matrix
    Sigma_sqz = [exp(-2*r),0;0,exp(2*r)];
    Sigma0 = kron(In, Sigma_sqz);
    V = Sr * Sigma0 * Sr';

    % 此处ZJ的操作实际上是把CM从x,p转换至xi,xi*
    % Vt = Zk*J*V*J'*Zk-diag(psv/nphotons)/2;
    Vt = Zk*J*V*J'*Zk - eye(2*K)/2;
    Matin = -0.5*X*Vt;

    % get normalization constant

    eta = response_Gaussian_linear(x, p, V);
    Hpoly = CFMatIn(xin,Vt,Matin,K,psv,nphotons)';

    normConstant = Hpoly((xsum+1)/2, (xsum+1)/2);

    eta_ps = real(Hpoly/normConstant.*eta);
    ratio_ps = real(Hpoly/normConstant);
    Fps(pindex) = sum(eta_ps .* CF_coh .* CF_coh','all')*dx(1)^2/pi;
    F(pindex) = sum(eta .* CF_coh .* CF_coh','all')*dx(1)^2/pi;
end
toc

%%
% data = readmatrix("C:\Users\PS\Documents\BaiduSyncdisk\0. 科研\0. 文章\2025\data.csv");
% 
% surf(data)
% figure
% surf(ratio_ps)
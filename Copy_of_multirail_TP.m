clc;clear
Nr = 1;
K = 5;
nphotons = 10;
psv = zeros(2*K, 1);
psv(3:4) = nphotons;
psv(end-1:end) = nphotons;
rdBbin = -4;


xplim = 4;
xnum = 51;
x = linspace(-xplim,xplim,xnum);
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

% Adjacent matrix for the cluster states
if K ~= 2
    Uadj = zeros(K);
    Uadj(1,2:K-1) = ones(1,K-2);
    Uadj(end,2:K-1) = ones(1,K-2);
else
    Uadj = [0,1;1,0];
end

for i = 2:K-1
    Uadj(i,1) = 1;
    Uadj(i,end) = 1;
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

[Mx, Mp] = meshgrid(x, p);
% 注意此处xin输入为复变量形式，所以会跟公式稍有不同
% 坐标变量x,p的矩阵化，与r值无关
xin = zeros(2*K, Nx, Np);
xin(1,:,:) = -Mp-1i*Mx;
xin(2,:,:) = -Mp+1i*Mx;
for i = 1:K-2
    xin(2*i+1,:,:) = -Mx/(K-2);
    xin(2*i+2,:,:) = -Mx/(K-2);
end
xin(2*K-1,:,:) = Mp-1i*Mx;
xin(2*K,:,:)   = Mp+1i*Mx;

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
    Vt = Zk*J*V*J'*Zk-diag(psv/nphotons)/2;
    Matin = -0.5*X*Vt;

    % get normalization constant
    normConstantp = CF(0,0,Vt,Matin,K,psv,nphotons);

    eta = response_Gaussian(x, p, V);

    Hpoly = CFMatIn(xin,Vt,Matin,K,psv,nphotons);
    normConstant = real(Hpoly((xnum+1)/2, (xnum+1)/2));

    eta_ps = real(Hpoly/normConstant.*eta);
    ratio_ps = real(Hpoly/normConstant);
    Fps(pindex) = sum(eta_ps .* CF_coh .* CF_coh','all')*dx(1)^2/pi
    F(pindex) = sum(eta .* CF_coh .* CF_coh','all')*dx(1)^2/pi
end
toc

%%
data = readmatrix("C:\Users\PS\Documents\data.csv");

surf(data)
figure
surf(eta_ps)
function [return_val, O, V, Sr] = OptFunction(theta,rdBbin,PS_flag)
% rdBbin = [-7, -6, -4, 0];
% theta = pi/4*ones(6,1);
% theta = zeros(6, 1);

rbin = rdBbin/8.69;

% Cluster state parameters
N_mode = 4;
K = N_mode;
I = eye(N_mode);

% nphotons: number of photons subtracted
nphotons = 1;
psv = zeros(2*K, 1);
psv(1:2) = nphotons;
psv(7:8) = nphotons;

% 用于计算Fidelity的grid
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

% adjacency matrix that defines a cluster state
Uadj = [0 1 1 0;...
        1 0 0 1;...
        1 0 0 1;...
        0 1 1 0];


% Rotation matrix - defined by theta
O = I;
k = 1;
for i = 1:N_mode-1
    for j = i+1:N_mode
        Rt = I;
        ag = theta(k);
        Rt(i,i) = cos(ag);
        Rt(i,j) = -sin(ag);
        Rt(j,i) = sin(ag);
        Rt(j,j) = cos(ag);
        O = O * Rt;
        k = k + 1;
    end
end


% 文献[1]中U的格式是按照x1,x2,...,p1,p2排列，需转换 [1]PHYSICAL REVIEW A 91, 032314 (2015)
U = (In + 1i*Uadj)/(sqrtm(Uadj*Uadj+In));

% 乘上旋转矩阵O
U = U * O;

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


% Initial covariance matrix for the squeezed modes
Sigma0_diag = zeros(1, 2*N_mode);
Sigma0_diag(1:2:end) = exp(-2*rbin);
Sigma0_diag(2:2:end) = exp(2*rbin);
Sigma0 = diag(Sigma0_diag);

V = Sr * Sigma0 * Sr';

% 此处ZJ的操作实际上是把CM从x,p转换至xi,xi*
Vt = Zk*J*V*J'*Zk-diag(psv/nphotons)/2;
Matin = -0.5*X*Vt;
xintiao = reshape(pagemtimes(X*Vt, reshape(xin, [2*K, 1, Nx, Np])), [2*K, Nx, Np]);

eta = response_Gaussian(x, p, V);

Hpoly = multiHermite_vectorized(psv, xintiao, Matin)';
multiHermite_vectorized([], [], [], true);

normConstant = Hpoly((xsum+1)/2, (xsum+1)/2);

eta_ps = real(Hpoly/normConstant.*eta);
% ratio_ps = real(Hpoly/normConstant);

Fps = sum(eta_ps .* CF_coh .* CF_coh','all')*dx(1)^2/pi;
F = sum(eta .* CF_coh .* CF_coh','all')*dx(1)^2/pi;

if PS_flag
    return_val = -Fps;
else
    return_val = -F;
end
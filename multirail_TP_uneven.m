% 采用递归的快速计算版本
% 多径簇态不均等压缩度
clc;clear

opt_flag = 1;

% number for photons subtracted for the two modes
nphotons = 1; 
% number of modes in the cluster state
N_mode = 4; 
K = N_mode;


% 此处不考虑rdB网格化
Nr = 1;
rdBbin = [-4,-4,-4,-4];

% 若设为参数[-7, -6, -4, -4]，可能优化效果还不如不优化，这可能是因为优化目标函数
% 设置为了平均方差，可能存在某些模的方差更重要的情况
% rdBbin = [-7, -6, -4, 0];


% 计算最优化平均归零算符方差旋转矩阵
OptOption = optimoptions(@fmincon,...
    'Algorithm', 'interior-point',...
    'Display','off');

theta_0 = pi/4*ones(6,1);
theta_min = zeros(6,1);
theta_max = pi*ones(6,1);
fun_no = @(theta)OptFunctionOrg(theta, rdBbin);
[theta_best, D_mean_no] = fmincon(fun_no,theta_0,[],[],[],[],theta_min,theta_max,[],OptOption);


% theta_best = [1.572015755981591
% 0.785408461806226
% 0.001337444738917
% 0.001712652426390
% 0.785405386204502
% 0.003719680588406];

[~, O] = OptFunction(theta_best, rdBbin, 1);
if ~opt_flag
    O = eye(4);
end

rbin = rdBbin/8.69;

psv = zeros(2*K, 1);
psv(1:2) = nphotons;
psv(7:8) = nphotons;

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


% 乘上旋转矩阵R
U = U*O;

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

pindex = 1;

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
ratio_ps = real(Hpoly/normConstant);
Fps(pindex) = sum(eta_ps .* CF_coh .* CF_coh','all')*dx(1)^2/pi;
F(pindex) = sum(eta .* CF_coh .* CF_coh','all')*dx(1)^2/pi;

toc


%% 测试各模之间的纠缠度（log-negativity）
mA = 2;
mB = 3;
V_in = zeros(4);
V_in(1:2,1:2) = V((2*mA)-1:(2*mA),(2*mA)-1:(2*mA));
V_in(1:2,3:4) = V((2*mA)-1:(2*mA),(2*mB)-1:(2*mB));
V_in(3:4,1:2) = V((2*mB)-1:(2*mB),(2*mA)-1:(2*mA));
V_in(3:4,3:4) = V((2*mB)-1:(2*mB),(2*mB)-1:(2*mB));
gaussian_entanglement(V_in);

%% 计算特定方差：应与Fidelity高度正相关
weight1 = [0 1 -2 0 -2 0 0 1];
weight2 = [1 0 0 -0.5 0 -0.5 1 0];
cov_sum1 = weight1*V*weight1'; % var(x_out)
cov_sum2 = weight2*V*weight2'; % var(p_out)


GetD(V, Uadj, 4)
display(D_mean_no)
display(F)
display(Fps)
% display(cov_sum1)
% display(cov_sum2)
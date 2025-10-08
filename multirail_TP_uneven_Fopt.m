% 采用递归的快速计算版本
% 多径簇态不均等压缩度
% 旋转因子针对保真度进行优化
clc;clear
N_mode = 4;


% rdBbin = [-8, -8, -8, -8];
% rdBbin = [-18, -18, 0, 0];
% rdBbin = [-18, 0, -0, -10];
% rdBbin = [0, -18, -18, 0];
rdBbin = [-18, -18, -18, -18];
% rdBbin = [-7, -6, -4, 0];
% rdBbin = [-7, -7, -7, -7];
% rdBbin = [-30, -30, -30, -30];


% 计算最优化平均归零算符方差旋转矩阵
OptOption = optimoptions(@fmincon,...
    'Algorithm', 'interior-point',...
    'Display','off');

theta_0 = pi/8*ones(6,1);
theta_min = zeros(6,1);
theta_max = pi*ones(6,1);

PS_Flag = 0;
fun = @(theta)OptFunction(theta, rdBbin, PS_Flag);
[theta_best, F_best] = fmincon(fun,theta_0,[],[],[],[],theta_min,theta_max,[],OptOption);
display(-F_best)


[~, O, V, Sr] = OptFunction(theta_best, rdBbin, PS_Flag);
[~, O0, V0, Sr0] = OptFunction(zeros(6,1), rdBbin, PS_Flag);


% Sr0' * V * Sr0


%%
% adjacency matrix that defines a cluster state
Uadj = [0 1 1 0;...
        1 0 0 1;...
        1 0 0 1;...
        0 1 1 0];
% Initial covariance matrix for the squeezed modes
rbin = rdBbin/8.69;
Sigma0_diag = zeros(1, 2*N_mode);
Sigma0_diag(1:2:end) = exp(-2*rbin);
Sigma0_diag(2:2:end) = exp(2*rbin);
Sigma0 = diag(Sigma0_diag);

% GetD(Sigma0,Uadj,4)
D = GetD(V,Uadj,4);
display(D)

%%
% V = V0;
mA = 1;
mB = 2;
V_in = zeros(4);
V_in(1:2,1:2) = V((2*mA)-1:(2*mA),(2*mA)-1:(2*mA));
V_in(1:2,3:4) = V((2*mA)-1:(2*mA),(2*mB)-1:(2*mB));
V_in(3:4,1:2) = V((2*mB)-1:(2*mB),(2*mA)-1:(2*mA));
V_in(3:4,3:4) = V((2*mB)-1:(2*mB),(2*mB)-1:(2*mB));

Eln = gaussian_entanglement(V_in);
display(Eln)

%% 计算特定方差：应与Fidelity高度正相关
weight1 = [0 1 0 0 0 0 0 -1];
weight2 = [1 0 0 -0.5 0 -0.5 1 0];
cov_sum1 = weight1*V*weight1'; % var(x_out)
cov_sum2 = weight2*V*weight2'; % var(p_out)

display(cov_sum1)
display(cov_sum2)
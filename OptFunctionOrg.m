% 根据mean variance优化旋转角的函数
function [ret, D, O] = OptFunctionOrg(theta,r_dB)
if nargin < 1
    theta = [pi/4 pi/5 pi/6 pi/4 pi/3 pi/6];
end

% Testing paramters
% err = 0.0; 
% r_dB = [-7, -6, -4, 0];
% % r_dB = [-10, -10, -10, -10];
% T = 1;
% % theta = [pi/4 pi/5 pi/6 pi/4 pi/3 pi/6];
% theta = [0 0 0 0 0 0];
% target = 2; % target mode to be sent through the channel

% Cluster state parameters
N_mode = 4;
I = eye(N_mode);

% adjacency matrix that defines a cluster state
% V = [0 1 0 0;...
%     1 0 1 0;...
%     0 1 0 1;...
%     0 0 1 0];

% V = [0 1 1 1;...
%     1 0 0 0;...
%     1 0 0 0;...
%     1 0 0 0];

V = [0 1 1 0;...
    1 0 0 1;...
    1 0 0 1;...
    0 1 1 0];

% Original squeezed state parameters
r = r_dB / 8.6859;

t = zeros(2 * N_mode, 1);

t(1:2:2*N_mode) = exp(-2 * r);
t((1:2:2*N_mode) + 1) = exp(2 * r);


Cov = diag(t);

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

% original transformation matrix for the linear network
U_i = (I+1i*V)/(sqrtm((V^2+I)));
U = U_i * O;

% symplectic transformation
% S = [real(U) -imag(U);imag(U) real(U)];
S = kron(real(U), [1,0;0,1]) + kron(imag(U), [0,-1;1,0]);

Cov_ini = S * Cov * S';

%% variance of the nullifiers
% var(p_i)
D = diag(Cov_ini((1:2:2*N_mode) + 1, (1:2:2*N_mode) + 1));

for n_d = 1:N_mode
    %     var(x_i) and cov(p_i,x_j)
    for n_i = 1:N_mode
        D(n_d) = D(n_d) + V(n_d, n_i) * Cov_ini(2 * n_i - 1, 2 * n_i - 1)...
            - 2 * V(n_d, n_i) * Cov_ini(2 * n_d, 2 * n_i - 1);
    end
    %     cov(x_i,x_j)
    for n_i = 1:N_mode
        for n_j = n_i+1:N_mode
            D(n_d) = D(n_d) + 2 * V(n_d, n_i) * V(n_d, n_j) *...
                Cov_ini(2 * n_i - 1, 2 * n_j - 1);
        end
    end
end


D = D./[3; 3; 3; 3];
for i = 1:2*N_mode
   if Cov_ini(i, i) < 0
      D = 2 * ones(1,N_mode);
      break
   end
end
D_mean = mean(D);

ret = D_mean;

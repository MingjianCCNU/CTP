clc;clear
% 设置压缩参数
r = 3/8.69;

% 生成TMSV态的协方差矩阵
cosh_2r = cosh(2 * r);
sinh_2r = sinh(2 * r);

cov_matrix = [cosh_2r, 0, sinh_2r, 0;
              0, cosh_2r, 0, -sinh_2r;
              sinh_2r, 0, cosh_2r, 0;
              0, -sinh_2r, 0, cosh_2r];

% 计算log-negativity
log_neg = gaussian_entanglement(cov_matrix);

fprintf('压缩参数 r = %.2f\n', r);
fprintf('协方差矩阵:\n');
disp(cov_matrix);
fprintf('Log-negativity E_N = %.4f\n', log_neg);

% 与理论值比较
theoretical_value = 2 * r / log(2);
fprintf('理论值 (2r/ln(2)) = %.4f\n', theoretical_value);
fprintf('差异: %.4f\n', abs(log_neg - theoretical_value));
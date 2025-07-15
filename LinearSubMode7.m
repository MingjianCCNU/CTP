clc;clear all
% 设置所有线条的默认粗细为2点
set(groot, 'DefaultLineLineWidth', 1.5);  
% number of modes in the cluster state
K = 7; 
Nr = 256;

figure
% 头尾两模
% number for photons subtracted for the two modes
nphotons = 1; 
psv = zeros(2*K, 1);
psv(1:2) = nphotons;
psv(end-1:end) = nphotons;
linear_TP_rec
plot(abs(rdBbin), Fps, '-.',Color="#ca4362", Marker="o",...
    MarkerIndices = 1:Nr/8:Nr)
hold on

% 头尾次两模
psv = zeros(2*K, 1);
psv(3:4) = nphotons;
psv(end-3:end-2) = nphotons;
linear_TP_rec
plot(abs(rdBbin), Fps, '-.', Color='#2c73d2', Marker="x",...
    MarkerIndices = 1:Nr/8:Nr)

% 头尾四模
psv = zeros(2*K, 1);
psv(1:4) = nphotons;
psv(end-3:end) = nphotons;
linear_TP_rec
plot(abs(rdBbin), Fps, Color="#8290bb", Marker="o",...
    MarkerIndices = 1:Nr/8:Nr)

% 头尾次四模
psv = zeros(2*K, 1);
psv(3:6) = nphotons;
psv(end-5:end-2) = nphotons;
linear_TP_rec
plot(abs(rdBbin), Fps, Color="#ffc572", Marker="x",...
    MarkerIndices = 1:Nr/8:Nr)

% %头尾六模
% nphotons = 1; 
% psv = zeros(2*K, 1);
% psv(1:6) = nphotons;
% psv(end-5:end) = nphotons;
% linear_TP_rec
% plot(abs(rdBbin), Fps, Color="#ffc572", Marker="x",...
%     MarkerIndices = 1:48:Nr)


% 计算顶部边界（例如曲线的最大值加一个偏移量）
y_top = 1.1;

% 构造填充区域的坐标（先绘制曲线，再连接到顶部边界，最后回到起点）
x_fill = [abs(rdBbin)-0.02, fliplr(abs(rdBbin))-0.02];
y_fill = [F', y_top*ones(size(rdBbin))];

% 填充区域
% plot(x, y, 'LineWidth', 2);  % 绘制原始曲线
% hold on;
fs = 14;
fill(x_fill, y_fill, 'r', 'FaceAlpha', 0.1, EdgeColor='k',LineStyle='--');  % 填充曲线以上区域（蓝色，30%透明度）
text(0.1,0.95,"Enhancement region")
ylim([0.2,1.01])
xlim([0, 14.5])
pbaspect([1, 0.618, 1])
ylabel("${F}$", Rotation = 0, FontSize = fs);
xlabel("$|r_\mathrm{dB}|$", "FontSize", fs)

% % for K = 5:
% legend("$\hat{a}_1\hat{a}_5$","$\hat{a}_1\hat{a}_2$",...
%     "$\hat{a}_2\hat{a}_4$",...
%     "$\hat{a}_1\hat{a}_2\hat{a}_4\hat{a}_5$",Location="southeast")

% for K = 7:
legend("$\hat{a}_1\hat{a}_7$","$\hat{a}_2\hat{a}_6$",...
    "$\hat{a}_1\hat{a}_2\hat{a}_6\hat{a}_7$",...
    "$\hat{a}_2\hat{a}_3\hat{a}_5\hat{a}_6$",Location="southeast")

if 0
    myfigure = gcf;
    figurename = 'fl7nps.svg';
    saveas(myfigure, figurename)
    exportgraphics(myfigure, figurename, ContentType="vector")
end
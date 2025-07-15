clc;clear
% 设置所有线条的默认粗细为2点
set(groot, 'DefaultLineLineWidth', 1.5);  

% number of modes in the cluster state
K = 5; 
Nr = 128;

% number for photons subtracted for the two modes
nphotons = 1; 
% 减单光子
psv = zeros(2*K, 1);
psv(1:2) = nphotons;
multirail_TP_rec
figure
plot(abs(rdBbin), Fps, Color="#2c73d2", Marker="o",...
    MarkerIndices = 1:16:Nr)
hold on

psv = zeros(2*K, 1);
psv(3:4) = nphotons;
multirail_TP_rec
plot(abs(rdBbin), Fps, '-.',Color="#ca4362", Marker="o",...
    MarkerIndices = 1:16:Nr)
hold on

% 减多光子
nphotons = 2; 
psv = zeros(2*K, 1);
psv(1:2) = nphotons;
multirail_TP_rec
plot(abs(rdBbin), Fps, Color="#ffc572", Marker="x",...
    MarkerIndices = 1:16:Nr)

psv = zeros(2*K, 1);
psv(3:4) = nphotons;
multirail_TP_rec
plot(abs(rdBbin), Fps, '-.', Color='#8290bb', Marker="x",...
    MarkerIndices = 1:16:Nr)


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
text(0.1,0.95,"Enhancement region",Position=[0.006711409567799,0.934065943746921,0],Units="normalized")
ylim([0.25,1.01])
xlim([0, 14.9])
pbaspect([1, 0.618, 1])
ylabel("${F}$", Rotation = 0, FontSize = fs);
xlabel("$|r_\mathrm{dB}|$", "FontSize", fs)
legend("$\hat{a}_1$","$\hat{a}_2$","$\hat{a}_1^2$","$\hat{a}_2^2$",Location="southeast")

if 0
    myfigure = gcf;
    figurename = 'f1ps.svg';
    saveas(myfigure, figurename)
    exportgraphics(myfigure, figurename, ContentType="vector")
end
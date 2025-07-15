clc;clear all
% 设置所有线条的默认粗细为1.5
set(groot, 'DefaultLineLineWidth', 1.5);  
% number of modes in the cluster state
K = 10; 
Nr = 512;

% number for photons subtracted for the two modes
nphotons = 1; 

psv = zeros(2*K, 1);
psv(1:2) = nphotons;
psv(end-1:end) = nphotons;
multirail_TP_rec

figure
plot(abs(rdBbin), Fps, '-.',Color="#ca4362", Marker="o",...
    MarkerIndices = 1:64:Nr)
hold on

nphotons = 2; 
psv = zeros(2*K, 1);
psv(1:2) = nphotons;
psv(end-1:end) = nphotons;
multirail_TP_rec
plot(abs(rdBbin), Fps, Color="#2c73d2", Marker="o",...
    MarkerIndices = 1:64:Nr)

nphotons = 3; 
psv = zeros(2*K, 1);
psv(1:2) = nphotons;
psv(end-1:end) = nphotons;
multirail_TP_rec
plot(abs(rdBbin), Fps, Color="#ffc572", Marker="x",...
    MarkerIndices = 1:64:Nr)

nphotons = 5; 
psv = zeros(2*K, 1);
psv(1:2) = nphotons;
psv(end-1:end) = nphotons;
multirail_TP_rec
plot(abs(rdBbin), Fps, '-.', Color="#8290bb", Marker="x",...
    MarkerIndices = 1:64:Nr)

nphotons = 6; 
psv = zeros(2*K, 1);
psv(1:2) = nphotons;
psv(end-1:end) = nphotons;
multirail_TP_rec
plot(abs(rdBbin), Fps, ':', Color="#8290bb",...
    MarkerIndices = 1:64:Nr)

% 计算顶部边界（例如曲线的最大值加一个偏移量）
y_top = 1.1;

% 构造填充区域的坐标（先绘制曲线，再连接到顶部边界，最后回到起点）
x_fill = [abs(rdBbin)-0.03, fliplr(abs(rdBbin))-0.02];
y_fill = [F', y_top*ones(size(rdBbin))];

% 填充区域
% plot(x, y, 'LineWidth', 2);  % 绘制原始曲线
% hold on;
fs = 14;
fill(x_fill, y_fill, 'r', 'FaceAlpha', 0.1, EdgeColor='k',LineStyle='--');  % 填充曲线以上区域（蓝色，30%透明度）
text(0.1,0.95,"Enhancement region",Position=[0.006711409567799,0.934065943746921,0],Units="normalized");
ylim([0.45,1.01])
xlim([0, 14.9])
pbaspect([1, 0.618, 1])
ylabel("${F}$", Rotation = 0, FontSize = fs);
xlabel("$|r_\mathrm{dB}|$", "FontSize", fs)
legend("$\hat{a}_1\hat{a}_5$","$(\hat{a}_1\hat{a}_5)^2$","$(\hat{a}_1\hat{a}_5)^3$",...
    "$(\hat{a}_1\hat{a}_5)^5$","$(\hat{a}_1\hat{a}_5)^{25}$",Location="southeast")

if 0
    myfigure = gcf;
    figurename = 'f215ps.svg';
    saveas(myfigure, figurename)
    exportgraphics(myfigure, figurename, ContentType="vector")
end


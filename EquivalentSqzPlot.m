% 多径簇态施加PS后的等效sqz
clc;clear
% 1245 sub:
Nr = 1e3;
% 设置所有线条的默认粗细为2点
set(groot, 'DefaultLineLineWidth', 1.5);  
r = linspace(-2.5, -0.0000001, Nr);
sqz5 = zeros(5, Nr);
sqz7 = zeros(5, Nr);
sqz9 = sqz5;

sqzini = [exp(1).^(2.*r)];
sqz5(1, :) = [3.*(35+(-46).*exp(1).^(2.*r)+59.*exp(1).^(4.*r)).*((-26)+170.* ...
  cosh(2.*r)).^(-1)];
sqz5(2, :) = [(89+(-82).*exp(1).^(2.*r)+137.*exp(1).^(4.*r)).*((-26)+170.*cosh( ...
  2.*r)).^(-1)];
sqz5(3, :) = [(89+(-82).*exp(1).^(2.*r)+137.*exp(1).^(4.*r)).*((-26)+170.*cosh( ...
  2.*r)).^(-1)];
sqz5(4, :) = [(89+(-82).*exp(1).^(2.*r)+137.*exp(1).^(4.*r)).*((-26)+170.*cosh( ...
  2.*r)).^(-1)];
sqz5(5, :) = [3.*(35+(-46).*exp(1).^(2.*r)+59.*exp(1).^(4.*r)).*((-26)+170.* ...
  cosh(2.*r)).^(-1)];

sqz7(1, :) = [(253+(-306).*exp(1).^(2.*r)+453.*exp(1).^(4.*r)).*((-42)+442.* ...
  cosh(2.*r)).^(-1)];
sqz7(2, :) = [5.*(45+(-26).*exp(1).^(2.*r)+61.*exp(1).^(4.*r)).*((-42)+442.* ...
  cosh(2.*r)).^(-1)];

sqz9(1, :) = [(465+(-538).*exp(1).^(2.*r)+857.*exp(1).^(4.*r)).*((-58)+842.* ...
  cosh(2.*r)).^(-1)];

sqz9(2, :) = [(425+(-178).*exp(1).^(2.*r)+537.*exp(1).^(4.*r)).*((-58)+842.* ...
  cosh(2.*r)).^(-1)];

ls = ["-","--","-."];
mk = ["x","o","+"];
figure;hold on
% plot(-r*8.69, sqzini./sqzini, 'k--')
% #ca4362 : 红色
for i = 1:2
    plot(-r*8.69, sqz5(i,:)./sqzini, Color="#ca4362", LineStyle=ls(i), Marker=mk(1), MarkerIndices=1:100:1e3)
end

% #2c73d2 : 蓝色
for i = 1:2
    plot(-r*8.69, sqz7(i,:)./sqzini, Color="#2c73d2", LineStyle=ls(i), Marker=mk(2), MarkerIndices=1:100:1e3)
end

% #ffc572 : 黄色

for i = 1:2
    plot(-r*8.69, sqz9(i,:)./sqzini, Color="#ffc572", LineStyle=ls(i), Marker=mk(3), MarkerIndices=1:100:1e3)
end

% 构造填充区域的坐标（先绘制曲线，再连接到顶部边界，最后回到起点）
x_fill = [-0.01, -0.01, 16, 16];
y_fill = [-0.01, 1, 1, 0];

% 填充区域
% plot(x, y, 'LineWidth', 2);  % 绘制原始曲线
% hold on;
fs = 15;
fill(x_fill, y_fill, 'r', 'FaceAlpha', 0.1, EdgeColor='k',LineStyle='--');  % 填充曲线以上区域（蓝色，30%透明度）
xlim([0,15])
ylim([0.8,1.2])
box on
fs = 14;
legend("End modes $(K=5)$","Mid modes $(K=5)\quad$",...
    "End modes $(K=7)$","Mid modes $(K=7)$",...
    "End modes $(K=9)$","Mid modes $(K=9)$",...
    Location="northwest", FontSize = 11)

text(0.1,0.95,"Squeezing enhancement region",Position=[0.412446338525327,0.072528928808196,0],Units="normalized", FontSize = fs);

ylabel("$\frac{\Delta^2p}{\mathrm{e}^{2r}}$", Rotation = 0, FontSize = fs);
xlabel("$|r_\mathrm{dB}|$", "FontSize", fs)
pbaspect([1, 0.618, 1])

if 0
    myfigure = gcf;
    figurename = 'sqz_equivalent.svg';
    saveas(myfigure, figurename)
    exportgraphics(myfigure, figurename, ContentType="vector")
end
% 单模真空压缩态施加PS后的压缩度
clc;clear
% 设置所有线条的默认粗细为2点
set(groot, 'DefaultLineLineWidth', 1.5);  

r = linspace(-2, -0.000001, 1e3);
V0 = exp(2*r);
V1 = 3.*exp(1).^(2.*r);
V2 = (7+(-18).*exp(1).^(2.*r)+15.*exp(1).^(4.*r)).*((-2)+6.*cosh(2.*r)) ...
  .^(-1);
V3 = (11+(-10).*exp(1).^(2.*r)+35.*exp(1).^(4.*r)).*(2+10.*cosh(2.*r)) ...
  .^(-1);
V4 = (132+(-15).*exp(1).^((-2).*r).*(5+7.*exp(1).^(4.*r).*(2+(-4).*exp( ...
  1).^(2.*r)+3.*exp(1).^(4.*r)))).*((-18)+40.*cosh(2.*r)+(-70).* ...
  cosh(4.*r)).^(-1);
V5 = ((-92)+7.*exp(1).^((-2).*r).*(19+34.*exp(1).^(4.*r)+(-36).*exp(1) ...
  .^(6.*r)+99.*exp(1).^(8.*r))).*(58+56.*cosh(2.*r)+126.*cosh(4.*r)) ...
  .^(-1);
V6 = (1125+7.*exp(1).^((-4).*r).*(69+(-114).*exp(1).^(2.*r)+(-220).* ...
  exp(1).^(6.*r)+315.*exp(1).^(8.*r)+(-594).*exp(1).^(10.*r)+429.* ...
  exp(1).^(12.*r))).*((-100)+210.*cosh(2.*r)+(-252).*cosh(4.*r)+ ...
  462.*cosh(6.*r)).^(-1);
V7 = 3.*(479+3.*exp(1).^((-4).*r).*(99+(-62).*exp(1).^(2.*r)+11.*exp(1) ...
  .^(6.*r).*((-12)+23.*exp(1).^(2.*r)+(-26).*exp(1).^(4.*r)+65.*exp( ...
  1).^(6.*r)))).*(212+774.*cosh(2.*r)+396.*cosh(4.*r)+858.*cosh(6.* ...
  r)).^(-1);
V8 = exp(1).^(2.*r).*((-1225)+2520.*cosh(2.*r)+(-2772).*cosh(4.*r)+ ...
  3432.*cosh(6.*r)+(-6435).*cosh(8.*r)).^(-1).*(49112.*cosh(2.*r)+( ...
  -56532).*cosh(4.*r)+87912.*cosh(6.*r)+(-61347).*cosh(8.*r)+(-7).*( ...
  3375+(-1696).*sinh(2.*r)+3936.*sinh(4.*r)+(-9504).*sinh(6.*r)+ ...
  6864.*sinh(8.*r)));


%%
fs = 14;
r = -r*8.69;
% plot(r,V0,'k--');

hold on;

plot(r,V2,Color="#004497");
plot(r,V4,Color="#4165bd");
plot(r,V6,Color="#6b87e4");
plot(r,V8,Color="#94acff");
plot(r,V1,'-.',Color="#970000");
plot(r,V3,'-.',Color="#d44649");
plot(r,V5,'-.',Color="#ff8090");
plot(r,V7,'-.',Color="#ffbcda");
box on


% legend("$\hat{a}_0$","$\hat{a}_1$","$\hat{a}_2$","$\hat{a}_3$","$\hat{a}_4$",...
    % "$\hat{a}_5$","$\hat{a}_6$","$\hat{a}_7$","$\hat{a}_8$")
xlabel("$|r_\mathrm{dB}|$", FontSize = fs)
pbaspect([1, 0.618, 1])
ylabel("$\Delta^2p$",Rotation=0, FontSize = fs)
xlim([0,16])
ylim([0,3])


% 构造填充区域的坐标（先绘制曲线，再连接到顶部边界，最后回到起点）
x_fill = [fliplr(abs(r))-0.02, abs(r)-0.03];
y_fill = [zeros(size(r))-0.02, V0];

% 填充区域
fill(x_fill, y_fill, 'r', 'FaceAlpha', 0.1, EdgeColor='k',LineStyle='--');  % 填充曲线以上区域（蓝色，30%透明度）
text(0.1,0.95,"Squeezing enhancement region",...
    Position=[-0.001172508108089,0.035122445017672,0],Units="normalized", FontSize = 12);

text(0.1,0.95,"1", FontSize = 12, Rotation=45, Position=[3.477584158415841,1.362128712871287,0]);
text(0.1,0.95,"3", FontSize = 12, Rotation=45, Position=[2.706613861386137,1.128217821782178,0]);
text(0.1,0.95,"5", FontSize = 12, Rotation=45, Position=[2.302772277227722,0.994554455445545,0]);
text(0.1,0.95,"7", FontSize = 12, Rotation=45, Position=[1.935643564356434,0.883168316831683,0]);

% 在主图右上角创建嵌入图
axes('Position', [0.427941176470588,0.471544715447154,0.427941176470588,0.43115058324496]);  % [左, 下, 宽, 高]，归一化坐标
hold on
plot(r,V2,Color="#004497");
plot(r,V4,Color="#4165bd");
plot(r,V6,Color="#6b87e4");
plot(r,V8,Color="#94acff");
plot(r,V1,'-.',Color="#970000");
plot(r,V3,'-.',Color="#d44649");
plot(r,V5,'-.',Color="#ff8090");
plot(r,V7,'-.',Color="#ffbcda");
xlim([0,3])
ylim([0.3,1])
fill(x_fill, y_fill, 'r', 'FaceAlpha', 0.1, EdgeColor='k',LineStyle='--');  % 填充曲线以上区域（蓝色，30%透明度）

text(0.1,0.95,"2", FontSize = 12, Rotation=45, Position=[0.748,0.6078125,0]);
text(0.1,0.95,"4", FontSize = 12, Rotation=45, Position=[0.532,0.5421875,0]);
text(0.1,0.95,"6", FontSize = 12, Rotation=45, Position=[0.412,0.5046875,0]);
text(0.1,0.95,"8", FontSize = 12, Rotation=45, Position=[0.316,0.4671875,0]);

box on


if 0
    myfigure = gcf;
    figurename = 'sqz_smode.pdf';
    saveas(myfigure, figurename)
    exportgraphics(myfigure, figurename, ContentType="vector")
end
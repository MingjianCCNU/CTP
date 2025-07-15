clc;clear

r = linspace(-1.25,0);
F14 = (0.112666E0+0.339043E-1.*exp(1).^((-2).*r)+0.339043E-1.*exp(1).^( ...
  2.*r)).^(-1/2).*pi.^(-1);

F23 = (0.106905E0+0.237846E-1.*exp(1).^((-2).*r)+0.237846E-1.*exp(1).^( ...
  2.*r)).^(-1/2).*pi.^(-1);

F24 = (0.105878E0+0.960194E-2.*exp(1).^((-2).*r)+0.48087E-1.*exp(1).^( ...
  2.*r)).^(-1/2).*pi.^(-1);


lw = 2;
fs = 14;
plot(abs(r),F14, Color = "#845ec2",LineWidth=lw)
hold on
plot(abs(r),F23,"--", Color = "#4b4453",LineWidth=lw)
plot(abs(r),F24,"-.", Color = "#b0a8b9",LineWidth=lw)
scatter(abs(r(20)), F24(20), 100, "pentagram", MarkerEdgeColor="#b0a8b9", MarkerFaceColor="#b0a8b9")
scatter(abs(r(20)), F23(20), 100, "pentagram", MarkerEdgeColor="#4b4453", MarkerFaceColor="#4b4453")
ylabel("${F}$", Rotation = 0, FontSize = fs);
xlabel("$|r_\mathrm{in}|$", "FontSize", fs)
legend("$r_1=r_4=r_\uparrow,\quad r_2=r_3=r_\downarrow\quad$",...
    "$r_2=r_3=r_\uparrow,\quad r_1=r_4=r_\downarrow$",...
    "$r_2=r_4=r_\uparrow,\quad r_1=r_3=r_\downarrow$",Location="southwest", FontSize = fs)
ylim([0.35,1])
xlim([0,1.25])
pbaspect([1, 0.618, 1])

%%
if 1
    myfigure = gcf;
    figurename = 'fsqz';
    saveas(myfigure, figurename)
    print(figurename,'-dsvg','-r300', '-vector',myfigure)
end
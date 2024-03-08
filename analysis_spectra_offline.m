T = readtable('shin_2013_data/shin2013_vertical.csv','HeaderLines',6);
outputpath = 'shin_2013_data/';
% 0.1:
titlepost = ' at z = 10% PBL';
fileout = 'spectra0.1.png';

kx = table2array(T(:,'x0_1hX_1'));
Fx = table2array(T(:,'x0_1hY_1'));

kz = table2array(T(:,'x0_1vX_1'));
Fz = table2array(T(:,'x0_1vY_1'));

% 
% % % 0.5:
% titlepost = ' at z = 50% PBL';
% fileout = 'spectra0.5.png';
% 
% kx = table2array(T(:,'x0_5hX_1'));
% Fx = table2array(T(:,'x0_5hY_1'));
% 
% kz = table2array(T(:,'x0_5vX_1'));
% Fz = table2array(T(:,'x0_5vY_1'));
% % 
% % % 0.8:
titlepost = ' at z = 80% PBL';
fileout = 'spectra0.8.png';

kx = table2array(T(:,'x0_8hX_1'));
Fx = table2array(T(:,'x0_8hY_1'));

kz = table2array(T(:,'x0_8vX_1'));
Fz = table2array(T(:,'x0_8vY_1'));

kx(isnan(kx)) = [];
Fx(isnan(Fx)) = [];
kz(isnan(kz)) = [];
Fz(isnan(Fz)) = [];
FFT_factor = Lx^2*Ly^2/(2*pi)^2/(nx*ny);
kx = kx/1000;
kz = kz/1000;
Fx = Fx*FFT_factor;
Fz = Fz*FFT_factor;
%%
int_scale_x = trapz(kx,(1./kx).*Fx)./trapz(kx,Fx);
int_scale_z = trapz(kz,(1./kz).*Fz)./trapz(kz,Fz);


[~,scalex] = min(abs(kx - 2*pi/min([int_scale_x,int_scale_z])));
[~,scalez] = min(abs(kz - 2*pi/min([int_scale_x,int_scale_z])));

anisotropy = trapz(kx(scalex:end),Fx(scalex:end))/trapz(kz(scalez:end),Fz(scalez:end));
%%
%close all;
h=figure('visible','on');
set(h,'Resize','on');
set(h,'PaperPosition',[0 0 16 9]);
set(h,'PaperUnits','inches');
set(h,'PaperSize',[16 9]); % IEEE columnwidth = 9cm
set(h,'Units','inches')
set(h,'Position',[0 0 16 9]);

ax = axes();

set(ax,'fontsize',13)

hold(ax,'on');
set(ax,'yscale','log');
set(ax,'xscale','log');
set(ax,'ylim',YL);%[10^(4),10^10]);
set(ax,'xlim',XL);%[min([kx;kz]),max([kx;kz])]);
xticks = [10^-4,5*10^-4,10^-3,5*10^-3,10^-2,5*10^-2,];
xticklabels = {'10^{-4}','5\times 10^{-4}','10^{-3}','5\times 10^{-3}','10^{-2}','5\times 10^{-2}'};
set(ax,'xtick',xticks)
set(ax,'xtickLabel',cellfun(@(c) [num2str(c,'%.0e'),'/m'],num2cell(get(ax,'xtick')),'uniformoutput',false))
grid(ax,'on');

ylabel(ax,'Average velocity spectrum (m^5/s^2)')
xlabel(ax,'radial-wavenumber (1/m)')
%text(ax,0,0,localname,'Position',[-0.1,1.05],'units','normalized','fontsize',16,'interpreter','none')
r_wl = 2*pi./kx;
ax_top = axes();
set(ax_top,'fontsize',13)

hold(ax_top,'on');
set(ax_top,'Position',get(ax,'Position'));
set(ax_top,'XAxisLocation','top');
set(ax_top,'Color','none');
%set(ax_top,'tickLength',5)
set(ax_top,'xdir','reverse')
set(ax_top,'xlim',[r_wl(end),r_wl(1)])
set(ax_top,'xtick',fliplr(2*pi./get(ax,'xtick')))
set(ax_top,'xscale','log')
set(ax_top,'yscale','log')
set(ax_top,'xtickLabel',cellfun(@(c) [num2str(c,2),'m'],num2cell(get(ax_top,'xtick')),'uniformoutput',false))
set(ax_top,'ylim',get(ax,'ylim'));
set(ax,'gridLineStyle','-')
xl = xlabel(ax_top,'radial-wavelength (m)');
set(xl,'units','inches')
xl_pos = get(xl,'position');
xl_pos(2) = [xl_pos(2)-0.6];
set(xl,'position',xl_pos);
set(xl,'backgroundcolor','w')
%set(ax_top,'grid','off')
plot(ax,kx,Fx,'r','linewidth',2);
%plot(ax,kz,Fz,'g','linewidth',2);
plot(ax,kz,Fz,'b','linewidth',2);
%plot(ax,r_wn,fr_co_norm,'m','linewidth',2);

%plot(ax,r_wn,fr_r,'m','linewidth',2);
%plot(ax,r_wn,ft_r,'c','linewidth',2);
lin = kx.^(-5/3);
ylim = get(ax,'ylim');

t = title(ax,['Radial average of velocity spectrum, from Shin and Hong 2013.',' Anisotropy= ',num2str(anisotropy,2),titlepost],'fontsize',16);
set(t,'units','inches')

t_pos = get(t,'position');
t_pos(2) = [t_pos(2)+0.25];
set(t,'position',t_pos);

lin = (10.^(-10:20))'*kx.^(-5/3)'/max(kx.^(-5/3));
plot(ax,kx,lin,'k--','linewidth',2);
set(ax,'ylim',ylim);

legend(ax,{['u^2 = ',num2str(trapz(kx,Fx),2),' m^4/s^2'],...
           %['v^2 = ',num2str(trapz(r_wn,fv_r),2),' m^4/s^2'],...
           ['w^2 = ',num2str(trapz(kz,Fz),2),' m^4/s^2'],...
           %['w''\theta'' = ',num2str(covar,2),' m^3K/s'],...
           '\kappa\^(-5/3)'},'location','southwest')
hgexport(h, fullfile(outputpath,fileout), hgexport('factorystyle'), 'Format', 'png');

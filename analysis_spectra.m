clear all
k = 37;
%k = 15;

datapath = 'blankpert/default/';
fileout2 = ['noturb_',num2str(k),'_2dspectrum.png'];
fileout3 = ['noturb_',num2str(k),'_3dspectrum.png'];
outputpath = datapath;
prefix = 'cm1out_noturb_';
% if ~exist('datapath','var')
%     datapath = 'example_data';
%     fileout = 'example_2dspectrum.png';
%     prefix = 'cm1out_';
%     outputpath = './';
% elseif exist('postfix','var');
%     [~,dp_header,~] = fileparts(datapath);
%     fileout = [dp_header,'_spectrum',postfix,'.png'];
% else
%     [~,dp_header,~] = fileparts(datapath);
%     fileout = [dp_header,'_spectrum','.png'];
% end

if ~exist('zfactor','var')
    zfactor = .5;
end

if ~exist('titlepost','var')
    titlepost = '';
end

e2d = 0;
e3d = 0;

time = zeros(1,max(k));
fu_r = 0;
fv_r = 0;
fw_r = 0;
fr_r = 0;
ft_r = 0;

fu_t = 0;
fv_t = 0;
fw_t = 0;
fr_t = 0;
ft_t = 0;

fr_co = 0;
num_files = 0;


for it=1:length(k);
    filename = fullfile(datapath,[prefix,num2str(k(it),'%06d'),'.nc']);
    if ~exist(filename,'file')
        continue
    end
    num_files = num_files + 1;

    time(k(it)) = double(ncread(filename,'time'));
    

    u = double(ncread(filename,'uinterp'));
    v = double(ncread(filename,'vinterp'));
    w = double(ncread(filename,'winterp'));
    th= double(ncread(filename,'th')); % potential temperature

    % boundary layer height:
    %hpbl = double(ncread(filename,'hpbl'));

    % grid size parameters:
    xh = double(ncread(filename,'xh'))*1000;
    dx = median(diff(xh));
    xh_mid = xh - max(xh)/2;
    yh = double(ncread(filename,'yh'))*1000;
    dy = median(diff(xh));
    yh_mid = yh - max(yh)/2;
    z = double(ncread(filename,'z'))*1000;  % for nodes
    zf = double(ncread(filename,'zf'))*1000; % for edges
    dz = median(diff(z));


    nx = numel(xh);
    ny = numel(yh);
    nz = numel(z);

    Lx = nx*dx;
    Ly = ny*dy;
    Lz = nz*dz;

    kx = 2*pi*(-nx/2:(nx/2-1))/Lx;
    ky = 2*pi*(-ny/2:(ny/2-1))/Ly;
    kz = 2*pi*(-nz/2:(nz/2-1))/Lz;

    % take FFT around pbl/2:
    iz = floor(nz/2);
    %[~,iz] = min(abs(z - mean(hpbl(:))*zfactor));
    Ninteg = 0; % radius of z-integration points around pbl/2
    integ = (iz-Ninteg):(iz+Ninteg);
    Ncirc = 250;

    [Xh_mid,Yh_mid] = meshgrid(xh_mid,yh_mid);
    R_rs = sqrt(Xh_mid.^2+Yh_mid.^2); % radius realspace
    Theta_rs = atan2(Yh_mid,Xh_mid);

    r_wn = linspace(0,min([max(kx),max(ky)]),Ncirc); % radius wavenumber
    rho_wn = linspace(0,min([max(kx),max(ky),max(kz)]),Ncirc); % radius wavenumber

    theta_wn = linspace(0,2*pi,Ncirc);
    theta_wn = theta_wn(1:(end-1));

    phi_wn = linspace(0,pi,Ncirc);
    % Factor so that FFT scales like /density/ kinetic energy
    FFT_factor_2d = FFT_factor('density',[nx,ny],[Lx,Ly]);
    FFT_factor_3d = FFT_factor('density',[nx,ny,nz],[Lx,Ly,Lz]);

    % Total kinetic energy (average over time):

    for ik = integ
        % kinetic energy:        
        for idir=1:3
            switch idir
                case 1
                    vel = squeeze(u(:,:,:,:));
                case 2
                    vel = squeeze(v(:,:,:,:));
                case 3
                    vel = squeeze(w(:,:,:,:));
                case 4
                    vel = (1./(R_rs)).*(Xh_mid.*squeeze(u(:,:,ik,:)) + Yh_mid.*squeeze(v(:,:,ik,:)));
                case 5
                    vel = (1./(R_rs)).*(Xh_mid.*squeeze(v(:,:,ik,:)) - Yh_mid.*squeeze(u(:,:,ik,:)));
            end
            %e2d = e2d + dx.*dy.*sum(sum(vel(:,:,ik).^2))/(Lx*Ly);
            %e3d = e3d + dx.*dy.*dz.*sum(sum(sum(vel.^2)))/(Lx*Ly*Lz);
            e3d = e3d + sum(sum(sum(vel.^2)))/(Lx*Ly*Lz);
            
            f2 = FFT_factor_2d*abs(fftshift(fft2(vel(:,:,ik)))).^2;
            f2_interp = griddedInterpolant({kx,ky},f2);
            f_2d_ik = zeros(size(r_wn));
            
            % 2d spectra:
            for ri = 1:length(r_wn)
                f_2d_ik(ri) = trapz(theta_wn,abs(r_wn(ri)).*f2_interp(r_wn(ri)*cos(theta_wn'),r_wn(ri)*sin(theta_wn')));
            end
            
            f3 = FFT_factor_3d*abs(fftshift(fftn(vel))).^2;
            f3_interp = griddedInterpolant({kx,ky,kz},f3);                        
            f_3d_ik = zeros(size(rho_wn));
            
            % 3d spectra:
            [T,P] = meshgrid(theta_wn,phi_wn);
            for ri = 1:length(rho_wn)
                r = rho_wn(ri);
                f_3d_ik(ri) = trapz(theta_wn,trapz(phi_wn,r.^2.*sin(P).*f3_interp(r.*cos(T).*sin(P),r.*sin(T).*sin(P),r.*cos(P)),1),2);
            end
%             for ti = 1:length(theta_wn)
%                 f_3d_ik(ti) = trapz(r_wn,f2_interp(r_wn'*cos(theta_wn(ti)),r_wn'*sin(theta_wn(ti))));
%             end

%             if idir == 3
%                 % calculate vertical cospectra of potential temperature
%                 % (and, when available, water vapour)
%                 % Base state:  q = qz(z,t) + pq
%                 % Perturbation potential temperature:
% %                    pw = w - repmat(wz,size(w,1),size(w,2),1);
%                 flxvar = squeeze(th(:,:,ik,:));
%                 flxz = sum(sum(flxvar,1),2).*dx.*dy / (Lx.*Ly);
%                 pertflx = flxvar - repmat(flxz,size(flxvar,1),size(flxvar,2),1);
% 
%                 fr_co_ik = zeros(size(r_wn));
%                 f_co_ik = FFT_factor_2d*abs(fftshift(fft2(vel)).*conj(fftshift(fft2(pertflx))));
%                 f_co_ik_interp = griddedInterpolant({kx,ky},f_co_ik);
%                 for ri = 1:length(r_wn);
%                     fr_co_ik(ri) = trapz(theta_wn,abs(r_wn(ri)).*f_co_ik_interp(r_wn(ri)*cos(theta_wn'),r_wn(ri)*sin(theta_wn')));
%                 end
%                 fr_co = fr_co + fr_co_ik;
%             end


            switch idir
                case 1
                    fu_r = fu_r+f_2d_ik;
                    fu_t = fu_t+f_3d_ik;
                case 2
                    fv_r = fv_r+f_2d_ik;
                    fv_t = fv_t+f_3d_ik;
                case 3
                    fw_r = fw_r+f_2d_ik;
                    fw_t = fw_t+f_3d_ik;
                case 4
                    fr_r = fr_r + f_2d_ik;
                    fr_t = fr_t + f_3d_ik;
                case 5
                    ft_r = ft_r + f_2d_ik;
                    ft_t = ft_t + f_3d_ik;
            end
        end
    end


end
e2d = e2d/num_files/length(integ);
e3d = e3d/num_files;


% after all sums are done, insure we normalize 
% via time integration length, and z-integration length

% r-spectra:
fu_r = abs(fu_r).^1/length(integ)/num_files;
fv_r = abs(fv_r).^1/length(integ)/num_files;
fw_r = abs(fw_r).^1/length(integ)/num_files;
fr_r = fr_r./length(integ)/num_files;
ft_r = ft_r./length(integ)/num_files;

% theta-spectra:
fu_t = abs(fu_t).^1/num_files;
fv_t = abs(fv_t).^1/num_files;
fw_t = abs(fw_t).^1/num_files;
fr_t = fr_t./num_files;
ft_t = ft_t./num_files;

e2d_f = trapz(r_wn,(fu_r+fv_r+fw_r));
e3d_f = trapz(rho_wn,(fu_t+fv_t+fw_t));
%%
int_scale_x = trapz(r_wn(2:end),(1./r_wn(2:end)).*fu_r(2:end))./trapz(r_wn(2:end),fu_r(2:end));
int_scale_z = trapz(r_wn(2:end),(1./r_wn(2:end)).*fw_r(2:end))./trapz(r_wn(2:end),fw_r(2:end));
[~,scale] = min(abs(r_wn - 2*pi/max([int_scale_x,int_scale_z])));
anisotropy = trapz(r_wn(scale:end),fu_r(scale:end))/trapz(r_wn(scale:end),fw_r(scale:end));
%%
% cospectra:
% fr_co = fr_co ./ length(integ) / num_files;
% 
% covar = trapz(r_wn,fr_co);
% fr_co_norm = fr_co *trapz(r_wn,fw_r)/covar;
% %fr_co_norm = fr_co;

%max(kx)
%%
close all;
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
set(ax,'ylim',[10^(-5),10^1]);
set(ax,'xlim',[r_wn(2),r_wn(end)]);
xticks = [10^-4,5*10^-4,10^-3,5*10^-3,10^-2,5*10^-2,];
xticklabels = {'10^{-4}','5\times 10^{-4}','10^{-3}','5\times 10^{-3}','10^{-2}','5\times 10^{-2}'};
set(ax,'xtick',xticks)
set(ax,'xtickLabel',cellfun(@(c) [num2str(c,'%.0e'),'/m'],num2cell(get(ax,'xtick')),'uniformoutput',false))
grid(ax,'on');

ylabel(ax,'Average velocity spectrum (m^2/s^2)')
xlabel(ax,'radial-wavenumber (1/m)')
%text(ax,0,0,localname,'Position',[-0.1,1.05],'units','normalized','fontsize',16,'interpreter','none')
r_wl = 2*pi./r_wn;
ax_top = axes();
set(ax_top,'fontsize',13)

hold(ax_top,'on');
set(ax_top,'Position',get(ax,'Position'));
set(ax_top,'XAxisLocation','top');
set(ax_top,'Color','none');
%set(ax_top,'tickLength',5)
set(ax_top,'xdir','reverse')
set(ax_top,'xlim',[r_wl(end),r_wl(2)])
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
%plot(ax,r_wn,fu_r,'r','linewidth',2);
%plot(ax,r_wn,fv_r,'g','linewidth',2);
%plot(ax,r_wn,fw_r,'b','linewidth',2);
%plot(ax,r_wn,fr_co_norm,'m','linewidth',2);
plot(ax,r_wn,fu_r,'r','linewidth',2);
plot(ax,r_wn,fv_r,'g','linewidth',2);
plot(ax,r_wn,fw_r,'b','linewidth',2);
%plot(ax,r_wn,fr_co_norm,'m','linewidth',2);
%plot(ax,r_wn,fr_r,'m','linewidth',2);
%plot(ax,r_wn,ft_r,'c','linewidth',2);
lin = r_wn.^(-5/3);

t = title(ax,['3d velocity spectrum, at t=',num2str(time(end)/60/60),'hr, z = 1.6km.',' Anisotropy= ',num2str(anisotropy,2),titlepost],'fontsize',16);
set(t,'units','inches')

t_pos = get(t,'position');
t_pos(2) = [t_pos(2)+0.25];
set(t,'position',t_pos);
YL = get(ax,'ylim');
XL = get(ax,'xlim');

lin = (10.^(-10:14))'*lin/max(lin(2:end));
plot(ax,r_wn,lin,'k--','linewidth',2);
set(ax,'ylim',YL);
legend(ax,{['u^2 = ',num2str(trapz(r_wn,fu_r),2),' m/s^2'],...
           ['v^2 = ',num2str(trapz(r_wn,fv_r),2),' m/s^2'],...
           ['w^2 = ',num2str(trapz(r_wn,fw_r),2),' m/s^2'],...
           '\kappa\^(-5/3)'},'location','southwest')%           ['w''\theta'' = ',num2str(covar,2),' m^3K/s'],...
hgexport(h, fullfile(outputpath,fileout2), hgexport('factorystyle'), 'Format', 'png');
%%
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
set(ax,'ylim',[10^(-5),10^1]);
set(ax,'xlim',[rho_wn(2),rho_wn(end)]);
xticks = [10^-4,5*10^-4,10^-3,5*10^-3,10^-2,5*10^-2,];
xticklabels = {'10^{-4}','5\times 10^{-4}','10^{-3}','5\times 10^{-3}','10^{-2}','5\times 10^{-2}'};
set(ax,'xtick',xticks)
set(ax,'xtickLabel',cellfun(@(c) [num2str(c,'%.0e'),'/m'],num2cell(get(ax,'xtick')),'uniformoutput',false))
grid(ax,'on');

ylabel(ax,'Average velocity spectrum (m^3/s^2)')
xlabel(ax,'radial-wavenumber (1/m)')
%text(ax,0,0,localname,'Position',[-0.1,1.05],'units','normalized','fontsize',16,'interpreter','none')
r_wl = 2*pi./rho_wn;
ax_top = axes();
set(ax_top,'fontsize',13)

hold(ax_top,'on');
set(ax_top,'Position',get(ax,'Position'));
set(ax_top,'XAxisLocation','top');
set(ax_top,'Color','none');
%set(ax_top,'tickLength',5)
set(ax_top,'xdir','reverse')
set(ax_top,'xlim',[r_wl(end),r_wl(2)])
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
%plot(ax,r_wn,fu_r,'r','linewidth',2);
%plot(ax,r_wn,fv_r,'g','linewidth',2);
%plot(ax,r_wn,fw_r,'b','linewidth',2);
%plot(ax,r_wn,fr_co_norm,'m','linewidth',2);
plot(ax,rho_wn,fu_t,'r','linewidth',2);
plot(ax,rho_wn,fv_t,'g','linewidth',2);
plot(ax,rho_wn,fw_t,'b','linewidth',2);
%plot(ax,r_wn,fr_co_norm,'m','linewidth',2);
%plot(ax,r_wn,fr_r,'m','linewidth',2);
%plot(ax,r_wn,ft_r,'c','linewidth',2);
lin = rho_wn.^(-5/3);

t = title(ax,['3d velocity spectrum, at t=',num2str(time(end)/60/60),'hr,',titlepost],'fontsize',16);
set(t,'units','inches')

t_pos = get(t,'position');
t_pos(2) = [t_pos(2)+0.25];
set(t,'position',t_pos);
YL = get(ax,'ylim');
XL = get(ax,'xlim');

lin = (10.^(-10:14))'*lin/max(lin(2:end));
plot(ax,rho_wn,lin,'k--','linewidth',2);
set(ax,'ylim',YL);
legend(ax,{['u^2 = ',num2str(trapz(rho_wn,fu_t),2),' m^2/s^2'],...
           ['v^2 = ',num2str(trapz(rho_wn,fv_t),2),' m^2/s^2'],...
           ['w^2 = ',num2str(trapz(rho_wn,fw_t),2),' m^2/s^2'],...
           '\kappa\^(-5/3)'},'location','southwest')%           ['w''\theta'' = ',num2str(covar,2),' m^3K/s'],...
hgexport(h, fullfile(outputpath,fileout3), hgexport('factorystyle'), 'Format', 'png');

% figure(3);clf(3); hold on;
% 
% plot(theta,fu_t,'r','linewidth',2);
% plot(theta,fv_t,'g','linewidth',2);
% plot(theta,fw_t,'b','linewidth',2);
% %lin = theta.^(-5/3);
% %lin = lin*1e7/max(lin(2:end));
% %plot(theta,lin,'k--','linewidth',2);
% 
% %set(gca,'xscale','log')
% %set(gca,'yscale','log')
% ylabel('Average velocity spectrum (m/s^2)')
% xlabel('theta-wavenumber (1/m)')
% %xlim([theta(1),theta(end)])
% %ylim([1,10^6])
% grid on;
% title(['Average 1-D velocity spectrum in mid-pbl'])
% legend('u','v','w')
% 
% figure(4);clf(4);hold on;
% plot(fu_t.*cos(theta_wn),fu_t.*sin(theta_wn),'r','linewidth',2);
% plot(fv_t.*cos(theta_wn),fv_t.*sin(theta_wn),'g','linewidth',2);
% plot(fw_t.*cos(theta_wn),fw_t.*sin(theta_wn),'b','linewidth',2);
% plot(fr_t.*cos(theta_wn),fr_t.*sin(theta_wn),'m','linewidth',2);
% plot(ft_t.*cos(theta_wn),ft_t.*sin(theta_wn),'c','linewidth',2);
% 
% xl = get(gca,'XLim');
% yl = get(gca,'YLim');
% axis([min([xl,yl]),max([xl,yl]),min([xl,yl]),max([xl,yl])]);
% tick = mean(diff(get(gca,'XTick')));
% for rplot=tick:tick:max([xl,yl]);
%     plot(rplot.*cos(theta_wn),rplot.*sin(theta_wn),':','linewidth',1,'color',[0.5,0.5,0.5]);
% end
% for thetaplot=0:30:360
%     rplot = [0,2*max([xl,yl])];
%     plot(rplot.*cosd(thetaplot),rplot.*sind(thetaplot),':','linewidth',1,'color',[0.5,0.5,0.5]);
% end
% legend('u','v','w','v_r','v_\theta')
% title('Axial average of velocity spectrum in mid-pbl','fontsize',13)
% ylabel('Average velocity spectrum (m/s^2)')
% xlabel('Average velocity spectrum (m/s^2)')

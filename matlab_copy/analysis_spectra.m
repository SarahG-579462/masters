k = k_spectra;
%k = 15;

if ~exist('dp_header','var')
       [~,dp_header,~] = fileparts(datapath);
end       
if ~exist('datapath','var')
    datapath = 'example_data';
    fileout = 'example_data_spectrum.png';
    outputpath = './';
    fileout2 = ['example_2dspectrum','.png'];
    fileout3 = ['example_3dspectrum','.png'];
    fileoutq = 'example_qspectrum.png'
    fileoutb = ['example_budget','.png'];
elseif exist('postfix','var');
    fileout = [dp_header,'_spectrum',postfix,'.png'];
    fileout2 = [dp_header,'_2dspectrum',postfix,'.png'];
    fileout3 = [dp_header,'_3dspectrum',postfix,'.png'];
    fileoutb = [dp_header,'_budget',postfix,'.png'];
    fileoutq = [dp_header,'_qspectrum',postfix,'.png'];
else
    fileout2 = [dp_header,'_2dspectrum','.png'];
    fileout3 = [dp_header,'_3dspectrum','.png'];
    fileoutb = [dp_header,'_budget','.png'];
    fileoutq = [dp_header,'_qspectrum','.png'];
end

if ~exist('zfactor','var')
    zfactor = .5;
end
if ~exist('z_spectra','var')
	z_spectra = 0;
end
if ~exist('titlepost','var')
    titlepost = '';
end

if ~exist('budget_spectra','var')
	budget_spectra = false;
end
if ~exist('compute_3d','var')
	compute_3d = false;
end
e2d = 0;
e3d = 0;

plotcloud = false;
time = zeros(1,max(k));
fu_2 = 0;
fv_2 = 0;
fw_2 = 0;
fq_2 = 0;
ft_r = 0;

fu_3 = 0;
fv_3 = 0;
fw_3 = 0;
fq_3 = 0;
ft_t = 0;

fr_co = 0;
num_files = 0;


for it=1:length(k);
  filename = fullfile(datapath,[dataprefix,num2str(k(it),'%06d'),'.nc']);
  if ~exist(filename,'file')
    continue
  end
  num_files = num_files + 1;

  time(k(it)) = double(ncread(filename,'time'));
  u = double(ncread(filename,'uinterp'));
  v = double(ncread(filename,'vinterp'));
  w = double(ncread(filename,'winterp'));
  th= double(ncread(filename,'th')); % potential temperature
  nc = ncinfo(filename);
  if ismember('qc',{nc.Variables.Name})
    plotcloud = true;
    q = double(ncread(filename,'qv'));
  else
    q = zeros(size(u));
  end
  % boundary layer height:

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

  dS = dx*dy/(Lx*Ly);

  uz = sum(sum(u,1),2)*dS;
  vz = sum(sum(v,1),2)*dS;
  wz = sum(sum(w,1),2)*dS;
  wq = sum(sum(q,1),2)*dS;

  pu = u - repmat(uz,size(u,1),size(u,2),1);
  pv = v - repmat(vz,size(v,1),size(v,2),1);
  pw = w - repmat(wz,size(w,1),size(w,2),1);
  pq = q - repmat(wq,size(q,1),size(q,2),1);

  kx = pi*(-nx/2:(nx/2-1))/Lx;
  ky = pi*(-ny/2:(ny/2-1))/Ly;
  kz = pi*(-nz/2:(nz/2-1))/Lz;

  % take FFT around pbl/2:
  if z_spectra == 0
    hpbl = double(ncread(filename,'hpbl'));
    [~,iz] = min(abs(z - mean(hpbl(:))*zfactor));
  else
    [~,iz] = min(abs(z - z_spectra));
  end
  Ninteg = 0; % radius of z-integration points around pbl/2
  integ = (iz-Ninteg):(iz+Ninteg);
  Ncirc = 5*numel(kx);

  [Xh_mid,Yh_mid] = meshgrid(xh_mid,yh_mid);
  R_rs = sqrt(Xh_mid.^2+Yh_mid.^2); % radius realspace
  Theta_rs = atan2(Yh_mid,Xh_mid);

  r_wn = linspace(0,min([max(kx),max(ky)]),Ncirc); % radius wavenumber
  rho_wn = linspace(0,min([max(kx),max(ky),max(kz)]),Ncirc); % radius wavenumber

  if (it == 1) && budget_spectra
    fb_r_ik = zeros(length(r_wn),3,8);
  end
  theta_wn = linspace(0,2*pi,Ncirc);
  theta_wn = theta_wn(1:(end-1));

  phi_wn = linspace(0,pi,Ncirc);

  % Factor so that FFT scales like density of kinetic energy
  FFT_factor_2d = Lx*Ly/(2*pi)^2/(nx^2*ny^2);
  FFT_factor_3d = Lx*Ly*Lz/(2*pi)^3/(nx*ny*nz)^2;
  % Total kinetic energy (average over time):

  
  for ik = integ
  % kinetic energy:
  idirmax = (plotcloud)*1 + 3;
  for idir=1:idirmax
    switch idir
      case 1
      vel = squeeze(pu(:,:,:,:));
      if budget_spectra
        sizebudg = [size(vel),8];
        sizebudg(1) = sizebudg(1)+1;
        budget = zeros(sizebudg);
        budget(:,:,:,1) = squeeze(double(ncread(filename,'ub_hadv'))); 
        budget(:,:,:,2) = squeeze(double(ncread(filename,'ub_vadv')));
        try
          budget(:,:,:,3) = squeeze(double(ncread(filename,'ub_hediff')));
        catch
        end
        budget(:,:,:,4) = squeeze(double(ncread(filename,'ub_hturb')));
        budget(:,:,:,5) = squeeze(double(ncread(filename,'ub_vturb')));
        budget(:,:,:,6) = squeeze(double(ncread(filename,'ub_pgrad')));
        budget(:,:,:,7) = squeeze(double(ncread(filename,'ub_rdamp')));

        budget = 0.5*(budget(1:(end-1),:,:,:) + budget(2:end,:,:,:));
      end
      case 2
      vel = squeeze(pv(:,:,:,:));
      if budget_spectra
        sizebudg = [size(vel),8];
        sizebudg(2) = sizebudg(2)+1;
        budget = zeros(sizebudg);
        budget(:,:,:,1) = squeeze(double(ncread(filename,'vb_hadv'))); 
        budget(:,:,:,2) = squeeze(double(ncread(filename,'vb_vadv')));
        try
          budget(:,:,:,3) = squeeze(double(ncread(filename,'vb_hediff')));
        catch
        end
        budget(:,:,:,4) = squeeze(double(ncread(filename,'vb_hturb')));
        budget(:,:,:,5) = squeeze(double(ncread(filename,'vb_vturb')));
        budget(:,:,:,6) = squeeze(double(ncread(filename,'vb_pgrad')));
        budget(:,:,:,7) = squeeze(double(ncread(filename,'vb_rdamp')));
        budget = 0.5*(budget(:,1:(end-1),:,:) + budget(:,2:end,:,:));
      end
      case 3
      vel = squeeze(pw(:,:,:,:));
      if budget_spectra
        sizebudg = [size(vel),8];
        sizebudg(3) = sizebudg(3)+1;
        budget = zeros(sizebudg);
        budget(:,:,:,1) = squeeze(double(ncread(filename,'wb_hadv'))); 
        budget(:,:,:,2) = squeeze(double(ncread(filename,'wb_vadv')));
        try
          budget(:,:,:,3) = squeeze(double(ncread(filename,'wb_hediff')));
        catch
        end
        budget(:,:,:,4) = squeeze(double(ncread(filename,'wb_hturb')));
        budget(:,:,:,5) = squeeze(double(ncread(filename,'wb_vturb')));
        budget(:,:,:,6) = squeeze(double(ncread(filename,'wb_pgrad')));
        budget(:,:,:,7) = squeeze(double(ncread(filename,'wb_rdamp')));
        budget(:,:,:,8) = squeeze(double(ncread(filename,'wb_buoy')));

        budget = 0.5*(budget(:,:,1:(end-1),:) + budget(:,:,2:end,:));
      end
      case 4
      vel = squeeze(pq(:,:,:,:));
      %vel = (1./(R_rs)).*(Xh_mid.*squeeze(u(:,:,:,:)) + Yh_mid.*squeeze(v(:,:,:,:)));
      case 5
      %vel = (1./(R_rs)).*(Xh_mid.*squeeze(v(:,:,:,:)) - Yh_mid.*squeeze(u(:,:,:,:)));
    end

    f = FFT_factor_2d*abs(fftshift(fftshift(fft2(vel(:,:,ik)),1),2)).^2;
    f_interp = griddedInterpolant({kx,ky},f);
    f_r_ik = zeros(size(r_wn));
  % f_t_ik = zeros(size(theta_wn));

    for ri = 1:length(r_wn)
      f_r_ik(ri) = trapz(theta_wn,abs(r_wn(ri)).*f_interp(r_wn(ri)*cos(theta_wn'),r_wn(ri)*sin(theta_wn')));
    end

    f_3d_ik = zeros(size(rho_wn));
    if compute_3d && (it == length(k))
      f3 = FFT_factor_3d*abs(fftshift(fftn(vel))).^2;
      f3_interp = griddedInterpolant({kx,ky,kz},f3);                        

      % 3d spectra:
      [T,P] = meshgrid(theta_wn,phi_wn);
      for ri = 1:length(rho_wn)
        r = rho_wn(ri);
        f_3d_ik(ri) = trapz(theta_wn,trapz(phi_wn,r.^2.*sin(P).*f3_interp(r.*cos(T).*sin(P),r.*sin(T).*sin(P),r.*cos(P)),1),2);
      end
    end
  %            for ti = 1:length(theta_wn)
  %                f_t_ik(ti) = trapz(r_wn,f_interp(r_wn'*cos(theta_wn(ti)),r_wn'*sin(theta_wn(ti))));
  %            end
    if budget_spectra && (idir <= 3)
      fb= squeeze(0.5*FFT_factor_2d*(fftshift(fftshift(fft2(budget(:,:,ik,:).*repmat(vel(:,:,ik),1,1,1,size(budget,4))),1),2)));
      for iB = 1:size(budget,4)
        fbi = real(squeeze(fb(:,:,iB)));
        fb_interp = griddedInterpolant({kx,ky},fbi);

        for ri = 1:length(r_wn)
          fb_r_ik(ri,idir,iB) = fb_r_ik(ri,idir,iB)+trapz(theta_wn,abs(r_wn(ri)).*fb_interp(r_wn(ri)*cos(theta_wn'),r_wn(ri)*sin(theta_wn')));
        end
      end
    end



    if idir == 3
    % calculate vertical cospectra of potential temperature
    % (and, when available, water vapour)
    % Base state:  q = qz(z,t) + pq
    % Perturbation potential temperature:
    %                    pw = w - repmat(wz,size(w,1),size(w,2),1);
    %                flxvar = squeeze(th(:,:,ik,:));
    %                flxz = sum(sum(flxvar,1),2).*dx.*dy / (Lx.*Ly);
    %                pertflx = flxvar - repmat(flxz,size(flxvar,1),size(flxvar,2),1);
    %
    %                fr_co_ik = zeros(size(r_wn));
    %                f_co_ik = FFT_factor*abs(fftshift(fft2(vel)).*conj(fftshift(fft2(pertflx))));
    %                f_co_ik_interp = griddedInterpolant({kx,ky},f_co_ik);
    %                for ri = 1:length(r_wn);
    %                    fr_co_ik(ri) = trapz(theta_wn,abs(r_wn(ri)).*f_co_ik_interp(r_wn(ri)*cos(theta_wn'),r_wn(ri)*sin(theta_wn')));
    %                end
    %                fr_co = fr_co + fr_co_ik;
    end


    switch idir
      case 1
      fu_2 = fu_2+f_r_ik;
      if compute_3d
        fu_3 = fu_3+f_3d_ik;
      end
      case 2
      fv_2 = fv_2+f_r_ik;
      if compute_3d
        fv_3 = fv_3+f_3d_ik;
      end              
      case 3
      fw_2 = fw_2+f_r_ik;
      if compute_3d
        fw_3 = fw_3+f_3d_ik;
      end              
      case 4
      fq_2 = fq_2 + f_r_ik;
      if compute_3d
        fq_3 = fq_3 + f_3d_ik;
      end  
      case 5
      ft_r = ft_r + f_r_ik;
      if compute_3d
        ft_t = ft_t + f_3d_ik;
      end  
  end
end
end


end
%Erh = Erh/num_files/length(integ);


% after all sums are done, insure we normalize 
% via time integration length, and z-integration length

% r-spectra:
fu_2 = abs(fu_2).^1/length(integ)/num_files;
fv_2 = abs(fv_2).^1/length(integ)/num_files;
fw_2 = abs(fw_2).^1/length(integ)/num_files;
fq_2 = fq_2./length(integ)/num_files;
ft_r = ft_r./length(integ)/num_files;
if budget_spectra
fb_r_ik = fb_r_ik ./ num_files;
end
% theta-spectra:
fu_3 = abs(fu_3).^1/length(integ)/num_files;
fv_3 = abs(fv_3).^1/length(integ)/num_files;
fw_3 = abs(fw_3).^1/length(integ)/num_files;
fq_3 = fq_3./length(integ)/num_files;
ft_t = ft_t./length(integ)/num_files;
%%
%int_scale_x = trapz(r_wn(2:end),(1./r_wn(2:end)).*fu_2(2:end))./trapz(r_wn(2:end),fu_2(2:end));
%int_scale_z = trapz(r_wn(2:end),(1./r_wn(2:end)).*fw_2(2:end))./trapz(r_wn(2:end),fw_2(2:end));
%[~,scale] = min(abs(r_wn - 2*pi/min([int_scale_x,int_scale_z])));
%anisotropy = trapz(r_wn(scale:end),fu_2(scale:end))/trapz(r_wn(scale:end),fw_2(scale:end));
%%
% cospectra:
%fr_co = fr_co ./ length(integ) / num_files;

%covar = trapz(r_wn,fr_co);
%fr_co_norm = fr_co *trapz(r_wn,fw_2)/covar;
%fr_co_norm = fr_co;
%Ekh = trapz(r_wn,(fu_2+fv_2));

%max(kx)
%%
if compute_3d
  close all;
  h=figure('visible','off');
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
  %set(ax,'ylim',[10^(4),10^10]);
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
  %set(ax_top,'tickLength',)
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
  plot(ax,rho_wn,fu_3,'r','linewidth',2);
  plot(ax,rho_wn,fv_3,'g','linewidth',2);
  plot(ax,rho_wn,fw_3,'b','linewidth',2);
  %plot(ax,r_wn,fr_co_norm,'m','linewidth',2);

  %plot(ax,r_wn,fq_2,'m','linewidth',2);
  %plot(ax,r_wn,ft_r,'c','linewidth',2);
  lin = rho_wn.^(-5/3);

  t = title(ax,['3D velocity spectrum, ',titlepost],'fontsize',16);
  set(t,'units','inches')

  t_pos = get(t,'position');
  t_pos(2) = [t_pos(2)+0.25];
  set(t,'position',t_pos);
  YL = get(ax,'ylim');
  lin = (10.^(-10:14))'*lin/max(lin(2:end));
  plot(ax,rho_wn,lin,'k--','linewidth',2);
  set(ax,'ylim',YL);
  legend(ax,{['u^2 = ',num2str(trapz(rho_wn,fu_3),2),' m^2/s^2'],...
  ['v^2 = ',num2str(trapz(rho_wn,fv_3),2),' m^2/s^2'],...
  ['w^2 = ',num2str(trapz(rho_wn,fw_3),2),' m^2/s^2'],...
  '\kappa\^(-5/3)'},'location','southwest')
  hgexport(h,fullfile(outputpath,fileout3),hgexport('factorystyle'), 'Format', 'png');
end
close all;

if plotcloud
  h=figure('visible','off');
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
  %set(ax,'ylim',[10^(4),10^10]);
  set(ax,'xlim',[r_wn(2),r_wn(end)]);
  xticks = [10^-4,5*10^-4,10^-3,5*10^-3,10^-2,5*10^-2,];
  xticklabels = {'10^{-4}','5\times 10^{-4}','10^{-3}','5\times 10^{-3}','10^{-2}','5\times 10^{-2}'};
  set(ax,'xtick',xticks)
  set(ax,'xtickLabel',cellfun(@(c) [num2str(c,'%.0e'),'/m'],num2cell(get(ax,'xtick')),'uniformoutput',false))
  grid(ax,'on');

  ylabel(ax,'Average wv spectrum')
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
  plot(ax,r_wn,fq_2,'r','linewidth',2);
  %plot(ax,r_wn,fr_co_norm,'m','linewidth',2);

  %plot(ax,r_wn,fq_2,'m','linewidth',2);
  %plot(ax,r_wn,ft_r,'c','linewidth',2);
  lin = r_wn.^(-5/3);

  t = title(ax,['2D water vapour spectrum, ',titlepost],'fontsize',16);
  set(t,'units','inches')

  t_pos = get(t,'position');
  t_pos(2) = [t_pos(2)+0.25];
  set(t,'position',t_pos);
  YL = get(ax,'ylim');

  lin = (10.^(-10:14))'*lin/max(lin(2:end));
  plot(ax,r_wn,lin,'k--','linewidth',2);
  set(ax,'ylim',[YL(2)*10^-4,YL(2)]);
  legend(ax,{['u^2 = ',num2str(trapz(r_wn,fu_2),2),' m^3/s^2'],...
              ['v^2 = ',num2str(trapz(r_wn,fv_2),2),' m^3/s^2'],...
              ['w^2 = ',num2str(trapz(r_wn,fw_2),2),' m^3/s^2'],...
              '\kappa\^(-5/3)'},'location','southwest')
  hgexport(h,fullfile(outputpath,fileoutq),hgexport('factorystyle'), 'Format', 'png');
end
close all;
h=figure('visible','off');
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
set(ax,'ylim',[10^(-3),10^1]);
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
plot(ax,r_wn,fu_2,'r','linewidth',2);
plot(ax,r_wn,fv_2,'g','linewidth',2);
plot(ax,r_wn,fw_2,'b','linewidth',2);
%plot(ax,r_wn,fr_co_norm,'m','linewidth',2);

%plot(ax,r_wn,fq_2,'m','linewidth',2);
%plot(ax,r_wn,ft_r,'c','linewidth',2);
lin = r_wn.^(-5/3);

t = title(ax,['2D velocity spectrum, ',titlepost],'fontsize',16);
set(t,'units','inches')

t_pos = get(t,'position');
t_pos(2) = [t_pos(2)+0.25];
set(t,'position',t_pos);
YL = get(ax,'ylim');
lin = (10.^(-10:14))'*lin/max(lin(2:end));
plot(ax,r_wn,lin,'k--','linewidth',2);
set(ax,'ylim',YL);
legend(ax,{['u^2 = ',num2str(trapz(r_wn,fu_2),2),' m^3/s^2'],...
['v^2 = ',num2str(trapz(r_wn,fv_2),2),' m^3/s^2'],...
['w^2 = ',num2str(trapz(r_wn,fw_2),2),' m^3/s^2'],...
'\kappa\^(-5/3)'},'location','southwest')
hgexport(h,fullfile(outputpath,fileout2),hgexport('factorystyle'), 'Format', 'png');
if budget_spectra
close all;
h=figure('visible','off');
set(h,'Resize','on');
set(h,'PaperPosition',[0 0 16 9]);
set(h,'PaperUnits','inches');
set(h,'PaperSize',[16 9]); % IEEE columnwidth = 9cm
set(h,'Units','inches')
set(h,'Position',[0 0 16 9]);

ax = axes();

set(ax,'fontsize',13)

hold(ax,'on');
set(ax,'yscale','lin');
set(ax,'xscale','log');
%set(ax,'ylim',[10^(4),10^10]);
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
set(ax_top,'yscale','lin')
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
fadv = sum(sum(fb_r_ik(:,:,1:2),2),3);
fturb = sum(sum(fb_r_ik(:,:,4:5),2),3);
fp = sum(fb_r_ik(:,:,6),2);
fbuoy = sum(fb_r_ik(:,:,8),2);
plot(ax,r_wn,fadv,'b','linewidth',3);
plot(ax,r_wn,fturb,'r','linewidth',3);
plot(ax,r_wn,fp + fbuoy,'g','linewidth',3);
%plot(ax,r_wn,fbuoy,'m','linewidth',3);
plot(ax,r_wn,sum(fb_r_ik(:,:,7),2),'y','linewidth',2);

plot(ax,r_wn,sum(fb_r_ik(:,:,3),2),'c','linewidth',3);

[~,kr] = min(abs(r_wl - 2000))
YL = [min(min(min(sum(fb_r_ik(kr:end,:,:),2)))),max(max(max(sum(fb_r_ik(kr:end,:,:),2))))]
if YL(1) >= YL(2)
YL = [0,1];
end
%plot(ax,r_wn,fr_co_norm,'m','linewidth',2);

%plot(ax,r_wn,fq_2,'m','linewidth',2);
%plot(ax,r_wn,ft_r,'c','linewidth',2);
lin = r_wn.^(-5/3);

t = title(ax,['2D energy budget spectrum, ',titlepost],'fontsize',16);
set(t,'units','inches')

t_pos = get(t,'position');
t_pos(2) = [t_pos(2)+0.25];
set(t,'position',t_pos);
lin = (10.^(-10:14))'*lin/max(lin(2:end));
plot(ax,r_wn,lin,':','linewidth',1,'color',[0.5,0.5,0.5]);
set(ax,'ylim',YL);
set(ax_top,'ylim',YL);
legend(ax,{'Advection',...
'Subgrid Turbulence',...
'Pressure + Buoyancy',...
'Damping',...
'Diffusion',...
'\kappa\^(-5/3)'},'location','southwest')
hgexport(h,fullfile(outputpath,fileoutb),hgexport('factorystyle'), 'Format', 'png');
end

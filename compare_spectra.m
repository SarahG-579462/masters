datapath1 = '/storage/sgammon/local_v19/run/les_ConvBoundLayer/run2_50m_lowershf';
datapath2 = '/storage/sgammon/local_v19/run/les_ConvBoundLayer/run4_anelastic_lowshf';


k = 37;

%pbl(37) = 1.1899;
%pbl(73) = 1.4322;
time = zeros(1,max(k));

fw1_r = 0;
fw2_r = 0;

num_files = 0;

nx = 256;
ny = 256;
nz = 60;
kx = 2*pi*(-nx/2:(nx/2-1))/(12800);
ky = 2*pi*(-ny/2:(ny/2-1))/(12800);
kz = 2*pi*(-nz/2:(nz/2-1))/3000;
for it=1:length(k);
    filename1 = fullfile(datapath1,['cm1out_',num2str(k(it),'%06d'),'.nc']);
    filename2 = fullfile(datapath2,['cm1out_',num2str(k(it),'%06d'),'.nc']);
    if ~exist(filename1,'file')
        continue
    end
    num_files = num_files + 1;
    time(k(it)) = double(ncread(filename1,'time'));
    w1 = double(ncread(filename1,'winterp'));
    w2 = double(ncread(filename2,'vinterp'));
    
    xh = double(ncread(filename1,'xh'));
    xh_mid = xh - max(xh)/2;
    yh = double(ncread(filename1,'yh'));
    yh_mid = yh - max(yh)/2;
    z = double(ncread(filename1,'z'));  % for nodes
    zf = double(ncread(filename1,'zf')); % for edges
    [~,iz] = min(abs(z - pbl(k(it))/2));
    

    Ninteg = 1;
    integ = (iz-Ninteg):(iz+Ninteg);
    Ncirc = 5*numel(kx);
    r_wn = linspace(0,max(kx),Ncirc); % radius wavenumber
    [Xh_mid,Yh_mid] = meshgrid(xh_mid,yh_mid);
    R_rs = sqrt(Xh_mid.^2+Yh_mid.^2); % radius realspace
    Theta_rs = atan2(Yh_mid,Xh_mid);
    %[R_rs,T_rs] = meshgrid(r_rs,theta_rs);
    theta = linspace(0,2*pi,Ncirc);
    theta = theta(1:(end-1));
    [R,Theta] = meshgrid(r_wn,theta);
    for ik = integ
        for idir=1:2
            switch idir
                case 1
                    vel = squeeze(w1(:,:,ik,:));
                case 2
                    vel = squeeze(w2(:,:,ik,:));
            end
            f = abs(fftshift(fft2(vel))).^2;
            f_interp = griddedInterpolant({kx,ky},f);
            f_r_ik = zeros(size(r_wn));
            f_t_ik = zeros(size(theta));

            for ri = 1:length(r_wn)
                f_r_ik(ri) = trapz(theta,abs(r_wn(ri)).*f_interp(r_wn(ri)*cos(theta),r_wn(ri)*sin(theta)));
            end

            for ti = 1:length(theta)
                f_t_ik(ti) = trapz(r_wn,f_interp(r_wn*cos(theta(ti)),r_wn*sin(theta(ti))));
            end
            
            switch idir
                case 1
                    fw1_r = fw1_r+f_r_ik;
                case 2
                    fw2_r = fw2_r+f_r_ik;
            end
        end
    end
    
end
%%
fw2_r = abs(fw2_r).^1/length(integ)/num_files;
fw1_r = abs(fw1_r).^1/length(integ)/num_files;


%%

%max(kx)
%%
figure(2);clf(2);
ax = axes();

hold(ax,'on');
set(ax,'yscale','log');
set(ax,'xscale','log');
set(ax,'ylim',[1,10^6]);
set(ax,'xlim',[r_wn(2),r_wn(end)]);
xticks = [10^-4,5*10^-4,10^-3,5*10^-3,10^-2,5*10^-2,];
xticklabels = {'10^{-4}','5\times 10^{-4}','10^{-3}','5\times 10^{-3}','10^{-2}','5\times 10^{-2}'};
set(ax,'xtick',xticks)
set(ax,'xtickLabel',cellfun(@(c) [num2str(c,'%.0e'),'km'],num2cell(get(ax,'xtick')),'uniformoutput',false))
grid(ax,'on');

ylabel(ax,'Average velocity spectrum (m/s^2)')
xlabel(ax,'radial-wavenumber (1/m)')

r_wl = 2*pi./r_wn;
ax_top = axes();
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
set(ax_top,'ylim',[1,10^6]);
set(ax,'gridLineStyle','-')
xlabel(ax_top,'radial-wavelength (m)')
%set(ax_top,'grid','off')
plot(ax,r_wn,fw1_r,'b','linewidth',2);
plot(ax,r_wn,fw2_r,'r','linewidth',2);

lin = r_wn.^(-5/3);

title(ax_top,['Radial average of velocity spectrum in mid-pbl'],'fontsize',13)

lin = lin*1e9/max(lin(2:end));
plot(ax,r_wn,lin,'k--','linewidth',2);
legend(ax,'w time-split','w anelastic','\kappa\^(-5/3)')


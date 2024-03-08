% example spectrum:
close all
plotting = false;

orgk = rho_wn;
orgk = orgk(2:end);
maxk = max(orgk);
E_k = fw_t;
E_k = E_k(2:end);

[~,i] = min(abs(orgk - 2*pi/50));
orgk = [orgk(1:i),fliplr(logspace(log10(2*pi/10),log10(orgk(i)+(orgk(i)-orgk(i-1))),100))];
E_k  = [E_k(1:i), E_k(i)*(orgk((i+1):end)/orgk(i)).^(-5)];
figure(2);clf(2);hold on;
plot(orgk,(orgk < (2*pi/200)).*orgk.^(-5/3)+(2*pi/200)^(-5/3)/((2*pi/200)^(-5))*(orgk >= (2*pi/200)).*orgk.^(-5),'-o')
plot(orgk,E_k,'-ro')
set(gca,'xscale','log');
set(gca,'yscale','log');

%pause
%    E_k(2*pi./orgk < 400) = 1e-1;

Lx = 3200;
Ly = 3200;
Lz = 3200;
Nx = 64;
Ny = 64;
Nz = 64;
FFT_factor_3d = FFT_factor('density',[Nx,Ny,Nz],[Lx,Ly,Lz]);
A = griddedInterpolant(orgk,sqrt(E_k/FFT_factor_3d./(4*pi*orgk.^2)),'linear','nearest');

% domain size:
n = 128;

Nz = n;
Ny = n;
Nx = n;
Lx = 3200;
Ly = Lx;
Lz = Lx;
FFT_factor_3d = FFT_factor('density',[Nx,Ny,Nz],[Lx,Ly,Lz]);

dx = Lx/Nx;
dy = Ly/Ny;
dz = Lz/Nz;

U = zeros(Nx,Ny,2*Nz);
V = zeros(Nx,Ny,2*Nz);
W = zeros(Nx,Ny,2*Nz);
dk = mean(diff(orgk));

%%
Kx = 2*pi*((-Nx/2 + 0):(Nx/2-1))/Lx';
Ky = 2*pi*((-Ny/2 + 0):(Ny/2-1))/Ly';
Kz1 = 2*pi*((-Nz   + 0):(Nz  -1))/Lz';

% halfspace code stolen shamelessly from fftshift
numDims = ndims(U);
idkpos = cell(1, numDims);
idkneg = cell(1, numDims);
idk0 = cell(1,numDims);
for k = 1:numDims
    m = size(U, k);
    p = ceil(m/2);
    idkpos{k} = p+2:m;
    idkneg{k} = p:-1:2;
    idk0{k} = [1, p+1];
end


for ix = idkpos{1}
    for iy = 1:Ny
        for iz = idkpos{3};
            kx = Kx(ix);
            ky = Ky(iy);
            kz = Kz1(iz);
            
            k = sqrt(kx.^2+ky.^2+kz.^2);
            kh = sqrt(kx.^2+ky.^2);
            khat = [kx,ky,kz]/k;
            %amp = (k < 4e-3).*k.^(2) + (4e-3)^2.*(k >= 4e-3).*(k < (2*pi/200)).*k.^(-5/3) + (4e-3)^2*(2*pi/200)^(-5/3)*(k >= (2*pi/200)).*k.^(-5);
            amp = (k < (2*pi/200)).*k.^(-5/3)+(2*pi/200)^(-5/3)/((2*pi/200)^(-5))*(k >= (2*pi/200)).*k.^(-5);
            % create random phase, defined amplitude signal:
            U(ix,iy,iz) = amp.*exp(1i*2*pi*rand(1,1));
            V(ix,iy,iz) = amp.*exp(1i*2*pi*rand(1,1));
            W(ix,iy,iz) = amp.*exp(1i*2*pi*rand(1,1));

            % remove solenoidal term:
%             if abs(k) > 1e-10
%             udot = dot([U(ix,iy,iz),V(ix,iy,iz),W(ix,iy,iz)],khat);
%             U(ix,iy,iz) = U(ix,iy,iz) - udot*kx/k;
%             V(ix,iy,iz) = V(ix,iy,iz) - udot*ky/k;
%             W(ix,iy,iz) = W(ix,iy,iz) - udot*kz/k;
%             end
        end
    end
end

%ensure conjugate symmetric and boundary conditions simultaneously:
%  there may be a better way than this, but this works for now...
[U1,V1,W1] = deal(zeros(size(U)),zeros(size(V)),zeros(size(W)));
% U1(idkpos{:}) = 0.5*(     U(idkpos{:})  + conj(U(idkneg{:})));
% U1(idkneg{:}) = 0.5*(conj(U(idkpos{:})) +      U(idkneg{:}));
% V1(idkpos{:}) = 0.5*(     V(idkpos{:})  + conj(V(idkneg{:})));
% V1(idkneg{:}) = 0.5*(conj(V(idkpos{:})) +      V(idkneg{:}));
U1(idkpos{1},idkpos{2},idkpos{3}) =      U(idkpos{1},idkpos{2},idkpos{3});
U1(idkneg{1},idkneg{2},idkneg{3}) = conj(U(idkpos{1},idkpos{2},idkpos{3}));
U1(idkneg{1},idkneg{2},idkpos{3}) = conj(U(idkpos{1},idkpos{2},idkpos{3}));
U1(idkpos{1},idkpos{2},idkneg{3}) =      U(idkpos{1},idkpos{2},idkpos{3});

U1(idkpos{1},idkneg{2},idkpos{3}) =      U(idkpos{1},idkneg{2},idkpos{3});
U1(idkneg{1},idkpos{2},idkneg{3}) = conj(U(idkpos{1},idkneg{2},idkpos{3}));
U1(idkneg{1},idkpos{2},idkpos{3}) = conj(U(idkpos{1},idkneg{2},idkpos{3}));
U1(idkpos{1},idkneg{2},idkneg{3}) =      U(idkpos{1},idkneg{2},idkpos{3});

V1(idkpos{1},idkpos{2},idkpos{3}) =      V(idkpos{1},idkpos{2},idkpos{3});
V1(idkneg{1},idkneg{2},idkneg{3}) = conj(V(idkpos{1},idkpos{2},idkpos{3}));
V1(idkneg{1},idkneg{2},idkpos{3}) = conj(V(idkpos{1},idkpos{2},idkpos{3}));
V1(idkpos{1},idkpos{2},idkneg{3}) =      V(idkpos{1},idkpos{2},idkpos{3});

V1(idkpos{1},idkneg{2},idkpos{3}) =      V(idkpos{1},idkneg{2},idkpos{3});
V1(idkneg{1},idkpos{2},idkneg{3}) = conj(V(idkpos{1},idkneg{2},idkpos{3}));
V1(idkneg{1},idkpos{2},idkpos{3}) = conj(V(idkpos{1},idkneg{2},idkpos{3}));
V1(idkpos{1},idkneg{2},idkneg{3}) =      V(idkpos{1},idkneg{2},idkpos{3});

W1(idkpos{1},idkpos{2},idkpos{3}) =      W(idkpos{1},idkpos{2},idkpos{3});
W1(idkneg{1},idkneg{2},idkneg{3}) = conj(W(idkpos{1},idkpos{2},idkpos{3}));
W1(idkneg{1},idkneg{2},idkpos{3}) =-conj(W(idkpos{1},idkpos{2},idkpos{3}));
W1(idkpos{1},idkpos{2},idkneg{3}) =     -W(idkpos{1},idkpos{2},idkpos{3});

W1(idkpos{1},idkneg{2},idkpos{3}) =      W(idkpos{1},idkneg{2},idkpos{3});
W1(idkneg{1},idkpos{2},idkneg{3}) = conj(W(idkpos{1},idkneg{2},idkpos{3}));
W1(idkneg{1},idkpos{2},idkpos{3}) =-conj(W(idkpos{1},idkneg{2},idkpos{3}));
W1(idkpos{1},idkneg{2},idkneg{3}) =     -W(idkpos{1},idkneg{2},idkpos{3});

U1(idk0{:}) = complex(real(U(idk0{:})),0);
V1(idk0{:}) = complex(real(V(idk0{:})),0);
W1(idk0{:}) = complex(real(W(idk0{:})),0);
[U,V,W] = deal(U1,V1,W1);

%% verify magnitudes:
rho_wn_new = linspace(0,min([max(Kx),max(Ky),max(Kz1)]),Ncirc); % radius wavenumber
rho_wn_new = rho_wn_new(2:end);
FU_new = griddedInterpolant({Kx,Ky,Kz1},FFT_factor_3d*abs(U).^2);
FV_new = griddedInterpolant({Kx,Ky,Kz1},FFT_factor_3d*abs(V).^2);
FW_new = griddedInterpolant({Kx,Ky,Kz1},FFT_factor_3d*abs(W).^2);

E_k_new_u = zeros(size(rho_wn_new));
E_k_new_v = zeros(size(rho_wn_new));
E_k_new_w = zeros(size(rho_wn_new));

for ri = 1:length(rho_wn_new)
    r = rho_wn_new(ri);
    E_k_new_u(ri) = trapz(theta_wn,trapz(phi_wn,r.^2.*sin(P).*FU_new(r.*cos(T).*sin(P),r.*sin(T).*sin(P),r.*cos(P)),1),2);
    E_k_new_v(ri) = trapz(theta_wn,trapz(phi_wn,r.^2.*sin(P).*FV_new(r.*cos(T).*sin(P),r.*sin(T).*sin(P),r.*cos(P)),1),2);
    E_k_new_w(ri) = trapz(theta_wn,trapz(phi_wn,r.^2.*sin(P).*FW_new(r.*cos(T).*sin(P),r.*sin(T).*sin(P),r.*cos(P)),1),2);

end

rmin = max(min(rho_wn_new),min(orgk));
rmax = min(max(rho_wn_new),max(orgk));

[~,i1] = min(abs(orgk - rmin));
[~,i2] = min(abs(orgk - rmax));
energy_org = trapz(orgk(i1:i2),E_k(i1:i2));
[~,i1] = min(abs(rho_wn_new - rmin));
[~,i2] = min(abs(rho_wn_new - rmax));
energy_new = trapz(rho_wn_new(i1:i2),E_k_new_u(i1:i2));
U = 2*U*sqrt(energy_org/energy_new);
energy_new = trapz(rho_wn_new(i1:i2),E_k_new_v(i1:i2));
V = 2*V*sqrt(energy_org/energy_new);
energy_new = trapz(rho_wn_new(i1:i2),E_k_new_w(i1:i2));
W = 2*W*sqrt(energy_org/energy_new);

if plotting
    FU_new = griddedInterpolant({Kx,Ky,Kz1},FFT_factor_3d*abs(U).^2);
    FV_new = griddedInterpolant({Kx,Ky,Kz1},FFT_factor_3d*abs(V).^2);
    FW_new = griddedInterpolant({Kx,Ky,Kz1},FFT_factor_3d*abs(W).^2);

    E_k_new_u = zeros(size(rho_wn_new));
    E_k_new_v = zeros(size(rho_wn_new));
    E_k_new_w = zeros(size(rho_wn_new));

    for ri = 1:length(rho_wn_new)
        r = rho_wn_new(ri);
        E_k_new_u(ri) = trapz(theta_wn,trapz(phi_wn,r.^2.*sin(P).*FU_new(r.*cos(T).*sin(P),r.*sin(T).*sin(P),r.*cos(P)),1),2);
        E_k_new_v(ri) = trapz(theta_wn,trapz(phi_wn,r.^2.*sin(P).*FV_new(r.*cos(T).*sin(P),r.*sin(T).*sin(P),r.*cos(P)),1),2);
        E_k_new_w(ri) = trapz(theta_wn,trapz(phi_wn,r.^2.*sin(P).*FW_new(r.*cos(T).*sin(P),r.*sin(T).*sin(P),r.*cos(P)),1),2);

    end
end
%%
u = zeros(Nx+1,Ny,Nz);
us = ifftn(ifftshift(U));
u(2:Nx,1:Ny,1:Nz) = 0.5*us(2:Nx,1:Ny,1:Nz) + 0.5*us(1:(Nx-1),1:Ny,1:Nz);
u(Nx+1,:,:) = 0.5*us(1,1:Ny,1:Nz)+0.5*us(Nx,1:Ny,1:Nz);
u(1,:,:) = 0.5*us(1,1:Ny,1:Nz)+0.5*us(Nx,1:Ny,1:Nz);

v = zeros(Nx,Ny+1,Nz);
vs = ifftn(ifftshift(V));
v(1:Nx,2:Ny,1:Nz) = 0.5*vs(1:Nx,2:Ny,1:Nz) + 0.5*vs(1:Nx,1:(Ny-1),1:Nz);
v(:,Ny+1,:) = 0.5*vs(1:Nx,1,1:Nz)+0.5*vs(1:Nx,Ny,1:Nz);
v(:,1,:) = 0.5*vs(1:Nx,1,1:Nz)+0.5*vs(1:Nx,Ny,1:Nz);

w = zeros(Nx,Ny,Nz+1);
ws = ifftn(ifftshift(W));
w(1:Nx,1:Ny,2:Nz) = 0.5*ws(1:Nx,1:Ny,2:Nz) + 0.5*ws(1:Nx,1:Ny,1:(Nz-1));
w(:,:,Nz+1) = 0;%ws(1:Nx,1:Ny,Nz);
w(:,:,1) = 0;
%w(:,:,1) = 0;
%u(:,:,1) = u(:,:,2);
%v(:,:,1) = v(:,:,2);
Kz = 2*pi*((-Nz/2   + 0):(Nz/2  -1))/Lz';

p = ((Nx+1)*(Ny+1)*(Nz+1));
[i,j,k] = ind2sub([Nx+1,Ny+1,Nz+1],1:p);
write_mat = zeros(p,6);
for ip = 1:p
    write_mat(ip,:) = [i(ip),j(ip),k(ip),0,0,0];
    try
        write_mat(ip,4) = u(i(ip),j(ip),k(ip));
    catch
    end
    try
        write_mat(ip,5) = v(i(ip),j(ip),k(ip));
    catch
    end
    try
        write_mat(ip,6) = w(i(ip),j(ip),k(ip));
    catch
    end
end

csvwrite(['init3d',num2str(n),'.csv'],write_mat);

%% Analyze obtained spectra:
if plotting
    fu_ri = 0;
    fv_ri = 0;
    fw_ri = 0;
    fr_ri = 0;
    ft_ri = 0;

    fu_ti = 0;
    fv_ti = 0;
    fw_ti = 0;
    fr_ti = 0;
    ft_ti = 0;
    iz = Nz/2;
    e2d = 0;
    e3d = 0;

    Ninteg = 0; % radius of z-integration points around pbl/2
    integ = (iz-Ninteg):(iz+Ninteg);
    Ncirc = 250;%5*numel(Kx);

    [Xh_mid,Yh_mid] = meshgrid(xh_mid,yh_mid);
    R_rs = sqrt(Xh_mid.^2+Yh_mid.^2); % radius realspace
    Theta_rs = atan2(Yh_mid,Xh_mid);

    r_wn = linspace(0,min([max(Kx),max(Ky)]),Ncirc); % radius wavenumber
    rho_wn_new = linspace(0,min([max(Kx),max(Ky),max(Kz)]),Ncirc); % radius wavenumber
    rho_wn_new = rho_wn_new(2:end);
    theta_wn = linspace(0,2*pi,Ncirc);
    theta_wn = theta_wn(1:(end-1));

    phi_wn = linspace(0,pi,Ncirc);
    % Factor so that FFT scales like kinetic energy density
    FFT_factor_2d = FFT_factor('density',[Nx,Ny],[Lx,Ly]);
    FFT_factor_3d = FFT_factor('density',[Nx,Ny,Nz],[Lx,Ly,Lz]);
    % Total kinetic energy (average over time):

    % kinetic energy:        
    for idir=1:3
        switch idir
            case 1
                vel = squeeze(.5*u(1:Nx,:,:,:)+.5*u(1+(1:Nx),:,:,:));
            case 2
                vel = squeeze(.5*v(:,1:Ny,:,:)+.5*v(:,1+(1:Ny),:,:));
            case 3
                vel = squeeze(.5*w(:,:,1:Nz,:)+.5*w(:,:,1+(1:Nz),:));
            case 4
                vel = (1./(R_rs)).*(Xh_mid.*squeeze(u(:,:,ik,:)) + Yh_mid.*squeeze(v(:,:,iz,:)));
            case 5
                vel = (1./(R_rs)).*(Xh_mid.*squeeze(v(:,:,ik,:)) - Yh_mid.*squeeze(u(:,:,iz,:)));
        end
        e2d = e2d + dx.*dy.*sum(sum(vel(:,:,iz).^2));
        e3d = e3d + dx.*dy.*dz.*sum(sum(sum(vel.^2)));

        f2 = FFT_factor_2d*abs(fftshift(fft2(vel(:,:,iz)))).^2;
        f2_interp = griddedInterpolant({Kx,Ky},f2);
        f_2d_ik = zeros(size(r_wn));

        % 2d spectra:
        for ri = 1:length(r_wn)
            f_2d_ik(ri) = trapz(theta_wn,abs(r_wn(ri)).*f2_interp(r_wn(ri)*cos(theta_wn'),r_wn(ri)*sin(theta_wn')));
        end

        f3 = FFT_factor_3d*abs(fftshift(fftn(vel))).^2;
        f3_interp = griddedInterpolant({Kx,Ky,Kz},f3);                        
        f_3d_ik = zeros(size(rho_wn_new));

        % 3d spectra:
        [T,P] = meshgrid(theta_wn,phi_wn);
        for ri = 1:length(rho_wn_new)
            r = rho_wn_new(ri);
            f_3d_ik(ri) = trapz(theta_wn,trapz(phi_wn,r.^2.*sin(P).*f3_interp(r.*cos(T).*sin(P),r.*sin(T).*sin(P),r.*cos(P)),1),2);
        end

        switch idir
            case 1
                fu_ri = fu_ri+f_2d_ik;
                fu_ti = fu_ti+f_3d_ik;
            case 2
                fv_ri = fv_ri+f_2d_ik;
                fv_ti = fv_ti+f_3d_ik;
            case 3
                fw_ri = fw_ri+f_2d_ik;
                fw_ti = fw_ti+f_3d_ik;
            case 4
                fr_ri = fr_ri + f_2d_ik;
                fr_ti = fr_ti + f_3d_ik;
            case 5
                ft_ri = ft_ri + f_2d_ik;
                ft_ti = ft_ti + f_3d_ik;
        end
    end


    fprintf('%i, %.0f, %.5f, %.5f, %.5f\n',n,Lx,energy_new/energy_org,energy_org,energy_new)


    e2d_f = trapz(r_wn,(fu_ri+fv_ri+fw_ri));
    e3d_f = trapz(rho_wn_new,(fu_ti+fv_ti+fw_ti));

    h = figure(1);clf(1)
    %h=figure('visible','on');
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
    %set(ax,'xlim',[r_wn(2),r_wn(end)]);
    xticks = [10^-4,5*10^-4,10^-3,5*10^-3,10^-2,5*10^-2,];
    xticklabels = {'10^{-4}','5\times 10^{-4}','10^{-3}','5\times 10^{-3}','10^{-2}','5\times 10^{-2}'};
    set(ax,'xtick',xticks)
    set(ax,'xtickLabel',cellfun(@(c) [num2str(c,'%.0e'),'/m'],num2cell(get(ax,'xtick')),'uniformoutput',false))
    grid(ax,'on');

    ylabel(ax,'Average velocity spectrum (m^5/s^2)')
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
    plot(ax,rho_wn_new,fu_ti,'r','linewidth',2);
    plot(ax,rho_wn_new,fv_ti,'g','linewidth',2);
    plot(ax,rho_wn_new,fw_ti,'b','linewidth',2);
    plot(ax,orgk,E_k,'m','linewidth',2);

    %plot(ax,r_wn,fr_r,'m','linewidth',2);
    %plot(ax,r_wn,ft_r,'c','linewidth',2);
    lin = r_wn.^(-5/3);

    t = title(ax,['Radial average of velocity spectrum, averaged over ' num2str(abs(min(time(time > 0)) - max(time(time > 0)))/3600,2), 'h. Anisotropy= ',num2str(anisotropy,2),titlepost],'fontsize',16);
    set(t,'units','inches')

    t_pos = get(t,'position');
    t_pos(2) = [t_pos(2)+0.25];
    set(t,'position',t_pos);
    YL = get(ax,'ylim');
    XL = get(ax,'xlim');

    lin = (10.^(-10:14))'*lin/max(lin(2:end));
    plot(ax,r_wn,lin,'k--','linewidth',2);
    set(ax,'ylim',YL);
    legend(ax,{['u^2 = ',num2str(trapz(rho_wn_new,fu_ti),2),' m^4/s^2'],...
               ['v^2 = ',num2str(trapz(rho_wn_new,fv_ti),2),' m^4/s^2'],...
               ['w^2 = ',num2str(trapz(rho_wn_new,fw_ti),2),' m^4/s^2'],...
               ['original = ',num2str(energy_org,2)],...
               '\kappa\^(-5/3)'},'location','southwest')
 end
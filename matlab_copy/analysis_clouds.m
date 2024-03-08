 colors = {'r','g','b','m','c','k'};
 res = 1;

 if ~exist('datapath','var')
    datapath = 'example_data';
    fileout = 'example_data_clouds';
    outputpath = './';
elseif exist('postfix','var');
    [~,dp_header,~] = fileparts(datapath);
    fileout = [dp_header,'_clouds_',postfix];
else
    [~,dp_header,~] = fileparts(datapath);
    fileout = [dp_header,'_clouds_'];
end

if ~exist('titlepost','var')
    titlepost = '';
end

init_u = 0;
init_v = 0;
init_q = 0;
init_T = 0;
cover_t = nan(numframes,1);
lwp_t   = nan(numframes,1);
tke_t   = nan(numframes,4);

iz = 0;
th_z = 0;
qv_z = 0;
u_z  = 0;
v_z  = 0;
ql_z = 0;

izz = 0;
upwp  = 0;
wpqlp = 0;
wpthvp= 0;
wpthlp= 0;
wpqvp = 0;
sigma_u = 0;
sigma_w = 0;
tke_z = 0;
ww_z  = 0;
uu_z  = 0;

time = nan(1,numframes);
num_files = 0;
shown = 0;
i_z_show = 0;
i_z_max = 10;

R = 287.04;
epsilon = R/461.5;
Lv = 2501000;
cpd = 1005.7;
for it = 1:numframes;
    % I/O:
    filename = fullfile(datapath,[dataprefix,num2str(it,'%06d'),'.nc']);
    if ~exist(filename,'file')
        continue
    end
    num_files = num_files + 1;
    time(it) = squeeze(double(ncread(filename,'time')))/60;
    try
            tke_all = squeeze(double(ncread(filename,'tke')));
    catch
	    tke_all = 0;
    end
    th = squeeze(double(ncread(filename,'th')));
    %hpbl = squeeze(double(ncread(filename,'hpbl')));

    u = squeeze(double(ncread(filename,'uinterp')));
    v = squeeze(double(ncread(filename,'vinterp')));
    w = squeeze(double(ncread(filename,'winterp')));
    xh = double(ncread(filename,'xh'));
    yh = double(ncread(filename,'yh'));
    z = double(ncread(filename,'z'));  % for nodes
    zf = double(ncread(filename,'zf')); % for edges

    qv  = squeeze(double(ncread(filename,'qv')));
    prs = squeeze(double(ncread(filename,'prs')));
    rho = squeeze(double(ncread(filename,'rho')));
    t =  prs./(rho.*R.*(1+qv/epsilon)); % temperature
    ql  = squeeze(double(ncread(filename,'qc')));
    thl = th - (th./t).*(Lv/cpd).*ql; % approximation for liquid water th
    thv = th.*(1 + 0.61*qv - ql); % virtual potential temp
    lwp = squeeze(double(ncread(filename,'lwp')));
    
    nx = numel(xh);
    ny = numel(yh);
    nz = numel(z);

    dx = mean(diff(xh));
    dy = mean(diff(yh));
    dz = mean(diff(zf));
    
    Lx = nx*dx;
    Ly = ny*dy;
    Lz = nz*dz;
    
    dS = dx*dy/(Lx*Ly);
    dV = dS*dz/Lz;
    %% initial state:
    if it == 1
    	init_u = dS.*squeeze(sum(sum(u,1),2));
	init_v = dS.*squeeze(sum(sum(v,1),2));
	init_q = dS.*squeeze(sum(sum(qv,1),2));
	init_T = dS* squeeze(sum(sum(thl,1),2));
    end 

    %% timeseries:
    % cloud cover:
    cover_t(it) = sum(sum(sum(ql > 1e-6,3),1),2)./(nx.*ny);
    % lwp::
    lwp_t(it) = dS.*sum(sum(lwp,1),2);
    % integrate TKE over all space, leaving only time as a variable.
    % (moved to below)

    %% Time Invariants:
    
    % Base state:  q = qz(z,t) + pq
    thz = sum(sum(th,1),2).*dS;
    wz  = sum(sum(w,1),2).*dS;
    uz  = sum(sum(u,1),2).*dS;
    vz  = sum(sum(v,1),2).*dS;

    qlz = dS*sum(sum(ql,1),2);
    qvz = dS*sum(sum(qv,1),2);
    thlz= dS*sum(sum(thl,1),2);
    thvz= dS*sum(sum(thv,1),2);
    % Perturbation potential temperature:
    pth = th - repmat(thz,size(th,1),size(th,2),1);
    pw = w - repmat(wz,size(w,1),size(w,2),1);
    pu = u - repmat(uz,size(u,1),size(u,2),1);
    pv = v - repmat(vz,size(v,1),size(v,2),1);
    
    pql = ql - repmat(qlz,size(ql,1),size(ql,2),1);
    pqv = qv - repmat(qvz,size(qv,1),size(qv,2),1);
    pthl=thl - repmat(thlz,size(thl,1),size(thl,2),1);
    pthv=thv - repmat(thvz,size(thv,1),size(thv,2),1);

    %tke_t(it,1) = sum(sum(sum(rho.*tke_all,1),2),3)*dV;
    tke_t(it,2) = sum(sum(sum(rho.*0.5.*(pu.^2 + pv.^2 + pw.^2),1),2),3).*dV;
    tke_t(it,3) = sum(sum(sum(rho.*0.5.*(u.^2 + v.^2 + w.^2),1),2),3)*dV;
    tke_t(it,4) = tke_t(it,1)+tke_t(it,2);
    
    
    if time(it) >= 5*60
	    iz = iz + 1;
	    th_z = ((iz-1)*th_z + squeeze(thz))/iz;
	    qv_z = ((iz-1)*qv_z + squeeze(qvz))/iz;
	    u_z  = ((iz-1)*u_z  + squeeze(uz))/iz;
	    v_z  = ((iz-1)*v_z  + squeeze(vz))/iz;
	    ql_z = ((iz-1)*ql_z + squeeze(qlz))/iz;
    end

    if time(it) >= 3*60
	    izz = izz+1;
	    upwp  = ((izz-1)*upwp + squeeze(sum(sum(pu.*pw,1),2).*dS))./izz;
	    wpqlp = ((izz-1)*wpqlp+ squeeze(sum(sum(pw.*pql,1),2).*dS))/izz;
	    wpthvp= ((izz-1)*wpthvp+squeeze(sum(sum(pw.*pthv,1),2).*dS))/izz;
	    wpthlp= ((izz-1)*wpthlp+squeeze(sum(sum(pw.*pthl,1),2).*dS))/izz;
	    wpqvp = ((izz-1)*wpqvp +squeeze(sum(sum(pw.*pqv,1),2).*dS))/izz;
	    sigma_u = squeeze(dS*sum(sum(pu.*pu + pv.*pv,1),2));
	    sigma_w = squeeze(dS*sum(sum(pw.*pw,1),2));
	    tke_z = ((iz-1)*tke_z + sigma_u + sigma_w)/izz;
	    ww_z  = ((iz-1)*ww_z  + sigma_w)/izz;
	    uu_z  = ((iz-1)*uu_z  + sigma_u)/izz;
    end
end

close all;
h=figure('visible','off'); hold on;
set(h,'Resize','on');
%set(h,'PaperPositionMode','auto');
set(h,'PaperPosition',[0 0 16 9]);
set(h,'PaperUnits','inches');
set(h,'PaperSize',[16 9]); % IEEE columnwidth = 9cm
set(h,'Units','inches')
set(h,'Position',[0 0 16 9]);
subplot(3,1,1)
set(gca,'fontsize',13)
plot(time, cover_t*100,'b-','linewidth',2)
title('Cloud cover percentage')
ylabel('Cloud Cover (%)');
grid on;
subplot(3,1,2);
set(gca,'fontsize',13)
plot(time, lwp_t,'b-','linewidth',2);
title('Liquid water path')
ylabel('Liquid water path (gm^{-2})')
grid on;
subplot(3,1,3);
set(gca,'fontsize',13); hold on;
plot(time,tke_t(:,1),'b-','linewidth',2)
plot(time,tke_t(:,2),'g-','linewidth',2)
plot(time,tke_t(:,4),'r-','linewidth',2)
legend('subgrid','resolved','total');
xlabel('time (min)');
ylabel('Average TKE (kgm^2/s^2)');
title(['Average TKE vs time ',titlepost])
grid on;
label = 'timeseries';
filename = fullfile(outputpath,[fileout,label,'.png']);
hgexport(h, filename, hgexport('factorystyle'), 'Format', 'png');

close all;
h=figure('visible','off'); hold on;
set(h,'Resize','on');
%set(h,'PaperPositionMode','auto');
set(h,'PaperPosition',[0 0 16 9]);
set(h,'PaperUnits','inches');
set(h,'PaperSize',[16 9]); % IEEE columnwidth = 9cm
set(h,'Units','inches')
set(h,'Position',[0 0 16 9]);
subplot(1,3,1);
set(gca,'fontsize',13); hold on;
plot(init_u,z,'bo-','linewidth',2)
plot(init_v,z,'ro-','linewidth',2)
ylabel('z (km)');
xlabel('velocity (m/s)');
title(['Initial wind profile',titlepost])
grid on;
subplot(1,3,2)
set(gca,'fontsize',13)
plot(init_q,z,'bo-','linewidth',2);
xlabel('q_v (g/g)')
title('Initial moisture profile');
grid on;
subplot(1,3,3)
set(gca,'fontsize',13)
plot(init_T,z,'bo-','linewidth',2)
grid on;
xlabel('\theta_l (K)')
title('Initial liquid water potential temp')
grid on;
label = 'init';
filename = fullfile(outputpath,[fileout,label,'.png']);
hgexport(h, filename, hgexport('factorystyle'), 'Format', 'png');

close all;
h=figure('visible','off'); hold on;
set(h,'Resize','on');
%set(h,'PaperPositionMode','auto');
set(h,'PaperPosition',[0 0 16 9]);
set(h,'PaperUnits','inches');
set(h,'PaperSize',[16 9]); % IEEE columnwidth = 9cm
set(h,'Units','inches')
set(h,'Position',[0 0 16 9]);
subplot(2,2,1);
set(gca,'fontsize',13); hold on;
plot(init_u,z,'b:','linewidth',2)
plot(u_z,z,'b-','linewidth',2)
plot(init_v,z,'r:','linewidth',2)
plot(v_z,z,'r-','linewidth',2)
ylabel('z (km)');
xlabel('velocity (m/s)');
legend('u_{init}','u_{ss}','v_{init}','v_{ss}')
title(['Steady-state wind profile',titlepost])
grid on;
subplot(2,2,2); hold on;
set(gca,'fontsize',13); hold on;
plot(init_q,z,'b:','linewidth',2);
plot(qv_z,z,'b-','linewidth',2);
xlabel('q_v (g/g)')
title('Steady-state vapour profile');
grid on;
subplot(2,2,3); hold on;
set(gca,'fontsize',13)
plot(init_T,z,'b:','linewidth',2)
plot(th_z,z,'b-','linewidth',2)
grid on;
xlabel('\theta_l (K)')
title('Steady-state potential temp')
grid on;
subplot(2,2,4); hold on;
set(gca,'fontsize',13);
plot(ql_z,z,'b-','linewidth',2);
grid on;
xlabel('q_l (g/g)');
title('Steady state liquid profile');
label = 'steady';
filename = fullfile(outputpath,[fileout,label,'.png']);
hgexport(h, filename, hgexport('factorystyle'), 'Format', 'png');

close all;
h=figure('visible','off'); hold on;
set(h,'Resize','on');
%set(h,'PaperPositionMode','auto');
set(h,'PaperPosition',[0 0 16 9]);
set(h,'PaperUnits','inches');
set(h,'PaperSize',[16 9]); % IEEE columnwidth = 9cm
set(h,'Units','inches')
set(h,'Position',[0 0 16 9]);
subplot(2,3,1);
set(gca,'fontsize',13); hold on;
plot(wpqvp,z,'r-','linewidth',2)
ylabel('z (km)');
xlabel('\langle w''q_v''\rangle');
title(['vapour flux',titlepost])
grid on;
subplot(2,3,2);
set(gca,'fontsize',13); hold on;
plot(wpthvp,z,'r-','linewidth',2)
ylabel('z (km)');
xlabel('\langle w''th_v''\rangle');
title(['virtual potential temp flux',titlepost])
grid on;
subplot(2,3,3);
set(gca,'fontsize',13); hold on;
plot(upwp,z,'r-','linewidth',2)
ylabel('z (km)');
xlabel('\langle w''u''\rangle');
title(['momentum flux',titlepost])
grid on;
subplot(2,3,4);
set(gca,'fontsize',13); hold on;
plot(wpqlp,z,'r-','linewidth',2)
ylabel('z (km)');
xlabel('\langle w''q_l''\rangle');
title(['liquid flux',titlepost])
grid on;
subplot(2,3,5);
set(gca,'fontsize',13); hold on;
plot(wpthlp,z,'r-','linewidth',2)
ylabel('z (km)');
xlabel('\langle w''th_v''\rangle');
title(['liquid pot.temp. flux',titlepost])
grid on;
label = 'flux';
filename = fullfile(outputpath,[fileout,label,'.png']);
hgexport(h, filename, hgexport('factorystyle'), 'Format', 'png');

close all;
h=figure('visible','off'); hold on;
set(h,'Resize','on');
%set(h,'PaperPositionMode','auto');
set(h,'PaperPosition',[0 0 16 9]);
set(h,'PaperUnits','inches');
set(h,'PaperSize',[16 9]); % IEEE columnwidth = 9cm
set(h,'Units','inches')
set(h,'Position',[0 0 16 9]);
set(gca,'fontsize',13)
hold on;
plot(tke_z,z,'b-','linewidth',2)
plot(uu_z,z,'r-','linewidth',2)
plot(ww_z,z,'g-','linewidth',2)
ylabel('z (km)');
xlabel('TKE components');
legend('TKE','u''^2','w''^2')
title(['TKE vs height, last 3 hrs',titlepost])
grid on;
label = 'variance';
filename = fullfile(outputpath,[fileout,label,'.png']);
hgexport(h, filename, hgexport('factorystyle'), 'Format', 'png');

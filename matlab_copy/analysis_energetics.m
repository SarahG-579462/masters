
%datapath = '/storage/sgammon/local_v19/run/les_ConvBoundLayer/run2_50m_lowershf';
% 
% datapaths = {...
%     '/storage/sgammon/local_v19/run/les_ConvBoundLayer/run7_50m_3k'
%     '/storage/sgammon/local_v19/run/les_ConvBoundLayer/run10_50m_anelastic'
%     '/storage/sgammon/local_v19/run/les_ConvBoundLayer/run11_50m_weno_off'
%     '/storage/sgammon/local_v19/run/les_ConvBoundLayer/run12_50m_weno_on'
%     '/storage/sgammon/local_v19/run/les_ConvBoundLayer/run13_50m_sgs_off'
%     '/storage/sgammon/local_v19/run/les_ConvBoundLayer/run14_50m_sgs_alt'
%     '/storage/sgammon/local_v19/run/les_ConvBoundLayer/run15_50m_vert'
%     '/storage/sgammon/local_v19/run/les_ConvBoundLayer/run16_50m_incomp'
%     '/storage/sgammon/local_v19/run/les_ConvBoundLayer/run17_50m_hiordweno'
%     '/storage/sgammon/local_v19/run/les_ConvBoundLayer/run18_50m_hiordall'
%     '/storage/sgammon/local_v19/run/les_ConvBoundLayer/run19_50m_noadaptdt'
% };



% prefix = '';
 colors = {'r','g','b','m','c','k'};
 res = 1;
% labels = {
%     '50m'
%     '25m'
%     '12.5m'};
% figure(1);clf(1);
% figure(2);clf(2);
% figure(3);clf(3);
% figure(4);clf(4);
% figure(5);clf(5);
% figure(6);clf(6);
% figure(7);clf(7);
% 
% for res = 1:length(datapaths);
% datapath = datapaths{res};

if ~exist('datapath','var')
    datapath = 'example_data';
    fileout = 'example_data_energetics';
    outputpath = './';
elseif exist('postfix','var');
    [~,dp_header,~] = fileparts(datapath);
    fileout = [dp_header,'_energetics_',postfix];
else
    [~,dp_header,~] = fileparts(datapath);
    fileout = [dp_header,'_energetics_'];
end

if ~exist('titlepost','var')
    titlepost = '';
end

tke = nan(numframes,3);
sigmath_t = nan(1,numframes);
sigmath_z = 0;
sigmaw_t = nan(1,numframes);
sigmau_t = nan(1,numframes);
sigmav_t = nan(1,numframes);
pbl = nan(1,numframes);
sigmaw_z = 0;
sigmau_z = 0;
sigmav_z = 0;
wth_t = nan(1,numframes);
wth_z = 0;
time = nan(1,numframes);
num_files = 0;
shown = 0;
i_z_show = 0;
i_z_max = 10;
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
    
    nx = numel(xh);
    ny = numel(yh);
    nz = numel(z);

    dx = mean(diff(xh));
    dy = mean(diff(yh));
    dz = mean(diff(zf));
    
    Lx = nx*dx;
    Ly = ny*dy;
    Lz = nz*dz;
    
    
    % integrate TKE over all space, leaving only time as a variable.
    tke(it,1) = sum(sum(sum(tke_all,1),2),3).*dx.*dy.*dz./(Lx.*Ly.*Lz);
    tke(it,2) = sum(sum(sum(0.5*(u.^2 + v.^2 + w.^2),1),2),3).*dx.*dy.*dz./(Lx.*Ly.*Lz);
    tke(it,3) = tke(it,1)+tke(it,2);
    % Base state:  q = qz(z,t) + pq
    thz = sum(sum(th,1),2).*dx.*dy / (Lx.*Ly);
    wz  = sum(sum(w,1),2).*dx.*dy / (Lx.*Ly);
    uz  = sum(sum(u,1),2).*dx.*dy / (Lx.*Ly);
    vz  = sum(sum(v,1),2).*dx.*dy / (Lx*Ly);
    % Perturbation potential temperature:
    pth = th - repmat(thz,size(th,1),size(th,2),1);
    pw = w - repmat(wz,size(w,1),size(w,2),1);
    pu = u - repmat(uz,size(u,1),size(u,2),1);
    pv = v - repmat(vz,size(v,1),size(v,2),1);
    
    % <pth^2> = 
    sigmath_t(it) = sum(sum(sum(pth.*pth,1),2),3).*dx.*dy.*dz/ (Lx.*Ly.*Lz);
    sigmath_z = sigmath_z + squeeze(sum(sum(pth.*pth,1),2).*dx.*dy)/ (Lx.*Ly);
    
    % <w^2> = 
    sigmaw_t(it) = sum(sum(sum(pw.*pw,1),2),3).*dx.*dy.*dz/ (Lx.*Ly.*Lz);
    sigmau_t(it) = sum(sum(sum(pu.*pu,1),2),3).*dx.*dy.*dz/ (Lx.*Ly.*Lz);
    sigmav_t(it) = sum(sum(sum(pv.*pv,1),2),3).*dx.*dy.*dz/ (Lx.*Ly.*Lz);
    sigmaw_z = sigmaw_z + squeeze(sum(sum(pw.*pw,1),2).*dx.*dy)/ (Lx.*Ly);
    sigmau_z = sigmau_z + squeeze(sum(sum(pu.*pu,1),2).*dx.*dy)/ (Lx.*Ly);
    sigmav_z = sigmav_z + squeeze(sum(sum(pv.*pv,1),2).*dx.*dy)/ (Lx.*Ly);
    % <pw.*pth> = 
    wth_t(it) = sum(sum(sum(pth.*pw,1),2),3).*dx.*dy.*dz/ (Lx.*Ly.*Lz);
    wth_z = wth_z + squeeze(sum(sum(pth.*pw,1),2).*dx.*dy)/ (Lx.*Ly);
    
    brunt_vaisalla = double(ncread(filename,'nm'));


    [~,iz] = max(brunt_vaisalla,[],3);
    pbl(it) = mean(zf(iz(:)));

    if ~shown && (num_files/numframes >= i_z_show/i_z_max);
        if i_z_show == 0
                theta_z_waterfall = zeros(i_z_max,numel(thz));
	        waterfall_label = cell(i_z_max,1);
	end
        shown = 1;
        i_z_show = i_z_show + 1;
        theta_z_waterfall(i_z_show,:) = thz - 300 + (i_z_show-1)*10;
        waterfall_label{i_z_show} = sprintf('%.1f min',time(it));
    end
    if ((num_files+1)/numframes >= i_z_show/i_z_max)
        shown = 0;
    end
end

wth_z = wth_z ./num_files;
sigmaw_z = sigmaw_z./num_files;
sigmath_z = sigmath_z./num_files;
sigmau_z = sigmau_z./num_files;
sigmav_z = sigmav_z./num_files;

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
plot(theta_z_waterfall',repmat(z,1,size(theta_z_waterfall,1)),'b','linewidth',3)
set(gca,'xtick',0:10:10*(size(theta_z_waterfall,1)-1));
set(gca,'xtickLabel',waterfall_label);
title('Waterfall plot for \theta(z) over time, spaced 10K apart.');
grid on;
label = 'waterfall';
filename = fullfile(outputpath,[fileout,label,'.png']);
hgexport(h, filename, hgexport('factorystyle'), 'Format', 'png');

%%
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

plot(time,tke(:,1),'bo-','linewidth',2)
plot(time,tke(:,2),'go-','linewidth',2)
plot(time,tke(:,3),'ro-','linewidth',2)
legend('subgrid','resolved','total');
xlabel('time (min)');
ylabel('Average TKE (m^2/s^2)');
title(['Average TKE vs time ',titlepost])
grid on;
label = 'tke';
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

plot(time,sigmath_t,'o-','linewidth',2,'color',colors{res})
xlabel('time (min)');
ylabel('\langle \theta''^2 \rangle (K^2)');
title(['Resolved \langle \theta''^2 \rangle vs time ',titlepost])
grid on;
label = 'ttt';
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

plot(time,sigmaw_t,'ro-','linewidth',2)
plot(time,sigmau_t,'go-','linewidth',2)
plot(time,sigmav_t,'bo-','linewidth',2)
xlabel('time (min)');
ylabel('\langle w''^2 \rangle (m^2/s^2)');
title(['Resolved \langle w''^2 \rangle vs time ',titlepost])
grid on;
legend('w^2','u^2','v^2');
label = 'wwt';
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

plot(time,wth_t,'o-','linewidth',2,'color',colors{res})
xlabel('time (min)');
ylabel('\langle w''\theta'' \rangle (mK/s)');
title(['Resolved \langle w''\theta'' \rangle vs time ',titlepost])
grid on;
label = 'wtt';
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

plot(wth_z,z,'o-','linewidth',2,'color',colors{res})
ylabel('z (km)');
xlabel('\langle w''\theta'' \rangle (mK/s)');
title(['Resolved \langle w''\theta'' \rangle_t vs height ',titlepost])
grid on;
label = 'wtz';
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

plot(sigmath_z,z,'o-','linewidth',2,'color',colors{res})
ylabel('z (km)');
xlabel('\langle \theta''^2 \rangle (K^2)');
title(['Resolved \langle \theta''^2 \rangle_t vs height ',titlepost])
grid on;
label = 'ttz';
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

plot(sigmaw_z,z,'ro-','linewidth',2)
plot(sigmau_z,z,'go-','linewidth',2)
plot(sigmav_z,z,'bo-','linewidth',2)
ylabel('z (km)');
xlabel('\langle w''^2 \rangle (m^2/s^2)');
title(['Resolved \langle w''^2 \rangle_t vs height ',titlepost])
grid on;
label = 'wwz';
legend('w^2','u^2','v^2');
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

plot(time,pbl,'o-','linewidth',2,'color',colors{res})
xlabel('time (min)');
ylabel('PBL Height (km)');
title(['PBL Height vs time ',titlepost])
grid on;

label = 'pbl';
filename = fullfile(outputpath,[fileout,label,'.png']);
hgexport(h, filename, hgexport('factorystyle'), 'Format', 'png');



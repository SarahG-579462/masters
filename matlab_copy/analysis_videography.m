close all

if ~exist('datapath','var')
    datapath = 'example_data';
    fileout = 'example_data.avi';
    outputpath = './';
elseif ~exist('fileout','var') && exist('postfix','var');
    [~,dp_header,~] = fileparts(datapath);
    fileout = [dp_header,postfix,'.avi'];
else
    [~,dp_header,~] = fileparts(datapath);
    fileout = [dp_header,'.avi'];
end
k = 1:numframes;
time = zeros(1,max(k));
num_files = 0;

F = figure('visible','off');

set(F,'Resize','on');
set(F,'PaperPosition',[0 0 16 9]);
set(F,'PaperUnits','inches');
set(F,'PaperSize',[16 9]); % IEEE columnwidth = 9cm
set(F,'Units','inches')
set(F,'Position',[0 0 16 9]);
set(F,'nextplot','replacechildren')

up_le =    subplot(1,2,1);
%up_ri =    subplot(2,2,2);set(up_ri,'FontSize',13)
dw_le =    subplot(1,2,2);
%dw_ri =    subplot(2,2,4);set(dw_ri,'FontSize',13)

%writerObj = VideoWriter(fullfile(outputpath,fileout),'Motion JPEG AVI');
%set(writerObj,'FrameRate',5);
%open(writerObj)
min_bv = -0.05;
max_bv =  0.05;
min_bw = -0.05;
max_bw =  0.05
M = 250;
for it=1:min(10,length(k));
    filename = fullfile(datapath,[dataprefix,num2str(k(it),'%06d'),'.nc']);
    if ~exist(filename,'file')
        continue
    end
    num_files = num_files + 1;
    framename = fullfile(outputpath,['frame_',num2str(k(it),'%02d'),'.png']);
    time(k(it)) = double(ncread(filename,'time'));
    
    u = double(ncread(filename,'uinterp'));
    v = double(ncread(filename,'vinterp'));
    w = double(ncread(filename,'winterp'));
    th= double(ncread(filename,'th')); % potential temperature

    % boundary layer height:
    % hpbl = double(ncread(filename,'hpbl'));

    % grid size parameters:
    xh = double(ncread(filename,'xh'))*1000;
    idx = median(diff(xh));
    xh_mid = xh - max(xh)/2;
    yh = double(ncread(filename,'yh'))*1000;
    idy = median(diff(xh));
    yh_mid = yh - max(yh)/2;
    z = double(ncread(filename,'z'))*1000;  % for nodes
    zf = double(ncread(filename,'zf'))*1000; % for edges
    idz = median(diff(z));
    
    [~,iz] = min(abs(z - 1000));
    [~,iy] = min(abs(yh - max(yh/2)));
    dz  = z((2:end-1)-1)-z((2:end-1)+1);
    dz = reshape(dz,[1,1,numel(dz)]);
    dz = repmat(dz,size(u,1),size(u,2),1);

    dx = xh((2:end-1)-1)-xh((2:end-1)+1);
    dx = reshape(dx,[numel(dx),1,1]);
    dx = repmat(dx,1,size(w,2),size(w,3));

    dy = yh((2:end-1)-1)-yh((2:end-1)+1);
    dy = reshape(dy,[1,numel(dy),1]);
    dy = repmat(dy,size(w,1),1,size(w,3));
    
    dzu = zeros(size(u));
    dzu(:,:,2:(end-1)) = 0.5*(u(:,:,(2:(end-1))-1)-u(:,:,(2:end-1)+1))./dz;
    dxw = zeros(size(w));
    dxw(2:(end-1),:,:) = 0.5*(w((2:(end-1))-1,:,:)-w((2:end-1)+1,:,:))./dx;
    dxv = zeros(size(v));
    dxv(2:(end-1),:,:) = 0.5*(v((2:(end-1))-1,:,:)-v((2:end-1)+1,:,:))./dx; 
    dyu = zeros(size(u));
    dyu(:,2:(end-1),:) = 0.5*(u(:,(2:(end-1))-1,:)-u(:,(2:end-1)+1,:))./dy;
    eta = dzu - dxw;
    nu  = dxv - dyu;
    min_bv = min([min_bv,min(eta(:)),min(nu(:))]);
    max_bv = max([max_bv,max(eta(:)),max(nu(:))]);
   
end

min_bv = -max(abs(min_bv),abs(max_bv));
max_bv =  max(abs(min_bv),abs(max_bv));
negCMN = ceil((abs(min_bv)/(abs(min_bv) + abs(max_bv))*M));
negCM = [linspace(0,1,negCMN)',linspace(0,1,negCMN)',linspace(1,1,negCMN)'];
posCMN = floor((abs(max_bv)/(abs(min_bv) + abs(max_bv))*M));
posCM = [linspace(1,1,posCMN)',linspace(1,0,posCMN)',linspace(1,0,posCMN)'];
colorMap = [negCM;posCM];
colormap(colorMap);
caxis([min_bv,max_bv])

for it=1:length(k);
    filename = fullfile(datapath,[dataprefix,num2str(k(it),'%06d'),'.nc']);
    if ~exist(filename,'file')
        continue
    end
    num_files = num_files + 1;
    framename = fullfile(outputpath,['frame_',num2str(k(it),'%02d'),'.png']);
    time(k(it)) = double(ncread(filename,'time'));
    
    u = double(ncread(filename,'uinterp'));
    v = double(ncread(filename,'vinterp'));
    w = double(ncread(filename,'winterp'));
    th= double(ncread(filename,'th')); % potential temperature

    % boundary layer height:
    % hpbl = double(ncread(filename,'hpbl'));

    % grid size parameters:
    xh = double(ncread(filename,'xh'))*1000;
    idx = median(diff(xh));
    xh_mid = xh - max(xh)/2;
    yh = double(ncread(filename,'yh'))*1000;
    idy = median(diff(xh));
    yh_mid = yh - max(yh)/2;
    z = double(ncread(filename,'z'))*1000;  % for nodes
    zf = double(ncread(filename,'zf'))*1000; % for edges
    idz = median(diff(z));
    
    [~,iz] = min(abs(z - 1000));
    [~,iy] = min(abs(yh - max(yh/2)));
    dz  = z((2:end-1)-1)-z((2:end-1)+1);
    dz = reshape(dz,[1,1,numel(dz)]);
    dz = repmat(dz,size(u,1),size(u,2),1);

    dx = xh((2:end-1)-1)-xh((2:end-1)+1);
    dx = reshape(dx,[numel(dx),1,1]);
    dx = repmat(dx,1,size(w,2),size(w,3));

    dy = yh((2:end-1)-1)-yh((2:end-1)+1);
    dy = reshape(dy,[1,numel(dy),1]);
    dy = repmat(dy,size(w,1),1,size(w,3));
    
    dzu = zeros(size(u));
    dzu(:,:,2:(end-1)) = 0.5*(u(:,:,(2:(end-1))-1)-u(:,:,(2:end-1)+1))./dz;
    dxw = zeros(size(w));
    dxw(2:(end-1),:,:) = 0.5*(w((2:(end-1))-1,:,:)-w((2:end-1)+1,:,:))./dx;
    dxv = zeros(size(v));
    dxv(2:(end-1),:,:) = 0.5*(v((2:(end-1))-1,:,:)-v((2:end-1)+1,:,:))./dx; 
    dyu = zeros(size(u));
    dyu(:,2:(end-1),:) = 0.5*(u(:,(2:(end-1))-1,:)-u(:,(2:end-1)+1,:))./dy;
    eta = dzu - dxw;
    nu  = dxv - dyu;
    plotclouds = false;
    nc = ncinfo(filename);
    if ismember('qc',{nc.Variables.Name})
	    plotclouds = true;
	    qc = double(ncread(filename,'qc'));
    end
    
    
    set(F,'CurrentAxes',up_le)
    imagesc(xh/1000,z/1000,squeeze(eta(:,iy,:,:))',[min_bv,max_bv])
    xlabel('x (km)');
    ylabel('z (km)');
    set(up_le,'ydir','normal')
    set(up_le,'xdir','normal')
    title('y-vorticity: x-z plane at y=L/2')
    colormap(colorMap);
    caxis([min_bv,max_bv]);
    set(up_le,'DataAspectRatio',[idx,idz,1]);
    if plotclouds
    hold on;
    if sum(sum(squeeze(qc(:,iy,:,:))')) > 1e-5
      contour(xh/1000,z/1000,squeeze(qc(:,iy,:,:))',[1,1]*1e-5,'color','k');
    end
    hold off; 
    end
    %set(F,'CurrentAxes',up_ri)
    %imagesc(xh/1000,z/1000,squeeze(w(:,iy,:,:))',[-2,2]);
    %xlabel('x (km)');
    %ylabel('z (km)');
    %set(up_ri,'ydir','normal')
    %set(up_ri,'xdir','normal')
    %title('w: x-z plane at y=L/2')
    %colormap(colorMap);
    %caxis([min_bv,max_bv]);   
    
    set(F,'CurrentAxes',dw_le)
    imagesc(xh/1000,yh/1000,squeeze(nu(:,:,iz,:))',[min_bw,max_bw]);
    set(dw_le,'ydir','normal')
    set(dw_le,'xdir','normal')
    xlabel('x (km)');
    ylabel('y (km)');
    title('z-vorticity: x-y plane at z=1km')
    colormap(colorMap);
    caxis([min_bv,max_bv])   
    set(dw_le,'DataAspectRatio',[idx,idy,1]);
    if plotclouds
    hold on;
    if sum(sum(squeeze(qc(:,:,iz,:))')) > 1e-5
    contour(xh/1000,yh/1000,squeeze(qc(:,:,iz,:))',[1,1]*1e-5,'color','k');
    end
    hold off;
    end
    %set(F,'CurrentAxes',dw_ri)
    %imagesc(xh/1000,yh/1000,squeeze(w(:,:,iz,:))',[-2,2])
    %xlabel('x (km)');
    %ylabel('y (km)');
    %set(dw_ri,'ydir','normal')
    %set(dw_ri,'xdir','normal')
    %title('w: x-y plane at z=pbl/2')
    %colormap(colorMap);
    %caxis([min_bv,max_bv])    
    hp4 = get(dw_le,'position');
    colorbar('Position', [hp4(1)+hp4(3)+0.01  hp4(2)  0.05  hp4(2)+hp4(4)*0.9],'AxisLocation','in')
    T = text(0,0,['time = ',sprintf('%i',floor(time(k(it))/60)),' min']);
    set(T,'units','normalized')
    Tp = get(T,'position');
    set(T,'FontSize',14)
    set(T,'position',[0,1.1])

    %writeVideo(writerObj, im2frame(hardcopy(F,'-dzbuffer','-r0')));
    hgexport(F, [framename], hgexport('factorystyle'), 'Format', 'png');
    
end
%close(writerObj)

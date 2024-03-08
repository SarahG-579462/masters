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

%up_ri =    subplot(2,2,2);set(up_ri,'FontSize',13)
up_le =    subplot(2,2,1);
up_ri =    subplot(2,2,2);
dw_le =    subplot(2,2,3);
dw_ri =    subplot(2,2,4);

%dw_ri =    subplot(2,2,4);set(dw_ri,'FontSize',13)

%writerObj = VideoWriter(fullfile(outputpath,fileout),'Motion JPEG AVI');
%set(writerObj,'FrameRate',5);
%open(writerObj)
min_bv = -0.00;
max_bv =  0;
min_bw = -0.00;
max_bw =  0;
M = 250;
plotclouds = false;
for it=1:min(10,length(k));
    filename = fullfile(datapath,[dataprefix,num2str(k(it),'%06d'),'.nc']);
    if ~exist(filename,'file')
        continue
    end

    nc = ncinfo(filename);
    if ismember('qc',{nc.Variables.Name})
	    plotclouds = true;
	    qc = double(ncread(filename,'qc'));
    else
	    return
    end
    
    num_files = num_files + 1;
    framename = fullfile(outputpath,['cloud_frame_',num2str(k(it),'%02d'),'.png']);
    time(k(it)) = double(ncread(filename,'time'));
    
    qc = double(ncread(filename,'qc'));
    max_bv = max(max(qc(:)),max_bw);
end
negCMN = ceil((abs(min_bv)/(abs(min_bv) + abs(max_bv))*M));
negCM = [linspace(0,1,negCMN)',linspace(0,1,negCMN)',linspace(1,1,negCMN)'];
posCMN = floor((abs(max_bv)/(abs(min_bv) + abs(max_bv))*M));
posCM = [linspace(1,1,posCMN)',linspace(1,0,posCMN)',linspace(1,0,posCMN)'];
colorMap = [posCM];
colormap(colorMap);
caxis([min_bv,max_bv])
for it=1:length(k);
    filename = fullfile(datapath,[dataprefix,num2str(k(it),'%06d'),'.nc']);
    if ~exist(filename,'file')
        continue
    end
    num_files = num_files + 1;
    framename = fullfile(outputpath,['cloud_frame_',num2str(k(it),'%02d'),'.png']);
    time(k(it)) = double(ncread(filename,'time'));
    
    qc = double(ncread(filename,'qc'));

 
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
    [~,ix] = min(abs(xh - max(xh/2)));
    [~,iy] = min(abs(yh - max(yh/2)));
    set(F,'CurrentAxes',up_le)
    imagesc(xh/1000,yh/1000,squeeze(qc(:,:,iz,:))',[min_bv,max_bv])
    set(up_le,'positionconstraint','outerposition')
    ylabel('y (km)');
    set(up_le,'ydir','normal')
    set(up_le,'xdir','normal')
    set(up_le,'units','centimeters');
    set(up_le,'xaxislocation','top')
    title('Cloud Water, kg/kg: x-y plane at z = 1km')
    colormap(colorMap);
    caxis([min_bv,max_bv]);
    set(up_le,'DataAspectRatio',[idx,idy,1]);
    pos = get(up_le,'innerposition');
    if sum(sum(squeeze(qc(:,:,iz,:))')) > 1e-5
    hold on
    contour(up_le,xh/1000,yh/1000,squeeze(qc(:,:,iz,:))',[1,1]*1e-5,'color','k');
    hold off
    end
    %
    set(F,'CurrentAxes',dw_le)
    imagesc(xh/1000,z/1000,squeeze(qc(:,iy,:,:))',[min_bv,max_bv])
    set(dw_le,'positionconstraint','outerposition')
    xlabel('x (km)');
    ylabel('z (km)');
    set(dw_le,'units','centimeters');
    set(dw_le,'ydir','normal')
    set(dw_le,'xdir','normal')
    set(dw_le,'DataAspectRatio',[idx,idz,1]);
    pos_dw_le = get(dw_le,'outerposition');
    set(dw_le,'innerposition',[pos(1) pos_dw_le pos(3) pos(3)*idz/idx]);
    
    if sum(sum(squeeze(qc(:,iy,:,:))')) > 1e-5
    set(dw_le,'DataAspectRatioMode','manual');
    set(dw_le,'PlotBoxAspectRatioMode','manual');
    set(dw_le,'CameraViewAngleMode','manual');
    hold on
    contour(dw_le,xh/1000,z/1000,squeeze(qc(:,iy,:,:))',[1,1]*1e-5,'color','k');
    hold off
    end
    colormap(colorMap);
    caxis([min_bv,max_bv]);
    
    set(F,'CurrentAxes',up_ri)
    imagesc(z/1000,yh/1000,squeeze(qc(ix,:,:,:)),[min_bv,max_bv])
    xlabel('z (km)');
    set(up_ri,'positionconstraint','outerposition')
    set(up_ri,'ydir','normal')
    set(up_ri,'xdir','reverse')
    set(up_ri,'units','centimeters')
    set(up_ri,'yaxislocation','right');
    pos_up_ri = get(up_ri,'outerposition');
    set(up_ri,'DataAspectRatio',[idz,idy,1]);
    set(up_ri,'innerposition',[pos_up_ri(1), pos(2), pos(4)*idy/idz, pos(4)]);    
    set(up_ri,'DataAspectRatioMode','manual');
    set(up_ri,'PlotBoxAspectRatioMode','manual');
    set(up_ri,'CameraViewAngleMode','manual');
    colormap(colorMap);
    caxis([min_bv,max_bv]);
    if sum(sum(squeeze(qc(ix,:,:,:))')) > 1e-5
    hold on
    contour(z/1000,yh/1000,squeeze(qc(ix,:,:,:)),[1,1]*1e-5,'color','k');
    hold off
    end
    set(dw_ri,'visible','off')

    %set(F,'CurrentAxes',up_ri)
    %imagesc(xh/1000,z/1000,squeeze(w(:,iy,:,:))',[-2,2]);
    %xlabel('x (km)');
    %ylabel('z (km)');
    %set(up_ri,'ydir','normal')
    %set(up_ri,'xdir','normal')
    %title('w: x-z plane at y=L/2')
    %colormap(colorMap);
    %caxis([min_bv,max_bv]);   
    
    %set(F,'CurrentAxes',dw_ri)
    %imagesc(xh/1000,yh/1000,squeeze(w(:,:,iz,:))',[-2,2])
    %xlabel('x (km)');
    %ylabel('y (km)');
    %set(dw_ri,'ydir','normal')
    %set(dw_ri,'xdir','normal')
    %title('w: x-y plane at z=pbl/2')
    %colormap(colorMap);
    %caxis([min_bv,max_bv])    
    colorbar('Position', [pos_up_ri(1)+pos_up_ri(3)+0.01,pos_dw_le(2)  0.05  pos_dw_le(4)+pos_up_ri(4)],'AxisLocation','in')
    T = text(0,0,['time = ',sprintf('%i',floor(time(k(it))/60)),' min']);
    set(T,'units','normalized')
    Tp = get(T,'position');
    set(T,'FontSize',14)
    set(T,'position',[0,1.1])

    set(dw_le,'DataAspectRatio',[idx,idy,1]);
    %writeVideo(writerObj, im2frame(hardcopy(F,'-dzbuffer','-r0')));
    hgexport(F, [framename], hgexport('factorystyle'), 'Format', 'png');
    
end
%close(writerObj)

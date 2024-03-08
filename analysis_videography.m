
clear all
close all

if ~exist('datapath','var')
    datapath = 'example_data';
    fileout = 'example_data.mp4';
    outputpath = './';
    numframes = 37;
elseif ~exist('fileout','var') && exist('postfix','var');
    [~,dp_header,~] = fileparts(datapath);
    fileout = [dp_header,postfix,'.mp4'];
else
    [~,dp_header,~] = fileparts(datapath);
    fileout = [dp_header,'.mp4'];
end
k = 1:numframes;
time = zeros(1,max(k));
num_files = 0;

F = figure('visible','on');

set(F,'Resize','on');
set(F,'PaperPosition',[0 0 9 9]);
set(F,'PaperUnits','inches');
set(F,'PaperSize',[9 9]); % IEEE columnwidth = 9cm
set(F,'Units','inches')
set(F,'Position',[0 0 9 9]);
set(F,'nextplot','replacechildren')

up_le =    subplot(2,1,1);set(up_le,'FontName','Arial')
%up_ri =    subplot(2,2,2);set(up_ri,'FontSize',13)
dw_le =    subplot(2,1,2);set(dw_le,'FontName','Arial')
%dw_ri =    subplot(2,2,4);set(dw_ri,'FontSize',13)

writerObj = VideoWriter(fullfile(outputpath,fileout),'MPEG-4');
set(writerObj,'FrameRate',5);
open(writerObj)
min_bv = -2;
max_bv =  2;
M = 250;
negCMN = ceil((abs(min_bv)/(abs(min_bv) + abs(max_bv))*M));
negCM = [linspace(0,1,negCMN)',linspace(0,1,negCMN)',linspace(1,1,negCMN)'];
posCMN = floor((abs(max_bv)/(abs(min_bv) + abs(max_bv))*M));
posCM = [linspace(1,1,posCMN)',linspace(1,0,posCMN)',linspace(1,0,posCMN)'];
colorMap = [negCM;posCM];
colormap(colorMap);
caxis([min_bv,max_bv])
for it=1:length(k);
    filename = fullfile(datapath,['cm1out_',num2str(k(it),'%06d'),'.nc']);
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
    hpbl = double(ncread(filename,'hpbl'));

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
    
    [~,iz] = min(abs(z - mean(hpbl(:))/2));
    [~,iy] = min(abs(yh - max(yh/2)));
    
    set(F,'CurrentAxes',up_le)
    imagesc(xh/1000,z/1000,squeeze(u(:,iy,:,:))',[-2,2])
    xlabel('x (km)');
    ylabel('z (km)');
    set(up_le,'ydir','normal')
    set(up_le,'xdir','normal')
    title('u: x-z plane at y=L/2')
    colormap(colorMap);
    caxis([min_bv,max_bv]);

%     set(F,'CurrentAxes',up_ri)
%     imagesc(xh/1000,z/1000,squeeze(w(:,iy,:,:))',[-2,2]);
%     xlabel('x (km)');
%     ylabel('z (km)');
%     set(up_ri,'ydir','normal')
%     set(up_ri,'xdir','normal')
%     title('w: x-z plane at y=L/2')
%     colormap(colorMap);
%     caxis([min_bv,max_bv]);   
    
    set(F,'CurrentAxes',dw_le)
    imagesc(xh/1000,yh/1000,squeeze(u(:,:,iz,:))',[-2,2]);
    set(dw_le,'ydir','normal')
    set(dw_le,'xdir','normal')
    xlabel('x (km)');
    ylabel('y (km)');
    title('u: x-y plane at z=pbl/2')
    colormap(colorMap);
    caxis([min_bv,max_bv])   
    
%     set(F,'CurrentAxes',dw_ri)
%     imagesc(xh/1000,yh/1000,squeeze(w(:,:,iz,:))',[-2,2])
%     xlabel('x (km)');
%     ylabel('y (km)');
%     set(dw_ri,'ydir','normal')
%     set(dw_ri,'xdir','normal')
%     title('w: x-y plane at z=pbl/2')
%     colormap(colorMap);
%     caxis([min_bv,max_bv])    
    hp4 = get(dw_le,'position');
    colorbar('Position', [hp4(1)+hp4(3)+0.01  hp4(2)  0.05  hp4(2)+hp4(3)*0.9])
    set(F,'CurrentAxes',up_le)
    T = text(0,0,['time = ',sprintf('%i',floor(time(k(it))/60)),' min']);
    set(T,'units','normalized')
    Tp = get(T,'position');
    set(T,'FontSize',14)
    set(T,'position',[0,1.1])
    writeVideo(writerObj, im2frame(hardcopy(F,'-dzbuffer','-r0')));
    break
end
close(writerObj)

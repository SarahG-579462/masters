%%
if ~exist('datapath','var')
    datapath = 'example_data';
    fileout = 'example_data_pbl.png';
    outputpath = './';
elseif exist('postfix','var');
    [~,dp_header,~] = fileparts(datapath);
    fileout = [dp_header,'_pbl',postfix,'.png'];
else
    [~,dp_header,~] = fileparts(datapath);
    fileout = [dp_header,'_pbl.png'];
end

fileoutfull = fullfile(outputpath,fileout);
numframes = 37;
time = nan(1,numframes);
pbl = nan(1,numframes);
num_files = 0;
for it = 1:numframes;
filename = fullfile(datapath,['cm1out_',num2str(it,'%06d'),'.nc']);
if ~exist(filename,'file')
    continue
end
num_files = num_files + 1;
% load data:
time(it) = double(ncread(filename,'time'));
zf = double(ncread(filename,'zf')); % for edges
brunt_vaisalla = double(ncread(filename,'nm'));
% get max vaisalla & average:
[~,iz] = max(brunt_vaisalla,[],3);
pbl(it) = mean(zf(iz(:)));

end
%%
close all;
h=figure('visible','off');
set(h,'Resize','on');
%set(h,'PaperPositionMode','auto');
set(h,'PaperPosition',[0 0 16 9]);
set(h,'PaperUnits','inches');
set(h,'PaperSize',[16 9]); % IEEE columnwidth = 9cm
set(h,'Units','inches')
set(h,'Position',[0 0 16 9]);
plot(time,pbl,'o--')
set(gca,'fontsize',13)
grid on;
xlabel('time (s)');
ylabel('PBL height (km)');
title(['averaged PBL height vs time',titlepost],'fontsize',15)

hgexport(h, fileoutfull, hgexport('factorystyle'), 'Format', 'png');

it = 1;
filename = ['C:\Users\Sarah\Documents\MATLAB\mesoscale_project\matlab_copy\Modularized\Plotters\moistwind\d_off\out\','cm1out_',num2str(it,'%06d'),'.nc'];

th = squeeze(double(ncread(filename,'th')));
hpbl = squeeze(double(ncread(filename,'hpbl')));
prs = squeeze(double(ncread(filename,'prs')));

xh = double(ncread(filename,'xh'));
yh = double(ncread(filename,'yh'));
z = double(ncread(filename,'z'));  % for nodes
zf = double(ncread(filename,'zf')); % for edges

th_z = zeros(size(z));
t_z = zeros(size(z));
for iz = 1:length(z);
    t = th(:,:,iz,:);
    th_z(iz) = mean(t(:));
    T = th(:,:,iz,:).*(prs(:,:,iz,:)./100000).^(0.286);
    t_z(iz) = mean(T(:));
end
figure(1);clf(1);hold on;

plot(th_z,z,'.-');
xlim([295,320])
grid on;
xlabel('Potential Temperature (K)','fontsize',11)
ylabel('z (km)','fontsize',11)
title('Potential Temperature Initial Profile','fontsize',13)


figure(2);clf(2);hold on;
plot(t_z,z,'.-');
%xlim([295,320])
grid on;
xlabel('Temperature (K)','fontsize',11)
ylabel('z (km)','fontsize',11)
title('Temperature Initial Profile','fontsize',13)
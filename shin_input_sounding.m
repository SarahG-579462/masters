z = linspace(0,3600,3601);
th = 300.*(z <= 925) + (300 + (z - 925).*0.0536).*(z > 925).*(z <= 1075) + (308.05 + (z - 1075).*0.003).*(z > 1075);
th_def = 300.*(z <= 974) + (300 + (z - 974).*0.08).*(z > 974).*(z <= 1074) + (308.05 + (z - 1074).*0.003).*(z > 1074);
qv = zeros(size(z));
u = zeros(size(z));
v = zeros(size(z));

for i = 1:length(z);
    fprintf('%i %.5f %.5f %.5f %.5f\n',z(i),th(i),qv(i),u(i),v(i));
end
%%
close all;
figure(1);clf(1);hold on;
plot(th_def,z,'r','linewidth',3);
plot(th,z,'b','linewidth',3);

xlim([295,320])
grid on;
xlabel('Potential Temperature (K)','fontsize',11)
ylabel('z (km)','fontsize',11)
legend('Sullivan, Patton, 2011','Shin, Hong, 2013')
title('Potential Temperature Initial Profile','fontsize',13)
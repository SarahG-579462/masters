function out = config_finder();

% config_finder
cores = 40;
dx = 12.5;
dy = 12.5;
t_max = 48;

Nz = 240;
dt = 0.75/2;
t_min = 1;
max_nodes = 9;

%% find possible domain sizes:
Nx = [];
Ny = [];
domain_minx = 6000;
domain_maxx = 7000;
domain_miny = 6000;
domain_maxy = 7000;
domain_sq = true; %square domain?
kx = 1; ky = 1;

for i = 1:100
  if (dx*i*cores >= domain_minx) && (dx*i*cores <= domain_maxx)
    Nx(kx) = cores*i;
    kx = kx+1;
  end
  
  if (dy*i*cores >= domain_miny) && (dy*i*cores <= domain_maxy)
    Ny(ky) = cores*i;
    ky = ky+1;
  end
end
if domain_sq
  Ny = Nx;
end
%%
n = 1;
max_rows = 1000;
out = table(zeros(max_rows,1),zeros(max_rows,1),zeros(max_rows,1),zeros(max_rows,1),zeros(max_rows,1),zeros(max_rows,1),zeros(max_rows,1),zeros(max_rows,1),zeros(max_rows,1),...
  'variableNames',{'Nx','Ny','nodex','nodey','difference','est_time','num_nodes','num_cores','mem_per_file'});

for kx = 1:length(Nx)
  factors_x = get_divisors(Nx(kx));
  for ky = 1:length(Ny);
    if domain_sq
      factors_y = factors_x;
      ky = kx;
    else
      factors_y = get_divisors(Ny(ky));
    end
    mem = 160/(60*128*128)*(Nz*Nx(kx)*Ny(ky)); %estimate from 50m simulation
    for i=1:length(factors_x)
      for j = 1:length(factors_y)
        px = factors_x(i);
        py = factors_y(j);
        time = Nx(kx)*Ny(ky)*Nz*21600/dt/(4*10^8) / (px*py);% estimate from Bryan for core-hours.
        if (mod(px*py,cores) == 0) && (time <= t_max) && (time >= t_min) && (px*py/cores <= max_nodes)
          out(n,'Nx') = {Nx(kx)};
          out(n,'Ny') = {Ny(ky)};
          out(n,'nodex') = {px};
          out(n,'nodey') = {py};
          out(n,'difference') = {abs(px-py)};
          out(n,'est_time') = {time};
          out(n,'num_nodes') = {px*py/cores};
          out(n,'num_cores') = {px*py};
          out(n,'mem_per_file') = {mem};
          
          %out(n,:) = {Nx(kx),Ny(ky),px,py,abs(px-py),time,px*py/cores,px*py,mem};
          n = n+1;

        end
      end
    end
    if domain_sq
      break
    end
  end
end
out = out((1:n-1),:);
out = sortrows(out,'difference');
out = sortrows(out,'est_time');
end

function div = get_divisors(N);
% quick and dirty, use divisors function if it exists...
if exist('divisors','builtin')
  div = div(N);
else
  px = factor(N);
  [ps] = unique(px,'last');
  counts = [];
  for ips=1:length(ps);
    counts(ips) = sum(px == ps(ips));
  end
  
  div = [];
  for ips = 1:length(ps)
    count = counts(ips);
    sm_factors = ps(ips).^(0:count);
    if ~isempty(div)
      prod_factors = div'*sm_factors;
    else
      prod_factors = [];
    end
    div = sort(unique([sm_factors,div,prod_factors(:)']));
  end
end
end
    
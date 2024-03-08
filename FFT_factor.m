function F = FFT_factor(scale,N,L);

switch scale
    case 'density'
        F = prod(L)/prod(N.^2)./(2*pi).^numel(N);
    case 'total'
        F = prod(L.^2)/prod(N.^2)./(2*pi).^numel(N);
    case 'none'
        F = 1;
    case 'discrete'
        % density, but divided by dx
        dx = prod(L)/prod(N);
        F = prod(L)/prod(N.^2)./(2*pi).^numel(N)/dx;
end
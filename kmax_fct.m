function [kmaxhet1,kmaxhet2,dhet1,dhet2] = kmax_fct(location)
%KMAX_FCT obtaines the growth and diffusion parameters based on the spatial landscape.

% Author: Lukas Eigentler (leigentler001@dundee.ac.uk)
% License: GNU GPL
% Last updated: 19/10/2021

global r dem 
X = linspace(-r(end),r(end),length(dem(:,1)));
L = interp2(X,X,dem,location.x,location.y); % interpolate fractal landscape onto PDE mesh
[kmaxhet1,kmaxhet2,dhet1,dhet2] = kmax_fct2(L);
end


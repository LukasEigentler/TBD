function [kmaxhet1,kmaxhet2,dhet1,dhet2] = kmax_fct2(L)
%KMAX_FCT2 calculates the growth and diffusion parameters based on the supplied landscape.

% Author: Lukas Eigentler (leigentler001@dundee.ac.uk)
% License: GNU GPL
% Last updated: 19/10/2021

global kmax1 kmax2 db1 db2 param_setting

if param_setting == 2
    % Param setting 2
    kmaxhet1 = 1*L+4; kmaxhet1(kmaxhet1<0) =0;
    kmaxhet2 = 1*L+4; kmaxhet2(kmaxhet2<0) =0;
    dhet1 = kmaxhet1/50;
    dhet2 = kmaxhet2/50;
    
elseif param_setting == 0
    kmaxhet1 = kmax1*ones(1,length(L));
    kmaxhet2 = kmax2*ones(1,length(L));
    dhet1 = db1*ones(1,length(L));
    dhet2 = db2*ones(1,length(L));
    
else
    % Param setting 1 or 3
    kmaxhet1 = kmax1*ones(1,length(L));
    kmaxhet2 = kmax2*ones(1,length(L));
    dhet1 = db1*ones(1,length(L));
    dhet2 = db2*ones(1,length(L));
    kmaxhet1(L<0) = kmax1/2;
    kmaxhet2(L<0) = kmax2/2;
    dhet1(L<0) = db1/2;
    dhet2(L<0) = db1/2;
end

end
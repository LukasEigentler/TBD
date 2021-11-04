function f = fcoeffunction_het(location,state)
%FOCOEFFFUNCTION_HET provides the f matrix for PDEToolbox.

% Author: Lukas Eigentler (leigentler001@dundee.ac.uk)
% License: GNU GPL
% Last updated: 19/10/2021

global kcap c12 c21
N = 2; % Number of equations
nr = length(location.x); % Number of columns
f = zeros(N,nr); % Allocate f
[kmaxhet1,kmaxhet2,~,~] = kmax_fct(location); % get parameter landscape

if length(location.x)>4
    kk1 = state.u(1,:) <0;
    state.u(1,kk1) = 0;
    kk2 = state.u(2,:) <0;
    state.u(2,kk2) = 0;
    kk = find(state.u(1,:)+state.u(2,:)>kcap);
    state.u(1,kk) = kcap*state.u(1,kk)./(state.u(1,kk)+state.u(2,kk));
    state.u(2,kk) = kcap*state.u(2,kk)./(state.u(1,kk)+state.u(2,kk));
end
f(1,:) = kmaxhet1.*state.u(1,:).*(1-(state.u(1,:)+state.u(2,:))/kcap)-c12*state.u(1,:).*state.u(2,:);
f(2,:) = kmaxhet2.*state.u(2,:).*(1-(state.u(1,:)+state.u(2,:))/kcap)-c21*state.u(1,:).*state.u(2,:);


end
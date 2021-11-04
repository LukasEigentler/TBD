function cmat = ccoeffunction_het(location,state)
%COCOEFFFUNCTION_HET provides the c matrix for PDEToolbox.

% Author: Lukas Eigentler (leigentler001@dundee.ac.uk)
% License: GNU GPL
% Last updated: 19/10/2021
global kcap N time_start

if toc(time_start) > 10000 && length(location.x)>1
    return
end
[~,~,dhet1,dhet2] = kmax_fct(location); % get parameter landscape
cmat = zeros(4*N^2,length(location.x));

if length(location.x)>4
    kk1 = state.u(1,:) <0;
    state.u(1,kk1) = 0;
    kk2 = state.u(2,:) <0;
    state.u(2,kk2) = 0;
    kk = find(state.u(1,:)+state.u(2,:)>kcap);
    state.u(1,kk) = kcap*state.u(1,kk)./(state.u(1,kk)+state.u(2,kk));
    state.u(2,kk) = kcap*state.u(2,kk)./(state.u(1,kk)+state.u(2,kk));
end
cmat(1,:) =  dhet1.*(1-(state.u(1,:)+state.u(2,:))/kcap); 
cmat(4,:) =   dhet1.*(1-(state.u(1,:)+state.u(2,:))/kcap); 
cmat(4*N+5,:) =  dhet2.*(1-(state.u(1,:)+state.u(2,:))/kcap); 
cmat(4*N+8,:) =  dhet2.*(1-(state.u(1,:)+state.u(2,:))/kcap); 



end
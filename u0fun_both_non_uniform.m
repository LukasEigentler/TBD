function uinit = u0fun_both_non_uniform(location)
% U0FUN_BOTH_NON_UNIFORM defines the model's initial condition on the
% spatial mesh 
%
% Author: Lukas Eigentler (leigentler001@dundee.ac.uk)
% License: GNU GPL
% Last updated: 19/10/2021

global kcap Cxloc_in_B1 Cyloc_in_B1 Cxloc_in_B2 Cyloc_in_B2
M = length(location.x);
uinit = zeros(2,M);
uinit(1,:) = 0;
uinit(2,:) = 0;
for cc = 1:length(Cxloc_in_B1)
    d = sqrt((location.x-Cxloc_in_B1(cc)).^2 + (location.y-Cyloc_in_B1(cc)).^2);
    [~,mind_ind] = min(d);
    uinit(1,mind_ind) = kcap;
end
for cc = 1:length(Cxloc_in_B2)
    d = sqrt((location.x-Cxloc_in_B2(cc)).^2 + (location.y-Cyloc_in_B2(cc)).^2);
    [~,mind_ind] = min(d);
    uinit(2,mind_ind) = kcap;
end
if length(location.x)>4
    for mm = 1:M
        if uinit(1,mm)+uinit(2,mm)>kcap
            uinit(1,mm) = kcap/2;
            uinit(2,mm) = kcap/2;
        end
    end
end
end
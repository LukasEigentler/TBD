%% Model simulation

%MODEL_SIMULATION performs a Voronoi tessellation and implements the 
% mathematical model using a finite element method. Note that the 
% PDEToolbox is required.

% Author: Lukas Eigentler (leigentler001@dundee.ac.uk)
% License: GNU GPL
% Last updated: 19/10/2021

clear; 
close all;
% rng(789) % optionally set random seed for reproducibility
%% parameters
global db1 db2 kmax1 kmax2 kcap c12 c21 Cxloc_in_B1 Cyloc_in_B1 Cxloc_in_B2 Cyloc_in_B2 N time_start dem r param_setting
db1 = 0.1; db2 = 0.1; kcap = 1; kmax1 = 5; kmax2 = 5; c12 = 0; c21 = 0;
r = 10; % domain radius 
n=6; % number of iterations for generation of fractal surface
H=0.5; % Hurst exponent, fractal dimension D = 3-H
param_setting = 1; % 1: fractal like structure with zero contour as separator between good and bad patches ("case (iii)"); 2: fractal like structure with continuous parameter landscape ("case (iv)") 3: exact chequerboard ("case (i)"); 4: same grid as exact chequerboard but randomly assigned good and bad patches ("case (ii)") 
tlist = linspace(0,5,10); % vector defining time outputs to be saved
%% create parameter landscape and plot

f9 = figure(9);
Axes(1) = axes('Position', [0.05 0.05 0.8 0.8]);

if param_setting < 3
    [row,col,dem] = fractalsurface(n,H); % create fractal like surface
    x = 1:length(dem); xm = length(dem)/2;
    xplot = 2*r*(x-xm)/length(dem);
    if param_setting == 1
        contourf(xplot,xplot,dem,[0,0],'linestyle', 'none') % plot landscape using the zero contour as separator between regions
    else
        contourf(xplot,xplot,dem,'linestyle', 'none') % plot landscape with more contours
    end
elseif param_setting == 3
    xplot = linspace(-r,r,1000); % vector to create chequerboard pattern - distance between elements should be a few order of magnitude smaller than sidelength of chequerboard tiles
    dem = NaN*ones(length(xplot));
    for xx = 1:length(xplot) % create chequerboard pattern in landscape 
        for yy = 1:length(xplot)
            if mod(ceil(xplot(xx)),2) + mod(ceil(xplot(yy)),2) == 1
                dem(xx,yy) = -1;
            else
                dem(xx,yy) = 1;
            end
        end
    end
    contourf(xplot,xplot,dem,[0,0],'linestyle', 'none') % plot landscape using the zero contour as separator between regions
    
elseif param_setting == 4
    xplot = linspace(0,2*r,1000); % vector to create chequerboard pattern - distance between elements should be a few order of magnitude smaller than sidelength of chequerboard tiles
    dem = NaN*ones(length(xplot));
    randmat = rand(21); % random square matrix used to randomly select values in grid
    for xx = 1:length(xplot)
        for yy = 1:length(xplot)
            xres = mod(ceil(xplot(xx)),100);
            yres = mod(ceil(xplot(yy)),100);
            draw = randmat(xres+1,yres+1);
            if draw <0.5
                dem(xx,yy) = -1;
            else
                dem(xx,yy) = 1;
            end
        end
    end
    contourf(xplot,xplot,dem,[0,0],'linestyle', 'none') % plot landscape using the zero contour as separator between regions
    xplot = linspace(-r,r,1000);
end
shading interp
pbaspect([1 1 1])
hold on
th = 0:pi/50:2*pi;
xunit = r * cos(th);
yunit = r * sin(th);
h = plot(xunit, yunit, 'k');
temp_ind1 = find(xunit>0); temp_ind2 = find(yunit>0); q1_ind = intersect(temp_ind1,temp_ind2);
temp_ind1 = find(xunit<0); temp_ind2 = find(yunit>0); q2_ind = intersect(temp_ind1,temp_ind2);
temp_ind1 = find(xunit<0); temp_ind2 = find(yunit<0); q3_ind = intersect(temp_ind1,temp_ind2);
temp_ind1 = find(xunit>0); temp_ind2 = find(yunit<0); q4_ind = intersect(temp_ind1,temp_ind2);
fill([-r,r, r, sort(xunit(q1_ind),'descend'), sort(xunit(q2_ind),'descend')],...
     [r, r, 0,      sort(yunit(q1_ind),'ascend'),  sort(yunit(q2_ind),'descend')], 'w')
fill([-r,  r, r, sort(xunit(q4_ind),'descend'),  sort(xunit(q3_ind),'descend')],...
     [-r, -r, 0,      sort(yunit(q4_ind),'descend'),  sort(yunit(q3_ind),'ascend')], 'w')
colormap(Axes(1),gray)
set(gca, 'XTick', []);
set(gca, 'YTick', []);
set(f9,'Windowstyle','normal')
set(findall(f9,'-property','FontSize'),'FontSize',11)
set(f9,'Units','centimeters')
set(f9,'Position',[2 10 2*17/4 17/4])

%% initialise pde
N = 2; %no of equations
model = createpde(N);

%% geometry
NC_in = 6; % number of initial population patches
NC_in_B1 = 0.5; % ratio of those allocated to B1
NC_in_B1 = ceil(NC_in_B1*NC_in);
R1 = [1;0;0;r]; % whole domain (disk)
sf = 'R001';
ns = char('R001');
gd = R1;
ns = ns';
[dl,bt] = decsg(gd,sf,ns);
geometryFromEdges(model,dl);

%% boundary conditions
applyBoundaryCondition(model,'neumann','edge',1:4,'g',zeros(N),'q',zeros(N));

%% PDE coefficients
specifyCoefficients(model,'f',@fcoeffunction_het,'m',0,'d',1,'c',@ccoeffunction_het,'a',0);

%% mesh
mesh = generateMesh(model, 'hmax',0.16,'hmin',0.16, 'GeometricOrder', 'Linear'); % note: choice of mesh element size determines possible choice of number of initial population patches

%% create graph
[graph, nodes_loc] = mesh2graph(mesh); % create graph used to calculate front propagation metric

%% initial condition

point_range = 2; %radius of centre region in which initial population is placed
% find mesh nodes in centre region
nrm = zeros(1,length(model.Mesh.Nodes(1,:)));
for kk = 1:length(model.Mesh.Nodes(1,:))
    nrm(kk) = norm(model.Mesh.Nodes(:,kk)); 
end
ic_ind = find(nrm<point_range);
ic_nodes =  model.Mesh.Nodes(:,ic_ind);
% randomly select locations of initial population patches from nodes in centre region
Cxloc_in = zeros(1,NC_in); Cyloc_in = Cxloc_in;
for ff = 1:NC_in
    select_ind = randi(length(ic_nodes(1,:)));
    entry = ic_nodes(:,select_ind);
    Cxloc_in(ff) = entry(1);
    Cyloc_in(ff) = entry(2);
    ic_nodes(:,select_ind) = [];
end
Cxloc_in_B1 = Cxloc_in(1:NC_in_B1); Cyloc_in_B1 = Cyloc_in(1:NC_in_B1);
Cxloc_in_B2 = Cxloc_in(NC_in_B1+1:end); Cyloc_in_B2 = Cyloc_in(NC_in_B1+1:end);

u0both_non_uniform = @u0fun_both_non_uniform;
setInitialConditions(model,u0both_non_uniform);

%% Voronoi tessellations

temp = solvepde(model,[0,1e-10]); % use pde solver only to bring initial condition into an easily readable format
B10 = temp.NodalSolution(:,1,1);
B20 = temp.NodalSolution(:,2,1);
B1_ind = find(B10 ~=0); % find B1 locations
B2_ind = find(B20 ~=0); % find B2 locations
B_ind = [B1_ind;B2_ind];
tm = length(B_ind); %total mass
com = sum(temp.Mesh.Nodes(:,B_ind),2)/tm; %centre of mass

% find distances between mesh nodes
dist_B1 = zeros(length(B1_ind),length(graph(:,1)));
prev_B1 = zeros(length(B1_ind),length(graph(:,1)));
dist_B2 = zeros(length(B1_ind),length(graph(:,1)));
prev_B2 = zeros(length(B1_ind),length(graph(:,1)));
for bb = 1:length(B1_ind) % loop over all initial B1 patches
    [~,node_ind] = min(vecnorm(nodes_loc - mesh.Nodes(:,B1_ind(bb)))); % find closest node
    [dist_B1(bb,:),prev_B1(bb,:)] = Dijkstra(graph,node_ind); % calculate shortest path from initial patch to all other nodes
end
for bb = 1:length(B2_ind) % loop over all initial B2 patches
    [~,node_ind] = min(vecnorm(nodes_loc - mesh.Nodes(:,B2_ind(bb))));
    [dist_B2(bb,:),prev_B2(bb,:)] = Dijkstra(graph,node_ind); % calculate shortest path from initial patch to all other nodes
end

% Voronoi tessellations - regular cartesian grid
xplot_measure = linspace(-r(end),r(end),1000); % grid where distances are measured
[xplot_measure_meshgrid,yplot_measure_meshgrid] = meshgrid(xplot_measure,xplot_measure);
dist_B1_measure = zeros(length(B1_ind),length(xplot_measure),length(xplot_measure));
initial_theta_b1 = 0;
initial_theta_b1_loc = zeros(length(xplot_measure),length(xplot_measure));
initial_theta_b2 = 0;
initial_theta_b2_loc = zeros(length(xplot_measure),length(xplot_measure));
for bb = 1:length(B1_ind)
    F = scatteredInterpolant(nodes_loc(1,:)',nodes_loc(2,:)',dist_B1(bb,:)'); % interpolate to a meshgrid 
    dist_B1_measure(bb,:,:) = F(xplot_measure_meshgrid,yplot_measure_meshgrid);
end
if length(B1_ind)>1
    dist_B1_measure_min = min(dist_B1_measure); % find distance of closest B1 patch 
else
    dist_B1_measure_min = dist_B1_measure;
end
dist_B2_measure = zeros(length(B2_ind),length(xplot_measure),length(xplot_measure));
for bb = 1:length(B2_ind)
    F = scatteredInterpolant(nodes_loc(1,:)',nodes_loc(2,:)',dist_B2(bb,:)'); % interpolate to a meshgrid 
    dist_B2_measure(bb,:,:) = F(xplot_measure_meshgrid,yplot_measure_meshgrid);
end
if length(B2_ind)>1
    dist_B2_measure_min = min(dist_B2_measure); % find distance of closest B2 patch 
else
    dist_B2_measure_min = dist_B2_measure;
end

footprint_pred_vec = zeros(1,length(dist_B1_measure_min(:)));
for bb = 1:length(dist_B1_measure_min(:)) % loop through all measurement points
    [dist_min,min_ind] = min([dist_B1_measure_min(bb),dist_B2_measure_min(bb)]); % find which strain is closer
    if dist_min < tlist(end) % only count if point can be reached in simulation time
        footprint_pred_vec(bb) = 1; % vector storing if prediction indicates that point is covered by range expansion
        if min_ind == 1 % if B1 is closer, increase count by area of circle sector
            initial_theta_b1_loc(bb) = 1;
            initial_theta_b1 = initial_theta_b1 +   1;
        else
            initial_theta_b2_loc(bb) = 1;
            initial_theta_b2 = initial_theta_b2 +   1;
        end
    end
end
initial_theta_b1 = initial_theta_b1/(initial_theta_b1+initial_theta_b2); % normalise to unity
initial_theta_b2 = 1-initial_theta_b1;
area_el = (xplot_measure(2)-xplot_measure(1))^2;
footprint_pred = area_el*sum(footprint_pred_vec);
disp("The prediction of the area covered by range expansion is " + num2str(footprint_pred, '%.3f'))

% intraspecies connectedness - regular polar grid
rplot_measure = linspace(0,2*r,1000); % grid where distances are measured
theta_plot_measure = linspace(0,2*pi,1000); % grid where distances are measured
dtheta = theta_plot_measure(2)-theta_plot_measure(1);
[rplot_measure_meshgrid,thetaplot_measure_meshgrid] = meshgrid(rplot_measure,theta_plot_measure);
[x_mg,y_mg] = pol2cart(thetaplot_measure_meshgrid,rplot_measure_meshgrid); x_mg = x_mg+com(1); y_mg = y_mg+com(2); % convert to cartesian and shift origin to center of mass
dist_B1_measure = zeros(length(B1_ind),length(rplot_measure),length(theta_plot_measure));
initial_theta_b1_loc_polar = zeros(1,length(theta_plot_measure));
initial_theta_b2_loc_polar = zeros(1,length(theta_plot_measure));
for bb = 1:length(B1_ind)
    F = scatteredInterpolant(nodes_loc(1,:)',nodes_loc(2,:)',dist_B1(bb,:)'); % interpolate to a meshgrid 
    dist_B1_measure(bb,:,:) = F(x_mg,y_mg);
end
dist_B1_measure_min = min(dist_B1_measure); % find distance of closest B1 patch 

dist_B2_measure = zeros(length(B2_ind),length(rplot_measure),length(theta_plot_measure));
for bb = 1:length(B2_ind)
    F = scatteredInterpolant(nodes_loc(1,:)',nodes_loc(2,:)',dist_B2(bb,:)'); % interpolate to a meshgrid 
    dist_B2_measure(bb,:,:) = F(x_mg,y_mg);
end
dist_B2_measure_min = min(dist_B2_measure); % find distance of closest B2 patch 
dist_all_min = min(dist_B1_measure_min,dist_B2_measure_min);

rmax = zeros(1,length(theta_plot_measure));
for bb = 1:length(theta_plot_measure) % loop through all angles
    rmax_ind = find(dist_all_min(1,bb,:)>tlist(end)); rmax_ind = rmax_ind(1);
    rmax(bb) = rplot_measure(rmax_ind); 
    [dist_min,min_ind] = min([dist_B1_measure_min(1,bb,rmax_ind),dist_B2_measure_min(1,bb,rmax_ind)]); % find which strain is closer
    if min_ind == 1 % if B1 is closer, increase count by area of circle sector
        initial_theta_b1_loc_polar(bb) = 1;
    else
        initial_theta_b2_loc_polar(bb) = 1;
    end
end
initial_theta_b1_loc_polar_sign_change = length(find(diff(initial_theta_b1_loc_polar) ~= 0)); % find initial intraspecies connectedness
if ~isempty(initial_theta_b1_loc_polar_sign_change)
    M = initial_theta_b1_loc_polar_sign_change/2;
elseif initial_theta_b1_loc_polar(1) == 1
    M = 1;
else 
    M = 0;
end
if initial_theta_b1_loc_polar(1) ~= initial_theta_b1_loc_polar
    M = M+0.5;
end

%% plot Voronoi tessellation
f2 =  figure(12);
hold on
plot(xplot_measure_meshgrid(initial_theta_b1_loc(:)==1), yplot_measure_meshgrid(initial_theta_b1_loc(:)==1), '.', 'color', [1, 0, 1])
plot(xplot_measure_meshgrid(initial_theta_b2_loc(:)==1), yplot_measure_meshgrid(initial_theta_b2_loc(:)==1), '.', 'color', [0, 1, 0])
title(['$V_1=',num2str(initial_theta_b1,'%.3f'), ', M = ', num2str(M), '$'], 'interpreter','latex')
xlim([-r,r])
ylim([-r,r])
pbaspect([1 1 1])
set(gca, 'XTick', []);
set(gca, 'YTick', []);
set(f2,'Windowstyle','normal')
set(findall(f2,'-property','FontSize'),'FontSize',11)
set(f2,'Units','centimeters')
set(f2,'Position',[24 10 5 5])


%% solve
time_start = tic;
try
    result = solvepde(model,tlist); % for time-dependent problems
catch
    disp('Solver failure!')
end
toc(time_start)   


%% check if solution is feasible 
x = linspace(-r,r,50);
[X,Y] = meshgrid(x,x);
interpol_sol = interpolateSolution(result,X,Y,[1,2],length(result.SolutionTimes));
interpol_sol_B1 = reshape(interpol_sol(:,1),size(X));
interpol_sol_B2 = reshape(interpol_sol(:,2),size(X));
if min(interpol_sol_B1(:))<-2e-2 || min(interpol_sol_B2(:)) <-2e-2
disp('Warning: infeasible solution detected')
end

%% area covered by range expansion in solution

[X,Y] = meshgrid(xplot_measure);
querypoints = [X(:),Y(:)]';
sol_interp_B1 = interpolateSolution(result,querypoints,1,10); % interpolate solution to a regular meshgrid
sol_interp_B2 = interpolateSolution(result,querypoints,2,10);
sol_interp_B1 = reshape(sol_interp_B1,size(X));
sol_interp_B2 = reshape(sol_interp_B2,size(X));
biofilm_ind = find(sol_interp_B1+sol_interp_B2 >0.1*kcap); % criterium for decision if solution occupies a point in space
footprint_sol = zeros(1,length(sol_interp_B1(:)));
footprint_sol(biofilm_ind) = 1; % vector storing if solution occupies a point in space
footprint_area = area_el*sum(footprint_sol);
disp("The area covered by range expansion is " + num2str(footprint_area, '%.3f'))

corr_temp = corrcoef(footprint_sol,footprint_pred_vec); % find correlation betweeen prediction and model simulation of area covered
corr = corr_temp(1,2);  
disp("The correlation between prediction and model outcome is " + num2str(corr, '%.3f'))

%% plot solution at t_final
f = figure(10);
B1 = result.NodalSolution(:,1,end);
B2 = result.NodalSolution(:,2,end);
mapsize = 100;
map2 = [1,1,1] - [0:1/mapsize:1; zeros(1,mapsize+1); 0:1/mapsize:1]'; % define colourmaps for plotting
map1 = [1,1,1] - [zeros(1,mapsize+1); 0:1/mapsize:1; zeros(1,mapsize+1)]';

Axes(3) = axes('Position', [0.8/3+0.1 0.05 0.8/3 0.8]); 
pdeplot(model,'XYData',B1,'Contour','off','ColorBar','off','FaceAlpha',1)
hold on
pdeplot(model,'XYData',0*B2,'Contour','off','FaceAlpha',0.5,'ColorBar','off')
xlabel('$B_1$','interpreter','latex')
pbaspect([1 1 1])
caxis manual
caxis([0 kcap]);

Axes(4) = axes('Position', [2*0.8/3+0.15 0.05 0.8/3 0.8]);
pdeplot(model,'XYData',0*B1,'Contour','off','FaceAlpha',1,'ColorBar','off')
hold on
pdeplot(model,'XYData',B2,'Contour','off','ColorBar','off','FaceAlpha',0.5)
xlabel('$B_2$','interpreter','latex')
pbaspect([1 1 1])
caxis manual
caxis([0 kcap]);
sgtitle(['$\overline{B_1}^\Omega = ', num2str(sum(B1)/(sum(B1)+sum(B2)), '%.3f'), '$'],'interpreter','latex')

Axes(1) = axes('Position', [0.05 0.05 0.8/3 0.8]);
pbaspect([1 1 1])
hold on
pdeplot(model,'XYData',B1,'Contour','off','FaceAlpha',1,'ColorBar','off')
xlabel('Merged','interpreter','latex')
caxis manual
caxis([0 kcap]);
Axes(2) = axes('Position', [0.05 0.05 0.8/3 0.8]);
pdeplot(model,'XYData',B2,'Contour','off','FaceAlpha',0.5,'ColorBar','off')
hold on
pbaspect([1 1 1])
caxis manual
caxis([0 kcap]);
plot([-0.5; 0.5], [-9; -9], '-k', 'LineWidth', 2)

set(Axes(2), 'XTick', []);
set(Axes(2), 'YTick', []);
set(Axes(2), 'visible', 'off');
set(Axes(1), 'XTick', []);
set(Axes(1), 'YTick', []);
set(Axes(3), 'XTick', []);
set(Axes(3), 'YTick', []);
set(Axes(4), 'XTick', []);
set(Axes(4), 'YTick', []);

colormap(Axes(1),map1)
colormap(Axes(2),map2)
colormap(Axes(3),map1)
colormap(Axes(4),map2)

set(f,'Windowstyle','normal')
set(findall(f,'-property','FontSize'),'FontSize',11)
set(f,'Units','centimeters')
set(f,'Position',[30 10 3*17/4 5])




%% plot initial condition
clear Axes
R_in = 3;

f1 = figure(11);
Axes(3) = axes('Position', [0.8/3+0.1 0.05 0.8/3 0.8]); 
pdeplot(model,'XYData',B10,'Contour','off','ColorBar','off','FaceAlpha',1)
hold on
pdeplot(model,'XYData',0*B20,'Contour','off','FaceAlpha',0.5,'ColorBar','off')
xlabel('$B_1$','interpreter','latex')
pbaspect([1 1 1])
ylim([com(2)-R_in,com(2)+R_in])
xlim([com(1)-R_in,com(1)+R_in])
caxis manual
caxis([0 kcap]);

Axes(4) = axes('Position', [2*0.8/3+0.15 0.05 0.8/3 0.8]);
pdeplot(model,'XYData',0*B10,'Contour','off','FaceAlpha',1,'ColorBar','off')
hold on
pdeplot(model,'XYData',B20,'Contour','off','ColorBar','off','FaceAlpha',0.5)
xlabel('$B_2$','interpreter','latex')
pbaspect([1 1 1])
ylim([com(2)-R_in,com(2)+R_in])
xlim([com(1)-R_in,com(1)+R_in])
caxis manual
caxis([0 kcap]);

Axes(1) = axes('Position', [0.05 0.05 0.8/3 0.8]);
pbaspect([1 1 1])
hold on
pdeplot(model,'XYData',B10,'Contour','off','FaceAlpha',1,'ColorBar','off')
xlabel('Merged','interpreter','latex')
ylim([com(2)-R_in,com(2)+R_in])
xlim([com(1)-R_in,com(1)+R_in])
caxis manual
caxis([0 kcap]);

Axes(2) = axes('Position', [0.05 0.05 0.8/3 0.8]);
pdeplot(model,'XYData',B20,'Contour','off','FaceAlpha',0.5,'ColorBar','off')
hold on
pbaspect([1 1 1])
ylim([com(2)-R_in,com(2)+R_in])
xlim([com(1)-R_in,com(1)+R_in])
caxis manual
caxis([0 kcap]);
plot([com(1)+R_in-1.1; com(1)+R_in-0.1], [com(2)-R_in+0.1; com(2)-R_in + 0.1], '-k', 'LineWidth', 2)

set(Axes(2), 'XTick', []);
set(Axes(2), 'YTick', []);
set(Axes(2), 'visible', 'off');
set(Axes(1), 'XTick', []);
set(Axes(1), 'YTick', []);
set(Axes(3), 'XTick', []);
set(Axes(3), 'YTick', []);
set(Axes(4), 'XTick', []);
set(Axes(4), 'YTick', []);

colormap(Axes(1),map1)
colormap(Axes(2),map2)
colormap(Axes(3),map1)
colormap(Axes(4),map2)

set(f1,'Windowstyle','normal')
set(findall(f1,'-property','FontSize'),'FontSize',11)
set(f1,'Units','centimeters')
set(f1,'Position',[11 10 3*17/4 5])


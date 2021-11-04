function [graph, nodes_loc] = mesh2graph(mesh)
%MESH2GRAPH converts the spatial mesh of the PDE solver into a weighted
%graph used in the calculation of the front propagation metric.

% Author: Lukas Eigentler (leigentler001@dundee.ac.uk)
% License: GNU GPL
% Last updated: 19/10/2021

global dem r
graph = sparse(length(mesh.Elements(1,:)),length(mesh.Elements(1,:))); % initialise graph as a sparse matrix
nodes_loc = zeros(2,length(mesh.Elements(1,:)));

% find locations of nodes
for ee = 1:length(mesh.Elements(1,:)) % loop through all graph nodes
    nd = mesh.Elements(:,ee); % nodes of element
    nodes_loc(:,ee) = sum(mesh.Nodes(:,nd),2)/length(nd); % find centre of mesh element
end

X = linspace(-r(end),r(end),length(dem(:,1)));
L = interp2(X,X,dem,nodes_loc(1,:),nodes_loc(2,:)); % interpolate fractal landscape onto graph nodes
% define kmax and d
[kmax,~,d,~] = kmax_fct2(L);

% assign weights to edges
for ee = 1:length(mesh.Elements(1,:)) % loop through all graph nodes
    % find neighbour graph nodes
    distvec = vecnorm(nodes_loc(:,ee)-nodes_loc); % distance to all other nodes
    [~, sort_ind] = sort(distvec); % sort by distance
    neighbour_nodes = sort_ind(2:5); % assign closest 3 nodes as neighbours (note the problem on the boundary - no problem here as cell popu
    for nn = 1:length(neighbour_nodes) % loop through all neighbouring graph nodes
        if graph(ee,neighbour_nodes(nn)) == 0
            dist = norm(nodes_loc(:,ee) - nodes_loc(:,neighbour_nodes(nn))); % distance between graph nodes
            speed = sqrt(kmax(ee)*d(ee)) + sqrt(kmax(neighbour_nodes(nn))*d(neighbour_nodes(nn))); % avg travelling wave speed between graph nodes
            graph(ee,neighbour_nodes(nn)) = dist/speed; % assing graph edge weight
            graph(neighbour_nodes(nn),ee) = dist/speed; % assing graph edge weight
        end
    end
end


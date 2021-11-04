function [dist,prev] = Dijkstra(graph,start)
%DIJKSTRA applies the Dijkstra algorithm to a graph
%
% [dist,prev] = Dijkstra(graph,start) applies finds the shortest path in the sense of the 
%front propagation metric from the source node to all other nodes in the graph.
%
% SYNTAX
%    [dist,prev] = Dijkstra(graph,start)
%
% INPUT
%        graph: a matrix representing a graph. Matrix entries are the
%        weights of the edges. Unconnected nodes need to be represented by
%        an edge of weight Inf or zero. Sparse input matrices are accepted.
%        start: index of the starting node
%     
%
% OUTPUT
%       dist: entry i of the vector is the shortest distance from start to node i
%       prev: entry i of the vector is the index of the previous node in
%       the shortest path to reach node i
% EXAMPLES
%     G = [Inf,Inf,1,Inf; Inf,Inf,2,3;1,2,Inf,10;Inf,3,10,Inf];  
%     [dist,prev] = Dijkstra(graph,1)
%     
%
%
%
% Author:       Lukas Eigentler
% email:        leigentler001 AT dundee.ac.uk
% Matlab ver.:  2020a
% Date:         05-Mar-2021
% Update:       

N = length(graph(:,1)); % determine number of nodes

if ~issymmetric(graph) % check if matrix is symmetric
    error('Input graph is not symmetric');
end
if mod(start,1) ~= 0 || start>N % check if start index is integer and a node of graph
    error('Input start needs to be an integer representing a node of the graph')
end

unvisited = 1:N; % create vector of unvisited nodes
dist = inf*ones(1,N); dist(start) = 0; % initialise distance vector
prev = NaN*ones(1,N); % initialise vector of pointers to previous node


while ~isempty(unvisited) % repeat until every node is visited
   [~,temp_ind] = min(dist(unvisited)); % make unvisited node with shortest distance the current node  
   current_node = unvisited(temp_ind);
   unvisited(unvisited == current_node) = []; % mark current node as visited
   new_dist = dist(current_node) + graph(current_node,:); % calculate new distances
   temp_ind1 = find(new_dist<dist == 1); % find indices for which new node is closer
   temp_ind2 = find(~(new_dist == dist(current_node))==1); % find indices for which weight was positive
   replace_ind = intersect(temp_ind1,temp_ind2);
   dist(replace_ind) = new_dist(replace_ind); % use new distance if smaller than current
   prev(replace_ind) = current_node; % assign current node as the previous node for neighbour
end
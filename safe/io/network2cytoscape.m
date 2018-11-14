function network_json = network2cytoscape(adjacency, nodePositions, varargin)

network.data.name = 'SAFE';
nnodes = size(adjacency,1);
npositions = nodePositions;

for i = 1 : nnodes

    network.elements.nodes(i).data.id = sprintf('n%d', i);
    network.elements.nodes(i).data.name = sprintf('n%d', i);
    network.elements.nodes(i).data.shared_name = sprintf('n%d', i);
    network.elements.nodes(i).position.x = npositions(i,1);
    network.elements.nodes(i).position.y = npositions(i,2);

end

[edgeS, edgeT] = find(triu(adjacency,1) > 0);

for i = 1 : length(edgeS)
    
    network.elements.edges(i).data.source = sprintf('n%d', edgeS(i));
    network.elements.edges(i).data.target = sprintf('n%d', edgeT(i));

end

network_json = savejson('',network);


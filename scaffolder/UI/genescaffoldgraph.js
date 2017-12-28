var scaffoldlibs = [new ScaffoldEdgeLib(0, '#0000ff', 'lib1', 'RNA_PAIR'), new ScaffoldEdgeLib(1, '#00ff00', 'lib2', 'SCAFF'), new ScaffoldEdgeLib(2, '#ff0000', 'lib3', 'SCAFF')];
var scaffoldnodes = [new ScaffoldNode(0, "node0", 10), new ScaffoldNode(1, "node1", 15), new ScaffoldNode(2, "node2", 20)];
var scaffoldedges = [new ScaffoldEdge(0, 0, 1, 0, 1), new ScaffoldEdge(1, 1, 2, 2, 5)];

var scaffoldgraph = new ScaffoldGraph(scaffoldlibs, scaffoldnodes, scaffoldedges);
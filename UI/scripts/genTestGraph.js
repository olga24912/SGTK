var scaffoldlibs = [new ScaffoldEdgeLib(0, '#0000ff', 'lib1', 'RNA_PAIR')];
var scaffoldnodes = [];
for (var i = 0; i < 10000000; ++i) {
    scaffoldnodes.push(new ScaffoldNode(i, "node" + i.toString(), 501));
}

var scaffoldedges = [];
for (i = 0; i < (10000000 - 1); ++i) {
    scaffoldedges.push(new ScaffoldEdge(i, i, i + 1, 0, 1));
}

var scaffoldgraph = new ScaffoldGraph(scaffoldlibs, scaffoldnodes, scaffoldedges);
var chromosomes = [];
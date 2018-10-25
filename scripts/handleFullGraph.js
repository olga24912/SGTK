/*
* Create graph for filtration "full graph" in free(dagre) layout
*/


function handleFullGraph() {
    for (i=0; i < scaffoldgraph.nodes.length; ++i) {
        if (scaffoldgraph.nodes[i].len >= min_contig_len) {
            nodes_to_draw.push(scaffoldgraph.nodes[i].id);
            special_nodes.add(scaffoldgraph.nodes[i].id);
        }
    }

    for (i=0; i < scaffoldgraph.edges.length; ++i) {
        if (isGoodEdge(i)) {
            edges_to_draw.push(scaffoldgraph.edges[i].id);
            special_edges.add(scaffoldgraph.edges[i].id);
        }
    }
    splitOnParts(nodes_to_draw, edges_to_draw);

    var notzero = [];
    for (i=0; i < edges_set.length; ++i) {
        if (edges_set[i].length > 0) {
            notzero.push(i);
        }
    }

    createComponentShowList(function(j) {
        var i = notzero[j];
        DrawGraph(nodes_set[i], edges_set[i]);
    }, function(j) {
        var i = notzero[j];
        return "comp " + i + "<br/> #nodes = " + nodes_set[i].length + "<br/>#edges = " + edges_set[i].length;
    }, function(i) {
        var compnum = notzero[i];
        return "Component #" + compnum;
    }, notzero.length);
}

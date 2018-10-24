/*
* Create graph for filtration "by edge id" in free(dagre) layout
*/

function handleEdgeLocalArea() {
    special_nodes.clear();
    special_edges.clear();
    area_size = document.getElementById("area_size").value;
    var edges = document.getElementById("vertext").value.replace(/\n/g, " ").split(" ");
    var edges_id = [];
    nodes_id = [];
    for (i = 0; i < edges.length; ++i) {
        edges_id.push(parseInt(edges[i]));
        nodes_id.push(scaffoldgraph.edges[edges_id[i]].from);
        nodes_id.push(scaffoldgraph.edges[edges_id[i]].to);
        special_edges.add(edges_id[edges_id.length - 1]);
    }

    var uniqueNode = nodes_id.filter(function (value, index, self) {
        return self.indexOf(value) === index;
    });

    if (area_size > 0) {
        findLocalArea(uniqueNode, area_size, min_contig_len, isGoodEdge);
    } else {
        for (i = 0; i < uniqueNode.length; ++i) {
            nodes_to_draw.push(uniqueNode[i]);
        }
        for (i = 0; i < edges_id.length; ++i) {
            edges_to_draw.push(edges_id[i]);
        }
    }
    splitOnParts(nodes_to_draw, edges_to_draw);

    createComponentShowList(function (i) {
        DrawGraph(nodes_set[i], edges_set[i]);
        if (edges_id.length > 0) {
            cy.fit(cy.$('#e' + edges_id[0]));
        }
    }, function (i) {
        var nm = "";
        for (var j = 0; j < edges_set[i].length; ++j) {
            if (special_edges.has(edges_set[i][j])) {
                nm += "edge id=" + edges_set[i][j].toString() + "<br/>";
            }
        }
        return nm;
    }, function (i) {
        var nm = "";
        for (var j = 0; j < edges_set[i].length; ++j) {
            if (special_edges.has(edges_set[i][j])) {
                nm += "edge id: " + edges_set[i][j].toString() + "<br/>";
                break;
            }
        }
        return nm;
    }, nodes_set.length);
}

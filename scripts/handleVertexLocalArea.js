/*
* Create graph for filtration "by vertex id" in free(dagre) layout
*/

function handleVertexLocalArea() {
    special_nodes.clear();
    special_edges.clear();
    area_size = document.getElementById("area_size").value;
    var nodes = document.getElementById("vertext").value.replace(/\n/g, " ").split(" ");
    var nodes_id = [];
    for (i = 0; i < nodes.length; ++i) {
        if (scaffoldgraph.id_by_name.has(nodes[i])) {
            nodes_id.push(scaffoldgraph.id_by_name.get(nodes[i]));
        } else {
            nodes_id.push(parseInt(nodes[i]));
        }
        special_nodes.add(nodes_id[nodes_id.length - 1]);
    }

    findLocalArea(nodes_id, area_size, min_contig_len, isGoodEdge);

    splitOnParts(nodes_to_draw, edges_to_draw);

    createComponentShowList(function (i) {
        DrawGraph(nodes_set[i], edges_set[i]);
        if (nodes_id.length > 0) {
            cy.fit(cy.$('#' + nodes_id[0]));
        }
    }, function (i) {
        var nm = "";
        for (var j = 0; j < nodes_set[i].length; ++j) {
            if (special_nodes.has(nodes_set[i][j])) {
                nm += scaffoldgraph.nodes[nodes_set[i][j]].name + " id=" + nodes_set[i][j].toString() + "<br/>";
            }
        }
        return nm;
    }, function (compnum) {
        var nm = "";
        for (var j = 0; j < nodes_set[compnum].length; ++j) {
            if (special_nodes.has(nodes_set[compnum][j])) {
                nm += scaffoldgraph.nodes[nodes_set[compnum][j]].name + " id=" + nodes_set[compnum][j].toString() + "<br/>";
                break;
            }
        }
        return nm;
    }, nodes_set.length);
}

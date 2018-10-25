/*
* Create graph for filtration "difference between sources" in free(dagre) layout
*/

function handleDiffInLibsFilter(areasize, min_contig_len, isGoodEdge) {
    special_nodes.clear();
    special_edges.clear();

    var presentLibMask = [];
    for (var i=0; i < scaffoldgraph.libs.length; ++i) {
        presentLibMask.push(2);
        if (document.getElementById("checkbox_present_"+i.toString()).checked) {
            presentLibMask[i] = 1;
        } else if (document.getElementById("checkbox_not_present_"+i.toString()).checked) {
            presentLibMask[i] = 0;
        }
    }

    var correct = document.getElementById("checkbox_correct").checked;
    var wrong = document.getElementById("checkbox_wrong").checked;

    var nodes_id = [];

    for (i = 0; i < scaffoldgraph.nodes.length; ++i) {
        var edgesSortByTo = [];
        for (var j = 0; j < scaffoldgraph.g[i].length; ++j) {
            if (isGoodEdge(scaffoldgraph.g[i][j].id)) {
                edgesSortByTo.push({to: scaffoldgraph.g[i][j].to, lib: scaffoldgraph.g[i][j].lib, id: scaffoldgraph.g[i][j].id});
            }
        }

        edgesSortByTo.sort(function (a, b) {
            return a.to - b.to;
        });

        edgesSortByTo.push({to: -2, lib: -1, id: -1});

        var lstval = -1;
        var curMask = [];
        var curEdges = [];
        for (j = 0; j < scaffoldgraph.libs.length; ++j) {
            curMask.push(0);
        }

        for (j = 0; j < edgesSortByTo.length; ++j) {
            if (!(lstval === -1) && !(lstval === edgesSortByTo[j].to)) {
                console.log("Process " + i.toString() + " " + lstval.toString());

                var isGoodConnection = false;
                var crcOrder = isCorrectOrder(i, lstval);
                if (crcOrder === true && correct === true) {
                    isGoodConnection = true;
                }
                if (crcOrder === false && wrong === true) {
                    isGoodConnection = true;
                }

                console.log(curMask[0].toString() + " " + curMask[1].toString());
                for (var g = 0; g < scaffoldgraph.libs.length && isGoodConnection; ++g) {
                    if ((presentLibMask[g] === 1 && curMask[g] === 0) || (presentLibMask[g] === 0 && curMask[g] === 1)) {
                        isGoodConnection = false;
                    }
                }

                for (g = 0; g < scaffoldgraph.libs.length; ++g) {
                    curMask[g] = 0;
                }

                if (isGoodConnection) {
                    console.log("connection " + i.toString() + " " + lstval.toString());
                    nodes_id.push(i);
                    nodes_id.push(lstval);
                    special_nodes.add(i);
                    special_nodes.add(lstval);

                    for (g = 0; g < curEdges.length; ++g) {
                        if (presentLibMask[curEdges[g].lib] === 1) {
                            special_edges.add(curEdges[g].id);
                        }
                    }
                }
                curEdges = [];
            }
            curEdges.push(edgesSortByTo[j]);
            lstval = edgesSortByTo[j].to;
            curMask[edgesSortByTo[j].lib] = 1;
        }
    }

    var uniqueNode = nodes_id.filter(function (value, index, self) {
        return self.indexOf(value) === index;
    });

    findLocalArea(uniqueNode, areasize, min_contig_len, isGoodEdge);

    splitOnParts(nodes_to_draw, edges_to_draw);
    createComponentShowList(function(i) {
        DrawGraph(nodes_set[i], edges_set[i]);
    }, function(i) {
        return "comp " + i + "<br/> #nodes = " + nodes_set[i].length + "<br/>#edges = " + edges_set[i].length;
    }, function(compnum) {
        return "Component #" + compnum;
    }, nodes_set.length);
}
/*
* Handle opening vertex
*/

//find coordinate for new vertex
function getYforNewVert(v, u, evt, isGoodEdge) {
    for (var i = 0; i < scaffoldgraph.g[v].length; ++i) {
        if (scaffoldgraph.g[v][i].to === u && isGoodEdge(scaffoldgraph.g[v][i].id)) {
            return evt.target.position().x + Math.random() * Math.floor(200);
        }
    }

    return evt.target.position().x - Math.random() * Math.floor(200);
}

function createAddNewNode(cy, curNodeSet) {
    cy.on('tap', 'node', function (evt) {
        var v = evt.target.id();
        var needAddVert = [];
        var needAddEdge = [];
        var newNode = new Set();
        findConnectedVertex(v, needAddVert, needAddEdge, newNode, area_size, min_contig_len, isGoodEdge, curNodeSet);

        for (g = 0; g < needAddVert.length; ++g) {
            u = needAddVert[g];
            var yc = getYforNewVert(v, u, evt, isGoodEdge);

            var spe = 0;
            var opt = document.getElementById("select_show_type").value;
            if (opt === "full graph") {
                spe = 1;
            }

            var nnode = {
                group: "nodes",
                data: {
                    id: u,
                    label: createLabelForNode(u),
                    len: 2*Math.log2(scaffoldgraph.nodes[nodes_to_draw[g]].len)/Math.log2(1.5),
                    /*shape: 'ellipse',*/
                    notALL: 0,
                    color: genColorNode(u),

                    color1: "#2A4986",
                    color2: "#2A4986",
                    color3: "#2A4986",
                    cnt1: 100,
                    cnt2: 0,
                    cnt3: 0,
                    special: spe
                },
                position: {
                    y: evt.target.position().y + Math.floor(Math.random() * Math.floor(200) - 100),
                    x: yc
                }
            };
            updateColorsNode(u, nnode);

            cy.add(nnode);
        }

        for (g = 0; g < needAddEdge.length; ++g) {
            var eid = needAddEdge[g].id;

            spe = 0;
            opt = document.getElementById("select_show_type").value;
            if (opt === "full graph") {
                spe = 1;
            }

            cy.add({
                group: "edges",
                data: {
                    id: "e" + eid.toString(),
                    source: getEdgeFrom(eid),
                    target: getEdgeTo(eid),
                    label: createLabelForEdge(eid),
                    faveColor: scaffoldgraph.libs[scaffoldgraph.edges[eid].lib].color,
                    weight: generateEdgeWeight(eid),
                    lstyle: 'dotted',
                    special: spe
                }
            });
        }

        for (g = 0; g < nodes_to_draw.length; ++g) {
            if (hasOtherEdges(nodes_to_draw[g], curNodeSet)) {
                console.log(nodes_to_draw[g].toString() + "  notAll");
                cy.$('#' + nodes_to_draw[g]).data('notALL', 1);
            } else {
                console.log(nodes_to_draw[g].toString() + " thats all");
                cy.$('#' + nodes_to_draw[g]).data('notALL', 0);
            }
        }
        //createTapInfo(cy);
    });
}

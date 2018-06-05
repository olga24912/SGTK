function hasOtherEdges(v, curNodeSet) {
    for (var h = 0; h < scaffoldgraph.g[v].length; ++h) {
        if (isGoodEdge(scaffoldgraph.g[v][h].id)) {
            if (scaffoldgraph.nodes[scaffoldgraph.g[v][h].to].len >= min_contig_len) {
                if (!curNodeSet.has(scaffoldgraph.g[v][h].to)) {
                    return true;
                }
            }
        }
    }

    for (h = 0; h < scaffoldgraph.gr[v].length; ++h) {
        if (isGoodEdge(scaffoldgraph.gr[v][h].id)) {
            if (scaffoldgraph.nodes[scaffoldgraph.gr[v][h].from].len >= min_contig_len) {
                if (!curNodeSet.has(scaffoldgraph.gr[v][h].from)) {
                    return true;
                }
            }
        }
    }

    return false;
}


function getYforNewVert(v, u, evt, isGoodEdge) {
    for (var i = 0; i < scaffoldgraph.g[v].length; ++i) {
        if (scaffoldgraph.g[v][i].to === u && isGoodEdge(scaffoldgraph.g[v][i].id)) {
            return evt.target.position().y + Math.random() * Math.floor(200);
        }
    }

    return evt.target.position().y - Math.random() * Math.floor(200);
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

            var nnode = {
                group: "nodes",
                data: {
                    id: u,
                    label: createLabelForNode(u),
                    len: 2*Math.log2(scaffoldgraph.nodes[nodes_to_draw[g]].len)/Math.log2(1.5),
                    shape: 'ellipse',
                    color: genColorNode(u),

                    color1: "#2A4986",
                    color2: "#2A4986",
                    color3: "#2A4986",
                    cnt1: 100,
                    cnt2: 0,
                    cnt3: 0,
                    special: 0
                },
                position: {
                    x: evt.target.position().x + Math.floor(Math.random() * Math.floor(200) - 100),
                    y: yc
                }
            };
            updateColorsNode(u, nnode);

            cy.add(nnode);
        }

        for (g = 0; g < needAddEdge.length; ++g) {
            var eid = needAddEdge[g].id;

            cy.add({
                group: "edges",
                data: {
                    id: "e" + eid.toString(),
                    source: getEdgeFrom(eid),
                    target: getEdgeTo(eid),
                    label: createLabelForEdge(eid),
                    faveColor: scaffoldgraph.libs[scaffoldgraph.edges[eid].lib].color,
                    weight: Math.min(5, Math.log(scaffoldgraph.edges[eid].weight) + 1),
                    lstyle: 'dotted'
                }
            });
        }

        for (g = 0; g < nodes_to_draw.length; ++g) {
            if (hasOtherEdges(nodes_to_draw[g], curNodeSet)) {
                cy.$('#' + nodes_to_draw[g]).data('shape', 'vee');
            } else {
                cy.$('#' + nodes_to_draw[g]).data('shape', 'ellipse');
            }
        }
        createTapInfo(cy);
    });
}

function DrawGraphCytoscapeWithPresetNode(dnodes, dedges, curNodeSet) {
    cy = cytoscape({
        container: document.getElementById('mainpanel'),

        boxSelectionEnabled: false,
        autounselectify: true,
        maxZoom: 20,
        minZoom: 0.02,

        elements: {
            nodes: dnodes,
            edges: dedges
        },

        layout: {
            name: 'dagre',
            rankDir: 'UD',
            ranker: 'tight-tree', //'longest-path',//'network-simplex',//
            edgeWeight: function( edge ){
                return edge.width();
                console.log(edge.id());
                console.log(100*special_edges.has(parseInt(edge.id().substring(1, edge.id().length))));
                return 100*special_edges.has(edge.id().substring(1, edge.id().length));
            }
        },

        ready: function () {
            window.cy = this;
        },

        style: cytoscape.stylesheet()
            .selector('node')
            .css({
                'content': 'data(label)',
                'color': '#2A4986',
                'background-color': '#2A4986',
                'width': 'data(len)',
                'height': 'data(len)',
                'border-width': 'mapData(special, 0, 1, 0px, 5px)',
                'pie-size': '100%',
                'pie-1-background-color': 'data(color1)',
                'pie-1-background-size': 'data(cnt1)',
                'pie-2-background-color': 'data(color2)',
                'pie-2-background-size': 'data(cnt2)',
                'pie-3-background-color': 'data(color3)',
                'pie-3-background-size': 'data(cnt3)',
                'shape' : 'data(shape)'
            })
            .selector('edge')
            .css({
                'curve-style': 'bezier',
                'target-arrow-shape': 'triangle',
                'line-color': 'data(faveColor)',
                'target-arrow-color': 'data(faveColor)',
                'width': 'data(weight)',
                'line-style': 'data(lstyle)',
                'content': 'data(label)'
            })
    });

    var element =  document.getElementById("cynav");
    if (typeof(element) != 'undefined' && element != null) {
        document.getElementById("cynav").remove();
    }

    cy.navigator({
        // options...
    });


    cy.on('zoom', function () {
        (document.getElementById("zoomInput")).innerText = Math.floor(cy.zoom() * 100).toString() +  "%";
    });

    createTapInfo(cy);
    createAddNewNode(cy, curNodeSet);
}

function genColorNode(vid) {
    cur_max = 0;
    cur_color = '#2a4986';
    for (var i = 0; i < scaffoldgraph.nodes[vid].alignments.length; ++i) {
        alig = scaffoldgraph.nodes[vid].alignments[i];
        if (alig.coorde - alig.coordb > cur_max) {
            cur_max = alig.coorde - alig.coordb;
            cur_color = chromosomes[alig.chr_id].color;
        }
    }
    return cur_color;
}

function updateColorsNode(vid, node) {
    var chrmCnt = [];
    for(var i = 0; i < chromosomes.length; ++i) {
        chrmCnt.push(0);
    }

    for (i = 0; i < scaffoldgraph.nodes[vid].alignments.length; ++i) {
        alig = scaffoldgraph.nodes[vid].alignments[i];
        chrmCnt[alig.chr_id] += alig.coorde - alig.coordb;
    }

    var sortChrVal = [];
    for(i = 0; i < chromosomes.length; ++i) {
        sortChrVal.push({
            id : i,
            val : chrmCnt[i]
        });
    }
    for (i = chromosomes.length; i < 4; ++i) {
        sortChrVal.push({
            id : 0,
            val : 0
        });
    }

    sortChrVal.sort(function (a, b) {
        return b.val - a.val;
    });

    var sum = sortChrVal[0].val + sortChrVal[1].val + sortChrVal[2].val;
    if (sum === 0) {
        return;
    }
    node.data.color1 = chromosomes[sortChrVal[0].id].color;
    node.data.color2 = chromosomes[sortChrVal[1].id].color;
    node.data.color3 = chromosomes[sortChrVal[2].id].color;
    node.data.cnt1 = sortChrVal[0].val*100/sum;
    node.data.cnt2 = sortChrVal[1].val*100/sum;
    node.data.cnt3 = sortChrVal[2].val*100/sum;
//    alert(node.data.color1);
//    alert(node.data.cnt1);
//    alert(node.data.color2);
//    alert(node.data.cnt2);
//    alert(node.data.color3);
//    alert(node.data.cnt3);
}

function DrawGraphCytoscape(nodes_to_draw, edges_to_draw) {
    var curNodeSet = new Set();
    for (var g = 0; g < nodes_to_draw.length; ++g) {
        curNodeSet.add(nodes_to_draw[g]);
    }

    var dnodes = [];
    var dedges = [];
    for (g = 0; g < nodes_to_draw.length; ++g) {
        var nall = 'ellipse';
        if (hasOtherEdges(nodes_to_draw[g], curNodeSet)) {
            nall = 'vee';
        }

        var sp = 0;
        if (special_nodes.has(nodes_to_draw[g])) {
            sp = 1;
        }

        dnodes.push({
            data: {
                id: nodes_to_draw[g],
                label: createLabelForNode(nodes_to_draw[g]),
                len: 2*Math.log2(scaffoldgraph.nodes[nodes_to_draw[g]].len)/Math.log2(1.5),
                shape: nall,
                color: genColorNode(nodes_to_draw[g]),
                color1: "#2A4986",
                color2: "#2A4986",
                color3: "#2A4986",
                cnt1: 100,
                cnt2: 0,
                cnt3: 0,
                special: sp
            }
        });

        updateColorsNode(nodes_to_draw[g], dnodes[dnodes.length - 1]);

    }

    for (g = 0; g < edges_to_draw.length; ++g) {
        sp = 'dotted';
        if (special_edges.has(edges_to_draw[g])) {
            sp = 'solid';
        }

        //TODO: small edge
        dedges.push({ data: {
                id: "e" + edges_to_draw[g].toString(),
                source: getEdgeFrom(edges_to_draw[g]),
                target: getEdgeTo(edges_to_draw[g]),
                label: createLabelForEdge(edges_to_draw[g]),
                faveColor: scaffoldgraph.libs[scaffoldgraph.edges[edges_to_draw[g]].lib].color,
                weight: Math.min(5, Math.log(scaffoldgraph.edges[edges_to_draw[g]].weight) + 1),
                lstyle: sp
            }});
    }

    DrawGraphCytoscapeWithPresetNode(dnodes, dedges, curNodeSet);
}

function DrawGraph(nodes_to_draw, edges_to_draw) {
    DrawGraphCytoscape(nodes_to_draw, edges_to_draw);
}

function findComponent(v, g, color, curc) {
    color.set(v, curc);
    var nb = g.get(v);
    for (var i=0; i < nb.length; ++i) {
        var u = nb[i];
        if (color.get(u) === -1) {
            findComponent(u, g, color, curc);
        }
    }
}

function findArea(ev, g, color, curc) {
    var edge_cnt = 0;
    var queue = [];
    queue.push(ev);
    color.set(ev, curc);
    edge_cnt += 1;
    var bg = 0;
    while (bg < queue.length) {
        ev = queue[bg];
        bg += 1;

        var nb = g.get(scaffoldgraph.edges[ev].from);
        for (i=0; i < nb.length; ++i) {
            var u = nb[i];
            if (color.get(u) === -1) {
                queue.push(u);
                color.set(u, curc);
                edge_cnt += 1;
            }
        }

        nb = g.get(scaffoldgraph.edges[ev].to);

        for (var i=0; i < nb.length; ++i) {
            u = nb[i];
            if (color.get(u) === -1) {
                queue.push(u);
                color.set(u, curc);
                edge_cnt += 1;
            }
        }


        if (edge_cnt > 200) {
            return;
        }
    }
}

function elemInList(elem, lst) {
   for (var i = 0; i < lst.length; ++i) {
       if (elem == lst[i]) {
           return true;
       }
   }

   return false;
}

function splitOnParts(nodes_to_draw, edges_to_draw) {
    var g = new Map();
    var color = new Map();
    for (var i = 0; i < edges_to_draw.length; ++i) {
        color.set(edges_to_draw[i], -1);
    }
    for (i = 0; i < nodes_to_draw.length; ++i) {
        g.set(nodes_to_draw[i], []);
    }

    for (i=0; i < edges_to_draw.length; ++i) {
        g.get(scaffoldgraph.edges[edges_to_draw[i]].from).push(edges_to_draw[i]);
        g.get(scaffoldgraph.edges[edges_to_draw[i]].to).push(edges_to_draw[i]);
    }

    var curc = 0;
    for (i=0; i < edges_to_draw.length; ++i) {
        if (color.get(edges_to_draw[i]) === -1) {
            findArea(edges_to_draw[i], g, color, curc);
            ++curc;
        }
    }

    nodes_set = [];
    edges_set = [];
    for (i=0; i < curc; ++i) {
        nodes_set.push([]);
        edges_set.push([]);
    }

    for(i=0; i < edges_to_draw.length; ++i) {
        edges_set[color.get(edges_to_draw[i])].push(edges_to_draw[i]);

        if (!elemInList(scaffoldgraph.edges[edges_to_draw[i]].from, nodes_set[color.get(edges_to_draw[i])])) {
            nodes_set[color.get(edges_to_draw[i])].push(scaffoldgraph.edges[edges_to_draw[i]].from);
        }
        if (!elemInList(scaffoldgraph.edges[edges_to_draw[i]].to, nodes_set[color.get(edges_to_draw[i])])) {
            nodes_set[color.get(edges_to_draw[i])].push(scaffoldgraph.edges[edges_to_draw[i]].to);
        }
    }
}

function isCorrectOrder(v, u) {
    var nv = scaffoldgraph.nodes[v];
    var nu = scaffoldgraph.nodes[u];

    for (var i = 0; i < nv.alignments.length; ++i) {
        //if (nv.alignments[i].coordne + 1000 >= nv.len) {
        for (var j = 0; j < nu.alignments.length; ++j) {
            //if (nu.alignments[j].coordnb - 1000 <= 0) {
            if (nv.alignments[i].chr_id === nu.alignments[j].chr_id) {
                if (nv.alignments[i].coorde - 100 < nu.alignments[j].coordb) {
                    if (nu.alignments[j].coordb - nv.alignments[i].coorde < 5000) {
                        return true;
                    }
                }
            }
            //}
        }
        //}
    }

    return false;
}

function updateText() {
    var nodes = cy.nodes();
    for (var i = 0; i < nodes.length; ++i) {
        cy.$('#' + nodes[i].id().toString()).data('label', createLabelForNode(nodes[i].id()));
    }
    var edges = cy.edges();
    for (var i = 0; i < edges.length; ++i) {
        cy.$('#' + edges[i].id().toString()).data('label', createLabelForEdge(parseInt(edges[i].id().substring(1))));
    }
}
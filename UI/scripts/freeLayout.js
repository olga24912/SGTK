/*
* Generate cytoscape graph for free(dagre) layout
*/


//get cytoscape edge width(int) by EdgeId(int)
function generateEdgeWeight(eid) {
    var edge = scaffoldgraph.edges[eid];

    var opt = document.getElementById("select_show_type").value;
    if (opt !== "full graph") {
        if (special_edges.has(eid)) {
            return 5;
        }
    }

    if (scaffoldgraph.libs[edge.lib].type === "FASTG" ||
        scaffoldgraph.libs[edge.lib].type === "GFA" || 
        scaffoldgraph.libs[edge.lib].type === "GFA2") {
        return 4;
    }

    if (scaffoldgraph.libs[edge.lib].type === "SCAFF") {
        return 3;
    }

    if (scaffoldgraph.libs[edge.lib].type === "DNA_PAIR" ||
        scaffoldgraph.libs[edge.lib].type === "LONG" ||
        scaffoldgraph.libs[edge.lib].type === "RNA_PAIR" ||
        scaffoldgraph.libs[edge.lib].type === "RNA_SPLIT_50" ||
        scaffoldgraph.libs[edge.lib].type === "RNA_SPLIT_30" ||
        scaffoldgraph.libs[edge.lib].type === "CONNECTION") {
        return 1;
    }

    return Math.min(5, Math.log(edge.weight) + 1);
}

//get vertex color by vertexId(int)
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

//Set pie color to cytoscape node by top free alignment of vid(int)
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

//Create cytoscape graph for cytoscape elements
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
            rankDir: 'LR',
            ranker: 'tight-tree', //'longest-path',//'network-simplex',//
            edgeWeight: function( edge ){
                return edge.width();
                //console.log(edge.id());
                //console.log(100*special_edges.has(parseInt(edge.id().substring(1, edge.id().length))));
                //return 100*special_edges.has(edge.id().substring(1, edge.id().length));
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
                'background-color': '#ffffff',
                'width': 'data(len)',
                'height': 'data(len)',
                'border-width': 'mapData(notALL, 0, 1, 0px, 5px)',
                'opacity':  'mapData(special, 0, 1, 0.6, 1)',
                'pie-size': '100%',
                'pie-1-background-color': 'data(color1)',
                'pie-1-background-size': 'data(cnt1)',
                'pie-2-background-color': 'data(color2)',
                'pie-2-background-size': 'data(cnt2)',
                'pie-3-background-color': 'data(color3)',
                'pie-3-background-size': 'data(cnt3)'//,
                //'pie-1-background-opacity': 'mapData(special, 0, 1, 1, 0.7)',
                //'pie-2-background-opacity': 'mapData(special, 0, 1, 1, 0.7)',
                //'pie-3-background-opacity': 'mapData(special, 0, 1, 1, 0.7)'
                //'shape' : 'data(shape)'
            })
            .selector('edge')
            .css({
                'curve-style': 'bezier',
                'target-arrow-shape': 'triangle',
                'line-color': 'data(faveColor)',
                'target-arrow-color': 'data(faveColor)',
                'width': 'data(weight)',
                'opacity':  'mapData(special, 0, 1, 0.6, 1)',
                //'line-style': 'data(lstyle)',
                'content': 'data(label)',
                //'target-endpoint': '-25% 0px',
                'target-distance-from-node': '1px'
            })
            .selector('.found')
            .css({
                'opacity': 0.3
            })
            .selector('.highlight')
            .css({
                'overlay-color': '#2A4986',
                'overlay-opacity': "0.5"
            })
            .selector('.fake')
            .css({
                'overlay-color': "#FF3503",
                'overlay-opacity': "0.5"
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

    cy.on('cxttap', 'node', function (evt) {
        var v = evt.target.id();
        cy.remove(cy.$("#" + v.toString()));
    });

    cy.on('cxttap', 'edge', function (evt) {
        var v = evt.target.id();
        cy.remove(cy.$("#" + v.toString()));
    });

    createInformationShown(cy);
    //createTapInfo(cy);
    createAddNewNode(cy, curNodeSet);
    highlightOnTap(cy)
}


//Crate cytoscape graph for list of nodesId(int[]) and edgesId(int[])
function DrawGraphCytoscape(nodes_to_draw, edges_to_draw) {
    var curNodeSet = new Set();
    for (var g = 0; g < nodes_to_draw.length; ++g) {
        curNodeSet.add(nodes_to_draw[g]);
    }

    var dnodes = [];
    var dedges = [];
    for (g = 0; g < nodes_to_draw.length; ++g) {
        var nall = 0; /*'ellipse';*/
        if (hasOtherEdges(nodes_to_draw[g], curNodeSet)) {
            nall = 1; /*'vee';*/
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
                //shape: nall,
                notALL: nall,
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
        var spe = 0;
        if (special_edges.has(edges_to_draw[g])) {
            sp = 'solid';
            spe = 1;
        }

        //TODO: small edge
        dedges.push({ data: {
                id: "e" + edges_to_draw[g].toString(),
                source: getEdgeFrom(edges_to_draw[g]),
                target: getEdgeTo(edges_to_draw[g]),
                label: createLabelForEdge(edges_to_draw[g]),
                faveColor: scaffoldgraph.libs[scaffoldgraph.edges[edges_to_draw[g]].lib].color,
                weight: generateEdgeWeight(edges_to_draw[g]),
                special: spe,
                lstyle: sp
            }});
    }

    DrawGraphCytoscapeWithPresetNode(dnodes, dedges, curNodeSet);
}

//Create graph for list of nodesId(int[]) and edgesId(int[])
function DrawGraph(nodes_to_draw, edges_to_draw) {
    DrawGraphCytoscape(nodes_to_draw, edges_to_draw);
}

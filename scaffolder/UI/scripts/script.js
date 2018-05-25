function createLabelForNode(node) {
    if (defZoom > 1000) {
        return "";
    }

    var label = "";
    if (document.getElementById("vert_checkbox_id").checked) {
        label += "id: " + scaffoldgraph.nodes[node].id + "\n";
    }
    if (document.getElementById("vert_checkbox_name").checked) {
        label += scaffoldgraph.nodes[node].name + "\n";
    }
    if (document.getElementById("vert_checkbox_len").checked) {
        label += "len: " + scaffoldgraph.nodes[node].len + "\n";
    }

    if (document.getElementById("vert_checkbox_align").checked) {
        if (scaffoldgraph.nodes[node].alignments.length > 0) {
            label += "Alignment: ";
            scaffoldgraph.nodes[node].alignments.sort(function (a, b) {
                return (b.coorde - b.coordb) - (a.coorde - a.coordb);
            });
            for (var i = 0; i < scaffoldgraph.nodes[node].alignments.length; ++i) {
                var cura = scaffoldgraph.nodes[node].alignments[i];
                label += chromosomes[cura.chr_id].name + " " + cura.coordb + " " + cura.coorde + " (" + ((cura.coorde - cura.coordb + 1) * 100/scaffoldgraph.nodes[node].len) +"%)\n";
            }
        }
    }
    if (document.getElementById("vert_checkbox_info").checked) {
        if (scaffoldgraph.nodes[node].info !== "") {
            label += scaffoldgraph.nodes[node].info + "\n";
        }
    }

    return label;
}

function createFullLabelForNode(node) {
    var label = "";
    label += "id: " + scaffoldgraph.nodes[node].id + "<br/>";
    label += scaffoldgraph.nodes[node].name + "<br/>";
    label += "len: " + scaffoldgraph.nodes[node].len + "<br/>";
    if (scaffoldgraph.nodes[node].alignments.length > 0) {
        label += "Alignment: ";
        scaffoldgraph.nodes[node].alignments.sort(function (a, b) {
            return (b.coorde - b.coordb) - (a.coorde - a.coordb);
        });
        for (var i = 0; i < scaffoldgraph.nodes[node].alignments.length; ++i) {
            var cura = scaffoldgraph.nodes[node].alignments[i];
            label += chromosomes[cura.chr_id].name + " " + cura.coordb + " " + cura.coorde + " (" + ((cura.coorde - cura.coordb + 1) * 100/scaffoldgraph.nodes[node].len) +"%)<br/>";
        }
    }
    if (scaffoldgraph.nodes[node].info !== "") {
        label += scaffoldgraph.nodes[node].info + "<br/>";
    }
    return label;
}

function createLabelForEdge(edge) {
    if (defZoom > 1000) {
        return "";
    }

    var label = "";
    if (document.getElementById("edge_checkbox_id").checked) {
        label += "id: " + scaffoldgraph.edges[edge].id + "\n";
    }
    if (document.getElementById("edge_checkbox_name").checked) {
        label += scaffoldgraph.libs[scaffoldgraph.edges[edge].lib].name + "\n";
        if (scaffoldgraph.libs[scaffoldgraph.edges[edge].lib].type === "SCAFF") {
            label += scaffoldgraph.edges[edge].name + "\n";
        }
    }
    if (document.getElementById("edge_checkbox_weight").checked) {
        label += "w: " + scaffoldgraph.edges[edge].weight + "\n";
    }
    if (document.getElementById("edge_checkbox_type").checked) {
        label += scaffoldgraph.libs[scaffoldgraph.edges[edge].lib].type + "\n";
    }
    if (document.getElementById("edge_checkbox_len").checked) {
        if (scaffoldgraph.edges[edge].len >= 0) {
            label += "len: " + scaffoldgraph.edges[edge].len + "\n";
        }
    }
    if (document.getElementById("edge_checkbox_info").checked) {
        if (scaffoldgraph.edges[edge].info !== "") {
            label += scaffoldgraph.edges[edge].info + "\n";
        }
    }
    return label;
}

function createFullLabelForEdge(edge) {
    var label = "";
    label += "id: " + scaffoldgraph.edges[edge].id + "<br/>";
    label += scaffoldgraph.libs[scaffoldgraph.edges[edge].lib].name + "<br/>";
    if (scaffoldgraph.libs[scaffoldgraph.edges[edge].lib].type === "SCAFF") {
        label += scaffoldgraph.edges[edge].name + "<br/>";
    }
    label += "w: " + scaffoldgraph.edges[edge].weight + "<br/>";
    label += scaffoldgraph.libs[scaffoldgraph.edges[edge].lib].type + "<br/>";
    if (scaffoldgraph.edges[edge].len >= 0) {
        label += "len: " + scaffoldgraph.edges[edge].len + "<br/>";
    }
    if (scaffoldgraph.edges[edge].info !== "") {
        label += scaffoldgraph.edges[edge].info + "<br/>";
    }
    return label;
}

function createTapInfo(cy) {
    cy.nodes().qtip({
        content: {
            text: function () {
                return createFullLabelForNode(this.id());
            }
        },
        show: {
            event: 'mouseover'
        },
        hide: {
            event: 'mouseout'
        },
        style: {
            classes: 'qtip-bootstrap',
            tip: {
                width: 16,
                height: 8
            }
        }
    });

    cy.edges().qtip({
        content: function(){
            return createFullLabelForEdge(parseInt(this.id().substring(1)));
        },
        position: {
            my: 'bottom left',
            at: 'bottom center'
        },
        show: {
            event: 'mouseover'
        },
        hide: {
            event: 'mouseout'
        },
        style: {
            classes: 'qtip-bootstrap',
            tip: {
                width: 10,
                height: 8
            }
        }
    });
}

var edges_set = [];
var nodes_set = [];
var special_nodes = new Set();
var special_edges = new Set();
var nodes_to_draw = [];
var edges_to_draw = [];
var isGoodEdge;
var min_contig_len = 0;
var area_size = 0;
var cur_show_id = 0;
var graph = null;
var cy = null;

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

function createAddNewNode(cy, curNodeSet) {
    cy.on('tap', 'node', function (evt) {
        var v = evt.target.id();
        var needAddVert = [];
        var needAddEdge = [];
        var newNode = new Set();
        findConnectedVertex(v, needAddVert, needAddEdge, newNode, area_size, min_contig_len, isGoodEdge, curNodeSet);

        for (g = 0; g < needAddVert.length; ++g) {
            u = needAddVert[g];
            var yc = calcYforV(u, area_size, min_contig_len, isGoodEdge, newNode, curNodeSet, cy, 500);

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
                    weight: Math.log(scaffoldgraph.edges[eid].weight) + 1,
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

function DrawGraphVis(nodes_to_draw, edges_to_draw) {
    var nodeslist = [];
    var edgeslist = [];
    var i = 0;
    var j = 0;
    for (i=0; i < nodes_to_draw.length; i++) {
        j = nodes_to_draw[i];
        var label = createLabelForNode(j);
        nodeslist.push({id: scaffoldgraph.nodes[j].id, label: label});
    }

    // create an array with nodes
    var nodes = new vis.DataSet(nodeslist);

    for (i=0; i < edges_to_draw.length; i++) {
        j = edges_to_draw[i];
        label = createLabelForEdge(j);

        edgeslist.push({from: scaffoldgraph.edges[j].from, to: scaffoldgraph.edges[j].to, label : label, arrows: 'to', color:{color: scaffoldgraph.libs[scaffoldgraph.edges[j].lib].color}});
    }

    // create an array with edges
    var edges = new vis.DataSet(edgeslist);

    // create a network
    var container = document.getElementById('mainpanel');

    // provide the data in the vis format
    var data = {
        nodes: nodes,
        edges: edges
    };
    var options = {
        nodes : {
            shape: 'dot',
            size: 7
        },
        layout:{
            randomSeed:34
        },
        physics: {
            forceAtlas2Based: {
                gravitationalConstant: -26,
                centralGravity: 0.005,
                springLength: 230,
                springConstant: 0.18
            },
            maxVelocity: 146,
            solver: 'forceAtlas2Based',
            timestep: 0.35,
            stabilization: {
                enabled:true,
                iterations:2000,
                updateInterval:25
            }
        }

        //,
        // interaction: {
        //     hideEdgesOnDrag: true,
        //     tooltipDelay: 200
        // },
        // physics: false
    };

    // initialize your network!
    var network = new vis.Network(container, data, options);
}

function DrawGraphViva(nodes_to_draw, edges_to_draw) {
    if (graph === null) {
        graph = Viva.Graph.graph();
    }
    graph.clear();

    for (g = 0; g < nodes_to_draw.length; ++g) {
        graph.addNode(nodes_to_draw[g]);
    }

    var graphics = Viva.Graph.View.svgGraphics(),
        nodeSize = 5;


    graphics.node(function(node) {
        return Viva.Graph.svg('image')
            .attr('width', nodeSize)
            .attr('height', nodeSize)
            .link('https://secure.gravatar.com/avatar/' + node.data);
    }).placeNode(function(nodeUI, pos) {
        nodeUI.attr('x', pos.x - nodeSize / 2).attr('y', pos.y - nodeSize / 2);
    });

    // To render an arrow we have to address two problems:
    //  1. Links should start/stop at node's bounding box, not at the node center.
    //  2. Render an arrow shape at the end of the link.

    // Rendering arrow shape is achieved by using SVG markers, part of the SVG
    // standard: http://www.w3.org/TR/SVG/painting.html#Markers
    var createMarker = function(id) {
            return Viva.Graph.svg('marker')
                .attr('id', id)
                .attr('viewBox', "0 0 10 10")
                .attr('refX', "10")
                .attr('refY', "5")
                .attr('markerUnits', "strokeWidth")
                .attr('markerWidth', "10")
                .attr('markerHeight', "5")
                .attr('orient', "auto");
        },

        marker = createMarker('Triangle');
    marker.append('path').attr('d', 'M 0 0 L 10 5 L 0 10 z');

    // Marker should be defined only once in <defs> child element of root <svg> element:
    var defs = graphics.getSvgRoot().append('defs');
    defs.append(marker);

    var geom = Viva.Graph.geom();

    graphics.link(function(link){
        // Notice the Triangle marker-end attribe:
        return Viva.Graph.svg('path')
            .attr('stroke', 'gray')
            .attr('marker-end', 'url(#Triangle)');
    }).placeLink(function(linkUI, fromPos, toPos) {
        // Here we should take care about
        //  "Links should start/stop at node's bounding box, not at the node center."

        // For rectangular nodes Viva.Graph.geom() provides efficient way to find
        // an intersection point between segment and rectangle
        var toNodeSize = nodeSize,
            fromNodeSize = nodeSize;

        var from = geom.intersectRect(
            // rectangle:
            fromPos.x - fromNodeSize / 2, // left
            fromPos.y - fromNodeSize / 2, // top
            fromPos.x + fromNodeSize / 2, // right
            fromPos.y + fromNodeSize / 2, // bottom
            // segment:
            fromPos.x, fromPos.y, toPos.x, toPos.y)
            || fromPos; // if no intersection found - return center of the node

        var to = geom.intersectRect(
            // rectangle:
            toPos.x - toNodeSize / 2, // left
            toPos.y - toNodeSize / 2, // top
            toPos.x + toNodeSize / 2, // right
            toPos.y + toNodeSize / 2, // bottom
            // segment:
            toPos.x, toPos.y, fromPos.x, fromPos.y)
            || toPos; // if no intersection found - return center of the node

        var data = 'M' + from.x + ',' + from.y +
            'L' + to.x + ',' + to.y;

        linkUI.attr("d", data);
    });

    for (g = 0; g < edges_to_draw.length; ++g) {
        graph.addLink(scaffoldgraph.edges[edges_to_draw[g]].from, scaffoldgraph.edges[edges_to_draw[g]].to);
    }

    var layout = Viva.Graph.Layout.forceDirected(graph, {
        springLength : 20,
        springCoeff : 0.0005,
        dragCoeff : 0.02,
        gravity : -1.2
    });

    var renderer = Viva.Graph.View.renderer(graph, {
        container: document.getElementById('mainpanel'),
        layout : layout,
        graphics: graphics
    });
    renderer.run();
}

function DrawGraphDagre(nodes_to_draw, edges_to_draw) {
    var dagre = require("dagre");

    // Create a new directed graph
    var g = new dagre.graphlib.Graph();

// Set an object for the graph label
    g.setGraph({});

// Default to assigning a new object as a label for each new edge.
    g.setDefaultEdgeLabel(function() { return {}; });

// Add nodes to the graph. The first argument is the node id. The second is
// metadata about the node. In this case we're going to add labels to each of
// our nodes.
    g.setNode("kspacey",    { label: "Kevin Spacey",  width: 144, height: 100 });
    g.setNode("swilliams",  { label: "Saul Williams", width: 160, height: 100 });
    g.setNode("bpitt",      { label: "Brad Pitt",     width: 108, height: 100 });
    g.setNode("hford",      { label: "Harrison Ford", width: 168, height: 100 });
    g.setNode("lwilson",    { label: "Luke Wilson",   width: 144, height: 100 });
    g.setNode("kbacon",     { label: "Kevin Bacon",   width: 121, height: 100 });

// Add edges to the graph.
    g.setEdge("kspacey",   "swilliams");
    g.setEdge("swilliams", "kbacon");
    g.setEdge("bpitt",     "kbacon");
    g.setEdge("hford",     "lwilson");
    g.setEdge("lwilson",   "kbacon");

    dagre.layout(g);
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

function splitOnParts(nodes_to_draw, edges_to_draw) {
    var g = new Map();
    var color = new Map();
    for (var i=0; i < nodes_to_draw.length; ++i) {
        g.set(nodes_to_draw[i], []);
        color.set(nodes_to_draw[i], -1);
    }

    for (i=0; i < edges_to_draw.length; ++i) {
        g.get(scaffoldgraph.edges[edges_to_draw[i]].from).push(scaffoldgraph.edges[edges_to_draw[i]].to);
        g.get(scaffoldgraph.edges[edges_to_draw[i]].to).push(scaffoldgraph.edges[edges_to_draw[i]].from);
    }

    var curc = 0;
    for (i=0; i < nodes_to_draw.length; ++i) {
        if (color.get(nodes_to_draw[i]) === -1) {
            findComponent(nodes_to_draw[i], g, color, curc);
            ++curc;
        }
    }

    nodes_set = [];
    edges_set = [];
    for (i=0; i < curc; ++i) {
        nodes_set.push([]);
        edges_set.push([]);
    }

    for(i=0; i < nodes_to_draw.length; ++i) {
        nodes_set[color.get(nodes_to_draw[i])].push(nodes_to_draw[i]);
    }

    for(i=0; i < edges_to_draw.length; ++i) {
        edges_set[color.get(scaffoldgraph.edges[edges_to_draw[i]].from)].push(edges_to_draw[i]);
    }
}

function handleFilterButton() {
    defZoom = 100;
    special_nodes.clear();
    special_edges.clear();

    var opt = document.getElementById("select_show_type").value;
    var min_edge_weight = [];
    min_contig_len = document.getElementById("min_contig_len").value;

    for (var i = 0; i < scaffoldgraph.libs.length; ++i) {
        min_edge_weight.push(document.getElementById("min_w_" + scaffoldgraph.libs[i].name).value);
    }

    isGoodEdge = function (e) {
        if (min_edge_weight[scaffoldgraph.edges[e].lib] > scaffoldgraph.edges[e].weight) {
            return false;
        }

        if (scaffoldgraph.nodes[scaffoldgraph.edges[e].from].len < min_contig_len) {
            return false;
        }

        if (scaffoldgraph.nodes[scaffoldgraph.edges[e].to].len < min_contig_len) {
            return false;
        }

        return true;
    };

    nodes_to_draw = [];
    edges_to_draw = [];

    if (opt == "vertices_local_area") {
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
    } else if (opt=="edges_local_area") {
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
    } else if (opt=="diff in libs") {
        area_size = document.getElementById("area_size").value;
        handleDiffInLibsFilter(area_size, min_contig_len, isGoodEdge);
    } else if (opt=="full graph") {
        for (i=0; i < scaffoldgraph.nodes.length; ++i) {
            if (scaffoldgraph.nodes[i].len >= min_contig_len) {
                nodes_to_draw.push(scaffoldgraph.nodes[i].id);
            }
        }

        for (i=0; i < scaffoldgraph.edges.length; ++i) {
            if (isGoodEdge(i)) {
                edges_to_draw.push(scaffoldgraph.edges[i].id);
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
    } else if (opt=="along chromosoms") {
        area_size = document.getElementById("area_size").value;
        handleAlongChromosomesFilter();
    } else if (opt=="scaffolds") {
        area_size = document.getElementById("area_size").value;
        handleScaffoldsFilter(document.getElementById("select_scaff_lib").value, area_size, min_contig_len, isGoodEdge);
    }
}

function InitLibTable() {
    var table = document.getElementById("lib_table");

    for (var i=0; i < scaffoldgraph.libs.length; ++i) {
        var tr = document.createElement("tr");

        var td_id = document.createElement("td");
        var lib_id = document.createElement("p");
        lib_id.appendChild(document.createTextNode("l" + scaffoldgraph.libs[i].id));
        td_id.align="center";
        td_id.appendChild(lib_id);

        var td_type = document.createElement("td");
        var lib_type = document.createElement("p");
        lib_type.appendChild(document.createTextNode((scaffoldgraph.libs[i].type).replace(/_/g, " ")));
        td_type.appendChild(lib_type);
        td_type.align="center";


        var td_name = document.createElement("td");
        var lib_name = document.createElement("p");
        lib_name.style.color = scaffoldgraph.libs[i].color;
        lib_name.appendChild(document.createTextNode(scaffoldgraph.libs[i].name));
        lib_name.id = "color" + scaffoldgraph.libs[i].name;
        td_name.align="center";
        td_name.appendChild(lib_name);

        var td_min_edge_weight = document.createElement("td");
        var input_weight = document.createElement("input");
        input_weight.type = "number";
        input_weight.min = 0;
        input_weight.size = 1;
        input_weight.value = getDefaultWeight(i);
        input_weight.id = "min_w_" + scaffoldgraph.libs[i].name;
        td_min_edge_weight.align="center";
        td_min_edge_weight.appendChild(input_weight);

        tr.appendChild(td_id);
        tr.appendChild(td_type);
        tr.appendChild(td_name);
        tr.appendChild(td_min_edge_weight);
        table.appendChild(tr);
    }
}

function InitAlignmentsForNodes() {
    for (var i=0; i < chromosomes.length; ++i) {
        for (var j=0; j < chromosomes[i].alignments.length; ++j) {
            scaffoldgraph.nodes[chromosomes[i].alignments[j].node_id].alignments.push(chromosomes[i].alignments[j]);
        }
    }
}

function disableNot(checkBox) {
    if (checkBox.checked) {
        document.getElementById("checkbox_not_present_" + checkBox.id.substring(17)).disabled = true;
        document.getElementById("notPresName" + checkBox.id.substring(17)).style.color = "#aaa";
    } else {
        document.getElementById("checkbox_not_present_" + checkBox.id.substring(17)).disabled = false;
        document.getElementById("notPresName" + checkBox.id.substring(17)).style.color = "#4F4F4F";
    }
}


function disable(checkBox) {
    if (checkBox.checked) {
        document.getElementById("checkbox_present_" + checkBox.id.substring(21)).disabled = true;
        document.getElementById("presName" + checkBox.id.substring(21)).style.color = "#aaa";
    } else {
        document.getElementById("checkbox_present_" + checkBox.id.substring(21)).disabled = false;
        document.getElementById("presName" + checkBox.id.substring(21)).style.color = "#4F4F4F";
    }
}

InitLibTable();
InitAlignmentsForNodes();
setupAutocompleteSearch();

document.getElementById("filter_button").addEventListener("click", handleFilterButton);
document.getElementById("select_show_type").addEventListener("change", function() {
    if(document.getElementById("select_show_type").value == "full graph") {
        document.getElementById("change_block").innerHTML = "";
    } else if (document.getElementById("select_show_type").value == "vertices_local_area" || document.getElementById("select_show_type").value == "edges_local_area") {
        document.getElementById("change_block").innerHTML = "<div class=\"block\">\n" +
            "                    <p>Area size:<br/>\n" +
            "                        <input type=\"number\" min=\"0\" id=\"area_size\" value=2>\n" +
            "                    </p>\n" +
            "                </div>\n" +
            "                <div class=\"block\">\n" +
            "                    <textarea rows=\"6\" id=\"vertext\"></textarea>\n" +
            "                </div>";
    } else if (document.getElementById("select_show_type").value == "diff in libs" ) {
        var html_code = "<div class=\"block\">\n" +
            "                    <p>Area size:<br/>\n" +
            "                        <input type=\"number\" min=\"0\" id=\"area_size\" value=2>\n" +
            "                    </p>\n" +
            "                </div>\n" +
            "<div class=\"one_line_block\">\n" +
            "                        <label class=\"container\">\n" +
            "                            <p>Wrong</p>\n" +
            "                            <input type=\"checkbox\" checked=\"true\" value=\"Wrong\" id=\"checkbox_wrong\">\n" +
            "                            <span class=\"checkmark\"></span>\n" +
            "                        </label>\n" +
            "                    </div>\n" +
            "\n" +
            "                    <div class=\"one_line_block\">\n" +
            "                        <label class=\"container\">\n" +
            "                            <p>Correct</p>\n" +
            "                            <input type=\"checkbox\" checked=\"true\" value=\"Correct\" id=\"checkbox_correct\">\n" +
            "                            <span class=\"checkmark\"></span>\n" +
            "                        </label>\n" +
            "                    </div>\n" +
            "                    <br/>\n";
        var present_block = "<div class=\"one_line_block\" id=\"present_block\">\n" +
            "                        Presnt:\n";

        var not_present_block = "<div class=\"one_line_block\" id=\"not_present_block\">\n" +
            "                        Not Presnt:\n";

        for (var i=0; i < scaffoldgraph.libs.length; ++i) {
            var lib_name = scaffoldgraph.libs[i].name;
            present_block +=
                "                        <label class=\"container\">\n" +
                "                            <p id=\"presName" + i.toString() +"\">" + lib_name + "</p>\n" +
                "                            <input type=\"checkbox\" value=\"" + lib_name + "\" id=\"checkbox_present_"+i.toString() + "\" onchange=\"disableNot(this)\">\n" +
                "                            <span class=\"checkmark\"></span>\n" +
                "                        </label>\n";

            not_present_block +=
                "                        <label class=\"container\">\n" +
                "                            <p id=\"notPresName" + i.toString() +"\">" + lib_name + "</p>\n" +
                "                            <input type=\"checkbox\" value=\"" + lib_name + "\" id=\"checkbox_not_present_"+i.toString() + "\" onchange=\"disable(this)\" >\n" +
                "                            <span class=\"checkmark\"></span>\n" +
                "                        </label>\n";
        }

        present_block += "</div>\n";
        not_present_block += "</div>\n";
        document.getElementById("change_block").innerHTML = html_code + present_block + not_present_block;
    } else if (document.getElementById("select_show_type").value == "along chromosoms") {
        document.getElementById("change_block").innerHTML = "<div class=\"block\">\n" +
            "                    <p>Area size:<br/>\n" +
            "                        <input type=\"number\" min=\"0\" id=\"area_size\" value=2>\n" +
            "                    </p>\n" +
            "                </div>";
    } else if (document.getElementById("select_show_type").value == "scaffolds") {
        document.getElementById("change_block").innerHTML = "<div class=\"block\">\n" +
            "                    <p>Area size:<br/>\n" +
            "                        <input type=\"number\" min=\"0\" id=\"area_size\" value=2>\n" +
            "                    </p>\n" +
            "                    <p>Scaffold lib:<br/>\n" +
            "                        <div class=\"styled-select\">\n" +
            "                           <select id=\"select_scaff_lib\">\n" +
            "                           </select>\n" +
            "                       </div>" +
            "                    </p>\n" +
            "                </div>";

        for (i=0; i < scaffoldgraph.libs.length; ++i) {
            if (scaffoldgraph.libs[i].type == 'SCAFF') {
                document.getElementById("select_scaff_lib").innerHTML += "<option value=\"" + scaffoldgraph.libs[i].name + "\">" + scaffoldgraph.libs[i].name + "</option>\n";
            }
        }
    }
});

window.addEventListener("keyup", function (evt) {
    if (evt.altKey) {
        var code = (evt.keyCode || evt.which);
        if (code === 187) {
            if (typeof cy !== 'undefined') {
                cy.zoom(cy.zoom() * 1.5);
            }
        } else if (code === 189) {
            if (typeof cy !== 'undefined') {
                cy.zoom(cy.zoom() / 1.5);
            }
        }
    }
});
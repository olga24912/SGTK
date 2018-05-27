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
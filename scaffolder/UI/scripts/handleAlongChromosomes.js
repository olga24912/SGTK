var defZoom = 100;

function createCoordinates(chr, cy) {
    var cur_zoom = cy.zoom();
    var cur_coord = cy.extent();

    var delta = 1;
    while (cur_zoom * delta < 50) {
        delta *= 10;
    }

    var start_pos = Math.max(delta, (Math.floor(cur_coord.y1 / delta) - 5) * delta);

    cy.remove(cy.$('#start'));
    cy.add({
        group: "nodes",
        data: {
            id: "start",
            label: "Start: 0",
            len: 10 / cy.zoom(),
            color: '#ff0000',
            width: 10 / cy.zoom(),
            faveShape: 'rectangle'
        },
        position: {
            x: 0,
            y: 0
        }
    });
    var ypos = cy.$('#start').renderedPosition().y;
    cy.$('#start').renderedPosition({
        x: 10,
        y: ypos
    });
    cy.$('#start').style({"font-size": 20 / cy.zoom(), "text-valign": "center", "text-halign": "right"});
    cy.$('#start').lock();

    cy.remove(cy.$('#end'));
    cy.add({
        group: "nodes",
        data: {
            id: "end",
            label: "End: " + chromosomes[chr].len,
            len: 10 / cy.zoom(),
            color: '#ff0000',
            width: 10 / cy.zoom(),
            faveShape: 'rectangle'
        },
        position: {
            x: 0,
            y: chromosomes[chr].len/defZoom
        }
    });
    ypos = cy.$('#end').renderedPosition().y;
    cy.$('#end').renderedPosition({
        x: 10,
        y: ypos
    });
    cy.$('#end').style({"font-size": 20 / cy.zoom(), "text-valign": "center", "text-halign": "right"});
    cy.$('#end').lock();


    for (var i = 0;  i < 50; ++i) {
        cy.remove(cy.$('#chrcoord' + i));
        if (start_pos + delta * i < chromosomes[chr].len/defZoom) {
            cy.add({
                group: "nodes",
                data: {
                    id: "chrcoord" + i,
                    label: (start_pos + delta * i)*defZoom,
                    len: 1 / cy.zoom(),
                    color: '#ffa500',
                    width: 1 / cy.zoom(),
                    faveShape: 'rectangle'
                },
                position: {
                    x: 0,
                    y: start_pos + delta * i
                }
            });
            ypos = cy.$('#chrcoord' + i).renderedPosition().y;
            cy.$('#chrcoord' + i).renderedPosition({
                x: 10,
                y: ypos
            });
            cy.$('#chrcoord' + i).style({"font-size": 20 / cy.zoom(), "text-valign": "center", "text-halign": "right"});
            cy.$('#chrcoord' + i).lock();
        }
    }
}

function calcY(curv, ypos, sumw) {
    return ypos/sumw;
}

function findNodeAroundChr(inode, area_size, min_contig_len, isGoodEdge) {
    var ypos = {};
    var newvert = [];
    var sumw = {};
    var rank = {};

    var used_id = new Set();
    var que = [];
    for (var i = 0; i < inode.length; ++i) {
        var v = inode[i];

        used_id.add(v.id);
        rank[v.id] = 0;
        que.push(v);
    }

    var bg = 0;
    while (bg < que.length) {
        v = que[bg];
        bg += 1;

        var curd = rank[v.id];
        var curv = v.id;
        for (i = 0; i < scaffoldgraph.g[curv].length; ++i) {
            var curedge = scaffoldgraph.g[curv][i];
            if (isGoodEdge(curedge.id)) {
                var curu = curedge.to;
                if (scaffoldgraph.nodes[curu].len >= min_contig_len) {
                    if ((!used_id.has(curu)) && curd < area_size) {
                        rank[curu] = curd + 1;
                        used_id.add(curu);
                        newvert.push(curu);
                        que.push({id: curu});

                        ypos[curu] = 0;
                        sumw[curu] = 0;
                    }

                    if (rank[curu] >= rank[curv]) {
                        var yc = 0;
                        if (curd === 0) {
                            yc = v.ce;
                        } else {
                            yc = calcY(curv, ypos[curv], sumw[curv]);
                        }
                        yc = yc + Math.random() * 1000;

                        ypos[curu] += curedge.weight * yc;
                        sumw[curu] += curedge.weight;
                    }

                    edges_to_draw.push(curedge.id);
                }
            }
        }


        for (i = 0; i < scaffoldgraph.gr[curv].length; ++i) {
            curedge = scaffoldgraph.gr[curv][i];
            if (isGoodEdge(curedge.id)) {
                curu = curedge.from;
                if (scaffoldgraph.nodes[curu].len >= min_contig_len) {
                    if ((!used_id.has(curu)) && curd < area_size) {
                        rank[curu] = curd + 1;
                        used_id.add(curu);
                        newvert.push(curu);
                        que.push({id: curu});

                        ypos[curu] = 0;
                        sumw[curu] = 0;
                    }

                    if (rank[curu] > rank[curv]) {
                        yc = 0;
                        if (curd === 0) {
                            yc = v.cb;
                        } else {
                            yc = calcY(curv, ypos[curv], sumw[curv]);
                        }
                        yc = yc  - Math.random() * 1000;

                        ypos[curu] += curedge.weight * yc;
                        sumw[curu] += curedge.weight;
                    }
                }
            }
        }
    }

    var res = [];

    for (i = 0; i < newvert.length; ++i) {
        res.push({id: newvert[i], rank: rank[newvert[i]], y: calcY(newvert[i], ypos[newvert[i]], sumw[newvert[i]])});
    }

    return res;
}

function getWeight(e) {
    var curW = scaffoldgraph.edges[e].weight;

    if (scaffoldgraph.libs[scaffoldgraph.edges[e].lib].type == 'SCAFF') {
        curW = 1024;
    }
    return curW
}

function calcYforV(u, area_size, min_contig_len, isGoodEdge, newNode, curNodeSet, cy, rcoef) {
    var ypos = 0;
    var sumw = 0;
    for (var h = 0; h < scaffoldgraph.g[u].length; ++h) {
        if (isGoodEdge(scaffoldgraph.g[u][h].id)) {
            if (scaffoldgraph.nodes[scaffoldgraph.g[u][h].to].len >= min_contig_len) {
                if (curNodeSet.has(scaffoldgraph.g[u][h].to)) {
                    if (!newNode.has(scaffoldgraph.g[u][h].to)) {
                        var v = scaffoldgraph.g[u][h].to;
                        var curedge = scaffoldgraph.g[u][h];
                        var yc = cy.$('#' + v).position().y - rcoef * Math.random();
                        ypos += curedge.weight * yc;
                        sumw += curedge.weight;
                    }
                }
            }
        }
    }


    for (h = 0; h < scaffoldgraph.gr[u].length; ++h) {
        if (isGoodEdge(scaffoldgraph.gr[u][h].id)) {
            if (scaffoldgraph.nodes[scaffoldgraph.gr[u][h].from].len >= min_contig_len) {
                if (curNodeSet.has(scaffoldgraph.gr[u][h].from)) {
                    if (!newNode.has(scaffoldgraph.gr[u][h].from)) {
                        v = scaffoldgraph.gr[u][h].from;
                        curedge = scaffoldgraph.gr[u][h];
                        yc = cy.$('#' + v).position().y + rcoef * Math.random();
                        ypos += curedge.weight * yc;
                        sumw += curedge.weight;
                    }
                }
            }
        }
    }

    return calcY(u, ypos, sumw);
}

function createNewVerAlongChr(cy, area_size, min_contig_len, isGoodEdge, curNodeSet) {
    cy.on('tap', 'node', function (evt) {
        var newNode = new Set();
        var v = evt.target.id();
        var needAddVert = [];
        var needAddEdge = [];

        findConnectedVertex(v, needAddVert, needAddEdge, newNode, area_size, min_contig_len, isGoodEdge, curNodeSet);

        for (g = 0; g < needAddVert.length; ++g) {
            u = needAddVert[g];
            var xc = evt.target.data('rank') + 1;
            var yc = calcYforV(u, area_size, min_contig_len, isGoodEdge, newNode, curNodeSet, cy, 1000);

            cy.add({
                group: "nodes",
                data: {
                    id: u,
                    label: createLabelForNode(u),
                    len: 20 * Math.log(scaffoldgraph.nodes[u].len),
                    width: 20 * Math.log(scaffoldgraph.nodes[u].len),
                    color: genColorNode(u),
                    rank: xc,
                    faveShape: 'ellipse'
                },
                position: {
                    x: 2000*xc + Math.random()*100,
                    y: yc
                }
            });
        }

        for (g = 0; g < needAddEdge.length; ++g) {
            var eid = needAddEdge[g].id;
            edges_to_draw.push(eid);
            cy.add({
                group: "edges",
                data: {
                    id: "e" + eid.toString(),
                    source: scaffoldgraph.edges[eid].from,
                    target: scaffoldgraph.edges[eid].to,
                    label: createLabelForEdge(eid),
                    faveColor: scaffoldgraph.libs[scaffoldgraph.edges[eid].lib].color,
                    weight: Math.log(getWeight(eid)) + 1,
                    curveStyle: "bezier",
                    controlPointDistances: 100
                }
            });
        }

        for (g = 0; g < nodes_to_draw.length; ++g) {
            if (hasOtherEdges(nodes_to_draw[g], curNodeSet)) {
                cy.$('#' + nodes_to_draw[g]).data('faveShape', 'triangle');
            } else {
                cy.$('#' + nodes_to_draw[g]).data('faveShape', 'ellipse');
            }
        }
        createTapInfo(cy);
    });
}


function updateZooming() {

}

function findContigs(cy, chr, inode, posx, posmin, posmax, curNodeSet) {
    for (i = 0; i < chromosomes[chr].alignments.length; ++i) {
        var curalig = chromosomes[chr].alignments[i];
        var vid = curalig.node_id;
        posx[vid] = Math.random() * 100;
        posmin[vid] = curalig.coordb / defZoom;
        posmax[vid] = curalig.coorde / defZoom;
        inode.push({id: vid, cb: curalig.coordb / defZoom, ce: curalig.coorde / defZoom});
        special_nodes.add(curalig.node_id);
        curNodeSet.add(curalig.node_id);
    }
}

function addContigs(cy, inode, posx, posmin, posmax) {
    for (i = 0; i < inode.length; ++i) {
        var vid = inode[i].id;
        cy.add({
            group: "nodes",
            data: {
                id: vid,
                label: createLabelForNode(vid),
                len: inode[i].ce - inode[i].cb,
                color: genColorNode(vid),
                width: 10,
                rank: 0,
                ymin: inode[i].cb,
                ymax: inode[i].ce,
                faveShape: 'rectangle'
            },
            position: {
                x: posx[vid],
                y: (posmin[vid] + posmax[vid])/2
            }
        });
    }
}

function addOtherNodes(cy, curNodeSet, vert_to_draw) {
    for (var g = 0; g < vert_to_draw.length; ++g) {
        curNodeSet.add(vert_to_draw[g].id);
    }
    for (g = 0; g < vert_to_draw.length; ++g) {
        nodes_to_draw.push(vert_to_draw[g].id);
        var nall = 'ellipse';
        if (hasOtherEdges(vert_to_draw[g].id, curNodeSet)) {
            nall = 'triangle';
        }

        cy.add({
            group: "nodes",
            data: {
                id: vert_to_draw[g].id,
                label: createLabelForNode(vert_to_draw[g].id),
                len: 20 * Math.log(scaffoldgraph.nodes[vert_to_draw[g].id].len),
                width: 20 * Math.log(scaffoldgraph.nodes[vert_to_draw[g].id].len),
                color: genColorNode(vert_to_draw[g].id),
                rank: vert_to_draw[g].rank,
                faveShape: nall
            },
            position: {
                x: 2000 * vert_to_draw[g].rank + Math.random() * 100,
                y: vert_to_draw[g].y
            }
        });
    }
}

function addEdges(cy) {
    alert(special_nodes.size);
    alert(special_nodes.has(71));
    alert('71' in special_nodes);

    for (g = 0; g < edges_to_draw.length; ++g) {
        if (!(special_nodes.has(scaffoldgraph.edges[edges_to_draw[g]].from)) ||
            !(special_nodes.has(scaffoldgraph.edges[edges_to_draw[g]].to))) {
            cy.add({
                group: "edges",
                data: {
                    id: "e" + edges_to_draw[g].toString(),
                    source: scaffoldgraph.edges[edges_to_draw[g]].from,
                    target: scaffoldgraph.edges[edges_to_draw[g]].to,
                    label: createLabelForEdge(edges_to_draw[g]),
                    faveColor: scaffoldgraph.libs[scaffoldgraph.edges[edges_to_draw[g]].lib].color,
                    weight: Math.log(getWeight(edges_to_draw[g])) + 1,
                    curveStyle: "bezier",
                    controlPointDistances: 100
                }
            });
        } else {
            cy.add({
                group: "edges",
                data: {
                    id: "e" + edges_to_draw[g].toString(),
                    source: scaffoldgraph.edges[edges_to_draw[g]].from,
                    target: scaffoldgraph.edges[edges_to_draw[g]].to,
                    label: createLabelForEdge(edges_to_draw[g]),
                    faveColor: scaffoldgraph.libs[scaffoldgraph.edges[edges_to_draw[g]].lib].color,
                    weight: Math.log(getWeight(edges_to_draw[g])) + 1,
                    curveStyle: "unbundled-bezier",
                    controlPointDistances: 100 + Math.random() * 10
                }
            });
        }
    }
}

function zoomAndFit() {

}

function createGraph(chr, cy, curNodeSet) {
    alert("draw");
    cy.elements().remove();
    special_nodes.clear();

    updateZooming();

    var inode = [];
    var posx = {};
    var posmin = {};
    var posmax = {};

    nodes_to_draw = [];
    edges_to_draw = [];

    findContigs(cy, chr, inode, posx, posmin, posmax, curNodeSet);
    alert("find Contig");
    addContigs(cy, inode, posx, posmin, posmax);
    alert("add Contigs");

    var vert_to_draw = findNodeAroundChr(inode, area_size, min_contig_len, isGoodEdge);
    addOtherNodes(cy, curNodeSet, vert_to_draw);
    addEdges(cy);
    zoomAndFit();
    createCoordinates(chr, cy);
}



function drawAlongChromosome(chr) {
    var curNodeSet = new Set();
    cy = cytoscape({
        container: document.getElementById('mainpanel'),

        boxSelectionEnabled: false,
        autounselectify: true,
        maxZoom: 25,
        minZoom: 0.005,

        layout: {
            name: 'preset'
        },


        ready: function () {
            window.cy = this;
        },

        style: cytoscape.stylesheet()
            .selector('node')
            .css({
                'shape': 'data(faveShape)',
                'content': 'data(label)',
                'color': '#2A4986',
                'height': 'data(len)',
                'width': 'data(width)',
                'background-color': 'data(color)'
            })
            .selector('edge')
            .css({
                'curve-style': 'data(curveStyle)',
                "control-point-distances": 'data(controlPointDistances)',
                "control-point-weights": 0.5,
                "arrow-scale": 7,
                'target-arrow-shape': 'triangle',
                'line-color': 'data(faveColor)',
                'target-arrow-color': 'data(faveColor)',
                'width': 'data(weight)',
                'content': 'data(label)'
            })
    });

    createTapInfo(cy);
    createNewVerAlongChr(cy, area_size, min_contig_len, isGoodEdge, curNodeSet);

    cy.on('zoom', function () {
        createCoordinates(chr, cy);
        var nodeWidth = 10 / cy.zoom();

        var y1 = cy.extent().y1;
        var y2 = cy.extent().y2;

        var selectstr = "[faveShape='rectangle'][ymin <= " + y2 + "][ymax >= " + y1 + "]";

        cy.nodes(selectstr).data('width', nodeWidth);
        for (var i = 0; i < edges_to_draw.length; ++i) {
            var e = edges_to_draw[i];
            var w = (Math.log(getWeight(e)) + 1) * 4/ Math.log((cy.zoom() * defZoom));
            var cpd = cy.getElementById("e" + e.toString()).data("curveStyle");
            if (cpd === "unbundled-bezier") {
                bg = scaffoldgraph.edges[e].from;
                ed = scaffoldgraph.edges[e].to;

                yminbg = cy.getElementById(bg).data('ymin');
                ymaxbg = cy.getElementById(bg).data('ymax');
                ymined = cy.getElementById(ed).data('ymin');
                ymaxed = cy.getElementById(ed).data('ymax');

                minDif = Math.min(Math.min(Math.abs(yminbg - ymaxed), Math.abs(yminbg - ymined)), Math.min(Math.abs(ymaxbg - ymaxed), Math.abs(ymaxbg - ymined)));

                minDif = Math.log(minDif);
                cpd = minDif * 10/ cy.zoom() + (Math.random() - 0.5) * minDif / cy.zoom();
                cy.getElementById("e" + e.toString()).data("controlPointDistances", cpd);
            }
            cy.getElementById("e" + e.toString()).data('weight', w);
        }
    });
    cy.on('pan', function() {
        createCoordinates(chr, cy);

        var nodeWidth = 10 / cy.zoom();
        var y1 = cy.extent().y1;
        var y2 = cy.extent().y2;

        var selectstr = "[faveShape='rectangle'][ymin <= " + y2 + "][ymax >= " + y1 + "]";
        cy.nodes(selectstr).data('width', nodeWidth);
    });


    cy.zoom(1);
    //if (inode.length > 0) {
    //    cy.fit(cy.$('#' + inode[0].id));
    //}
    createGraph(chr, cy, curNodeSet);
}

function handleAlongChromosomesFilter() {
    createComponentShowList(drawAlongChromosome, function(i) {
        return chromosomes[i].name;
    }, function(i) {
        return "Chromosome " + chromosomes[i].name;
    }, chromosomes.length);
}
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
            len: 10 / cur_zoom,
            color: '#ff0000',
            width: 10 / cur_zoom,
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
    cy.$('#start').style({"font-size": 20 / cur_zoom, "text-valign": "center", "text-halign": "right"});
    cy.$('#start').lock();

    cy.remove(cy.$('#end'));
    cy.add({
        group: "nodes",
        data: {
            id: "end",
            label: "End: " + chromosomes[chr].len,
            len: 10 / cur_zoom,
            color: '#ff0000',
            width: 10 / cur_zoom,
            faveShape: 'rectangle'
        },
        position: {
            x: 0,
            y: chromosomes[chr].len
        }
    });
    ypos = cy.$('#end').renderedPosition().y;
    cy.$('#end').renderedPosition({
        x: 10,
        y: ypos
    });
    cy.$('#end').style({"font-size": 20 / cur_zoom, "text-valign": "center", "text-halign": "right"});
    cy.$('#end').lock();


    for (var i = 0; i < 20; ++i) {
        cy.remove(cy.$('#chrcoord' + i));
        if (start_pos + delta * i < chromosomes[chr].len) {
            cy.add({
                group: "nodes",
                data: {
                    id: "chrcoord" + i,
                    label: start_pos + delta * i,
                    len: 1 / cur_zoom,
                    color: '#ffa500',
                    width: 1 / cur_zoom,
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
            cy.$('#chrcoord' + i).style({"font-size": 20 / cur_zoom, "text-valign": "center", "text-halign": "right"});
            cy.$('#chrcoord' + i).lock();
        }
    }
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

    cy.on('tap', 'edge', function (evt) {
        this.qtip({
            content: function () {
                return 'Expression: '
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
    });

    cy.on('tap', 'node', function (evt) {
        var v = evt.target.id();
        var needAddVert = [];

        for (var h = 0; h < scaffoldgraph.g[v].length; ++h) {
            if (isGoodEdge(scaffoldgraph.g[v][h].id)) {
                if (scaffoldgraph.nodes[scaffoldgraph.g[v][h].to].len >= min_contig_len) {
                    if (!curNodeSet.has(scaffoldgraph.g[v][h].to)) {
                        curNodeSet.add(scaffoldgraph.g[v][h].to);
                        needAddVert.push(scaffoldgraph.g[v][h].to);
                        nodes_to_draw.push(scaffoldgraph.g[v][h].to);
                    }
                }
            }
        }

        for (h = 0; h < scaffoldgraph.gr[v].length; ++h) {
            if (isGoodEdge(scaffoldgraph.gr[v][h].id)) {
                if (scaffoldgraph.nodes[scaffoldgraph.gr[v][h].from].len >= min_contig_len) {
                    if (!curNodeSet.has(scaffoldgraph.gr[v][h].from)) {
                        curNodeSet.add(scaffoldgraph.gr[v][h].from);
                        needAddVert.push(scaffoldgraph.gr[v][h].from);
                        nodes_to_draw.push(scaffoldgraph.gr[v][h].from);
                    }
                }
            }
        }

        var needAddEdge = [];


        for (var g = 0; g < needAddVert.length; ++g) {
            var u = needAddVert[g];

            for (h = 0; h < scaffoldgraph.g[u].length; ++h) {
                if (isGoodEdge(scaffoldgraph.g[u][h].id)) {
                    if (curNodeSet.has(scaffoldgraph.g[u][h].to)) {
                        needAddEdge.push(scaffoldgraph.g[u][h]);
                    }
                }
            }

            for (h = 0; h < scaffoldgraph.gr[u].length; ++h) {
                if (isGoodEdge(scaffoldgraph.gr[u][h].id)) {
                    if (curNodeSet.has(scaffoldgraph.gr[u][h].from)) {
                        needAddEdge.push(scaffoldgraph.gr[u][h]);
                    }
                }
            }
        }

        needAddEdge = needAddEdge.filter(function (value, index, self) {
            return self.indexOf(value) === index;
        });


        for (g = 0; g < needAddVert.length; ++g) {
            u = needAddVert[g];

            cy.add({
                group: "nodes",
                data: {
                    id: u,
                    label: createLabelForNode(u),
                    len: scaffoldgraph.nodes[u].len,
                    notAll: 0
                },
                position: {
                    x: evt.target.position().x + Math.floor(Math.random() * Math.floor(200) - 100),
                    y: evt.target.position().y + Math.floor(Math.random() * Math.floor(200) - 100)
                }
            });
        }

        for (g = 0; g < needAddEdge.length; ++g) {
            var eid = needAddEdge[g].id;

            cy.add({
                group: "edges",
                data: {
                    source: scaffoldgraph.edges[eid].from,
                    target: scaffoldgraph.edges[eid].to,
                    label: createLabelForEdge(eid),
                    faveColor: scaffoldgraph.libs[scaffoldgraph.edges[eid].lib].color,
                    weight: Math.log(scaffoldgraph.edges[eid].weight) + 1
                }
            });
        }

        for (g = 0; g < nodes_to_draw.length; ++g) {
            if (hasOtherEdges(nodes_to_draw[g], curNodeSet)) {
                cy.$('#' + nodes_to_draw[g]).data('notAll', 1);
            } else {
                cy.$('#' + nodes_to_draw[g]).data('notAll', 0);
            }
        }

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
    });
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

function drawAlongChromosome(chr) {
    var dnodes = [];
    var dedges = [];
    var inode = [];
    var pos = {};

    special_nodes.clear();
    var curNodeSet = new Set();

    for (i = 0; i < chromosomes[chr].alignments.length; ++i) {
        var curalig = chromosomes[chr].alignments[i];
        if (1.2 * (curalig.coorde - curalig.coordb) > scaffoldgraph.nodes[curalig.node_id].len) {
            inode.push({id: curalig.node_id, cb: curalig.coordb, ce: curalig.coorde});
            special_nodes.add(curalig.node_id);
            curNodeSet.add(curalig.node_id);
            dnodes.push({
                data: {
                    id: curalig.node_id,
                    label: createLabelForNode(curalig.node_id),
                    len: scaffoldgraph.nodes[curalig.node_id].len,
                    color: '#2A4986',
                    width: 10,
                    faveShape: 'rectangle'
                }
            });
            pos[curalig.node_id] = {x: Math.random()*1000, y: (curalig.coorde + curalig.coordb)/2};
        }
    }

    var vert_to_draw = findNodeAroundChr(inode, area_size, min_contig_len, isGoodEdge);
    alert(vert_to_draw.length);

    for (var g = 0; g < vert_to_draw.length; ++g) {
        curNodeSet.add(vert_to_draw[g].id);
    }

    for (g = 0; g < vert_to_draw.length; ++g) {
        var nall = '#2A4986';
        if (hasOtherEdges(vert_to_draw[g].id, curNodeSet)) {
            nall = '#FFD700';
        }

        dnodes.push({
            data: {
                id: vert_to_draw[g].id,
                label: createLabelForNode(vert_to_draw[g].id),
                len: 20 * Math.log(scaffoldgraph.nodes[vert_to_draw[g].id].len),
                width: 20 * Math.log(scaffoldgraph.nodes[vert_to_draw[g].id].len),
                color: nall,
                faveShape: 'ellipse'
            }
        });

        pos[vert_to_draw[g].id] = {x: 1000*vert_to_draw[g].rank + Math.random()*1000, y: vert_to_draw[g].y};
    }

    for (g = 0; g < edges_to_draw.length; ++g) {
        dedges.push({
            data: {
                id: "e" + edges_to_draw[g].toString(),
                source: scaffoldgraph.edges[edges_to_draw[g]].from,
                target: scaffoldgraph.edges[edges_to_draw[g]].to,
                label: createLabelForEdge(edges_to_draw[g]),
                faveColor: scaffoldgraph.libs[scaffoldgraph.edges[edges_to_draw[g]].lib].color,
                weight: Math.log(getWeight(edges_to_draw[g])) + 1
            }
        });
    }

    var cy = cytoscape({
        container: document.getElementById('mainpanel'),

        boxSelectionEnabled: false,
        autounselectify: true,
        maxZoom: 3,
        minZoom: 0.001,

        elements: {
            nodes: dnodes,
            edges: dedges
        },

        layout: {
            name: 'preset',

            positions: pos
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
                'curve-style': 'bezier',
                'target-arrow-shape': 'triangle',
                'line-color': 'data(faveColor)',
                'target-arrow-color': 'data(faveColor)',
                'width': 'data(weight)',
                'content': 'data(label)'
            })
    });

    cy.zoom(1);
    createCoordinates(chr, cy);
    createTapInfo(cy);

    cy.on('zoom', function() {
        createCoordinates(chr, cy);

        var v = inode[i];
        var nodeWidth = 10 / cy.zoom();
        cy.nodes("[faveShape='rectangle']").data('width', nodeWidth);
        for (var i = 0; i < edges_to_draw.length; ++i) {
            var e = edges_to_draw[i];
            var w = (Math.log(getWeight(e)) + 1)/cy.zoom();
            cy.getElementById("e" + e.toString()).data('weight', w);
        }
    });

    cy.on('pan', function () {
        createCoordinates(chr, cy);
    })
}

function handleAlongChromosomesFilter() {
    createComponentShowList(drawAlongChromosome, function(i) {
        return chromosomes[i].name;
    }, function(i) {
        return "Chromosome " + chromosomes[i].name;
    }, chromosomes.length);
}
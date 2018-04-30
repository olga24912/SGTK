var defZoom = 100;
var maxZoom = 10000000;
var IntervalTree = {};


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
                    x: 10,
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

function getYC_D() {
    return 100/defZoom;
}

function findNodeAroundChr(inode, area_size, min_contig_len, isGoodEdge, curNodeSet) {
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
                    if ((!used_id.has(curu)) && (curd < area_size || curNodeSet.has(curu))) {
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
                        yc = yc + Math.random() * getYC_D();

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
                    if ((!used_id.has(curu)) && (curd < area_size || curNodeSet.has(curu))) {
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
                        yc = yc  - Math.random() * getYC_D();

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


function getWidth(cy) {
    return 10/cy.zoom();
}

function getDispersion() {
    return 10;
}

function getRankDist() {
    return 25;
}

function getEdgeWeight(cy, e) {
    return Math.log(getWeight(e))/cy.zoom();
}

function getScala(cy) {
    return 1.5/cy.zoom();
}

function geOtherNodeWidth(id) {
    return Math.log(scaffoldgraph.nodes[id].len)/cy.zoom();
}

function createNewVerAlongChr(cy, area_size, min_contig_len, isGoodEdge, curNodeSet) {
    cy.on('tap', 'node', function (evt) {
        var newNode = new Set();
        var v = evt.target.id();
        var needAddVert = [];
        var needAddEdge = [];

        findConnectedVertex(v, needAddVert, needAddEdge, newNode, area_size, min_contig_len, isGoodEdge, curNodeSet);

        for (g = 0; g < needAddVert.length; ++g) {
            var u = needAddVert[g];
            var xc = evt.target.data('rank') + 1;
            var yc = calcYforV(u, area_size, min_contig_len, isGoodEdge, newNode, curNodeSet, cy, getDispersion());

            nodes_to_draw.push(u);
            cy.add({
                group: "nodes",
                data: {
                    id: u,
                    label: createLabelForNode(u),
                    len: geOtherNodeWidth(u),
                    width: geOtherNodeWidth(u),
                    color: genColorNode(u),
                    rank: xc,
                    faveShape: 'ellipse'
                },
                position: {
                    x: getRankDist() * xc + Math.random() * getDispersion(),
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
                    weight: getEdgeWeight(cy, eid),
                    curveStyle: "bezier",
                    controlPointDistances: 1,
                    scala: getScala(cy)
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

function getPointDistances(cy, e) {
    var bg = scaffoldgraph.edges[e].from;
    var ed = scaffoldgraph.edges[e].to;

    var order1 = cy.getElementById(bg).data('order');
    var order2 = cy.getElementById(ed).data('order');
    //yminbg = cy.getElementById(bg).data('ymin');
    //ymaxbg = cy.getElementById(bg).data('ymax');
    //ymined = cy.getElementById(ed).data('ymin');
    //ymaxed = cy.getElementById(ed).data('ymax');

    //minDif = Math.min(Math.min(Math.abs(yminbg - ymaxed), Math.abs(yminbg - ymined)), Math.min(Math.abs(ymaxbg - ymaxed), Math.abs(ymaxbg - ymined)));

    //minDif = Math.log(minDif);

    var oneStepDistant = 75;
    var randFree = 30;

    return Math.abs(order1 - order2)*oneStepDistant/cy.zoom() + randFree*Math.random()/cy.zoom();
}

function updateZooming(cy, posx, posmin, posmax, oldPosition) {
    //alert(cy.zoom());
    //alert(x / cy.zoom());
    //alert(cy.extent().x1);
    //alert(cy.extent().x2);
    var mul = 1;
    while (cy.zoom() > 10 && defZoom > 1) {
        cy.zoom(cy.zoom()/10);
        mul = mul * 10;
    }
    while (cy.zoom() < 1) {
        cy.zoom(cy.zoom()*10);
        mul = mul / 10;
    }

    //alert(cy.extent().x1);
    //alert(cy.extent().x2);

    cy.nodes().forEach(function (ele) {
        oldPosition.set(ele.id(), {x: ele.position("x"), y: ele.position("y")});
    });

    if (mul !== 1) {
        defZoom /= mul;
        oldPosition.clear();
        var ks = Array.from(posx.keys());
        for (var k = 0; k < ks.length; ++k) {
            posx.set(ks[k], posx.get(ks[k]) * mul);
        }
        posmin.clear();
        posmax.clear();
    }
}

function processFoundContig(elem, inode, posx, posmin, posmax, curNodeSet, order) {
    var vid = elem.id;
    if (!(posx.has(vid))) {
        posx.set(vid, Math.random() * getDispersion());
    }
    if (!(posmin.has(vid))) {
        posmin.set(vid, elem.cb);
        posmax.set(vid, elem.ce);
    }
    inode.push({id: elem.id, cb: elem.cb, ce: elem.ce, order: order});
    special_nodes.add(vid);
    curNodeSet.add(vid);
}

function findContigsByTree(tr, inode, posx, posmin, posmax, curNodeSet, ymin, ymax, lsm) {
    if (tr["lstL"].length === 0) {
        return;
    }

    if (ymin <= tr["md"] && ymax >= tr["md"]) {
        for (var i = 0; i < tr["lstL"].length; ++i) {
            processFoundContig(tr["lstL"][i], inode, posx, posmin, posmax, curNodeSet, tr["lt"]["size"] + i + lsm);
        }
    } else if (ymax < tr["md"]) {
        for (i = 0; i < tr["lstL"].length; ++i) {
            if (tr["lstL"][i].cb > ymax) {
                break;
            }
            processFoundContig(tr["lstL"][i], inode, posx, posmin, posmax, curNodeSet, tr["lt"]["size"] + i + lsm);
        }
    } else if (ymin > tr["md"]) {
        for (i = 0; i < tr["lstR"].length; ++i) {
            if (tr["lstR"][i].ce < ymax) {
                break;
            }
            processFoundContig(tr["lstR"][i], inode, posx, posmin, posmax, curNodeSet, tr["lt"]["size"] + tr["lstR"].length - i - 1 + lsm);
        }
    }

    if (ymin < tr["md"]) {
        findContigsByTree(tr["lt"], inode, posx, posmin, posmax, curNodeSet, ymin, ymax, lsm);
    }
    if (ymax > tr["md"]) {
        findContigsByTree(tr["rt"], inode, posx, posmin, posmax, curNodeSet, ymin, ymax, lsm + tr["lt"]["size"] + tr["lstR"].length);
    }
}

function findContigs(cy, chr, inode, posx, posmin, posmax, curNodeSet) {
    findContigsByTree(IntervalTree[defZoom], inode, posx, posmin, posmax, curNodeSet, cy.extent().y1, cy.extent().y2, 0);
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
                width: getWidth(cy),
                rank: 0,
                ymin: inode[i].cb,
                ymax: inode[i].ce,
                faveShape: 'rectangle',
                order: inode[i].order
            },
            position: {
                x: posx.get(vid),
                y: (posmin.get(vid) + posmax.get(vid))/2
            }
        });
    }
}

function addOtherNodes(cy, curNodeSet, vert_to_draw, oldPosition) {
    for (var g = 0; g < vert_to_draw.length; ++g) {
        curNodeSet.add(vert_to_draw[g].id);
    }
    for (g = 0; g < vert_to_draw.length; ++g) {
        nodes_to_draw.push(vert_to_draw[g].id);
        var nall = 'ellipse';
        if (hasOtherEdges(vert_to_draw[g].id, curNodeSet)) {
            nall = 'triangle';
        }

        if (!(oldPosition.has(vert_to_draw[g].id))) {
            oldPosition.set(vert_to_draw[g].id, {x: getRankDist() * vert_to_draw[g].rank + Math.random() * getDispersion(),
                y: vert_to_draw[g].y + Math.random() * getDispersion()});
        }

        cy.add({
            group: "nodes",
            data: {
                id: vert_to_draw[g].id,
                label: createLabelForNode(vert_to_draw[g].id),
                len: geOtherNodeWidth(vert_to_draw[g].id),
                width: geOtherNodeWidth(vert_to_draw[g].id),
                color: genColorNode(vert_to_draw[g].id),
                rank: vert_to_draw[g].rank,
                faveShape: nall
            },
            position: {
                x: oldPosition.get(vert_to_draw[g].id).x,
                y: oldPosition.get(vert_to_draw[g].id).y
            }
        });
    }
}

function addEdges(cy) {
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
                    weight: getEdgeWeight(cy, edges_to_draw[g]),
                    curveStyle: "bezier",
                    controlPointDistances: 1,
                    scala: getScala(cy)
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
                    weight: getEdgeWeight(cy, edges_to_draw[g]),
                    curveStyle: "unbundled-bezier",
                    controlPointDistances: getPointDistances(cy, edges_to_draw[g]),
                    scala: getScala(cy)
                }
            });
        }
    }
}

function createGraph(chr, cy, curNodeSet, posx, posmin, posmax, oldPosition) {
    updateZooming(cy, posx, posmin, posmax, oldPosition);
    cy.elements().remove();
    special_nodes.clear();

    var inode = [];

    nodes_to_draw = [];
    edges_to_draw = [];

    findContigs(cy, chr, inode, posx, posmin, posmax, curNodeSet);
    addContigs(cy, inode, posx, posmin, posmax);

    var vert_to_draw = findNodeAroundChr(inode, area_size, min_contig_len, isGoodEdge, curNodeSet);
    addOtherNodes(cy, curNodeSet, vert_to_draw, oldPosition);
    addEdges(cy);

    createCoordinates(chr, cy);
}


//var tmp = [];
function buildITree(lst) {
    if (lst.length === 0) {
        return {lt: -1, rt: -1, lstL: [], lstR: [], md: -1, size: 0};
    }

    /*tmp = [];
    for (var i = 0; i < lst.length; ++i) {
        tmp.push(lst[i].cb);
        tmp.push(lst[i].ce);
    }
    tmp.sort();*/

    var tr = {lt: -1, rt: -1, lstL: [], lstR: [], md: lst[Math.floor(lst.length/2)].cb/*tmp[Math.floor(tmp.length/2)]*/, size: lst.length};
    var lstR = [];
    var lstL = [];
    for (i = 0; i < lst.length; ++i) {
        if (lst[i].ce < tr.md) {
            lstL.push(lst[i]);
        } else if (lst[i].cb > tr.md) {
            lstR.push(lst[i]);
        } else {
            tr["lstL"].push(lst[i]);
            tr["lstR"].push(lst[i]);
        }
    }

    tr["lstL"].sort(function (a, b) { return a.cb - b.cb });
    tr["lstR"].sort(function (a, b) { return b.ce - a.ce });

    tr["lt"] = buildITree(lstL);
    tr["rt"] = buildITree(lstR);
    return tr;
}

function buildIT(chr, dz) {
    var lst = [];
    for (var i = 0; i < chromosomes[chr].alignments.length; ++i) {
        var curalig = chromosomes[chr].alignments[i];
        lst.push({id: curalig.node_id, cb: curalig.coordb / dz, ce: curalig.coorde / dz})
    }

    return buildITree(lst);
}

function drawAlongChromosome(chr) {
    defZoom = 100;
    for (var i = 1; i <= maxZoom; i *= 10) {
        IntervalTree[i] = buildIT(chr, i);
    }

    var curNodeSet = new Set();
    var posx = new Map();
    var posmin = new Map();
    var posmax = new Map();
    var oldPosition = new Map();
    cy = cytoscape({
        container: document.getElementById('mainpanel'),

        boxSelectionEnabled: false,
        autounselectify: true,
        maxZoom: 11,
        minZoom: 0.1,

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
                "arrow-scale": 'data(scala)',
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
        createGraph(chr, cy, curNodeSet, posx, posmin, posmax, oldPosition);
    });
    cy.on('pan', function() {
        createGraph(chr, cy, curNodeSet, posx, posmin, posmax, oldPosition);
    });

    createGraph(chr, cy, curNodeSet, posx, posmin, posmax, oldPosition);
}

function handleAlongChromosomesFilter() {
    createComponentShowList(drawAlongChromosome, function(i) {
        return chromosomes[i].name;
    }, function(i) {
        return "Chromosome " + chromosomes[i].name;
    }, chromosomes.length);
}
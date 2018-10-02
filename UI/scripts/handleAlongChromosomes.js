var defZoom = 100;
var maxZoom = 10000000;
var IntervalTree = {};
var lastMinX = 0;
var lastMaxX = 0;

function generateCoordinateLabel(x, delta) {
    if (defZoom*delta >= 1000000) {
        return (x/1000000).toString() + "M";
    } else if (defZoom*delta >= 1000) {
        return (x/1000).toString() + "K";
    } else {
        return x
    }
}

function createCoordinates(chr, cy) {
    var cur_zoom = cy.zoom();
    var cur_coord = cy.extent();

    var delta = 1;
    while (cur_zoom * delta < 100) {
        delta *= 10;
    }

    var start_pos = Math.max(delta, (Math.floor(cur_coord.x1 / delta) - 5) * delta);

    cy.remove(cy.$('#start'));
    cy.add({
        group: "nodes",
        data: {
            id: "start",
            label: "Start: 0",
            len: 2 / cy.zoom(),
            color: '#ff0000',
            width: 20 / cy.zoom(),
            faveShape: 'rectangle'
        },
        position: {
            y: 0,
            x: 0
        }
    });
    var xpos = cy.$('#start').renderedPosition().x;
    cy.$('#start').renderedPosition({
        y: 32,
        x: xpos
    });
    cy.$('#start').style({"font-size": 20 / cy.zoom(),"text-valign": "top", "text-halign": "left"});
    cy.$('#start').lock();

    cy.remove(cy.$('#end'));
    cy.add({
        group: "nodes",
        data: {
            id: "end",
            label: "End: " + chromosomes[chr].len,
            len: 2 / cy.zoom(),
            color: '#ff0000',
            width: 20 / cy.zoom(),
            faveShape: 'rectangle'
        },
        position: {
            y: 0,
            x: chromosomes[chr].len/defZoom
        }
    });
    xpos = cy.$('#end').renderedPosition().x;
    cy.$('#end').renderedPosition({
        y: 32,
        x: xpos
    });
    cy.$('#end').style({"font-size": 20 / cy.zoom(), "text-valign": "top", "text-halign": "right"});
    cy.$('#end').lock();


    for (var i = 0;  i < 50; ++i) {
        cy.remove(cy.$('#chrcoord' + i));
        if (start_pos + delta * i < chromosomes[chr].len/defZoom) {
            cy.add({
                group: "nodes",
                data: {
                    id: "chrcoord" + i,
                    label: generateCoordinateLabel((start_pos + delta * i)*defZoom, delta),
                    len: 2 / cy.zoom(),
                    color: '#2A4986',
                    width: 20 / cy.zoom(),
                    faveShape: 'rectangle'
                },
                position: {
                    y: 30,
                    x: start_pos + delta * i
                }
            });
            xpos = cy.$('#chrcoord' + i).renderedPosition().x;
            cy.$('#chrcoord' + i).renderedPosition({
                y: 32,
                x: xpos
            });

            cy.$('#chrcoord' + i).style({"font-size": 20 / cy.zoom(), "text-valign": "top", "text-halign": "center"});
            cy.$('#chrcoord' + i).lock();
        }
    }
}

function calcY(curv, ypos, sumw) {
    return ypos/sumw;
}

function getYC_D() {
    return 1000/defZoom;
}


function isBigContig(cb, ce, dz) {
    return (ce - cb > dz/10 && ce - cb > min_contig_len);
}

function getWeight(e) {
    var curW = scaffoldgraph.edges[e].weight;

    if (scaffoldgraph.libs[scaffoldgraph.edges[e].lib].type == 'SCAFF') {
        curW = 16;
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
                        var yc = cy.$('#' + v).position().x - rcoef * Math.random();
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
                        yc = cy.$('#' + v).position().x + rcoef * Math.random();
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
    return 10;//50;
}

function getContigXPosD() {
    return 1000/defZoom;
}

function contigHasEdgesInThisScala(cy, v) {
    return cy.getElementById(v).data('len') >= 5;
}

function getRankDist() {
    return 50;
}

function getEdgeWeight(cy, e) {
    var edge = scaffoldgraph.edges[e];

    if (scaffoldgraph.libs[edge.lib].type === "FASTG" ||
        scaffoldgraph.libs[edge.lib].type === "GFA") {
        return 5/cy.zoom();
    }

    if (scaffoldgraph.libs[edge.lib].type === "SCAFF") {
        return 3/cy.zoom();
    }

    if (scaffoldgraph.libs[edge.lib].type === "DNA_PAIR" ||
        scaffoldgraph.libs[edge.lib].type === "LONG" ||
        scaffoldgraph.libs[edge.lib].type === "RNA_PAIR" ||
        scaffoldgraph.libs[edge.lib].type === "RNA_SPLIT_50" ||
        scaffoldgraph.libs[edge.lib].type === "RNA_SPLIT_30" ||
        scaffoldgraph.libs[edge.lib].type === "CONNECTION") {
        return 1/cy.zoom();
    }

    return Math.min(5, Math.log(getWeight(e)) + 1)/cy.zoom();
}

function getScala(cy) {
    return 1.5/cy.zoom();
}

function geOtherNodeWidth(id) {
    return Math.log(scaffoldgraph.nodes[id].len)*2/cy.zoom();
}

function createNewVerAlongChr(cy, area_size, min_contig_len, isGoodEdge, curNodeSet, openNode) {
    cy.on('tap', 'node', function (evt) {
        var newNode = new Set();
        var v = evt.target.id();
        openNode.add(v);
        var needAddVert = [];
        var needAddEdge = [];

        findConnectedVertex(v, needAddVert, needAddEdge, newNode, area_size, min_contig_len, isGoodEdge, curNodeSet);

        for (g = 0; g < needAddVert.length; ++g) {
            var u = needAddVert[g];
            openNode.add(u);
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
                    notALL: 0,
                    faveShape: 'ellipse'
                },
                position: {
                    y: getRankDist() * xc + Math.random() * getDispersion(),
                    x: yc + (Math.random() - 0.5) * getDispersion()
                }
            });
            cy.$("#" + u.toString()).style({"font-size": 10 / cy.zoom(), "text-valign": "center", "text-halign": "right"});
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
            cy.$("#" + "e" + eid.toString()).style({"font-size": 10 / cy.zoom(), "edge-text-rotation": "autorotate"});
        }

        for (g = 0; g < nodes_to_draw.length; ++g) {
            if (hasOtherEdges(nodes_to_draw[g], curNodeSet)) {
                cy.$('#' + nodes_to_draw[g]).data('notALL', 1);
            } else {
                cy.$('#' + nodes_to_draw[g]).data('notALL', 0);
            }
        }
        //createTapInfo(cy);
    });
}

function updateZooming(cy, posx, posmin, posmax, oldPosition) {
    var mul = 1;
    while (cy.zoom()/mul > 10 && defZoom/mul > 1) {
        mul = mul * 10;
    }
    while (cy.zoom()/mul < 1 && defZoom/mul < maxZoom) {
        mul = mul / 10;
    }
    cy.zoom(cy.zoom()/mul);

    cy.nodes().forEach(function (ele) {
        oldPosition.set(ele.id(), {x: ele.position("x"), y: ele.position("y")});
    });

    if (mul !== 1) {
        defZoom /= mul;
        cy.maxZoom(defZoom);
        cy.minZoom(defZoom/maxZoom);

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
        posx.set(vid, Math.random() * getContigXPosD());
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
            if (tr["lstR"][i].ce < ymin) {
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
    lastMinX = cy.extent().x1 - (cy.extent().x2 - cy.extent().x1);
    lastMaxX = cy.extent().x2 + (cy.extent().x2 - cy.extent().x1);
    findContigsByTree(IntervalTree[defZoom], inode, posx, posmin, posmax, curNodeSet, lastMinX, lastMaxX, 0);
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
                y: posx.get(vid),
                x: (posmin.get(vid) + posmax.get(vid))/2
            }
        });
        if (defZoom < 10000) {
            cy.$("#" + vid.toString()).style({
                "font-size": 15 / cy.zoom(),
                "text-valign": "top",
                "text-halign": "center"
            });
        }
    }
}

function addOtherNodes(cy, curNodeSet, vert_to_draw, oldPosition) {
    for (var g = 0; g < vert_to_draw.length; ++g) {
        curNodeSet.add(vert_to_draw[g].id);
    }
    for (g = 0; g < vert_to_draw.length; ++g) {
        nodes_to_draw.push(vert_to_draw[g].id);
        var nall = 0;
        if (hasOtherEdges(vert_to_draw[g].id, curNodeSet)) {
            nall = 1;
        }

        if (!(oldPosition.has(vert_to_draw[g].id.toString()))) {
            oldPosition.set(vert_to_draw[g].id.toString(), {y: getRankDist() * vert_to_draw[g].rank + Math.random() * getDispersion(),
                x: vert_to_draw[g].x + (Math.random() - 0.5) * getDispersion()});
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
                faveShape: 'ellipse',
                notALL: nall
            },
            position: {
                x: oldPosition.get(vert_to_draw[g].id.toString()).x,
                y: oldPosition.get(vert_to_draw[g].id.toString()).y
            }
        });
        if (defZoom < 10000) {
            cy.$("#" + vert_to_draw[g].id.toString()).style({
                "font-size": 15 / cy.zoom(),
                "text-valign": "bottom",
                "text-halign": "left"
            });
        }
    }
}


function getPointDistances(cy, e, deepsPOS, toSmallCoord) {
    var bg = scaffoldgraph.edges[e].from;
    var ed = scaffoldgraph.edges[e].to;

    var order1 = cy.getElementById(bg).data('order');
    var order2 = cy.getElementById(ed).data('order');

    var oneStepDistant = 75;
    var randFree = 100;
    return -(deepsPOS[toSmallCoord[Math.min(order1, order2)]][toSmallCoord[Math.max(order1, order2)] -
    toSmallCoord[Math.min(order1, order2)]] * oneStepDistant/cy.zoom() + randFree*Math.random()/cy.zoom());
    /*deepsPOS[toSmallCoord[Math.min(order1, order2)]][toSmallCoord[Math.max(order1, order2)] -
       toSmallCoord[Math.min(order1, order2)]]*/
    //return Math.min(5, (Math.max(order1, order2) - Math.min(order1, order2))) * oneStepDistant/cy.zoom() + randFree*Math.random()/cy.zoom();
}

function calculateDinamicForDistPoint(cy, deepsPOS, toSmallCoord) {
    var idslist = [];

    var mxLg = 0;

    for (var g = 0; g < edges_to_draw.length; ++g) {
        if (special_nodes.has(scaffoldgraph.edges[edges_to_draw[g]].from) &&
            special_nodes.has(scaffoldgraph.edges[edges_to_draw[g]].to)) {
            idslist.push(parseInt(cy.getElementById(scaffoldgraph.edges[edges_to_draw[g]].from).data('order')));
            idslist.push(parseInt(cy.getElementById(scaffoldgraph.edges[edges_to_draw[g]].to).data('order')));
        }
    }

    idslist.sort(function (a, b) {
        return a - b;
    });
    idslist = idslist.filter(function (v, i, a) {
        return a.indexOf(v) === i;
    });

    for (g = 0; g < idslist.length; ++g) {
        toSmallCoord[idslist[g]] = g;
    }

    for (g = 0; g < edges_to_draw.length; ++g) {
        var fv = scaffoldgraph.edges[edges_to_draw[g]].from;
        var tv = scaffoldgraph.edges[edges_to_draw[g]].to;
        if (special_nodes.has(fv) &&
            special_nodes.has(tv)) {
            mxLg = Math.max(mxLg, toSmallCoord[cy.getElementById(tv).data('order')] - toSmallCoord[cy.getElementById(fv).data('order')]);
            mxLg = Math.max(mxLg, toSmallCoord[cy.getElementById(fv).data('order')] - toSmallCoord[cy.getElementById(tv).data('order')]);
        }
    }

    for (g = 0; g < idslist.length; ++g) {
        deepsPOS.push([]);
        for (var j = 0; j <= mxLg; ++j) {
            deepsPOS[g].push(0);
        }
    }

    for (g = 0; g < edges_to_draw.length; ++g) {
        var fv = scaffoldgraph.edges[edges_to_draw[g]].from;
        var tv = scaffoldgraph.edges[edges_to_draw[g]].to;
        if (special_nodes.has(fv) &&
            special_nodes.has(tv)) {
            deepsPOS[toSmallCoord[cy.getElementById(fv).data('order')]]
                [toSmallCoord[cy.getElementById(tv).data('order')] - toSmallCoord[cy.getElementById(fv).data('order')]] = 1;
        }
    }

    for (var lg = 1; lg <= mxLg; ++lg) {
        for (g = 0; g < idslist.length - lg + 1; ++g) {
            if (g === idslist.length - 1) {
                deepsPOS[g][lg] += deepsPOS[g][lg - 1];
            } else {
                deepsPOS[g][lg] += Math.max(deepsPOS[g][lg - 1], deepsPOS[g + 1][lg - 1]);
            }
        }
    }
}


function isAlign(u) {
    return cy.getElementById(u).data('faveShape') === 'rectangle';
}

function addEdges(cy) {
    var deepsPOS = [];
    var toSmallCoord = {};
    calculateDinamicForDistPoint(cy, deepsPOS, toSmallCoord);

    for (var g = 0; g < edges_to_draw.length; ++g) {
        if (!(isAlign(scaffoldgraph.edges[edges_to_draw[g]].from)) ||
            !(isAlign(scaffoldgraph.edges[edges_to_draw[g]].to))) {
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
                    controlPointDistances: getPointDistances(cy, edges_to_draw[g], deepsPOS, toSmallCoord),
                    scala: getScala(cy)
                }
            });
        }
        if (defZoom < 10000) {
            cy.$("#" + "e" + edges_to_draw[g].toString()).style({
                "font-size": 15 / cy.zoom(),
                "edge-text-rotation": "autorotate"
            });
        }
    }
}

function createGraph(chr, cy, curNodeSet, posx, posmin, posmax, oldPosition, openNode) {
    updateZooming(cy, posx, posmin, posmax, oldPosition);
    cy.elements().remove();
    special_nodes.clear();
    curNodeSet.clear();
    oldPosition.clear();

    var inode = [];

    nodes_to_draw = [];
    edges_to_draw = [];

    findContigs(cy, chr, inode, posx, posmin, posmax, curNodeSet);
    addContigs(cy, inode, posx, posmin, posmax);

    var vert_to_draw = findNodeAroundChr(inode, area_size, min_contig_len, isGoodEdge, curNodeSet, openNode, cy);
    addOtherNodes(cy, curNodeSet, vert_to_draw, oldPosition);

    addEdges(cy);

    createInformationShown(cy);
    createCoordinates(chr, cy);
}

function updateGraph(chr, cy) {
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
        if (isBigContig(curalig.coordb, curalig.coorde, dz)) {
            lst.push({id: curalig.node_id, cb: curalig.coordb / dz, ce: curalig.coorde / dz})
        }
    }

    return buildITree(lst);
}

function drawAlongChromosome(chr) {
    defZoom = 100;
    lastMinX = 0;
    lastMaxX = 0;
    for (var i = 1; i <= maxZoom; i *= 10) {
        IntervalTree[i] = buildIT(chr, i);
    }

    var curNodeSet = new Set();
    var openNode = new Set();
    var posx = new Map();
    var posmin = new Map();
    var posmax = new Map();
    var oldPosition = new Map();

    cy = cytoscape({
        container: document.getElementById('mainpanel'),

        boxSelectionEnabled: false,
        autounselectify: true,
        maxZoom: 11000,
        minZoom: 0.0001,

        layout: {
            name: 'preset'
        },

        style: cytoscape.stylesheet()
            .selector('node')
            .css({
                'shape': 'data(faveShape)',
                'content': 'data(label)',
                'color': '#2A4986',
                'width': 'data(len)',
                'height': 'data(width)',
                'background-color': 'data(color)',
                'border-width': 'mapData(notALL, 0, 1, 0px, 1px)'
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

    var element =  document.getElementById("cynav");
    if (typeof(element) != 'undefined' && element != null) {
        document.getElementById("cynav").remove();
    }
    createNewVerAlongChr(cy, area_size, min_contig_len, isGoodEdge, curNodeSet, openNode);

    cy.on('zoom', function () {
        createGraph(chr, cy, curNodeSet, posx, posmin, posmax, oldPosition, openNode);
        (document.getElementById("zoomInput")).innerText = Math.floor(cy.zoom() * 100 * 100/ defZoom).toString() +  "%";
    });

    cy.on('pan', function() {
        if (cy.extent().x1 >= lastMinX && cy.extent().x2 <= lastMaxX) {
            updateGraph(chr, cy);
        } else {
            createGraph(chr, cy, curNodeSet, posx, posmin, posmax, oldPosition, openNode);
        }
    });


    cy.ready(function () {
        window.cy = this;
    });
    setTimeout(function() {
        createGraph(chr, cy, curNodeSet, posx, posmin, posmax, oldPosition, openNode);
    }, 1);
}

function handleAlongChromosomesFilter() {
    createComponentShowList(drawAlongChromosome, function(i) {
        return chromosomes[i].name;
    }, function(i) {
        return "Chromosome " + chromosomes[i].name;
    }, chromosomes.length);
}
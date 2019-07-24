/*
* Create graph for filtration in genome browser mode
*/

var defZoom = 10;
var maxZoom = 10000000;
var IntervalTree = {};
var lastMinX = 0;
var lastMaxX = 0;
var mulConst = 10;
var minZoomUpdate = 0.1;
var maxZoomUpdate = 1;
var widthAddConst = 1;
var EPS = 0.0001;
var graphUpdating = false;

function generateCoordinateLabel(x, delta) {
    if (defZoom*delta >= 1000000) {
        return (x/1000000).toString() + "M";
    } else if (defZoom*delta >= 1000) {
        return (x/1000).toString() + "K";
    } else {
        return x
    }
}


//Generate nodes for coordinate
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
            faveShape: 'rectangle',
            notALL: 0
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
            faveShape: 'rectangle',
            notALL: 0
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
                    faveShape: 'rectangle',
                    notALL: 0
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

function getWeight(e) {
    var curW = scaffoldgraph.edges[e].weight;

    if (scaffoldgraph.libs[scaffoldgraph.edges[e].lib].type == 'SCAFF') {
        curW = 16;
    }
    return curW
}

function getWidth(cy) {
    return 10/cy.zoom();
}

function getScala(cy) {
    return 2;
}

function getDispersion() {
    return 250;
}

function getContigXPosD() {
    return 11/cy.zoom();//;1000/defZoom;
}

function getRankDist() {
    return 500;
}

function getNotAllValue() {
    return 5/cy.zoom() + 'px';
}

function contigHasEdgesInThisScala(cy, v) {
    return cy.getElementById(v).data('len') >= 50;
}

function isBigContig(cb, ce, dz) {
    return (/*ce - cb > dz*5 &&*/ ce - cb > min_contig_len);
}

function printThisContig(cb, ce, dz) {
    if (isBigContig(cb, ce, dz)) {
        if (Math.random() < (ce - cb)/(10* dz)) {
            return true;
        } else {
            return false;
        }
    }
}

function geOtherNodeWidth(id) {
    return (Math.log(scaffoldgraph.nodes[id].len)*4)/cy.zoom();
}

function getEdgeWeight(cy, e) {
    var edge = scaffoldgraph.edges[e];

    if (scaffoldgraph.libs[edge.lib].type === "FASTG" ||
        scaffoldgraph.libs[edge.lib].type === "GFA" || 
        scaffoldgraph.libs[edge.lib].type === "GFA2") {
        return (widthAddConst + 5)/cy.zoom();
    }

    if (scaffoldgraph.libs[edge.lib].type === "SCAFF") {
        return (widthAddConst + 3)/cy.zoom();
    }

    if (scaffoldgraph.libs[edge.lib].type === "DNA_PAIR" ||
        scaffoldgraph.libs[edge.lib].type === "LONG" ||
        scaffoldgraph.libs[edge.lib].type === "RNA_PAIR" ||
        scaffoldgraph.libs[edge.lib].type === "RNA_SPLIT_50" ||
        scaffoldgraph.libs[edge.lib].type === "RNA_SPLIT_30" ||
        scaffoldgraph.libs[edge.lib].type === "CONNECTION") {
        return (widthAddConst + 1)/cy.zoom();
    }

    return (widthAddConst + Math.min(5, Math.log(getWeight(e)) + 1))/cy.zoom();
}

function getPointDistances(cy, e, deepsPOS, toSmallCoord) {
    var bg = scaffoldgraph.edges[e].from;
    var ed = scaffoldgraph.edges[e].to;

    var order1 = cy.getElementById(bg).data('order');
    var order2 = cy.getElementById(ed).data('order');

    var oneStepDistant = 75;
    var randFree = 100;
    return -5*(deepsPOS[toSmallCoord[Math.min(order1, order2)]][toSmallCoord[Math.max(order1, order2)] -
    toSmallCoord[Math.min(order1, order2)]] * oneStepDistant + randFree*Math.random());


    /*deepsPOS[toSmallCoord[Math.min(order1, order2)]][toSmallCoord[Math.max(order1, order2)] -
       toSmallCoord[Math.min(order1, order2)]]*/
    //return Math.min(5, (Math.max(order1, order2) - Math.min(order1, order2))) * oneStepDistant/cy.zoom() + randFree*Math.random()/cy.zoom();
}


function nodePositionChange(cy, posmin, posmax) {
    cy.on('position', 'node', function (evt) {
        var v = evt.target.id();
        if (graphUpdating === false && cy.getElementById(v).data('rank') === 0) {
            if (posmin.has(parseInt(v))) {
                if (cy.getElementById(v).position('x') != ((posmin.get(parseInt(v)) + posmax.get(parseInt(v))) / 2)) {
                    cy.getElementById(v).position('x', (posmin.get(parseInt(v)) + posmax.get(parseInt(v))) / 2);
                }
            }

            var order = cy.getElementById(v).data('order');
            var y1 = cy.nodes("[order = " + (order + 1).toString() + "]").position('y');
            var y2 = cy.nodes("[order = " + (order - 1).toString() + "]").position('y');
            var y = cy.getElementById(v).position('y');
            var direction = 1 - 2 * Math.round(Math.random());
            if (y1 !== undefined && Math.abs(y - y1) < getContigXPosD() - EPS) {
                cy.getElementById(v).position('y', y1 + direction * (getContigXPosD()));
            }

            if (y2 !== undefined && Math.abs(y - y2) < getContigXPosD() - EPS) {
                cy.getElementById(v).position('y', y2 + direction * (getContigXPosD()));
            }
        }
    });
}


//open vertex
function createNewVerAlongChr(cy, area_size, min_contig_len, isGoodEdge, curNodeSet, openNode) {
    cy.on('taphold', 'node', function (evt) {
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
                cy.$('#' + nodes_to_draw[g]).data('notALL', getNotAllValue());
            } else {
                cy.$('#' + nodes_to_draw[g]).data('notALL', 0);
            }
        }
        //createTapInfo(cy);
        cy.ignoreTap = true
    });
}


//Update current zooming
function updateZooming(cy, posx, posmin, posmax, oldPosition) {
    var mul = 1;
    while (cy.zoom()/mul > maxZoomUpdate && defZoom/mul > 1) {
        mul = mul * mulConst;
    }
    while (cy.zoom()/mul < minZoomUpdate && defZoom/mul < maxZoom) {
        mul = mul / mulConst;
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

function processFoundContig(elem, inode, posx, posmin, posmax, curNodeSet, order,  levelX) {
    var vid = elem.id;
    if (printThisContig(elem.cb, elem.ce, defZoom)) {
        posx.set(vid, -1 * levelX.get(vid) * getContigXPosD());
        if (!(posmin.has(vid))) {
            posmin.set(vid, elem.cb / defZoom);
            posmax.set(vid, elem.ce / defZoom);
        }
        inode.push({id: elem.id, cb: elem.cb / defZoom, ce: elem.ce / defZoom, order: order});
        special_nodes.add(vid);
        curNodeSet.add(vid);
    }
}


/*
*Find list of contigs in interval tree for current position and zoom
 */
function findContigsByTree(tr, inode, posx, posmin, posmax, curNodeSet, ymin, ymax, lsm,  levelX) {
    if (tr["lstL"].length === 0) {
        return;
    }

    if (ymin <= tr["md"] && ymax >= tr["md"]) {
        for (var i = 0; i < tr["lstL"].length; ++i) {
            processFoundContig(tr["lstL"][i], inode, posx, posmin, posmax, curNodeSet, tr["lt"]["size"] + i + lsm,  levelX);
        }
    } else if (ymax < tr["md"]) {
        for (i = 0; i < tr["lstL"].length; ++i) {
            if (tr["lstL"][i].cb > ymax) {
                break;
            }
            processFoundContig(tr["lstL"][i], inode, posx, posmin, posmax, curNodeSet, tr["lt"]["size"] + i + lsm,  levelX);
        }
    } else if (ymin > tr["md"]) {
        for (i = 0; i < tr["lstR"].length; ++i) {
            if (tr["lstR"][i].ce < ymin) {
                break;
            }
            processFoundContig(tr["lstR"][i], inode, posx, posmin, posmax, curNodeSet, tr["lt"]["size"] + tr["lstR"].length - i - 1 + lsm,  levelX);
        }
    }

    if (ymin < tr["md"]) {
        findContigsByTree(tr["lt"], inode, posx, posmin, posmax, curNodeSet, ymin, ymax, lsm,  levelX);
    }
    if (ymax > tr["md"]) {
        findContigsByTree(tr["rt"], inode, posx, posmin, posmax, curNodeSet, ymin, ymax, lsm + tr["lt"]["size"] + tr["lstR"].length,  levelX);
    }
}

function findContigs(cy, chr, inode, posx, posmin, posmax, curNodeSet,  levelX) {
    lastMinX = cy.extent().x1 - (cy.extent().x2 - cy.extent().x1);
    lastMaxX = cy.extent().x2 + (cy.extent().x2 - cy.extent().x1);
    findContigsByTree(IntervalTree, inode, posx, posmin, posmax, curNodeSet, lastMinX*defZoom, lastMaxX*defZoom, 0,  levelX);
}


//Add contigs to cytoscape
function addContigs(cy, inode, posx, posmin, posmax) {
    for (i = 0; i < inode.length; ++i) {
        var vid = inode[i].id;
        cy.add({
            group: "nodes",
            data: {
                id: vid,
                label: createLabelForNode(vid),
                len: Math.max(inode[i].ce - inode[i].cb, 2),
                color: genColorNode(vid),
                width: getWidth(cy),
                rank: 0,
                ymin: inode[i].cb,
                ymax: inode[i].ce,
                faveShape: 'rectangle',
                order: inode[i].order,
                notALL: 0
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


//Add not align vertex to cytoscape
function addOtherNodes(cy, curNodeSet, vert_to_draw, oldPosition) {
    for (var g = 0; g < vert_to_draw.length; ++g) {
        curNodeSet.add(vert_to_draw[g].id);
    }
    for (g = 0; g < vert_to_draw.length; ++g) {
        nodes_to_draw.push(vert_to_draw[g].id);
        var nall = 0;
        if (hasOtherEdges(vert_to_draw[g].id, curNodeSet)) {
            nall = getNotAllValue();
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


//calculate edges distance point
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


//Add edges to cytoscape
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

//Rebuilt graph
function createGraph(chr, cy, curNodeSet, posx, posmin, posmax, oldPosition, openNode, levelX) {
    updateZooming(cy, posx, posmin, posmax, oldPosition);
    cy.elements().remove();
    special_nodes.clear();
    curNodeSet.clear();
    oldPosition.clear();

    var inode = [];

    nodes_to_draw = [];
    edges_to_draw = [];

    findContigs(cy, chr, inode, posx, posmin, posmax, curNodeSet,  levelX);
    addContigs(cy, inode, posx, posmin, posmax);

    var vert_to_draw = findNodeAroundChr(inode, area_size, min_contig_len, isGoodEdge, curNodeSet, openNode, cy);
    addOtherNodes(cy, curNodeSet, vert_to_draw, oldPosition);

    addEdges(cy);

    createInformationShown(cy);
    createCoordinates(chr, cy);
}


//update graph on small pan or zoom
function updateGraph(chr, cy, levelX) {
    graphUpdating = true;

    createCoordinates(chr, cy);
    var wght = getWidth(cy);
    cy.nodes().filter(function (ele) {
        return ele.data('rank') === 0;
    }).forEach(function (node, index) {
        node.data('width', wght);
        node.position('y', -1*getContigXPosD()*levelX.get(parseInt(node.id())));
    });


    cy.nodes().filter(function (ele) {
        return ele.data('faveShape') === 'ellipse';
    }).forEach(function (node, i) {
        node.data('len', geOtherNodeWidth(node.id()));
        node.data('width', geOtherNodeWidth(node.id()));
        if (node.data('notALL') != 0) {
            node.data('notALL', getNotAllValue());
        }
    });

    cy.edges().forEach(function (edge) {
        edge.data('weight', getEdgeWeight(cy, parseInt(edge.id().substr(1))));
        edge.data('scala', getScala(cy));
    });

    graphUpdating = false;
}


/*
* Built interval Tree for fast search of contigs
*/

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


//Generate vertical position level for chromosomes node
function calculateContigLevelX(chr, levelX) {
    var lst = [];
    var j = 0;
    var levelList = new Array(20);
    for (var i = 0; i < chromosomes[chr].alignments.length; ++i) {
        var curalig = chromosomes[chr].alignments[i];
        for (j = 0; j < 20; ++j) {
            levelList[parseInt(j)] = 0;
        }
        for (j = 0; j < lst.length; ++j) {
            if ((curalig.coordb <= lst[j].cb && lst[j].cb <= curalig.coorde) ||
                (curalig.coordb <= lst[j].ce && lst[j].ce <= curalig.coorde) ||
                (lst[j].cb <= curalig.coordb && curalig.coordb <= lst[j].ce)) {
                levelList[parseInt(levelX.get(lst[j].id))] = 1;
            }
        }
        lst.push({id: curalig.node_id, cb: curalig.coordb, ce: curalig.coorde});
        for (j = 19; j >= 0; --j) {
            if (levelList[j] === 0) {
                levelX.set(curalig.node_id, j);
            }
        }
    }
}

//generate cytoscape graph for genome browser
function drawAlongChromosome(chr) {
    defZoom = 100;
    lastMinX = 0;
    lastMaxX = 0;
    IntervalTree = buildIT(chr, 1);

    var curNodeSet = new Set();
    var openNode = new Set();
    var posx = new Map();
    var posmin = new Map();
    var posmax = new Map();
    var oldPosition = new Map();
    var levelX = new Map();
    calculateContigLevelX(chr, levelX);
    
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
                'border-width': 'data(notALL)'
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
    createNewVerAlongChr(cy, area_size, min_contig_len, isGoodEdge, curNodeSet, openNode);
    nodePositionChange(cy, posmin, posmax);
    highlightOnTap(cy);

    cy.on('zoom', function () {
        if (cy.zoom() < maxZoomUpdate && cy.zoom() > minZoomUpdate && cy.extent().x1 >= lastMinX && cy.extent().x2 <= lastMaxX) {
            updateGraph(chr, cy, levelX);
        } else {
            document.getElementById("UpdateGraph").style.visibility = "visible";
            document.getElementById("updateGrpahButton").focus();
            updateGraph(chr, cy, levelX);
        }

        (document.getElementById("zoomInput")).innerText = Math.floor(cy.zoom() * 100 * 100/ defZoom).toString() +  "%";
    });

    cy.on('pan', function() {
        if (cy.extent().x1 >= lastMinX && cy.extent().x2 <= lastMaxX) {
            updateGraph(chr, cy, levelX);
        } else {
            document.getElementById("UpdateGraph").style.visibility = "visible";
            document.getElementById("updateGrpahButton").focus();
            updateGraph(chr, cy, levelX);
        }
    });

    document.getElementById("updateGrpahButton").onclick = function(){
        createGraph(chr, cy, curNodeSet, posx, posmin, posmax, oldPosition, openNode, levelX);
        document.getElementById("UpdateGraph").style.visibility = "hidden";
        (document.getElementById("zoomInput")).innerText = Math.floor(cy.zoom() * 100 * 100/ defZoom).toString() +  "%";
    };

    cy.ready(function () {
        window.cy = this;
    });
    setTimeout(function() {
        createGraph(chr, cy, curNodeSet, posx, posmin, posmax, oldPosition, openNode, levelX);
        var height =  document.getElementById('mainpanel').clientHeight;
        cy.pan({x: cy.pan().x, y: height/3});
    }, 1);
}

function handleAlongChromosomesFilter() {
    createComponentShowList(drawAlongChromosome, function(i) {
        return chromosomes[i].name;
    }, function(i) {
        return "Chromosome " + chromosomes[i].name;
    }, chromosomes.length);
}

function calcY(curv, ypos, sumw) {
    return ypos/sumw;
}

function getYC_D() {
    return 1000/defZoom;
}

function getLeftOrder(curu, area_size, isGoodEdge, curNodeSet, openNode, cy) {
    var mnOrder = -1;
    for (var i = 0; i < scaffoldgraph.g[curu].length; ++i) {
        var curedge = scaffoldgraph.g[curu][i];
        if (isGoodEdge(curedge.id)) {
            var u = curedge.to;
            if (curNodeSet.has(u)) {
                if (mnOrder === -1 || mnOrder > cy.getElementById(u).data('order')) {
                    mnOrder = cy.getElementById(u).data('order');
                }
            }
        }
    }

    for (i = 0; i < scaffoldgraph.gr[curu].length; ++i) {
        curedge = scaffoldgraph.gr[curu][i];
        if (isGoodEdge(curedge.id)) {
            u = curedge.from;
            if (curNodeSet.has(u)) {
                if (mnOrder === -1 || mnOrder > cy.getElementById(u).data('order')) {
                    mnOrder = cy.getElementById(u).data('order');
                }
            }
        }
    }
    return mnOrder;
}

function getRightOrder(curu, area_size, isGoodEdge, curNodeSet, openNode, cy) {
    var mxOrder = -1;
    for (var i = 0; i < scaffoldgraph.g[curu].length; ++i) {
        var curedge = scaffoldgraph.g[curu][i];
        if (isGoodEdge(curedge.id)) {
            var u = curedge.to;
            if (curNodeSet.has(u)) {
                if (mxOrder === -1 || mxOrder < cy.getElementById(u).data('order')) {
                    mxOrder = cy.getElementById(u).data('order');
                }
            }
        }
    }

    for (i = 0; i < scaffoldgraph.gr[curu].length; ++i) {
        curedge = scaffoldgraph.gr[curu][i];
        if (isGoodEdge(curedge.id)) {
            u = curedge.from;
            if (curNodeSet.has(u)) {
                if (mxOrder === -1 || mxOrder < cy.getElementById(u).data('order')) {
                    mxOrder = cy.getElementById(u).data('order');
                }
            }
        }
    }
    return mxOrder;
}

function handleVertexConnectedToAlignContig(curu, curv, curedge, curNodeSet, cy, used_id, isGoodEdge,
                                            area_size, newvert, que, ypos, sumw, process, curd, openNode, rank, maxRank) {
    if (!isBigContig(0, scaffoldgraph.nodes[curu].len, defZoom)) {
        return -1;
    }
    var rightPos = getRightOrder(curu, area_size, isGoodEdge, curNodeSet, openNode, cy);
    if (cy.getElementById(curv).data('order') !== rightPos) {
        return -1;
    }

    var leftPos = getLeftOrder(curu, area_size, isGoodEdge, curNodeSet, openNode, cy);
    var maxR = 0;
    for (var i = leftPos; i < rightPos; ++i) {
        if ((i in maxRank) && maxRank[i] > maxR) {
            maxR = maxRank[i]
        }
    }

    maxRank[rightPos] = maxR + 1;

    return maxR;
}

function handleEdge(curedge, curv, curNodeSet, cy, used_id, isGoodEdge, area_size, newvert, que, ypos, sumw, process, curd, openNode, rank, maxRank, v) {
    var isForward = (curedge.to !== curv);
    if (isGoodEdge(curedge.id)) {
        var curu = curedge.to;
        if (!isForward) {
            curu = curedge.from;
        }

        var curdd = curd;
        if (!curNodeSet.has(curu) && curNodeSet.has(curv)) {
            curdd = handleVertexConnectedToAlignContig(curu, curv, curedge, curNodeSet, cy, used_id, isGoodEdge,
                area_size, newvert, que, ypos, sumw, process, curd, openNode, rank, maxRank);
        }

        if (isBigContig(0, scaffoldgraph.nodes[curu].len, defZoom) && (!curNodeSet.has(curu) || contigHasEdgesInThisScala(cy, curu))) {
            if (curdd !== -1) {
                if ((!used_id.has(curu)) && (curd < area_size || curNodeSet.has(curu) || openNode.has(curu))) {
                    rank[curu] = curdd + 1;
                    used_id.add(curu);
                    newvert.push(curu);
                    que.push({id: curu});
                }
            }

            if (!(curu in ypos)) {
                ypos[curu] = 0;
                sumw[curu] = 0;
            }

            if (!process.has(curu)) {
                var yc = 0;
                if (curd === 0) {
                    yc = v.ce;
                    if (!isForward) {
                        yc = v.cb;
                    }
                } else {
                    yc = calcY(curv, ypos[curv], sumw[curv]);
                }
                if (isForward) {
                    yc = yc + getYC_D();
                } else {
                    yc = yc - getYC_D();
                }

                ypos[curu] += /*curedge.weight*/ yc;
                sumw[curu] += /*curedge.weight*/ 1;
            }

            if (curdd !== -1) {
                //TODO: mistake calculate real dist, not rank;
                if ((curd < area_size || curNodeSet.has(curu) || openNode.has(curu) || (curd === area_size && used_id.has(curu)))) {
                    edges_to_draw.push(curedge.id);
                }
            }
        }
    }
}

function findNodeAroundChr(inode, area_size, min_contig_len, isGoodEdge, curNodeSet, openNode, cy) {
    var maxRank = {};
    var ypos = {};
    var newvert = [];
    var sumw = {};
    var rank = {};

    var used_id = new Set();
    var que = [];
    var process = new Set();

    inode.sort(function (a, b) {
        return a.order - b.order;
    });

    for (var i = 0; i < inode.length; ++i) {
        var v = inode[i];

        if (contigHasEdgesInThisScala(cy, v.id)) {
            used_id.add(v.id);
            rank[v.id] = 0;
            que.push(v);
        }
    }

    var bg = 0;
    while (bg < que.length) {
        v = que[bg];
        bg += 1;

        var curd = rank[v.id];
        var curv = v.id;
        process.add(curv);
        for (i = 0; i < scaffoldgraph.g[curv].length; ++i) {
            var curedge = scaffoldgraph.g[curv][i];
            handleEdge(curedge, curv, curNodeSet, cy, used_id, isGoodEdge, area_size,
                newvert, que, ypos, sumw, process, curd, openNode, rank, maxRank, v);
        }

        for (i = 0; i < scaffoldgraph.gr[curv].length; ++i) {
            curedge = scaffoldgraph.gr[curv][i];
            handleEdge(curedge, curv, curNodeSet, cy, used_id, isGoodEdge, area_size,
                newvert, que, ypos, sumw, process, curd, openNode, rank, maxRank, v);
        }

        edges_to_draw = edges_to_draw.filter(function (value, index, self) {
            return self.indexOf(value) === index;
        });
    }

    var res = [];

    for (i = 0; i < newvert.length; ++i) {
        res.push({id: newvert[i], rank: rank[newvert[i]], x: calcY(newvert[i], ypos[newvert[i]], sumw[newvert[i]])});
    }

    return res;
}
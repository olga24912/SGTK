/*
* Find local area around list of nodes
*/

function findLocalArea(inodes, area_size, min_contig_len, isGoodEdge) {
    var dist = new Map();
    var queue = [];
    var bp = 0;
    var i = 0;
    for (i=0; i < inodes.length; ++i) {
        queue.push(inodes[i]);
        dist.set(inodes[i], 0);
    }

    nodes_to_draw = [];
    edges_to_draw = [];

    while (bp < queue.length) {
        var curv = queue[bp];
        ++bp;
        var curd = dist.get(curv);
        if ((curd <= area_size) && (scaffoldgraph.nodes[curv].len >= min_contig_len)) {
            nodes_to_draw.push(curv);
            for (i = 0; i < scaffoldgraph.g[curv].length; ++i) {
                var curedge = scaffoldgraph.g[curv][i];
                if (isGoodEdge(curedge.id)) {
                    var curu = curedge.to;

                    if (!dist.has(curu)) {
                        queue.push(curu);
                        dist.set(curu, curd + 1);
                    }

                    if (dist.get(curu) <= area_size) {
                        edges_to_draw.push(curedge.id);
                    }
                }
            }


            for (i = 0; i < scaffoldgraph.gr[curv].length; ++i) {
                curedge = scaffoldgraph.gr[curv][i];
                if (isGoodEdge(curedge.id)) {
                    curu = curedge.from;

                    if (!dist.has(curu)) {
                        queue.push(curu);
                        dist.set(curu, curd + 1);
                    }
                }
            }
        }
    }
}


function findConnectedVertex(v, needAddVert, needAddEdge, newNode, area_size, min_contig_len, isGoodEdge, curNodeSet) {
    for (var h = 0; h < scaffoldgraph.g[v].length; ++h) {
        if (isGoodEdge(scaffoldgraph.g[v][h].id)) {
            if (scaffoldgraph.nodes[scaffoldgraph.g[v][h].to].len >= min_contig_len) {
                if (!curNodeSet.has(scaffoldgraph.g[v][h].to)) {
                    curNodeSet.add(scaffoldgraph.g[v][h].to);
                    newNode.add(scaffoldgraph.g[v][h].to);
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
                    newNode.add(scaffoldgraph.gr[v][h].from);
                    needAddVert.push(scaffoldgraph.gr[v][h].from);
                    nodes_to_draw.push(scaffoldgraph.gr[v][h].from);
                }
            }
        }
    }

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
}
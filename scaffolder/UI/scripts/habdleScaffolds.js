var libnum = 0;
var global_areasize = 0;
var global_min_contig_len = 0;
var global_isGoodEdge;
var nodes_to_draw = [];
var edges_to_draw = [];

function drawScaffold(i) {
    special_edges.clear();
    special_nodes.clear();
    var inodes = [];

    inodes.push(scaffoldgraph.libs[libnum].scaffolds[i].edges[0].from);
    special_nodes.add(inodes[0]);
    for (var g = 0; g < scaffoldgraph.libs[libnum].scaffolds[i].edges.length; ++g) {
        var e = scaffoldgraph.libs[libnum].scaffolds[i].edges[g];
        inodes.push(e.to);
        special_nodes.add(e.to);
        special_edges.add(e.id);
    }

    findLocalArea(inodes, global_areasize, global_min_contig_len, global_isGoodEdge);

    for (g = 0; g < scaffoldgraph.libs[libnum].scaffolds[i].edges.length; ++g) {
        e = scaffoldgraph.libs[libnum].scaffolds[i].edges[g];
        if (scaffoldgraph.nodes[e.to].len >= global_min_contig_len && scaffoldgraph.nodes[e.from].len < global_min_contig_len) {
            edges_to_draw.push(e.id);
        }
    }


    DrawGraph(nodes_to_draw, edges_to_draw);
    cy.fit(cy.$('#' + scaffoldgraph.libs[libnum].scaffolds[i].edges[0].from));
}

function handleScaffoldsFilter(scafflibname, areasize, min_contig_len, isGoodEdge) {
    for (var j = 0; j < scaffoldgraph.libs.length; ++j) {
        if (scaffoldgraph.libs[j].name == scafflibname) {
            libnum = j;
        }
    }

    global_areasize = areasize;
    global_min_contig_len = min_contig_len;
    global_isGoodEdge = isGoodEdge;

    notzero = [];
    for (j = 0; j < scaffoldgraph.libs[libnum].scaffolds.length; ++j) {
        var vertCnt = 0;
        for (var g = 0; g < scaffoldgraph.libs[libnum].scaffolds[j].edges.length; ++g) {
            var e = scaffoldgraph.libs[libnum].scaffolds[j].edges[g];
            if (scaffoldgraph.nodes[e.from].len >= min_contig_len) {
                ++vertCnt;
            }
            if (scaffoldgraph.nodes[e.to].len >= min_contig_len) {
                ++vertCnt;
            }
        }
        if (scaffoldgraph.libs[libnum].scaffolds[j].edges.length > 0 && vertCnt > 0) {
            notzero.push(j);
        }
    }

    createComponentShowList(function(i) {
            drawScaffold(notzero[i]);
        }, function(i) {
            var vertCnt = 0;
            var j = notzero[i];
            var e = scaffoldgraph.libs[libnum].scaffolds[j].edges[0];
            if (scaffoldgraph.nodes[e.from].len >= min_contig_len) {
                ++vertCnt;
            }
            for (var g = 0; g < scaffoldgraph.libs[libnum].scaffolds[j].edges.length; ++g) {
                e = scaffoldgraph.libs[libnum].scaffolds[j].edges[g];
                if (scaffoldgraph.nodes[e.to].len >= min_contig_len) {
                    ++vertCnt;
                }
            }
            return scaffoldgraph.libs[libnum].scaffolds[notzero[i]].name + "<br/> scaffold len: " + vertCnt;
        }, function(i) {
            return "Scaffold " + scaffoldgraph.libs[libnum].scaffolds[notzero[i]].name;
        }, notzero.length
    );
}

function getEdgeFrom(e) {
    var v = scaffoldgraph.edges[e].from;
    var edges_name = scaffoldgraph.edges[e].name;
    var lib = scaffoldgraph.edges[e].lib;
    while (scaffoldgraph.nodes[v].len < min_contig_len) {
        var nxte = -1;

        for (var i = 0; i < scaffoldgraph.gr[v].length; ++i) {
            var ne = scaffoldgraph.gr[v][i].id;
            if (scaffoldgraph.edges[ne].name === edges_name && scaffoldgraph.edges[e].lib === lib) {
                nxte = ne;
            }
        }

        if (nxte === -1) {
            return v;
        } else {
            v = scaffoldgraph.edges[nxte].from;
        }
    }

    return v;
}

function getEdgeTo(e) {
    var v = scaffoldgraph.edges[e].to;
    var edges_name = scaffoldgraph.edges[e].name;
    var lib = scaffoldgraph.edges[e].lib;
    while (scaffoldgraph.nodes[v].len < min_contig_len) {
        var nxte = -1;

        for (var i = 0; i < scaffoldgraph.g[v].length; ++i) {
            var ne = scaffoldgraph.g[v][i].id;
            if (scaffoldgraph.edges[ne].name === edges_name && scaffoldgraph.edges[e].lib === lib) {
                nxte = ne;
            }
        }

        if (nxte === -1) {
            return v;
        } else {
            v = scaffoldgraph.edges[nxte].to;
        }
    }

    return v;
}
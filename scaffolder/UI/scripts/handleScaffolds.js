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

function containWrong(libnum, j) {
    for (var g = 0; g < scaffoldgraph.libs[libnum].scaffolds[j].edges.length; ++g) {
        var e = scaffoldgraph.libs[libnum].scaffolds[j].edges[g].id;
        var v = getEdgeFrom(e);
        var u = getEdgeTo(e);
        if (!isCorrectOrder(v, u) && scaffoldgraph.nodes[v].len >= min_contig_len && scaffoldgraph.nodes[u].len >= min_contig_len) {
            console.log("wrong " + v.toString() + " " + u.toString())
            return true;
        }
    }
    return false;
}


function containContinuation(libnum, j) {
    global_min_contig_len = min_contig_len;
    global_isGoodEdge = isGoodEdge;

    var fv = -1;
    for (var g = 0; g < scaffoldgraph.libs[libnum].scaffolds[j].edges.length; ++g) {
        var e = scaffoldgraph.libs[libnum].scaffolds[j].edges[g].id;
        var v = getEdgeFrom(e);
        if (scaffoldgraph.nodes[v].len >= min_contig_len) {
            fv = v;
            break;
        }
    }

    var lv = -1;
    for (g = scaffoldgraph.libs[libnum].scaffolds[j].edges.length - 1; g >= 0; --g) {
        e = scaffoldgraph.libs[libnum].scaffolds[j].edges[g].id;
        v = getEdgeTo(e);
        if (scaffoldgraph.nodes[v].len >= min_contig_len) {
            lv = v;
            break;
        }
    }

    if (fv != -1) {
        for (g = 0; g < scaffoldgraph.gr[fv].length; ++g) {
            if (isGoodEdge(scaffoldgraph.gr[fv][g].id) &&
            scaffoldgraph.nodes[scaffoldgraph.gr[fv][g].from].len >= min_contig_len) {
                if (isCorrectOrder(scaffoldgraph.gr[fv][g].from, fv)) {
                    return true;
                }
            }

        }
    }

    if (lv != -1) {
        for (g = 0; g < scaffoldgraph.g[lv].length; ++g) {
            if (isGoodEdge(scaffoldgraph.g[lv][g].id) &&
                scaffoldgraph.nodes[scaffoldgraph.g[lv][g].to].len >= min_contig_len) {
                if (isCorrectOrder(lv, scaffoldgraph.g[lv][g].to)) {
                    return true;
                }
            }
        }
    }
    return false;
}

function containAmbiguous(libnum, j) {
    for (var g = 0; g < scaffoldgraph.libs[libnum].scaffolds[j].edges.length; ++g) {
        var e = scaffoldgraph.libs[libnum].scaffolds[j].edges[g].id;
        var v = getEdgeFrom(e);
        var u = getEdgeTo(e);

        if (scaffoldgraph.nodes[v].len >= min_contig_len && scaffoldgraph.nodes[u].len >= min_contig_len) {
            for (var h = 0; h < scaffoldgraph.g[v].length; ++h) {
                if (isGoodEdge(scaffoldgraph.g[v][h].id) && scaffoldgraph.g[v][h].to != u) {
                    return true;
                }
            }

            for (h = 0; h < scaffoldgraph.gr[u].length; ++h) {
                if (isGoodEdge(scaffoldgraph.gr[u][h].id) && scaffoldgraph.gr[u][h].from != v) {
                    return true;
                }
            }
        }
    }
    return false;
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
            if (g === 0 && scaffoldgraph.nodes[e.from].len >= min_contig_len) {
                ++vertCnt;
            }
            if (scaffoldgraph.nodes[e.to].len >= min_contig_len) {
                ++vertCnt;
            }
        }
        if (scaffoldgraph.libs[libnum].scaffolds[j].edges.length > 0 && vertCnt >= document.getElementById("min_scaffold_len").value) {
            if (!document.getElementById("scaff_wrng").checked && !document.getElementById("scaff_cont").checked && !document.getElementById("scaff_ambig").checked) {
                notzero.push(j);
            } else if (document.getElementById("scaff_wrng").checked && containWrong(libnum, j)) {
                notzero.push(j);
            } else if (document.getElementById("scaff_cont").checked && containContinuation(libnum, j)) {
                notzero.push(j);
            } else if (document.getElementById("scaff_ambig").checked && containAmbiguous(libnum, j)) {
                notzero.push(j);
            }
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
    var nm = scaffoldgraph.edges[e].num;
    while (scaffoldgraph.nodes[v].len < min_contig_len) {
        var nxte = -1;

        for (var i = 0; i < scaffoldgraph.gr[v].length; ++i) {
            var ne = scaffoldgraph.gr[v][i].id;
            if (scaffoldgraph.edges[ne].name === edges_name &&
                scaffoldgraph.edges[ne].lib === lib &&
                scaffoldgraph.edges[ne].num === nm - 1) {
                nxte = ne;
            }
        }

        if (nxte === -1) {
            return v;
        } else {
            v = scaffoldgraph.edges[nxte].from;
            nm -= 1;
        }
    }

    return v;
}

function getEdgeTo(e) {
    var v = scaffoldgraph.edges[e].to;
    var edges_name = scaffoldgraph.edges[e].name;
    var lib = scaffoldgraph.edges[e].lib;
    var nm = scaffoldgraph.edges[e].num;
    while (scaffoldgraph.nodes[v].len < min_contig_len) {
        var nxte = -1;

        for (var i = 0; i < scaffoldgraph.g[v].length; ++i) {
            var ne = scaffoldgraph.g[v][i].id;
            if (scaffoldgraph.edges[ne].name === edges_name && scaffoldgraph.edges[ne].lib === lib &&
                scaffoldgraph.edges[ne].num === nm + 1) {
                nxte = ne;
            }
        }

        if (nxte === -1 || nm === -1) {
            return v;
        } else {
            v = scaffoldgraph.edges[nxte].to;
            nm += 1;
        }
    }

    return v;
}
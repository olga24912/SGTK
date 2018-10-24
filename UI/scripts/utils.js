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

function elemInList(elem, lst) {
    for (var i = 0; i < lst.length; ++i) {
        if (elem == lst[i]) {
            return true;
        }
    }

    return false;
}

function isCorrectOrder(v, u) {
    var nv = scaffoldgraph.nodes[v];
    var nu = scaffoldgraph.nodes[u];

    for (var i = 0; i < nv.alignments.length; ++i) {
        //if (nv.alignments[i].coordne + 1000 >= nv.len) {
        for (var j = 0; j < nu.alignments.length; ++j) {
            //if (nu.alignments[j].coordnb - 1000 <= 0) {
            if (nv.alignments[i].chr_id === nu.alignments[j].chr_id) {
                if (nv.alignments[i].coorde - 100 < nu.alignments[j].coordb) {
                    if (nu.alignments[j].coordb - nv.alignments[i].coorde < 5000) {
                        return true;
                    }
                }
            }
            //}
        }
        //}
    }

    return false;
}

function updateText() {
    var nodes = cy.nodes();
    for (var i = 0; i < nodes.length; ++i) {
        cy.$('#' + nodes[i].id().toString()).data('label', createLabelForNode(nodes[i].id()));
    }
    var edges = cy.edges();
    for (var i = 0; i < edges.length; ++i) {
        cy.$('#' + edges[i].id().toString()).data('label', createLabelForEdge(parseInt(edges[i].id().substring(1))));
    }
}
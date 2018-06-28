function getCntVertWithLibConnection(libNum) {
    var cnt = 0;
    var used = 0;
    for (var i = 0; i < scaffoldgraph.nodes.length; ++i) {
        var v = scaffoldgraph.nodes[i].id;
        used = 0;
        for (var j = 0; j < scaffoldgraph.g[v].length && used===0; ++j) {
            if (scaffoldgraph.g[v][j].lib === libNum) {
                cnt += 1;
                used = 1;
            }
        }
    }
    return cnt;
}

function getDefaultWeight(libNum) {
    var vertCnt = getCntVertWithLibConnection(libNum);
    var listWeightEdges = [];
    for (var i = 0; i < scaffoldgraph.edges.length; ++i) {
        if (scaffoldgraph.edges[i].lib === libNum) {
            listWeightEdges.push(scaffoldgraph.edges[i].weight);
        }
    }

    listWeightEdges = listWeightEdges.sort(function(a, b) {
        return b - a;
    });

    if (listWeightEdges.length === 0) {
        return 1;
    } else if (listWeightEdges.length <= vertCnt) {
        return 1;
    } else {
        var pos = parseInt(vertCnt);
        while (pos > 0 && listWeightEdges[pos] === listWeightEdges[parseInt(vertCnt)]) {
            --pos;
        }
        return listWeightEdges[pos];
    }
}
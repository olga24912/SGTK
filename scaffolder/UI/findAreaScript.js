function findLocalArea(inodes, area_size, min_contig_len, isGoodEdge) {
    alert("find local area");
    var dist = new Map();
    var queue = [];
    var bp = 0;
    var i = 0;
    for (i=0; i < inodes.length; ++i) {
        queue.push(inodes[i]);
        dist.set(inodes[i], 0);
        alert(inodes[i].toString() + " " + dist.has(inodes[i]).toString());
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

                    alert(curu.toString() + " " + dist.has(curu).toString());

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

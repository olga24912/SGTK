function handleAmbiguousFilter(areasize, min_contig_len, isGoodEdge) {
    special_nodes.clear();
    special_edges.clear();

    var nodes_id = [];

    for (var i = 0; i < scaffoldgraph.nodes.length; ++i) {
        if (scaffoldgraph.nodes[i].len < min_contig_len) {
            continue;
        }
        var edgesSortByTo = [];
        for (var j = 0; j < scaffoldgraph.g[i].length; ++j) {
            if (isGoodEdge(scaffoldgraph.g[i][j].id)) {
                edgesSortByTo.push(scaffoldgraph.g[i][j].to);
            }
        }

        edgesSortByTo.sort(function (a, b) {
            return a - b;
        });

        edgesSortByTo.push(-1);

        var lstval = -1;
        var cntTo = 0;
        var hasRight = 0;

        for (j = 0; j < edgesSortByTo.length; ++j) {
            if (!(lstval === -1) && !(lstval === edgesSortByTo[j])) {
                cntTo += 1;
                if (isCorrectOrder(i, lstval) || chromosomes.length === 0) {
                    hasRight += 1;
                }
            }

            lstval = edgesSortByTo[j];
        }
        if (cntTo > 1 && hasRight > 0) {
            nodes_id.push(i);
            special_nodes.add(i);

            for (j = 0; j < edgesSortByTo.length; ++j) {
                if (!(lstval === -1) && !(lstval === edgesSortByTo[j])) {
                    nodes_id.push(lstval);
                }

                lstval = edgesSortByTo[j];
            }
        }
    }

    var uniqueNode = nodes_id.filter(function (value, index, self) {
        return self.indexOf(value) === index;
    });

    findLocalArea(uniqueNode, areasize, min_contig_len, isGoodEdge);

    console.log(nodes_id.length);
    console.log(nodes_to_draw.length);
    splitOnParts(nodes_to_draw, edges_to_draw);

    createComponentShowList(function(i) {
        DrawGraph(nodes_set[i], edges_set[i]);
    }, function(i) {
        return "comp " + i + "<br/> #nodes = " + nodes_set[i].length + "<br/>#edges = " + edges_set[i].length;
    }, function(compnum) {
        return "Component #" + compnum;
    }, nodes_set.length);
}
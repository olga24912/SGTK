/*
* Find components to draw
*/

//Maximum edges in one component
MAX_EDGES_CNT_IN_ONE_PART = 200;

//Maximum nodes in one component
MAX_NODES_CNT_IN_ONE_PART = 100;


//Find component around edgeID ev(int) (MAX_EDGES_CNT_IN_ONE_PART nearest edges)
function findArea(ev, g, color, curc) {
    var edge_cnt = 0;
    var queue = [];
    queue.push(ev);
    color.set(ev, curc);
    edge_cnt += 1;
    var bg = 0;
    while (bg < queue.length) {
        ev = queue[bg];
        bg += 1;

        var nb = g.get(scaffoldgraph.edges[ev].from);
        for (i=0; i < nb.length; ++i) {
            var u = nb[i];
            if (color.get(u) === -1) {
                queue.push(u);
                color.set(u, curc);
                edge_cnt += 1;
            }
        }

        nb = g.get(scaffoldgraph.edges[ev].to);

        for (var i=0; i < nb.length; ++i) {
            u = nb[i];
            if (color.get(u) === -1) {
                queue.push(u);
                color.set(u, curc);
                edge_cnt += 1;
            }
        }


        if (edge_cnt > MAX_EDGES_CNT_IN_ONE_PART) {
            return;
        }
    }
}

//Split graph on components. nodes_to draw = NodesId(int[]), edges_to_draw = EdgesId(int[])
function splitOnParts(nodes_to_draw, edges_to_draw) {
    var g = new Map();
    var color = new Map();
    for (var i = 0; i < edges_to_draw.length; ++i) {
        color.set(edges_to_draw[i], -1);
    }
    for (i = 0; i < nodes_to_draw.length; ++i) {
        g.set(nodes_to_draw[i], []);
    }

    for (i=0; i < edges_to_draw.length; ++i) {
        g.get(scaffoldgraph.edges[edges_to_draw[i]].from).push(edges_to_draw[i]);
        g.get(scaffoldgraph.edges[edges_to_draw[i]].to).push(edges_to_draw[i]);
    }

    var curc = 0;
    for (i=0; i < edges_to_draw.length; ++i) {
        if (color.get(edges_to_draw[i]) === -1) {
            findArea(edges_to_draw[i], g, color, curc);
            ++curc;
        }
    }

    nodes_set = [];
    edges_set = [];
    for (i=0; i < curc; ++i) {
        nodes_set.push([]);
        edges_set.push([]);
    }

    for(i=0; i < edges_to_draw.length; ++i) {
        edges_set[color.get(edges_to_draw[i])].push(edges_to_draw[i]);

        if (!elemInList(scaffoldgraph.edges[edges_to_draw[i]].from, nodes_set[color.get(edges_to_draw[i])])) {
            nodes_set[color.get(edges_to_draw[i])].push(scaffoldgraph.edges[edges_to_draw[i]].from);
        }
        if (!elemInList(scaffoldgraph.edges[edges_to_draw[i]].to, nodes_set[color.get(edges_to_draw[i])])) {
            nodes_set[color.get(edges_to_draw[i])].push(scaffoldgraph.edges[edges_to_draw[i]].to);
        }
    }
}
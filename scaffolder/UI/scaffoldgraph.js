function ScaffoldGraph(libs, nodes, edges) {
    this.libs = libs;
    this.nodes = nodes;
    this.edges = edges;

    this.id_by_name = new Map();
    this.g = [];
    this.gr = [];

    var i = 0;
    for (i=0; i < this.nodes.length; ++i) {
        this.id_by_name.set(this.nodes[i].name, this.nodes[i].id);
        this.g.push([]);
        this.gr.push([]);
    }

    for (i = 0; i < this.edges.length; ++i) {
        this.g[this.edges[i].from].push(this.edges[i]);
        this.gr[this.edges[i].to].push(this.edges[i]);
    }
}

function ScaffoldEdge(id, from, to, lib, weight) {
    this.id = id;
    this.from = from;
    this.to = to;
    this.lib = lib;
    this.weight = weight;
}

function ScaffoldNode(id, name, len) {
    this.id = id;
    this.name = name;
    this.len = len;
}

function ScaffoldEdgeLib(id, color, name, type) {
    this.id = id;
    this.color = color;
    this.name = name;
    this.type = type;
}

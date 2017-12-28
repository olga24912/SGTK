function ScaffoldGraph(libs, nodes, edges) {
    this.libs = libs;
    this.nodes = nodes;
    this.edges = edges;
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
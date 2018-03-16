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
    this.name = "";
}

function Scaffold(name) {
    this.name = name;
    this.edges = [];
}

function ScaffoldNode(id, name, len) {
    this.id = id;
    this.name = name;
    this.len = len;
    this.alignments = [];
}

function ScaffoldEdgeLib(id, color, name, type) {
    this.id = id;
    this.color = color;
    this.name = name;
    this.type = type;
    this.scaffolds = [];
}

function Chromosome(id, name, len) {
    this.id = id;
    this.color = "#"+((1<<24)*Math.random()|0).toString(16);
    this.name = name;
    this.len = len;
    this.alignments = [];
}

function Alignment(coord_begin, coord_end, chr_id, node_id/*, coordnb, coordne*/) {
    this.coordb = coord_begin;
    this.coorde = coord_end;
    //this.coordnb = coordnb;
    //this.coordne = coordne;
    this.chr_id = chr_id;
    this.node_id = node_id;
}
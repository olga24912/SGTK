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

var scaffoldlibs = [new ScaffoldEdgeLib(0, '#0000ff', 'lib1', 'RNA_PAIR'), new ScaffoldEdgeLib(1, '#00ff00', 'lib2', 'SCAFF'), new ScaffoldEdgeLib(2, '#ff0000', 'lib3', 'SCAFF')];
var scaffoldnodes = [new ScaffoldNode(0, "node0", 10), new ScaffoldNode(1, "node1", 15), new ScaffoldNode(2, "node2", 20)];
var scaffoldedges = [new ScaffoldEdge(0, 0, 1, 0, 1), new ScaffoldEdge(1, 1, 2, 2, 5)];

var scaffoldgraph = new ScaffoldGraph(scaffoldlibs, scaffoldnodes, scaffoldedges);

function DrawGraph(nodes_to_draw, edges_to_draw) {
    var nodeslist = [];
    var edgeslist = [];
    var i = 0;
    var j = 0;
    for (i=0; i < nodes_to_draw.length; i++) {
        j = nodes_to_draw[i];
        nodeslist.push({id: scaffoldgraph.nodes[j].id, label: scaffoldgraph.nodes[j].name});
    }

    // create an array with nodes
    var nodes = new vis.DataSet(nodeslist);

    for (i=0; i < edges_to_draw.length; i++) {
        j = edges_to_draw[i];
        edgeslist.push({from: scaffoldgraph.edges[j].from, to: scaffoldgraph.edges[j].to, arrows: 'to', color:{color: scaffoldgraph.libs[scaffoldgraph.edges[j].lib].color}});
    }

    // create an array with edges
    var edges = new vis.DataSet(edgeslist);

    // create a network
    var container = document.getElementById('mainpanel');

    // provide the data in the vis format
    var data = {
        nodes: nodes,
        edges: edges
    };
    var options = {};

    // initialize your network!
    var network = new vis.Network(container, data, options);
}

function drawLocalArea(inodes, area_size) {
    var dist = new Map();
    var queue = [];
    var bp = 0;
    var i = 0;
    for (i=0; i < inodes.length; ++i) {
        queue.push(inodes[i]);
        dist.set(inodes[i], 0);
    }

    var nodes_to_draw = [];
    var edges_to_draw = [];

    while (bp < queue.length) {
        var curv = queue[bp];
        ++bp;
        var curd = dist.get(curv);
        //alert("curv " + curv + " curd " + curd);
        if (curd <= area_size) {
            nodes_to_draw.push(curv);
            for (i = 0; i < scaffoldgraph.g[curv].length; ++i) {
                var curedge = scaffoldgraph.g[curv][i];
                var curu = curedge.to;

                if (!dist.has(curu)) {
                    queue.push(curu);
                    dist.set(curu, curd + 1);
                }

                if (dist.get(curu) <= area_size) {
                    edges_to_draw.push(curedge.id);
                }
            }


            for (i = 0; i < scaffoldgraph.gr[curv].length; ++i) {
                curedge = scaffoldgraph.gr[curv][i];
                curu = curedge.from;

                if (!dist.has(curu)) {
                    queue.push(curu);
                    dist.set(curu, curd + 1);
                }
            }
        }
    }

    DrawGraph(nodes_to_draw, edges_to_draw);
}

function handleFilterButton() {
    var opt = document.getElementById("select_show_type").value;
    if (opt=="vertices_local_area") {
        var areasize = document.getElementById("area_size").value;
        var nodes = document.getElementById("vertext").value.replace(/\n/g, " ").split(" ");
        var nodes_id = [];
        for (var i=0; i < nodes.length; ++i) {
            nodes_id.push(scaffoldgraph.id_by_name.get(nodes[i]));
        }
        drawLocalArea(nodes_id, areasize);
    }
}

function InitLibTable() {
    var table = document.getElementById("lib_table");

    for (var i=0; i < scaffoldgraph.libs.length; ++i) {
        var tr = document.createElement("tr");

        var td_type = document.createElement("td");
        var lib_type = document.createElement("p");
        lib_type.appendChild(document.createTextNode(scaffoldgraph.libs[i].type));
        td_type.appendChild(lib_type);
        td_type.align="center";

        var td_id = document.createElement("td");
        var lib_id = document.createElement("p");
        lib_id.appendChild(document.createTextNode("l" + scaffoldgraph.libs[i].id));
        td_id.align="center";
        td_id.appendChild(lib_id);



        var td_name = document.createElement("td");
        var lib_name = document.createElement("p");
        lib_name.style.color = scaffoldgraph.libs[i].color;
        lib_name.appendChild(document.createTextNode(scaffoldgraph.libs[i].name));
        td_name.align="center";
        td_name.appendChild(lib_name);

        var td_min_edge_weight = document.createElement("td");
        var input_weight = document.createElement("input");
        input_weight.type = "number";
        input_weight.min = 0;
        td_min_edge_weight.align="center";
        td_min_edge_weight.appendChild(input_weight);

        tr.appendChild(td_type);
        tr.appendChild(td_id);
        tr.appendChild(td_name);
        tr.appendChild(td_min_edge_weight);
        table.appendChild(tr);
    }
}

DrawGraph([0, 1, 2], [0, 1]);
InitLibTable();
document.getElementById("filter_button").addEventListener("click", handleFilterButton);
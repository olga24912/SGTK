function fill_nodes(obj) {
    for (var i = 0; i < obj["nodes"].length; ++i) {
        scaffoldnodes.push(
            new ScaffoldNode(obj["nodes"][i]["id"], obj["nodes"][i]["name"], obj["nodes"][i]["len"]));
        if ("info" in obj["nodes"][i]) {
            scaffoldnodes[scaffoldnodes.length - 1].info = obj["nodes"][i]["info"]
        }
    }
}

function add_curedge(cur_edge, i) {
    scaffoldedges.push(new ScaffoldEdge(cur_edge["id"], cur_edge["from"], cur_edge["to"], i, cur_edge["weight"]));
    if ("len" in cur_edge) {
        scaffoldedges[scaffoldedges.length - 1].len = cur_edge["len"]
    }
    if ("info" in cur_edge) {
        scaffoldedges[scaffoldedges.length - 1].info = cur_edge["info"]
    }
}

function fill_libs(obj) {
    for (var i = 0; i < obj["libs"].length; ++i) {
        scaffoldlibs.push(
            new ScaffoldEdgeLib(obj["libs"][i]["id"], obj["libs"][i]["color"],
                obj["libs"][i]["name"], obj["libs"][i]["type"]));
        if ("edges" in obj["libs"][i]) {
            for (var j = 0; j < obj["libs"][i]["edges"].length; ++j) {
                var cur_edge = obj["libs"][i]["edges"][j];
                add_curedge(cur_edge, i)
            }
        }

        if ("scaffolds" in obj["libs"][i]) {
            for (var g = 0; g < obj["libs"][i]["scaffolds"].length; ++g) {
                var cur_scaf = obj["libs"][i]["scaffolds"][g];
                scaffoldlibs[i].scaffolds.push(new Scaffold(cur_scaf["name"]));
                for (var j = 0; j < cur_scaf["edges"].length; ++j) {
                    cur_edge = cur_scaf["edges"][j];
                    add_curedge(cur_edge, i);
                    scaffoldedges[scaffoldedges.length - 1].name = cur_scaf["name"];
                    scaffoldedges[scaffoldedges.length - 1].num = j;
                    scaffoldlibs[i].scaffolds[g].edges.push(scaffoldedges[scaffoldedges.length - 1])
                }
            }
        }
    }
}

function fill_chromosomes(obj) {
    for (var i = 0; i < obj["chromosomes"].length; ++i) {
        chromosomes.push(new Chromosome(obj["chromosomes"][i]["id"], obj["chromosomes"][i]["name"], obj["chromosomes"][i]["len"]));
    }
}

function fill_alignments(obj) {
    for (var i = 0; i < obj["alignments"].length; ++i) {
        var cur_algn = obj["alignments"][i];
        chromosomes[cur_algn["chr_id"]].alignments.push(
            new Alignment(cur_algn["coord_begin"], cur_algn["coord_end"], cur_algn["chr_id"], cur_algn["node_id"])
        );
    }
}

function fill_scaffoldGraph(obj) {
    scaffoldgraph = new ScaffoldGraph(scaffoldlibs, scaffoldnodes, scaffoldedges);
}

//var toJson = function(obj){ return obj.json(); };
//var json_graph = fetch('scripts/data.json').then(toJson);
fill_nodes(json_graph);
fill_libs(json_graph);
fill_chromosomes(json_graph);
fill_alignments(json_graph);
fill_scaffoldGraph(json_graph);
delete json_graph;
main();
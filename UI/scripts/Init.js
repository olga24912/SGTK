/*
* Init main page on start
*/


//Generate name for main page by TYPE name
function TypeToStr(s) {
    var typeToRep = {"DNA_PAIR": "Paired reads", "RNA_PAIR": "RNA-Seq (paired)", "RNA_SPLIT_50": "RNA-Seq (single)",
        "RNA_SPLIT_30": "RNA-Seq (single)", "SCAFF": "Scaffolds", "CONNECTION" : "Connection", "LONG" : "Long reads",
        "FASTG": "FASTG", "GFA": "GFA", "GFA2": "GFA2"};
    if (s in typeToRep) {
        return typeToRep[s]
    }
    return s.replace(/_/g, " ");
}

function getTdType(i) {
    var td_type = document.createElement("td");
    var lib_type = document.createElement("p");
    lib_type.appendChild(document.createTextNode(TypeToStr(scaffoldgraph.libs[i].type)));
    td_type.appendChild(lib_type);
    td_type.align="left";
    return td_type
}

function getTdName(i) {
    var td_name = document.createElement("td");
    var lib_name = document.createElement("p");
    lib_name.style.color = scaffoldgraph.libs[i].color;
    lib_name.appendChild(document.createTextNode(scaffoldgraph.libs[i].name));
    lib_name.id = "color" + scaffoldgraph.libs[i].name;
    td_name.align="center";
    td_name.appendChild(lib_name);
    return td_name;
}

//Init information about sources
function InitLibTable() {
    var table = document.getElementById("lib_table");
    var included_libs_tabel = document.getElementById("included_libs");

    for (var i=0; i < scaffoldgraph.libs.length; ++i) {
        var tr = document.createElement("tr");
        var tr_highlight = document.createElement("tr");



        var td_min_edge_weight = document.createElement("td");
        var input_weight = document.createElement("input");
        input_weight.type = "number";
        input_weight.min = 0;
        input_weight.size = 1;
        input_weight.value = getDefaultWeight(i);
        input_weight.id = "min_w_" + scaffoldgraph.libs[i].name;
        td_min_edge_weight.align="center";
        td_min_edge_weight.appendChild(input_weight);

        var td_check_box = document.createElement("td");
        var label = document.createElement("label");
        label.className = "container";
        var input_check_box = document.createElement("input");
        input_check_box.type = "checkbox";
        input_check_box.id = "include_lib_" + scaffoldgraph.libs[i].name;
        var span = document.createElement("span");
        span.className = "checkmark";
        span.style="left: 20px";
        label.appendChild(input_check_box);
        label.appendChild(span);
        td_check_box.appendChild(label);
        td_check_box.align="right";

        tr.appendChild(getTdType(i));
        tr.appendChild(getTdName(i));
        tr.appendChild(td_min_edge_weight);
        table.appendChild(tr);

        tr_highlight.appendChild(getTdType(i));
        tr_highlight.appendChild(getTdName(i));
        tr_highlight.appendChild(td_check_box);
        included_libs_tabel.appendChild(tr_highlight);
    }
}

//Generate alignments for nodes
function InitAlignmentsForNodes() {
    for (var i=0; i < chromosomes.length; ++i) {
        for (var j=0; j < chromosomes[i].alignments.length; ++j) {
            scaffoldgraph.nodes[chromosomes[i].alignments[j].node_id].alignments.push(chromosomes[i].alignments[j]);
        }
    }
}


//Genrate edges for scaffolds
function putEdgesNumForScaffolds() {
    for (var i = 0; i < scaffoldgraph.libs.length; ++i) {
        var lb = scaffoldgraph.libs[i];
        for (var j = 0; j < lb.scaffolds.length; ++j) {
            var scaffs = lb.scaffolds[j];
            for (var g = 0; g < scaffs.edges.length; ++g) {
                scaffs.edges[g].num = g;
                scaffs.edges[g].name = scaffs.name;
            }
        }
    }
}

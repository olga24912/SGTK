function createLabelForNode(node) {
    if (defZoom > 1000) {
        return "";
    }

    var label = "";
    if (document.getElementById("vert_checkbox_id").checked) {
        label += "id: " + scaffoldgraph.nodes[node].id + "\n";
    }
    if (document.getElementById("vert_checkbox_name").checked) {
        label += scaffoldgraph.nodes[node].name + "\n";
    }
    if (document.getElementById("vert_checkbox_len").checked) {
        label += "len: " + scaffoldgraph.nodes[node].len + "\n";
    }

    if (document.getElementById("vert_checkbox_align").checked) {
        if (scaffoldgraph.nodes[node].alignments.length > 0) {
            label += "Alignment: ";
            scaffoldgraph.nodes[node].alignments.sort(function (a, b) {
                return (b.coorde - b.coordb) - (a.coorde - a.coordb);
            });
            for (var i = 0; i < scaffoldgraph.nodes[node].alignments.length; ++i) {
                var cura = scaffoldgraph.nodes[node].alignments[i];
                label += chromosomes[cura.chr_id].name + " " + cura.coordb + " " + cura.coorde + " (" + ((cura.coorde - cura.coordb + 1) * 100/scaffoldgraph.nodes[node].len) +"%)\n";
            }
        }
    }
    if (document.getElementById("vert_checkbox_info").checked) {
        if (scaffoldgraph.nodes[node].info !== "") {
            label += scaffoldgraph.nodes[node].info + "\n";
        }
    }

    return label;
}

function createFullLabelForNode(node) {
    var label = "";
    label += "id: " + scaffoldgraph.nodes[node].id + "</br>";
    label += scaffoldgraph.nodes[node].name + "</br>";
    label += "len: " + scaffoldgraph.nodes[node].len + "</br>";
    if (scaffoldgraph.nodes[node].alignments.length > 0) {
        label += "Alignment: ";
        scaffoldgraph.nodes[node].alignments.sort(function (a, b) {
            return (b.coorde - b.coordb) - (a.coorde - a.coordb);
        });
        for (var i = 0; i < scaffoldgraph.nodes[node].alignments.length; ++i) {
            var cura = scaffoldgraph.nodes[node].alignments[i];
            label += chromosomes[cura.chr_id].name + " " + cura.coordb + " " + cura.coorde + " (" + ((cura.coorde - cura.coordb + 1) * 100/scaffoldgraph.nodes[node].len) +"%)</br>";
        }
    }
    if (scaffoldgraph.nodes[node].info !== "") {
        label += scaffoldgraph.nodes[node].info + "</br>";
    }
    return label;
}

function createLabelForEdge(edge) {
    if (defZoom > 1000) {
        return "";
    }

    var label = "";
    if (document.getElementById("edge_checkbox_id").checked) {
        label += "id: " + scaffoldgraph.edges[edge].id + "\n";
    }
    if (document.getElementById("edge_checkbox_name").checked) {
        label += scaffoldgraph.libs[scaffoldgraph.edges[edge].lib].name + "\n";
        if (scaffoldgraph.libs[scaffoldgraph.edges[edge].lib].type === "SCAFF") {
            label += scaffoldgraph.edges[edge].name + "\n";
        }
    }
    if (document.getElementById("edge_checkbox_weight").checked) {
        label += "w: " + scaffoldgraph.edges[edge].weight + "\n";
    }
    if (document.getElementById("edge_checkbox_type").checked) {
        label += scaffoldgraph.libs[scaffoldgraph.edges[edge].lib].type + "\n";
    }
    if (document.getElementById("edge_checkbox_len").checked) {
        if (scaffoldgraph.edges[edge].len >= 0) {
            label += "len: " + scaffoldgraph.edges[edge].len + "\n";
        }
    }
    if (document.getElementById("edge_checkbox_info").checked) {
        if (scaffoldgraph.edges[edge].info !== "") {
            label += scaffoldgraph.edges[edge].info + "\n";
        }
    }
    return label;
}

function createFullLabelForEdge(edge) {
    var label = "";
    label += "id: " + scaffoldgraph.edges[edge].id + "</br>";
    label += scaffoldgraph.libs[scaffoldgraph.edges[edge].lib].name + "</br>";
    if (scaffoldgraph.libs[scaffoldgraph.edges[edge].lib].type === "SCAFF") {
        label += scaffoldgraph.edges[edge].name + " " + scaffoldgraph.edges[edge].num + "</br>";
    }
    label += "w: " + scaffoldgraph.edges[edge].weight + "</br>";
    label += scaffoldgraph.libs[scaffoldgraph.edges[edge].lib].type + "</br>";
    if (scaffoldgraph.edges[edge].len >= 0) {
        label += "len: " + scaffoldgraph.edges[edge].len + "</br>";
    }
    if (scaffoldgraph.edges[edge].info !== "") {
        label += scaffoldgraph.edges[edge].info + "</br>";
    }
    return label;
}

function generateGeneralInfo() {
    return "Nodes: " + scaffoldgraph.nodes.length.toString() + "</br>" +
        "Edges: " + scaffoldgraph.edges.length.toString() + "</br>" +
        "Chromosomes: " + chromosomes.length.toString() + "</br>" +
        "Sources: " + scaffoldgraph.libs.length.toString() + "</br>";
}

function createInformationShown(cy) {
    var def_height = 140;

    cy.on('mouseover', 'node', function (evt) {
        var v = evt.target.id();
        var printInfo = createFullLabelForNode(v);
        document.getElementById("extra_info").style = "";
        document.getElementById("extra_info").innerHTML =
            "<p style='font-size: 14px; margin-top: 0px; margin-bottom: 0px;'>" + printInfo + "</p>";
        if (document.getElementById("extra_info").clientHeight <= def_height) {
            document.getElementById("extra_info").style.height = def_height + 'px';
        }
    });

    cy.on('mouseout', 'node', function (evt) {
        document.getElementById("extra_info").style.height = def_height + 'px';
        document.getElementById("extra_info").innerHTML = "<p style='margin-top: 0px; margin-bottom: 0px;'>" + generateGeneralInfo() + "</p>";
    });

    cy.on('mouseover', 'edge', function (evt) {
        var v = evt.target.id();
        var printInfo = createFullLabelForEdge(v.substring(1));
        document.getElementById("extra_info").style = "";
        document.getElementById("extra_info").innerHTML =
            "<p id='innerTextExtraInfo' style='font-size: 14px; margin-top: 0px; margin-bottom: 0px;'>" + printInfo + "</p>";
        console.log(document.getElementById("extra_info").clientHeight);
        if (document.getElementById("extra_info").clientHeight <= def_height) {
            document.getElementById("extra_info").style.height = def_height + 'px';
        }
    });

    cy.on('mouseout', 'edge', function (evt) {
        document.getElementById("extra_info").style.height = def_height + 'px';
        document.getElementById("extra_info").innerHTML = "<p style='margin-top: 0px; margin-bottom: 0px;'>" + generateGeneralInfo() + "</p>";
    });
}

var edges_set = [];
var nodes_set = [];
var special_nodes = new Set();
var special_edges = new Set();
var nodes_to_draw = [];
var edges_to_draw = [];
var isGoodEdge;
var min_contig_len = 0;
var area_size = 0;
var cur_show_id = 0;
var graph = null;
var cy = null;

function handleFullGraph() {
    for (i=0; i < scaffoldgraph.nodes.length; ++i) {
        if (scaffoldgraph.nodes[i].len >= min_contig_len) {
            nodes_to_draw.push(scaffoldgraph.nodes[i].id);
            special_nodes.add(scaffoldgraph.nodes[i].id);
        }
    }

    for (i=0; i < scaffoldgraph.edges.length; ++i) {
        if (isGoodEdge(i)) {
            edges_to_draw.push(scaffoldgraph.edges[i].id);
            special_edges.add(scaffoldgraph.edges[i].id);
        }
    }
    splitOnParts(nodes_to_draw, edges_to_draw);

    var notzero = [];
    for (i=0; i < edges_set.length; ++i) {
        if (edges_set[i].length > 0) {
            notzero.push(i);
        }
    }

    createComponentShowList(function(j) {
        var i = notzero[j];
        DrawGraph(nodes_set[i], edges_set[i]);
    }, function(j) {
        var i = notzero[j];
        return "comp " + i + "<br/> #nodes = " + nodes_set[i].length + "<br/>#edges = " + edges_set[i].length;
    }, function(i) {
        var compnum = notzero[i];
        return "Component #" + compnum;
    }, notzero.length);
}

function handleVertexLocalArea() {
    special_nodes.clear();
    special_edges.clear();
    area_size = document.getElementById("area_size").value;
    var nodes = document.getElementById("vertext").value.replace(/\n/g, " ").split(" ");
    var nodes_id = [];
    for (i = 0; i < nodes.length; ++i) {
        if (scaffoldgraph.id_by_name.has(nodes[i])) {
            nodes_id.push(scaffoldgraph.id_by_name.get(nodes[i]));
        } else {
            nodes_id.push(parseInt(nodes[i]));
        }
        special_nodes.add(nodes_id[nodes_id.length - 1]);
    }

    findLocalArea(nodes_id, area_size, min_contig_len, isGoodEdge);

    splitOnParts(nodes_to_draw, edges_to_draw);

    createComponentShowList(function (i) {
        DrawGraph(nodes_set[i], edges_set[i]);
        if (nodes_id.length > 0) {
            cy.fit(cy.$('#' + nodes_id[0]));
        }
    }, function (i) {
        var nm = "";
        for (var j = 0; j < nodes_set[i].length; ++j) {
            if (special_nodes.has(nodes_set[i][j])) {
                nm += scaffoldgraph.nodes[nodes_set[i][j]].name + " id=" + nodes_set[i][j].toString() + "<br/>";
            }
        }
        return nm;
    }, function (compnum) {
        var nm = "";
        for (var j = 0; j < nodes_set[compnum].length; ++j) {
            if (special_nodes.has(nodes_set[compnum][j])) {
                nm += scaffoldgraph.nodes[nodes_set[compnum][j]].name + " id=" + nodes_set[compnum][j].toString() + "<br/>";
                break;
            }
        }
        return nm;
    }, nodes_set.length);
}

function handleEdgeLocalArea() {
    special_nodes.clear();
    special_edges.clear();
    area_size = document.getElementById("area_size").value;
    var edges = document.getElementById("vertext").value.replace(/\n/g, " ").split(" ");
    var edges_id = [];
    nodes_id = [];
    for (i = 0; i < edges.length; ++i) {
        edges_id.push(parseInt(edges[i]));
        nodes_id.push(scaffoldgraph.edges[edges_id[i]].from);
        nodes_id.push(scaffoldgraph.edges[edges_id[i]].to);
        special_edges.add(edges_id[edges_id.length - 1]);
    }

    var uniqueNode = nodes_id.filter(function (value, index, self) {
        return self.indexOf(value) === index;
    });

    if (area_size > 0) {
        findLocalArea(uniqueNode, area_size, min_contig_len, isGoodEdge);
    } else {
        for (i = 0; i < uniqueNode.length; ++i) {
            nodes_to_draw.push(uniqueNode[i]);
        }
        for (i = 0; i < edges_id.length; ++i) {
            edges_to_draw.push(edges_id[i]);
        }
    }
    splitOnParts(nodes_to_draw, edges_to_draw);

    createComponentShowList(function (i) {
        DrawGraph(nodes_set[i], edges_set[i]);
        if (edges_id.length > 0) {
            cy.fit(cy.$('#e' + edges_id[0]));
        }
    }, function (i) {
        var nm = "";
        for (var j = 0; j < edges_set[i].length; ++j) {
            if (special_edges.has(edges_set[i][j])) {
                nm += "edge id=" + edges_set[i][j].toString() + "<br/>";
            }
        }
        return nm;
    }, function (i) {
        var nm = "";
        for (var j = 0; j < edges_set[i].length; ++j) {
            if (special_edges.has(edges_set[i][j])) {
                nm += "edge id: " + edges_set[i][j].toString() + "<br/>";
                break;
            }
        }
        return nm;
    }, nodes_set.length);
}

function handleFilterButton() {
    document.getElementById("UpdateGraph").style.visibility = "hidden";
    defZoom = 100;
    special_nodes.clear();
    special_edges.clear();

    var min_edge_weight = [];
    min_contig_len = document.getElementById("min_contig_len").value;

    for (var i = 0; i < scaffoldgraph.libs.length; ++i) {
        min_edge_weight.push(document.getElementById("min_w_" + scaffoldgraph.libs[i].name).value);
    }

    isGoodEdge = function (e) {
        if (min_edge_weight[scaffoldgraph.edges[e].lib] > scaffoldgraph.edges[e].weight) {
            return false;
        }

        if (scaffoldgraph.nodes[scaffoldgraph.edges[e].from].len < min_contig_len) {
            return false;
        }

        if (scaffoldgraph.nodes[scaffoldgraph.edges[e].to].len < min_contig_len) {
            return false;
        }

        return true;
    };

    nodes_to_draw = [];
    edges_to_draw = [];

    if ((document.getElementById("select_layout").value !== "free_layout")) {
        area_size = document.getElementById("area_size").value;
        handleAlongChromosomesFilter();
        return;
    }

    var opt = document.getElementById("select_show_type").value;
    if (opt == "vertices_local_area") {
        handleVertexLocalArea();
    } else if (opt=="edges_local_area") {
        handleEdgeLocalArea();
    } else if (opt=="diff in libs") {
        area_size = document.getElementById("area_size").value;
        handleDiffInLibsFilter(area_size, min_contig_len, isGoodEdge);
    } else if (opt=="full graph") {
        handleFullGraph();
    } else if (opt=="scaffolds") {
        area_size = document.getElementById("area_size").value;
        handleScaffoldsFilter(document.getElementById("select_scaff_lib").value, area_size, min_contig_len, isGoodEdge);
    } else if (opt=="ambiguous") {
        area_size = document.getElementById("area_size").value;
        handleAmbiguousFilter(area_size, min_contig_len, isGoodEdge);
    }
}

function InitLibTable() {
    var table = document.getElementById("lib_table");

    for (var i=0; i < scaffoldgraph.libs.length; ++i) {
        var tr = document.createElement("tr");

        var td_id = document.createElement("td");
        var lib_id = document.createElement("p");
        lib_id.appendChild(document.createTextNode("l" + scaffoldgraph.libs[i].id));
        td_id.align="center";
        td_id.appendChild(lib_id);

        var td_type = document.createElement("td");
        var lib_type = document.createElement("p");
        lib_type.appendChild(document.createTextNode((scaffoldgraph.libs[i].type).replace(/_/g, " ")));
        td_type.appendChild(lib_type);
        td_type.align="center";


        var td_name = document.createElement("td");
        var lib_name = document.createElement("p");
        lib_name.style.color = scaffoldgraph.libs[i].color;
        lib_name.appendChild(document.createTextNode(scaffoldgraph.libs[i].name));
        lib_name.id = "color" + scaffoldgraph.libs[i].name;
        td_name.align="center";
        td_name.appendChild(lib_name);

        var td_min_edge_weight = document.createElement("td");
        var input_weight = document.createElement("input");
        input_weight.type = "number";
        input_weight.min = 0;
        input_weight.size = 1;
        input_weight.value = getDefaultWeight(i);
        input_weight.id = "min_w_" + scaffoldgraph.libs[i].name;
        td_min_edge_weight.align="center";
        td_min_edge_weight.appendChild(input_weight);

        tr.appendChild(td_id);
        tr.appendChild(td_type);
        tr.appendChild(td_name);
        tr.appendChild(td_min_edge_weight);
        table.appendChild(tr);
    }
}

function InitAlignmentsForNodes() {
    for (var i=0; i < chromosomes.length; ++i) {
        for (var j=0; j < chromosomes[i].alignments.length; ++j) {
            scaffoldgraph.nodes[chromosomes[i].alignments[j].node_id].alignments.push(chromosomes[i].alignments[j]);
        }
    }
}

function disableNot(checkBox) {
    if (checkBox.checked) {
        document.getElementById("checkbox_not_present_" + checkBox.id.substring(17)).disabled = true;
        document.getElementById("notPresName" + checkBox.id.substring(17)).style.color = "#aaa";
    } else {
        document.getElementById("checkbox_not_present_" + checkBox.id.substring(17)).disabled = false;
        document.getElementById("notPresName" + checkBox.id.substring(17)).style.color = "#4F4F4F";
    }
}


function disable(checkBox) {
    if (checkBox.checked) {
        document.getElementById("checkbox_present_" + checkBox.id.substring(21)).disabled = true;
        document.getElementById("presName" + checkBox.id.substring(21)).style.color = "#aaa";
    } else {
        document.getElementById("checkbox_present_" + checkBox.id.substring(21)).disabled = false;
        document.getElementById("presName" + checkBox.id.substring(21)).style.color = "#4F4F4F";
    }
}

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

InitLibTable();
InitAlignmentsForNodes();
setupAutocompleteSearch();
putEdgesNumForScaffolds();
handleFilterButton();
document.getElementById("extra_info").innerHTML = "<p style='margin-top: 0px; margin-bottom: 0px;'>" + generateGeneralInfo() + "</p>";

function updeteChangeBlock() {
    if(document.getElementById("select_show_type").value == "full graph") {
        document.getElementById("change_block").innerHTML = "";
    } else if (document.getElementById("select_show_type").value == "vertices_local_area" || document.getElementById("select_show_type").value == "edges_local_area") {
        document.getElementById("change_block").innerHTML = "<div class=\"block\">\n" +
            "                    <p>Area size:<br/>\n" +
            "                        <input type=\"number\" min=\"0\" id=\"area_size\" value=1>\n" +
            "                    </p>\n" +
            "                </div>\n" +
            "                <div class=\"block\">\n" +
            "                    <textarea rows=\"6\" id=\"vertext\"></textarea>\n" +
            "                </div>";
    } else if (document.getElementById("select_show_type").value == "diff in libs" ) {
        var html_code = "<div class=\"block\">\n" +
            "                    <p>Area size:<br/>\n" +
            "                        <input type=\"number\" min=\"0\" id=\"area_size\" value=1>\n" +
            "                    </p>\n" +
            "                </div>\n" +
            "<div class=\"one_line_block\">\n" +
            "                        <label class=\"container\">\n" +
            "                            <p>Wrong</p>\n" +
            "                            <input type=\"checkbox\" checked=\"true\" value=\"Wrong\" id=\"checkbox_wrong\">\n" +
            "                            <span class=\"checkmark\"></span>\n" +
            "                        </label>\n" +
            "                    </div>\n" +
            "\n" +
            "                    <div class=\"one_line_block\">\n" +
            "                        <label class=\"container\">\n" +
            "                            <p>Correct</p>\n" +
            "                            <input type=\"checkbox\" checked=\"true\" value=\"Correct\" id=\"checkbox_correct\">\n" +
            "                            <span class=\"checkmark\"></span>\n" +
            "                        </label>\n" +
            "                    </div>\n" +
            "                    <br/>\n";
        var present_block = "<div class=\"one_line_block\" id=\"present_block\">\n" +
            "                        Presnt:\n";

        var not_present_block = "<div class=\"one_line_block\" id=\"not_present_block\">\n" +
            "                        Not Presnt:\n";

        for (var i=0; i < scaffoldgraph.libs.length; ++i) {
            var lib_name = scaffoldgraph.libs[i].name;
            present_block +=
                "                        <label class=\"container\">\n" +
                "                            <p id=\"presName" + i.toString() +"\">" + lib_name + "</p>\n" +
                "                            <input type=\"checkbox\" value=\"" + lib_name + "\" id=\"checkbox_present_"+i.toString() + "\" onchange=\"disableNot(this)\">\n" +
                "                            <span class=\"checkmark\"></span>\n" +
                "                        </label>\n";

            not_present_block +=
                "                        <label class=\"container\">\n" +
                "                            <p id=\"notPresName" + i.toString() +"\">" + lib_name + "</p>\n" +
                "                            <input type=\"checkbox\" value=\"" + lib_name + "\" id=\"checkbox_not_present_"+i.toString() + "\" onchange=\"disable(this)\" >\n" +
                "                            <span class=\"checkmark\"></span>\n" +
                "                        </label>\n";
        }

        present_block += "</div>\n";
        not_present_block += "</div>\n";
        document.getElementById("change_block").innerHTML = html_code + present_block + not_present_block;
    } else if (document.getElementById("select_show_type").value == "scaffolds") {
        document.getElementById("change_block").innerHTML = "<div class=\"block\">\n" +
            "                    <p>Area size: \n<br>" +
            "                        <input type=\"number\" min=\"0\" id=\"area_size\" value=1>\n" +
            "                    </p>\n" +
            "                    <p>Scaffold source:\n<br>" +
            "                        <div class=\"styled-select\">\n" +
            "                           <select id=\"select_scaff_lib\">\n" +
            "                           </select>\n" +
            "                       </div>" +
            "                    </p>\n" +
            "                    <p> Min scaffold len:\n" +
            "                    <input type=\"number\" min=\"2\" id=\"min_scaffold_len\" value=2> </p>\n" +
            "                    <label class=\"container\">\n" +
            "                    <p>Wrong connection(s)\n" +
            "                    <input type=\"checkbox\" id=\"scaff_wrng\">\n" +
            "                    <span class=\"checkmark\"></span>" +
            "                    </p></label>" +
            "                    <label class=\"container\">\n" +
            "                    <p>Possibly incomplete\n" +
            "                    <input type=\"checkbox\" id=\"scaff_cont\">\n" +
            "                    <span class=\"checkmark\"></span>" +
            "                    </p></label>" +
            "                    <label class=\"container\">\n" +
            "                    <p>Ambiguous connection(s)\n" +
            "                    <input type=\"checkbox\" id=\"scaff_ambig\" style='width: 13px'>\n" +
            "                    <span class=\"checkmark\"></span>" +
            "                    </p></label>" +
            "                    <p style='font-size: 12px;'>When multiple boxes are checked, components that satisfy at least one connection will be shown.</p>" +
            "                </div>";

        for (i=0; i < scaffoldgraph.libs.length; ++i) {
            if (scaffoldgraph.libs[i].type == 'SCAFF') {
                document.getElementById("select_scaff_lib").innerHTML += "<option value=\"" + scaffoldgraph.libs[i].name + "\">" + scaffoldgraph.libs[i].name + "</option>\n";
            }
        }
    } else if (document.getElementById("select_show_type").value == "ambiguous") {
        document.getElementById("change_block").innerHTML = "<div class=\"block\">\n" +
            "                    <p>Area size:<br/>\n" +
            "                        <input type=\"number\" min=\"0\" id=\"area_size\" value=1>\n" +
            "                    </p>\n" +
            "                </div>";
    }
}

document.getElementById("filter_button").addEventListener("click", handleFilterButton);

document.getElementById("select_layout").addEventListener("change", function(){
    if (document.getElementById("select_layout").value == "free_layout") {
        document.getElementById("filtration").innerHTML = "<div class=\"block\">\n" +
            "                <h2 class=\"section_name\">\n" +
            "                Filtration:\n" +
            "                </h2>\n" +
            "                <div class=\"styled-select\">\n" +
            "                <select id=\"select_show_type\">\n" +
            "                    <option value=\"full graph\">Full graph</option>\n" +
            "                    <option value=\"scaffolds\">Scaffolds</option>\n" +
            "                    <option value=\"diff in libs\">Difference in sources</option>\n" +
            "                    <option value=\"ambiguous\">Ambiguous</option>\n" +
            "                    <option value=\"vertices_local_area\">Vertices local area</option>\n" +
            "                    <option value=\"edges_local_area\">Edges local area</option>\n" +
            "                </select>\n" +
            "                </div>\n" +
            "                </div>";

        document.getElementById("change_block").innerHTML = "";
        document.getElementById("select_show_type").addEventListener("change", updeteChangeBlock);
    } else {
        document.getElementById("filtration").innerHTML = "";
        document.getElementById("change_block").innerHTML = "<div class=\"block\">\n" +
            "                    <p>Area size:<br/>\n" +
            "                        <input type=\"number\" min=\"0\" id=\"area_size\" value=1>\n" +
            "                    </p>\n" +
            "                </div>";
    }

});

document.getElementById("select_show_type").addEventListener("change", updeteChangeBlock);
/*
* generate labels for nodes and edges
*/


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
        label += "length: " + scaffoldgraph.nodes[node].len + "\n";
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
    label += "length: " + scaffoldgraph.nodes[node].len + "</br>";
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
            label += "length: " + scaffoldgraph.edges[edge].len + "\n";
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
        label += "length: " + scaffoldgraph.edges[edge].len + "</br>";
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
        if (document.getElementById("extra_info").clientHeight <= def_height) {
            document.getElementById("extra_info").style.height = def_height + 'px';
        }
    });

    cy.on('mouseout', 'edge', function (evt) {
        document.getElementById("extra_info").style.height = def_height + 'px';
        document.getElementById("extra_info").innerHTML = "<p style='margin-top: 0px; margin-bottom: 0px;'>" + generateGeneralInfo() + "</p>";
    });
}
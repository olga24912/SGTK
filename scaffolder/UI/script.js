function createLabelForNode(node) {
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
    return label;
}

function createLabelForEdge(edge) {
    var label = "";
    if (document.getElementById("edge_checkbox_id").checked) {
        label += "id: " + scaffoldgraph.edges[edge].id + "\n";
    }
    if (document.getElementById("edge_checkbox_name").checked) {
        label += scaffoldgraph.libs[scaffoldgraph.edges[edge].lib].name + "\n";
    }
    if (document.getElementById("edge_checkbox_weight").checked) {
        label += "w: " + scaffoldgraph.edges[edge].weight + "\n";
    }

    if (document.getElementById("edge_checkbox_type").checked) {
        label += scaffoldgraph.libs[scaffoldgraph.edges[edge].lib].type + "\n";
    }
    return label;
}

var edges_set = [];
var nodes_set = [];
var cur_show_id = 0;

function DrawGraph(nodes_to_draw, edges_to_draw) {
    var nodeslist = [];
    var edgeslist = [];
    var i = 0;
    var j = 0;
    for (i=0; i < nodes_to_draw.length; i++) {
        j = nodes_to_draw[i];
        var label = createLabelForNode(j);
        nodeslist.push({id: scaffoldgraph.nodes[j].id, label: label});
    }

    // create an array with nodes
    var nodes = new vis.DataSet(nodeslist);

    for (i=0; i < edges_to_draw.length; i++) {
        j = edges_to_draw[i];
        label = createLabelForEdge(j);

        edgeslist.push({from: scaffoldgraph.edges[j].from, to: scaffoldgraph.edges[j].to, label : label, arrows: 'to', color:{color: scaffoldgraph.libs[scaffoldgraph.edges[j].lib].color}});
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
    var options = {
        nodes : {
            shape: 'dot',
            size: 7
        },
        layout:{
            randomSeed:34
        },
        physics: {
            forceAtlas2Based: {
                gravitationalConstant: -26,
                centralGravity: 0.005,
                springLength: 230,
                springConstant: 0.18
            },
            maxVelocity: 146,
            solver: 'forceAtlas2Based',
            timestep: 0.35,
            stabilization: {
                enabled:true,
                iterations:2000,
                updateInterval:25
            }
        }

        //,
       // interaction: {
       //     hideEdgesOnDrag: true,
       //     tooltipDelay: 200
       // },
       // physics: false
    };

    // initialize your network!
    var network = new vis.Network(container, data, options);
}

function findComponent(v, g, color, curc) {
    color.set(v, curc);
    var nb = g.get(v);
    for (var i=0; i < nb.length; ++i) {
        var u = nb[i];
        if (color.get(u) === -1) {
            findComponent(u, g, color, curc);
        }
    }
}

function splitOnParts(nodes_to_draw, edges_to_draw) {
    var g = new Map();
    var color = new Map();
    for (var i=0; i < nodes_to_draw.length; ++i) {
        g.set(nodes_to_draw[i], []);
        color.set(nodes_to_draw[i], -1);
    }

    for (i=0; i < edges_to_draw.length; ++i) {
        g.get(scaffoldgraph.edges[edges_to_draw[i]].from).push(scaffoldgraph.edges[edges_to_draw[i]].to);
        g.get(scaffoldgraph.edges[edges_to_draw[i]].to).push(scaffoldgraph.edges[edges_to_draw[i]].from);
    }

    var curc = 0;
    for (i=0; i < nodes_to_draw.length; ++i) {
        if (color.get(nodes_to_draw[i]) === -1) {
            findComponent(nodes_to_draw[i], g, color, curc);
            ++curc;
        }
    }

    nodes_set = [];
    edges_set = [];
    for (i=0; i < curc; ++i) {
        nodes_set.push([]);
        edges_set.push([]);
    }

    for(i=0; i < nodes_to_draw.length; ++i) {
        nodes_set[color.get(nodes_to_draw[i])].push(nodes_to_draw[i]);
    }

    for(i=0; i < edges_to_draw.length; ++i) {
        edges_set[color.get(scaffoldgraph.edges[edges_to_draw[i]].from)].push(edges_to_draw[i]);
    }
}

function drawLocalArea(inodes, area_size, min_edge_weight) {
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
        if (curd <= area_size) {
            nodes_to_draw.push(curv);
            for (i = 0; i < scaffoldgraph.g[curv].length; ++i) {
                var curedge = scaffoldgraph.g[curv][i];
                if (min_edge_weight[curedge.lib] <= curedge.weight) {
                    var curu = curedge.to;

                    if (!dist.has(curu)) {
                        queue.push(curu);
                        dist.set(curu, curd + 1);
                    }

                    if (dist.get(curu) <= area_size) {
                        edges_to_draw.push(curedge.id);
                    }
                }
            }


            for (i = 0; i < scaffoldgraph.gr[curv].length; ++i) {
                curedge = scaffoldgraph.gr[curv][i];
                if (min_edge_weight[curedge.lib] <= curedge.weight) {
                    curu = curedge.from;

                    if (!dist.has(curu)) {
                        queue.push(curu);
                        dist.set(curu, curd + 1);
                    }
                }
            }
        }
    }

    DrawGraph(nodes_to_draw, edges_to_draw);
}

function handleFilterButton() {
    var opt = document.getElementById("select_show_type").value;
    var min_edge_weight = [];
    var min_contig_len = document.getElementById("min_contig_len").value;

    for (var i=0; i < scaffoldgraph.libs.length; ++i) {
        min_edge_weight.push(document.getElementById("min_w_" + scaffoldgraph.libs[i].name).value);
    }

    var isGoodEdge = function(e) {
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


    if (opt=="vertices_local_area") {
        var areasize = document.getElementById("area_size").value;
        var nodes = document.getElementById("vertext").value.replace(/\n/g, " ").split(" ");
        var nodes_id = [];
        for (i=0; i < nodes.length; ++i) {
            nodes_id.push(scaffoldgraph.id_by_name.get(nodes[i]));
        }
        drawLocalArea(nodes_id, areasize, min_edge_weight);
    } else if (opt=="full graph") {
        var nodes_to_draw = [];
        var edges_to_draw = [];
        for (i=0; i < scaffoldgraph.nodes.length; ++i) {
            if (scaffoldgraph.nodes[i].len >= min_contig_len) {
                nodes_to_draw.push(scaffoldgraph.nodes[i].id);
            }
        }

        for (i=0; i < scaffoldgraph.edges.length; ++i) {
            if (isGoodEdge(i)) {
                edges_to_draw.push(scaffoldgraph.edges[i].id);
            }
        }
        splitOnParts(nodes_to_draw, edges_to_draw);
        createComponentShowList(function(i) {
            DrawGraph(nodes_set[i], edges_set[i]);
        }, function(i) {
            return "comp " + i + "<br/> #nodes = " + nodes_set[i].length + "<br/>#edges = " + edges_set[i].length;
        }, function(compnum) {
            return "Component #" + compnum;
        }, nodes_set.length);
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
        input_weight.id = "min_w_" + scaffoldgraph.libs[i].name;
        td_min_edge_weight.align="center";
        td_min_edge_weight.appendChild(input_weight);

        tr.appendChild(td_type);
        tr.appendChild(td_id);
        tr.appendChild(td_name);
        tr.appendChild(td_min_edge_weight);
        table.appendChild(tr);
    }
}

InitLibTable();
document.getElementById("filter_button").addEventListener("click", handleFilterButton);
document.getElementById("select_show_type").addEventListener("change", function() {
    if(document.getElementById("select_show_type").value == "full graph") {
        document.getElementById("change_block").innerHTML = "";
    } else if (document.getElementById("select_show_type").value == "vertices_local_area" || document.getElementById("select_show_type").value == "edges_local_area") {
        document.getElementById("change_block").innerHTML = "<div class=\"block\">\n" +
            "                    <p>Area size:<br/>\n" +
            "                        <input type=\"number\" min=\"0\" id=\"area_size\">\n" +
            "                    </p>\n" +
            "                </div>\n" +
            "                <div class=\"block\">\n" +
            "                    <textarea rows=\"6\" id=\"vertext\"></textarea>\n" +
            "                </div>";
    } else if (document.getElementById("select_show_type").value == "diff in libs" ) {
        var html_code = "<div class=\"block\">\n" +
            "                    <p>Area size:<br/>\n" +
            "                        <input type=\"number\" min=\"0\" id=\"area_size\">\n" +
            "                    </p>\n" +
            "                </div>\n" +
            "<div class=\"one_line_block\">\n" +
            "                        <label class=\"container\">\n" +
            "                            <p>Wrong</p>\n" +
            "                            <input type=\"checkbox\" value=\"Wrong\" id=\"checkbox_wrong\">\n" +
            "                            <span class=\"checkmark\"></span>\n" +
            "                        </label>\n" +
            "                    </div>\n" +
            "\n" +
            "                    <div class=\"one_line_block\">\n" +
            "                        <label class=\"container\">\n" +
            "                            <p>Correct</p>\n" +
            "                            <input type=\"checkbox\" value=\"Correct\" id=\"checkbox_correct\">\n" +
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
                "                            <p>" + lib_name + "</p>\n" +
                "                            <input type=\"checkbox\" value=\"" + lib_name + "\" id=\"checkbox_present_"+lib_name + "\">\n" +
                "                            <span class=\"checkmark\"></span>\n" +
                "                        </label>\n";

            not_present_block +=
                "                        <label class=\"container\">\n" +
                "                            <p>" + lib_name + "</p>\n" +
                "                            <input type=\"checkbox\" value=\"" + lib_name + "\" id=\"checkbox_not_present_"+lib_name + "\">\n" +
                "                            <span class=\"checkmark\"></span>\n" +
                "                        </label>\n";
        }

        present_block += "</div>\n";
        not_present_block += "</div>\n";
        document.getElementById("change_block").innerHTML = html_code + present_block + not_present_block;
    }
});
/*
* Main executable script.
* Handle filtration updating.
*/

//int[][], for edges_set[i] list of edgeID in component #i
var edges_set = [];

//int[][], for nodes_set[i] list of nodesID in component #i
var nodes_set = [];

//List of nodes ID(int) to draw
var nodes_to_draw = [];

//List of edges ID(int) to draw
var edges_to_draw = [];

//Set of specific nodes which will be highlighted
var special_nodes = new Set();

//Set of specific edges which will be highlighted
var special_edges = new Set();

//function which get edge ID(int) return if this edge not filtered(bool)
var isGoodEdge;

//minimum contig length to draw
var min_contig_len = 0;

//Maximum distance to draw
var area_size = 0;

//Id of selected component
var cur_show_id = 0;

//Scaffold graph
var graph = null;

//Cytoscape graph
var cy = null;

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

//Updating disabled in Absent colon in "difference of sources" mode
function disableNot(checkBox) {
    if (checkBox.checked) {
        document.getElementById("checkbox_not_present_" + checkBox.id.substring(17)).disabled = true;
        document.getElementById("notPresName" + checkBox.id.substring(17)).style.color = "#aaa";
    } else {
        document.getElementById("checkbox_not_present_" + checkBox.id.substring(17)).disabled = false;
        document.getElementById("notPresName" + checkBox.id.substring(17)).style.color = "#4F4F4F";
    }
}


//Updating disabled in Present colon in "difference of sources" mode
function disable(checkBox) {
    if (checkBox.checked) {
        document.getElementById("checkbox_present_" + checkBox.id.substring(21)).disabled = true;
        document.getElementById("presName" + checkBox.id.substring(21)).style.color = "#aaa";
    } else {
        document.getElementById("checkbox_present_" + checkBox.id.substring(21)).disabled = false;
        document.getElementById("presName" + checkBox.id.substring(21)).style.color = "#4F4F4F";
    }
}

//Init page
InitLibTable();
InitAlignmentsForNodes();
setupAutocompleteSearch();
putEdgesNumForScaffolds();
handleFilterButton();
document.getElementById("extra_info").innerHTML = "<p style='margin-top: 0px; margin-bottom: 0px;'>" + generateGeneralInfo() + "</p>";

/*
* Updating filter panele
*/
function updateChangeBlock() {
    if(document.getElementById("select_show_type").value == "full graph") {
        document.getElementById("change_block").innerHTML = "";
    } else if (document.getElementById("select_show_type").value == "vertices_local_area" || document.getElementById("select_show_type").value == "edges_local_area") {
        document.getElementById("change_block").innerHTML = "<div class=\"block\">\n" +
            "                    <p>Distance:<br/>\n" +
            "                        <input type=\"number\" min=\"0\" id=\"area_size\" value=1>\n" +
            "                    </p>\n" +
            "                </div>\n" +
            "                <div class=\"block\">\n" +
            "                    <textarea rows=\"6\" id=\"vertext\" placeholder=\"ID1, ID2...\"></textarea>\n" +
            "                </div>";
    } else if (document.getElementById("select_show_type").value == "diff in libs" ) {
        var html_code_1 = "<div class=\"block\">\n" +
            "                    <p>Distance:<br/>\n" +
            "                        <input type=\"number\" min=\"0\" id=\"area_size\" value=1>\n" +
            "                    </p>\n" +
            "                </div>\n";

        var text_wrong_correct = "<p>Display only</p>";
        var wrong_correct = "<div class=\"one_line_block\">\n" +
            "                        <label class=\"container\">\n" +
            "                            <p>Wrong connection</p>\n" +
            "                            <input type=\"checkbox\" checked=\"true\" value=\"Wrong\" id=\"checkbox_wrong\">\n" +
            "                            <span class=\"checkmark\"></span>\n" +
            "                        </label>\n" +
            "                    </div>\n" +
            "\n" +
            "                    <div class=\"one_line_block\">\n" +
            "                        <label class=\"container\">\n" +
            "                            <p>Correct connection</p>\n" +
            "                            <input type=\"checkbox\" checked=\"true\" value=\"Correct\" id=\"checkbox_correct\">\n" +
            "                            <span class=\"checkmark\"></span>\n" +
            "                        </label>\n" +
            "                    </div>\n" +
            "                    <br/>\n";

        var text_after = "<p style='font-size: 12px;'>Wrong and correct connections are detected using reference genome</p>";
        var text_before_present = "<p style='padding-bottom: 0px; margin-bottom: 0px;'>Display connection where following sources are</p>";

        var present_block = "<div class=\"one_line_block\" id=\"present_block\">\n" +
            "                        <p style='padding-top: 0px; margin-top: 2px;'>present:</p>\n";

        var not_present_block = "<div class=\"one_line_block\" id=\"not_present_block\">\n" +
            "                        <p style='padding-top: 0px; margin-top: 2px;'>absent:</p>\n";

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
        document.getElementById("change_block").innerHTML = html_code_1 + text_before_present + present_block +
            not_present_block + text_wrong_correct + wrong_correct + text_after;
    } else if (document.getElementById("select_show_type").value == "scaffolds") {
        document.getElementById("change_block").innerHTML = "<div class=\"block\">\n" +
            "                    <p>Distance: \n<br>" +
            "                        <input type=\"number\" min=\"0\" id=\"area_size\" value=1>\n" +
            "                    </p>\n" +
            "                    <p>Scaffold source:\n<br>" +
            "                        <div class=\"styled-select\">\n" +
            "                           <select id=\"select_scaff_lib\">\n" +
            "                           </select>\n" +
            "                       </div>" +
            "                    </p>\n" +
            "                    <p> Minimum scaffold size:\n" +
            "                    <input type=\"number\" min=\"2\" id=\"min_scaffold_len\" value=2> </p>\n" +
            "                    <label class=\"container\">\n" +
            "                    <p id='scaff_wrng_p'>Wrong connection(s)\n" +
            "                    <input type=\"checkbox\" id=\"scaff_wrng\">\n" +
            "                    <span class=\"checkmark\"></span>" +
            "                    </p></label>" +
            "                    <label class=\"container\">\n" +
            "                    <p id='scaff_cont_p'>Possibly incomplete\n" +
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

        if (chromosomes.length === 0) {
            document.getElementById("scaff_wrng").disabled = true;
            document.getElementById("scaff_cont").disabled = true;

            document.getElementById("scaff_wrng_p").style.color = "#aaa";
            document.getElementById("scaff_cont_p").style.color = "#aaa";
        }

        for (i=0; i < scaffoldgraph.libs.length; ++i) {
            if (scaffoldgraph.libs[i].type == 'SCAFF') {
                document.getElementById("select_scaff_lib").innerHTML += "<option value=\"" + scaffoldgraph.libs[i].name + "\">" + scaffoldgraph.libs[i].name + "</option>\n";
            }
        }
    } else if (document.getElementById("select_show_type").value == "ambiguous") {
        document.getElementById("change_block").innerHTML = "<div class=\"block\">\n" +
            "                    <p>Distance:<br/>\n" +
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
            "                    <option value=\"full graph\">full graph</option>\n" +
            "                    <option value=\"scaffolds\">by scaffold name</option>\n" +
            "                    <option value=\"vertices_local_area\">by vertex id</option>\n" +
            "                    <option value=\"edges_local_area\">by edge id</option>\n" +
            "                    <option value=\"diff in libs\">difference between sources</option>\n" +
            "                    <option value=\"ambiguous\">ambiguous connection</option>\n" +
            "                </select>\n" +
            "                </div>\n" +
            "                </div>";

        document.getElementById("change_block").innerHTML = "";
        document.getElementById("select_show_type").addEventListener("change", updateChangeBlock);
    } else {
        document.getElementById("filtration").innerHTML = "";
        document.getElementById("change_block").innerHTML = "<div class=\"block\">\n" +
            "                    <p>Distance:<br/>\n" +
            "                        <input type=\"number\" min=\"0\" id=\"area_size\" value=1>\n" +
            "                    </p>\n" +
            "                </div>";
    }

});

document.getElementById("select_show_type").addEventListener("change", updateChangeBlock);
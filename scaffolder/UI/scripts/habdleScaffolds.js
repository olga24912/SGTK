var libnum = 0;
var global_areasize = 0;
var global_min_contig_len = 0;
var global_isGoodEdge;
var nodes_to_draw = [];
var edges_to_draw = [];

function drawScaffold(i) {
    var inodes = [];

    inodes.push(scaffoldgraph.libs[libnum].scaffolds[i].edges[0].from);
    for (var g = 0; g < scaffoldgraph.libs[libnum].scaffolds[i].edges.length; ++g) {
        var e = scaffoldgraph.libs[libnum].scaffolds[i].edges[g];
        inodes.push(e.to);
    }

    findLocalArea(inodes, global_areasize, global_min_contig_len, global_isGoodEdge);


    DrawGraph(nodes_to_draw, edges_to_draw);
}

function handleScaffoldsFilter(scafflibname, areasize, min_contig_len, isGoodEdge) {
    for (var j = 0; j < scaffoldgraph.libs.length; ++j) {
        if (scaffoldgraph.libs[j].name == scafflibname) {
            libnum = j;
        }
    }

    global_areasize = areasize;
    global_min_contig_len = min_contig_len;
    global_isGoodEdge = isGoodEdge;

    notzero = [];
    for (j = 0; j < scaffoldgraph.libs[libnum].scaffolds.length; ++j) {
        if (scaffoldgraph.libs[libnum].scaffolds[j].edges.length > 0) {
            notzero.push(j);
        }
    }

    createComponentShowList(function(i) {
            drawScaffold(notzero[i]);
        }, function(i) {
            return scaffoldgraph.libs[libnum].scaffolds[notzero[i]].name + "<br/> scaffold len: " + (scaffoldgraph.libs[libnum].scaffolds[notzero[i]].edges.length + 1);
        }, function(i) {
            return "Scaffold " + scaffoldgraph.libs[libnum].scaffolds[notzero[i]].name;
        }, notzero.length
    );
}
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

    var dnodes = [];
    var dedges = [];
    for (g = 0; g < nodes_to_draw.length; ++g) {
        dnodes.push({ data: {id: nodes_to_draw[g], label: createLabelForNode(nodes_to_draw[g])}});
    }

    for (g = 0; g < edges_to_draw.length; ++g) {
        dedges.push({ data: { source: scaffoldgraph.edges[edges_to_draw[g]].from,
                target: scaffoldgraph.edges[edges_to_draw[g]].to, label: createLabelForEdge(edges_to_draw[g]),
                faveColor: scaffoldgraph.edges[edges_to_draw[g]].color}});
    }

    var cy = cytoscape({
        container: document.getElementById('mainpanel'),

        boxSelectionEnabled: false,
        autounselectify: true,
        maxZoom: 2,
        minZoom: 0.25,

        elements: {
            nodes: dnodes,
            edges: dedges
        },

        layout: {
            name: 'breadthfirst',
            directed: true,
            padding: 10
        },

        ready: function(){
            window.cy = this;
        },

        style:  cytoscape.stylesheet()
            .selector('node')
            .css({
                'content': 'data(name)',
                'text-valign': 'center',
                'color': 'white',
                'text-outline-width': 2,
                'background-color': '#999',
                'text-outline-color': '#999'
            })
            .selector('edge')
            .css({
                'curve-style': 'bezier',
                'target-arrow-shape': 'triangle',
                'width': 1,
                'source-arrow-color': 'data(faveColor)',
                'target-arrow-color': 'data(faveColor)'
            })
            .selector(':selected')
            .css({
                'background-color': 'black',
                'line-color': 'black',
                'target-arrow-color': 'black',
                'source-arrow-color': 'black'
            })
            .selector('.faded')
            .css({
                'opacity': 0.25,
                'text-opacity': 0
            })
    });
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

    createComponentShowList(drawScaffold, function(i) {
        return scaffoldgraph.libs[libnum].scaffolds[i].name + "<br/> scaffold len: " + (scaffoldgraph.libs[libnum].scaffolds[i].edges.length + 1);
    }, function(i) {
        return "Scaffold " + scaffoldgraph.libs[libnum].scaffolds[i].name;
    }, scaffoldgraph.libs[libnum].scaffolds.length);
}
function getALLNamesFrom(intTree, result) {
    if (intTree.lt !== -1) {
        getALLNamesFrom(intTree.lt, result);
    }
    if (intTree.rt !== -1) {
        getALLNamesFrom(intTree.rt, result);
    }

    for (var i = 0; i < intTree.lstL.length; ++i) {
        result.push("V " + "id: " + intTree.lstL[i].id.toString() + " " + scaffoldgraph.nodes[intTree.lstL[i].id].name);
    }
}

function getNodesId(intTree, nodes) {
    if (intTree.lt !== -1) {
        getNodesId(intTree.lt, nodes);
    }
    if (intTree.rt !== -1) {
        getNodesId(intTree.rt, nodes);
    }

    for (var i = 0; i < intTree.lstL.length; ++i) {
        nodes.push({id: intTree.lstL[i].id, cb: intTree.lstL[i].cb, ce: intTree.lstL[i].ce});
    }
}


function setupAutocompleteSearch(){
    $("#search_text").autocomplete({
        minLength: 1,
        maxHeight: 200,
        deferRequestBy: 50,
        source:  function(request, response) {
            var opt = document.getElementById("select_layout").value;
            var result = [];
            if (opt === "genome_browser") {
                getALLNamesFrom(IntervalTree[1], result);
            } else {
                var nodes = cy.filter('node');
                for (var i = 0; i < nodes.length; ++i) {
                    result.push("V " + "id: " + nodes[i].id().toString() + " " + scaffoldgraph.nodes[nodes[i].id()].name);
                }

                var edges = cy.filter('edge');
                for (i = 0; i < edges.length; ++i) {
                    var eid = edges[i].id().substring(1);
                    var ieid = parseInt(eid);
                    result.push("E " + "id: " + eid + " " + scaffoldgraph.edges[ieid].from + "->" + scaffoldgraph.edges[ieid].to);
                }
            }

            var matcher = new RegExp(request.term.toString());
            var results = $.grep(result, function (value) {
                return matcher.test(value);
            });

            results.length = Math.min(results.length, 20);

            response(results);
        }
    });
}

function parseNode(nodes, scr, vid, request) {
    for (var i = 0; i < nodes.length; ++i) {
        var cur_scr = 0;
        var matcher = new RegExp("V");
        if (matcher.test(request)) {
            cur_scr += 1;
        }
        matcher = new RegExp(nodes[i].id);
        if (matcher.test(request)) {
            cur_scr += 2;
        }
        console.log(request);
        console.log(nodes[i].id);
        matcher = new RegExp(scaffoldgraph.nodes[nodes[i].id].name);
        if (matcher.test(request)) {
            cur_scr += 3;
        }
        console.log(cur_scr);
        if (cur_scr > scr) {
            scr = cur_scr;
            vid = nodes[i];
        }
        console.log(vid.id);
    }
    return [scr, vid];
}

function search() {
    var request = document.getElementById("search_text").value;
    var opt = document.getElementById("select_layout").value;
    var vid = 0;
    var scr = 0;
    var nodes = [];

    var highlightFoundNode = function() {
        if (cntChanges % 2 === 0) {
            cy.getElementById(idd).addClass('found');
        } else {
            cy.getElementById(idd).removeClass('found');
        }
        ++cntChanges;
        if (cntChanges < 6) {
            setTimeout(highlightFoundNode, 300);
        }
    };

    if (opt === "genome_browser") {
        getNodesId(IntervalTree[1], nodes);
        var res = parseNode(nodes, scr, vid, request);
        vid = res[1];

        cy.zoom((1000/(vid.ce - vid.cb))*defZoom);

        var width =  document.getElementById('mainpanel').clientWidth;
        var height =  document.getElementById('mainpanel').clientHeight;

        cy.pan({y: height/3, x: width/2 -cy.zoom()*(vid.cb + vid.ce)/(2*defZoom)});

        var cntChanges = 0;
        var idd = vid.id;

        highlightFoundNode();
    } else {
        nodes = cy.filter('node');
        for (var i = 0; i < nodes.length; ++i) {
            nodes[i] = {id: nodes[i].id()};
        }
        res = parseNode(nodes, scr, vid, request);
        vid = res[1].id;
        scr = res[0];

        var eid = 0;
        var scre = 0;

        var edges = cy.filter('edge');
        for (i = 0; i < edges.length; ++i) {
            cur_scr = 0;
            matcher = new RegExp("E");
            if (matcher.test(request)) {
                cur_scr += 1;
            }
            matcher = new RegExp(edges[i].id().substring(1));
            if (matcher.test(request)) {
                cur_scr += 2;
            }

            var ieid = parseInt(edges[i].id().substring(1));

            matcher = new RegExp(scaffoldgraph.edges[ieid].from + "->" + scaffoldgraph.edges[ieid].to);
            if (matcher.test(request)) {
                cur_scr += 3;
            }

            if (cur_scr > scr) {
                scre = cur_scr;
                eid = edges[i].id();
            }
        }


        if (scr > 0 && scr >= scre) {
            cy.fit(cy.$('#' + vid));

            cntChanges = 0;
            idd = vid;

            highlightFoundNode();
        }


        if (scre > 0 && scre > scr) {
            cy.fit(cy.$('#' + eid));

            cntChanges = 0;
            idd = eid;

            highlightFoundNode();
        }
    }
}

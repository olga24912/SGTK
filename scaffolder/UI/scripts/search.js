

function setupAutocompleteSearch(){
    $( "#search_text" ).autocomplete({
        minLength: 1,
        maxHeight: 200,
        deferRequestBy: 50,
        source:  function(request, response) {
            var nodes = cy.filter('node');
            var result = [];
            for (var i = 0; i < nodes.length; ++i) {
                result.push("V " + "id: " + nodes[i].id().toString() + " " + scaffoldgraph.nodes[nodes[i].id()].name);
            }

            var edges = cy.filter('edge');
            for (i = 0; i < edges.length; ++i) {
                var eid = edges[i].id().substring(1);
                var ieid = parseInt(eid);
                result.push("E " + "id: " + eid + " " + scaffoldgraph.edges[ieid].from + "->" + scaffoldgraph.edges[ieid].to);
            }

            var matcher = new RegExp(request.term.toString());
            var results = $.grep(result, function (value) {
                return matcher.test(value);
            });

            response(results);
        }
    });
}

function search() {
    var request = document.getElementById("search_text").value;
    var vid = 0;
    var scr = 0;

    var nodes = cy.filter('node');
    for (var i = 0; i < nodes.length; ++i) {
        var cur_scr = 0;
        var matcher = new RegExp("V");
        if (matcher.test(request)) {
            cur_scr += 1;
        }
        matcher = new RegExp(nodes[i].id());
        if (matcher.test(request)) {
            cur_scr += 2;
        }
        matcher = new RegExp(scaffoldgraph.nodes[nodes[i].id()].name);
        if (matcher.test(request)) {
            cur_scr += 3;
        }
        if (cur_scr > scr) {
            scr = cur_scr;
            vid = nodes[i].id();
        }
    }

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
        cy.fit(cy.$('#' + vid))
    }


    if (scre > 0 && scre > scr) {
        cy.fit(cy.$('#' + eid))
    }
}

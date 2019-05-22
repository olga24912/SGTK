function isHighlightEdge(eid) {
    if (!isGoodEdge(eid)) {
        return false;
    }

    var edge = scaffoldgraph.edges[eid];
    var lib = scaffoldgraph.libs[edge.lib];
    var checkmark_inbox = document.getElementById("include_lib_" + lib.id);
    return checkmark_inbox.checked;
}

function getShortestPath(sv, fv) {
    var dist = {};
    dist[sv] = 0;
    var edges = {};
    var queue = [];
    queue.push(sv);
    var bg = 0;
    while (bg < queue.length) {
        var cv = queue[bg];
        if (dist.size > 200) {
            break;
        }
        bg += 1;

        for (var i = 0; i < scaffoldgraph.g[cv].length; ++i) {
            var cedge = scaffoldgraph.g[cv][i];
            var tov = cedge.to;
            if (isHighlightEdge(cedge.id)) {
                if (!(tov in dist) || (dist[tov] > dist[cv] + 1)) {
                    dist[tov] = dist[cv] + 1;
                    edges[tov] = scaffoldgraph.g[cv][i];
                    queue.push(tov);
                }
            }
        }
    }

    if (!(fv in dist)) {
        return [];
    } else {
        var connectionList = [];
        cv = fv;
        while (!(cv == sv)) {
            var pv = edges[cv].from;
            connectionList.push({
                "type": "edge",
                "from": pv.toString(),
                "to": cv.toString()
            });
            cv = pv;
        }

        connectionList = connectionList.reverse();
        return connectionList
    }
}

function addNoneConnection(connectionList, from, to) {
    connectionList.push({
        "type": "none",
        "from": from,
        "to": to
    });
}

function addEdgeConnection(connectionList, from, to) {
    connectionList.push({
        "type": "edge",
        "from": from,
        "to": to
    });
}

function addPathConnection(connectionList, from, to) {
    var clist = getShortestPath(from, to);
    for (var j = 0; j < clist.length; ++j) {
        connectionList.push(clist[j]);
    }
}

function getTokens(text_elements) {
    text_elements = text_elements.replace(/;/g, ' $1 ').trim();
    text_elements = text_elements.replace(/([^-])(-)(>)/g, '$1 $2$3 ').trim();
    text_elements = text_elements.replace(/(-)(-)(>)/g, ' $1$2$3 ').trim();
    return text_elements.split(/[\s\n\t,]+/);
}

function parseTextWithElements(text_elements) {
    var tokens = getTokens(text_elements);
    var connectionList = [];
    var vertexList = [];
    var res = [];
    for (var i = 0; i < tokens.length; ++i) {
        if (scaffoldgraph.id_by_name.has(tokens[i])) {
            res.push(scaffoldgraph.id_by_name.get(tokens[i]).toString())
        } else if (scaffoldgraph.scaffold_by_name.has(tokens[i])) {
            scaffols = scaffoldgraph.scaffold_by_name.get(tokens[i]).edges;
            res.push(scaffols[0].from.toString());
            for (var j = 0; j < scaffols.length; ++j) {
                res.push(scaffols[j].to.toString())
            }
        } else {
                res.push(tokens[i])
        }
    }

    for (i = 0; i < res.length; ++i) {
        if (res[i].length > 1 && res[i][0] === '-' && res[i][res[i].length - 1] === '>') {
            continue
        }

        if (res[i] !== "->" && res[i] !== ';' && res[i] !== "-->") {
            vertexList.push(res[i])
        }
    }

    for (i = 1; i < res.length; ++i) {
        if (res[i] === ";") {
            addNoneConnection(connectionList, res[i - 1], res[i + 1]);
            i += 1;
        } else if (res[i] === "->") {
            addEdgeConnection(connectionList, res[i - 1], res[i + 1]);
            i += 1;
        } else if (res[i] === "-->") {
            addPathConnection(connectionList, res[i - 1], res[i + 1]);
            i += 1;
        } else if (res[i].length > 1 && res[i][0] === '-' && res[i][res[i].length - 1] === '>') {
            connectionList.push({
                "type": "edge",
                "from": res[i - 1],
                "to": res[i + 1],
                "edge_id": res[i].substring(1, res[i].length - 1)
            });
            i += 1;
        } else {
            var edge_type = document.getElementById("default_connection");
            edge_type = edge_type.options[edge_type.selectedIndex].value;
            if (edge_type === "one_edge") {
                addEdgeConnection(connectionList, res[i - 1], res[i]);
            } else if (edge_type === "shortest_path") {
                addPathConnection(connectionList, res[i - 1], res[i]);
            } else {
                addNoneConnection(connectionList, res[i - 1], res[i]);
            }
        }
    }
    return [vertexList, connectionList];
}

function updateHighlight() {
    var text_elements = document.getElementById("highlight_elements").value;
    var res = parseTextWithElements(text_elements);
    var elements_id = res[0], connectionList = res[1];
    cy.elements().forEach(function (element, index) {
        element.removeClass("highlight");
    });

    for (var i = 0; i < connectionList.length; ++i) {
        if (connectionList[i]["type"] === "none") {
            continue;
        }
        cy.edges().filter(function (ele) {
            if ("edge_id" in connectionList[i]) {
                return ele.data("id") === ("e" + connectionList[i]["edge_id"]);
            }
            if (isHighlightEdge(ele.data("id").substring(1))) {
                return ele.data('source') === connectionList[i]["from"] &&
                    ele.data('target') === connectionList[i]["to"]
            } else {
                return false;
            }
        }).forEach(function(edge, index){
           edge.addClass("highlight");
        });
    }

    cy.nodes().filter(function (ele) {
        return elements_id.indexOf(ele.data('id')) >= 0;
    }).forEach(function (node, index) {
        node.addClass("highlight");
    });
}

function highlightAutocompleteSetUp() {
    function getVertexNameList() {
        var result = [];
        for (var i = 0; i < scaffoldgraph.nodes.length; ++i) {
            result.push(scaffoldgraph.nodes[i].name);
        }
        for (i = 0; i < scaffoldgraph.nodes.length; ++i) {
            result.push(scaffoldgraph.nodes[i].id.toString());
        }

        return result;
    }

    function getPathNamesList() {
        var result = [];
        let keys = Array.from(scaffoldgraph.scaffold_by_name.keys());
        for (var i = 0; i < keys.length; ++i) {
            result.push(keys[i]);
        }
        return result;
    }

    function valueInCurrentCy(value) {
        var inCy = 0;
        cy.nodes().forEach(function (node, index) {
            if (node.data('id') == value) {
                inCy = 1;
            }

            /*if (value in scaffoldgraph.id_by_name) {
                if (ele.data('id') == scaffoldgraph.id_by_name[value]) {
                    inCy = 1;
                }
            }*/
        });
        return inCy;
    }

    $("#highlight_elements").autocomplete({
        minLength : 0,
        maxLength : 200,
        autoFocus : true,
        source : function(request, response) {
            var tokens = getTokens(request.term.toString());
            var lastToken = tokens[tokens.length - 1];
            var res = getPathNamesList().concat(getVertexNameList());

            var matcher = new RegExp(lastToken);
            res = $.grep(res, function (value) {
                return matcher.test(value);
            });

            res = res.sort(function(a, b) {
                var via = valueInCurrentCy(a);
                var vib = valueInCurrentCy(b);
                return vib - via;
            });

            res.length = Math.min(res.length, 20);
            response(res);
        },
        select : function(event, ui) {
            var tokens = getTokens(this.value);
            tokens.pop();
            tokens.push(ui.item.value);
            this.value = tokens.join(" ");
            return false;
        }
    });
}
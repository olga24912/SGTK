function isHighlightEdge(eid) {
    if (eid.toString().includes("_")) {
        return false;
    }

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
        return [{
            "type": "edge",
            "from": sv.toString(),
            "to": fv.toString()
        }];
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

function allVertexInCy(connectionList) {
    for (var i = 0; i < connectionList.length; ++i) {
        var elem = cy.getElementById(connectionList[i]["from"]);
        if (elem.length === 0) {
            return false;
        }
        var elem = cy.getElementById(connectionList[i]["to"]);
        if (elem.length === 0) {
            return false;
        }
    }

    return true;
}

function add_vertex(vid, prev_vid) {
    var xc = 0;
    var yc = 0;

    if (prev_vid !== undefined) {
        xc = cy.getElementById(prev_vid).position('x') + 100;
        yc = cy.getElementById(prev_vid).position('y');
    }

    var spe = 0;
    var opt = document.getElementById("select_show_type").value;
    if (opt === "full graph") {
        spe = 1;
    }

    var nnode = {
        group: "nodes",
        data: {
            id: vid,
            label: createLabelForNode(vid),
            len: 2*Math.log2(scaffoldgraph.nodes[vid].len)/Math.log2(1.5),
            notALL: 0,
            color: genColorNode(vid),

            color1: "#2A4986",
            color2: "#2A4986",
            color3: "#2A4986",
            cnt1: 100,
            cnt2: 0,
            cnt3: 0,
            special: spe
        },
        position: {
            y: yc,
            x: xc
        }
    };
    updateColorsNode(vid, nnode);
    cy.add(nnode);

    function inCy(vid) {
        var elem = cy.getElementById(vid);
        if (elem.length === 0) {
            return false;
        }
        return true
    }

    function addEdgeWithId(eid) {
        if (!inCy(getEdgeTo(eid)) || !inCy(getEdgeFrom(eid))) {
            return;
        }
        if (inCy("e" + eid.toString())) {
            return;
        }
        if (!isGoodEdge(eid)) {
            return;
        }

        spe = 0;
        opt = document.getElementById("select_show_type").value;
        if (opt === "full graph") {
            spe = 1;
        }

        cy.add({
            group: "edges",
            data: {
                id: "e" + eid.toString(),
                source: getEdgeFrom(eid),
                target: getEdgeTo(eid),
                label: createLabelForEdge(eid),
                faveColor: scaffoldgraph.libs[scaffoldgraph.edges[eid].lib].color,
                weight: generateEdgeWeight(eid),
                lstyle: 'dotted',
                special: spe
            }
        });
    }

    for (g = 0; g < scaffoldgraph.g[vid].length; ++g) {
        var eid = scaffoldgraph.g[vid][g].id;
        addEdgeWithId(eid);
    }

    for (g = 0; g < scaffoldgraph.gr[vid].length; ++g) {
        var eid = scaffoldgraph.gr[vid][g].id;
        addEdgeWithId(eid);
    }
}

function add_miss_vertex(connectionList) {
    for (var i = 0; i < connectionList.length; ++i) {
        var elem = cy.getElementById(connectionList[i]["from"]);
        if (elem.length === 0) {
            add_vertex(connectionList[i]["from"]);
        }
        
        var elem = cy.getElementById(connectionList[i]["to"]);
        if (elem.length === 0) {
            add_vertex(connectionList[i]["to"], connectionList[i]["from"]);
        }
    }
}

function check_miss_vertex(connectionList) {
    if (allVertexInCy(connectionList) === false) {
        document.getElementById("AlertBoxNotAllVertexPresent").style.visibility = "visible";
    } else {
        document.getElementById("AlertBoxNotAllVertexPresent").style.visibility = "hidden";
    }
}

function updateHighlight(with_update=false) {
    var text_elements = document.getElementById("highlight_elements").value;
    var res = parseTextWithElements(text_elements);
    var elements_id = res[0], connectionList = res[1];

    if (with_update === false){
        check_miss_vertex(connectionList);
    } else {
        add_miss_vertex(connectionList);
    }

    cy.elements().forEach(function (element, index) {
        element.removeClass("highlight");
    });

    cy.edges().filter(function (value) {
        return value.hasClass("fake");
    }).remove();

    for (var i = 0; i < connectionList.length; ++i) {
        if (connectionList[i]["type"] === "none") {
            continue;
        }

        var cntHighlight = 0;

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
            cntHighlight += 1;
           edge.addClass("highlight");
        });

        edge_id = "f" + connectionList[i]["from"] + "_" + connectionList[i]["to"];
        if (cntHighlight === 0) {
            cy.add({
                group: 'edges',
                data: {
                    id: edge_id,
                    source: connectionList[i]["from"],
                    target: connectionList[i]["to"],
                    label: "",
                    faveColor: "#000",
                    weight: 1,
                    special: 0,
                    lstyle: 'solid'}
            });
        }
    }

    cy.edges().filter(function (ele) {
        return ele.data("id")[0] === 'f';
    }).forEach(function (value) { value.addClass("fake"); });

    cy.nodes().filter(function (ele) {
        return elements_id.indexOf(ele.data('id')) >= 0;
    }).forEach(function (node, index) {
        node.addClass("highlight");
    });

    document.getElementById("DontRedrawButton").onclick = function(){
        document.getElementById("AlertBoxNotAllVertexPresent").style.visibility = "hidden";
    };

    document.getElementById("RedrawButton").onclick = function(){
        document.getElementById("AlertBoxNotAllVertexPresent").style.visibility = "hidden";
        updateHighlight(true);
    };
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

            if (scaffoldgraph.id_by_name.has(value)) {
                if (node.data('id') == scaffoldgraph.id_by_name.get(value)) {
                    inCy = 1;
                }
            }
        });
        return inCy;
    }

    function isNewToken(reqStr, tokens) {
        return tokens.length == 0 || (" \n\t,->".indexOf(reqStr[reqStr.length - 1]) >= 0);
    }

    $("#highlight_elements").autocomplete({
        minLength : 0,
        maxLength : 200,
        autoFocus : true,
        source : function(request, response) {
            var reqStr = request.term.toString();
            var tokens = getTokens(reqStr);
            var lastToken = "";
            if (!isNewToken(reqStr, tokens)) {
                lastToken = tokens[tokens.length - 1];
            }
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
            if (!isNewToken(this.value, tokens)) {
                tokens.pop();
            }
            tokens.push(ui.item.value);
            this.value = tokens.join(" ");
            return false;
        },
        focus: function (event, ui) {
            return false;
        }
    });
}

function highlightOnTap(cy) {
    cy.on('tap', 'node', function (evt) {
        if (cy.ignoreTap) {
            delete cy.ignoreTap;
            return
        }
        var v = evt.target.id();
        var highlight_text = document.getElementById('highlight_elements');
        highlight_text.value += " " + v;
        updateHighlight();
    });
}
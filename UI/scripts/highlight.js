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
            if (isGoodEdge(cedge.id)) {
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

function parseTextWithElements(text_elements) {
    text_elements = text_elements.replace(/;/g, ' $1 ').trim();
    text_elements = text_elements.replace(/([^-])(-)(>)/g, '$1 $2$3 ').trim();
    text_elements = text_elements.replace(/(-)(-)(>)/g, ' $1$2$3 ').trim();
    var res = text_elements.split(/[\s\n\t,]+/);
    var connectionList = [];
    var vertexList = [];
    for (var i = 0; i < res.length; ++i) {
        if (scaffoldgraph.id_by_name.has(res[i])) {
            res[i] = scaffoldgraph.id_by_name.get(res[i]).toString()
        }

        if (res[i] !== "->" && res[i] !== ';' && res[i] !== "-->") {
            vertexList.push(res[i])
        }
    }

    for (i = 1; i < res.length; ++i) {
        if (res[i] === ";") {
            connectionList.push({
                "type": "none",
                "from": res[i - 1],
                "to": res[i + 1]
            });
            i += 1;
        } else if (res[i] === "->") {
            connectionList.push({
                "type": "edge",
                "from": res[i - 1],
                "to": res[i + 1]
            });
            i += 1;
        } else if (res[i] === "-->") {
            var clist = getShortestPath(res[i - 1], res[i + 1]);
            for (var j = 0; j < clist.length; ++j) {
                connectionList.push(clist[j]);
            }
        } else {
            connectionList.push({
                "type": "edge",
                "from": res[i - 1],
                "to": res[i]
            });
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
            return ele.data('source') === connectionList[i]["from"] &&
                ele.data('target') === connectionList[i]["to"]
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
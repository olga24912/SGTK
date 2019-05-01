function parseTextWithElements(text_elements) {
    var res = text_elements.split(/[\s\n\t;,]+/);
    var connectionList = [];
    var vertexList = [];
    for (var i = 0; i < res.length; ++i) {
        if (scaffoldgraph.id_by_name.has(res[i])) {
            res[i] = scaffoldgraph.id_by_name.get(res[i]).toString()
        }

        if (res[i] !== "->") {
            vertexList.push(res[i])
        }
    }

    for (var i = 1; i < res.length; ++i) {
        if (res[i] === "->") {
            connectionList.push({
                "type": "edge",
                "from": res[i - 1],
                "to": res[i + 1]
            });
            i += 1;
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
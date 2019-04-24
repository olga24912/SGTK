function parseTextWithElements(text_elements) {
    var res = text_elements.split(/[\s\n\t;,]+/);
    /*for (var i = 0; i < res.length; ++i) {
        res[i] = parseInt(res[i])
    }*/
    return res;
}

function updateHighlight() {
    var text_elements = document.getElementById("highlight_elements").value;
    var elements_id = parseTextWithElements(text_elements);

    cy.elements().forEach(function (element, index) {
        element.removeClass("highlight");
    });

    for (var i = 0; i < elements_id.length - 1; ++i) {
        cy.edges().filter(function (ele) {
            return ele.data('source') === elements_id[i] &&
                ele.data('target') === elements_id[i + 1]
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
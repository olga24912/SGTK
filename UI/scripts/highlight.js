function parseTextWithElements(text_elements) {
    var res = text_elements.split(" ");
    /*for (var i = 0; i < res.length; ++i) {
        res[i] = parseInt(res[i])
    }*/
    return res;
}

function updateHighlight() {
    var text_elements = document.getElementById("highlight_elements").value;
    var elements_id = parseTextWithElements(text_elements);

    cy.nodes().forEach(function (node, index) {
        node.removeClass("highlight");
    });

    cy.nodes().filter(function (ele) {
        return elements_id.indexOf(ele.data('id')) >= 0;
    }).forEach(function (node, index) {
        node.addClass("highlight");
    });
}
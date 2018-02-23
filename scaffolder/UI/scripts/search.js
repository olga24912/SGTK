

function setupAutocompleteSearch(){
    $( "#search_text" ).autocomplete({
        minLength: 1,
        maxHeight: 200,
        deferRequestBy: 50,
        source:  function(request, response) {
            var nodes = cy.filter('node');
            var result = [];
            for (var i = 0; i < nodes.length; ++i) {
                result.push(nodes[i].id().toString());
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
    var vid = document.getElementById("search_text").value;
    cy.fit(cy.$('#' + vid))
}

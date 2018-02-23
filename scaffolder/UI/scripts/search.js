function search() {
    var vid = document.getElementById("search_text").value;

    cy.$('#' + vid).style({'border-width': '5px'});
    cy.fit(cy.$('#' + vid))
}
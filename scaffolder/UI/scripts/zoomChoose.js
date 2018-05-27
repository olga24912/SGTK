function zoomDropDown() {
    document.getElementById("zoomDropdown").classList.toggle("show");
}

function zoomPlus() {
    if (cy !== null) {
        cy.zoom(cy.zoom() * 1.5);
    }
}

function zoomMinus() {
    if (cy !== null) {
        cy.zoom(cy.zoom() / 1.5);
    }
}

// Close the dropdown menu if the user clicks outside of it
window.onclick = function(event) {
    if (!event.target.matches('.zoomSelect')) {
        var dropdowns = document.getElementsByClassName("dropdown-content");
        var i;
        for (i = 0; i < dropdowns.length; i++) {
            var openDropdown = dropdowns[i];
            if (openDropdown.classList.contains('show')) {
                openDropdown.classList.remove('show');
            }
        }
    }
};

window.addEventListener("keyup", function (evt) {
    if (evt.altKey) {
        var code = (evt.keyCode || evt.which);
        if (code === 187) {
            zoomPlus();
        } else if (code === 189) {
            zoomMinus();
        }
    }
});

function updateZoomFromInput() {
    var inputVal = parseInt(document.getElementById("zoomInput").innerText);

    if (cy !== null) {
        if (inputVal > 0 && inputVal < 10000) {
            cy.zoom(inputVal/100);
        }
        (document.getElementById("zoomInput")).innerText = (cy.zoom() * 100).toString() +  "%";
    } else {
        (document.getElementById("zoomInput")).innerText = "100%";
    }
}

document.getElementById("zoomInput").addEventListener("focusout", updateZoomFromInput);

document.getElementById("zoomInput").addEventListener("keypress", function (e) {
    //if enter press
    if (e.keyCode === 13) {
        document.getElementById("zoomInput").blur();
        updateZoomFromInput()
    }
});
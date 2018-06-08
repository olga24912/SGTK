function zoomDropDown() {
    document.getElementById("zoomDropdown").classList.toggle("show");
}

function zoomPlus() {
    if (cy !== null) {
        var width =  document.getElementById('mainpanel').clientWidth;
        var height =  document.getElementById('mainpanel').clientHeight;
        cy.panBy({x: -(width - width/1.5)/2, y: -(height - height/1.5)/2});
        cy.pan({x: cy.pan().x * 1.5, y: cy.pan().y * 1.5});
        cy.zoom(cy.zoom() * 1.5);
        if ((document.getElementById("select_layout").value !== "free_layout")) {
            cy.pan({x:cy.pan().x, y: height/3});
        }
    }
}

function zoomMinus() {
    if (cy !== null) {
        var width =  document.getElementById('mainpanel').clientWidth;
        var height =  document.getElementById('mainpanel').clientHeight;
        cy.panBy({x: -(width - width*1.5)/2, y: -(height - height*1.5)/2});
        cy.pan({x: cy.pan().x / 1.5, y: cy.pan().y / 1.5});
        cy.zoom(cy.zoom() / 1.5);
        if ((document.getElementById("select_layout").value !== "free_layout")) {
            cy.pan({x: cy.pan().x, y: height/3});
        }
        //cy.pan({x: 0, y: -100});
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

window.addEventListener("keydown", function (evt) {
    var delta = 20;
    if (evt.shiftKey) {
        delta = 100;
    }
    var code = (evt.keyCode || evt.which);
    if (code === 37) { //left arrow
        cy.panBy({
            x: delta,
            y:0
        });
    } else if (code === 39) { //right arrow
        cy.panBy({
            x:-delta,
            y:0
        });
    } else if (code === 38) { //up arrow
        cy.panBy({
            x:0,
            y:delta
        });
    } else if (code === 40) { //down arrow
        cy.panBy({
            x:0,
            y:-delta
        });
    }
});


function updateZoomFromInput() {
    var inputVal = parseInt(document.getElementById("zoomInput").innerText);

    if (cy !== null) {
        if (inputVal > 0 && inputVal < 10000) {
            cy.zoom((inputVal*defZoom/(100 * 100)));
        }
        (document.getElementById("zoomInput")).innerText = (cy.zoom() * 100 * 100/defZoom).toString() +  "%";
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
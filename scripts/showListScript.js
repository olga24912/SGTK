/*
* Create list of components
*/

function createComponentShowList(draw, write, get_title, cnt_comp) {
    cur_show_id = 0;
    var show_block = document.getElementById("show_block");
    show_block.innerHTML = "";
    for (var i = 0; i < cnt_comp; ++i) {
        var divb = document.createElement("div");
        divb.className = "block";
        var dive = document.createElement("div");
        dive.className = "pic";
        dive.id = "pic_" + i;
        var pe = document.createElement("p");
        pe.innerHTML = write(i);
        dive.appendChild(pe);
        dive.addEventListener("click", function (){
            (document.getElementById("pic_" + cur_show_id)).className = "pic";
            this.className = "choosePic";
            cur_show_id = parseInt(this.id.substring(4));

            (document.getElementById("choose_page")).value = cur_show_id;
            (document.getElementById("titel_text")).innerHTML = get_title(cur_show_id);

            (document.getElementById("mainpanel")).innerHTML = '<div class="loader" id="loader"></div>';
            setTimeout(function() {
                draw(cur_show_id);
                (document.getElementById("loader")).style.display = 'none';
            }, 100);
        });
        divb.appendChild(dive);
        show_block.appendChild(divb);
    }

    document.getElementById("cnt_page").innerHTML = "<input type=\"number\" min=\"0\" max=" + cnt_comp + " id = \"choose_page\">\n" +
        "                of "+ cnt_comp;

    document.getElementById("choose_page").addEventListener("keypress", function (e) {
        //if enter press
        if (e.keyCode === 13) {
            if ((document.getElementById("choose_page")).value < cnt_comp &&
                (document.getElementById("choose_page")).value >= 0) {
                (document.getElementById("pic_" + cur_show_id)).className = "pic";
                cur_show_id = (document.getElementById("choose_page")).value;
                (document.getElementById("pic_" + cur_show_id)).className = "choosePic";
                (document.getElementById("titel_text")).innerHTML = get_title(cur_show_id);
                (document.getElementById("show_block")).scrollTop = (document.getElementById("pic_" + cur_show_id)).offsetTop - 50;
                (document.getElementById("mainpanel")).innerHTML = '<div class="loader" id="loader"></div>';
                setTimeout(function() {
                    draw(cur_show_id);
                    (document.getElementById("loader")).style.display = 'none';
                }, 100);
            }
        }
    });

    document.getElementById("choose_page").addEventListener("focusout", function() {
        (document.getElementById("choose_page")).value = cur_show_id;
    });

    if (cnt_comp > 0) {
        cur_show_id = 0;
        setTimeout(function() {
            draw(cur_show_id);
            (document.getElementById("loader")).style.display = 'none';
        }, 100);
    }
}


$(document).ready(function() {

    var myVar = setInterval(update, 3000);

    function update() {
        $.getJSON("/status/get_update")
            .done(function (msg) {
                $.each(msg, function (k, v) {

                    if (document.getElementById(k) !== null) {
                        var element = document.getElementById(k)
                        if ((v.includes('not') || v.includes('stopped') || v.includes('closed')) && element.tagName == 'TD') {
                            element.innerHTML = v;
                            element.style.background = '#ff6666'
                        }
                        else if (element.tagName == 'TD') {
                            element.innerHTML = v;
                            element.style.background = 'white'
                        }
                        else {
                            element.innerHTML = v;
                        }

                    }

                });


            });
    }


});
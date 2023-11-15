document.getElementById("fileIn").addEventListener("change", function(e) {

            let files = e.target.files;
            var arr = new Array(files.length*1);
            var fname = new Array(1);
            for (let i=0; i<files.length; i++) {

            //console.log(files[i].webkitRelativePath);
            //console.log(files[i].name);
           
            fname = files[1].webkitRelativePath;
            arr[i] = files[i].name;


            }

            Shiny.onInputChange("mydata", arr);
            Shiny.onInputChange("mypath", fname);
    });

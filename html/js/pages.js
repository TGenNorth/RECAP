// default display of 20 contigs on a single page
var maxContig = 30;
var currentPage = 0;

// for imported data
var samples;
var contigs;
var display = "snps";

// handle multiple pages and navigation between them
jQuery(document).ready(function($) {

  // on window load (only runs once)
  window.history.replaceState(null, null, '#1') // push to default window state (home page)
  var speed = 0;                                // animation speed for window state change
  readData();                                   // read JSONP file and execute function

  // on window state change (runs every time a navigation bar link is clicked)
  if (window.history && window.history.replaceState) {
    $(window).on('popstate', function() {

      $('.loading').css('display', 'block');
      var url = window.location.hash;
      $(".active").removeClass("active");
      $('.page.current').css('opacity', 0);
      $('.page.current').removeClass('current');

      if(url == '#1'){
        var i=1;
        $(".button1").addClass("active");
      }
      else if(url == '#2'){
        var i=2;
        $(".button2").addClass("active");
      }
      else if(url == '#3'){
        var i=3;
        $(".button3").addClass("active");
      }
      else if(url == '#4'){
        var i=4;
        $(".button4").addClass("active");
      }
      else if(url == '#5'){
        var i=5;
        $(".button5").addClass("active");
      }
      else if(url == ""){
        var i=1;
        $(".button1").addClass("active");
      }
      else{
        var i=1;
        $(".button1").addClass("active");
        window.history.replaceState(null, null, '#1')
      }
      $('.page').eq(i).addClass('current');

      // build graphs before changing CSS
      var timer = setTimeout(function() {
        if (i === 2) {
          drawContigs();
        }
        $('.loading').css('display', 'none');
        $('.page.current').css('opacity', 1);
      }, 0);
    });
  }
});

/*
TODO:
  1. make custom loading image
  3. add search option for samples & contigs
  4. fix contigs with single position (cannot see on graph)
*/

// create and build all graph objects
function drawContigs() {
  checkedSamples = CHECKBOXreturnChecked();
  var out = [];

  var max;
  if ((currentPage*maxContig)+maxContig > contigs.length) {
    max = contigs.length;
  } else {
    max = (currentPage*maxContig)+maxContig;
  }

  document.getElementById("graphs").innerHTML = ''
  for (var i = currentPage*maxContig; i < max; i++) {
    document.getElementById("graphs").innerHTML +=[
     "  <div class='graph'>",
     "  <div class='container-fluid'>",
     "  <div class='row'>",
     "    <div class='col-xs-12 col-sm-12 col-md-12'>",
     "      <h2><font size='5'>Contig "+parseInt(i+1)+":</font><font size='4'> "+contigs[i].name+"</font></h2>",
     "    </div>",
     //"    <div class='col-xs-2 col-sm-2 col-md-2'><button type='button' class='btn btn-primary btn-block btn-sm'>SNPs</button></div>",
     //"    <div class='col-xs-2 col-sm-2 col-md-2'><button type='button' class='btn btn-outline-primary btn-block btn-sm'>Depth</button></div>",
     //"    <div class='col-xs-2 col-sm-2 col-md-2'><button type='button' class='btn btn-outline-primary btn-block btn-sm'>Phi</button></div>",
     "    <div class='col-xs-12 col-sm-12 col-md-12'>",
     "      <div id='graph"+i+"' style='margin:auto; max-width:100%; height:250px;'></div>",
     "    </div>",
     "    <div class='col-xs-12 col-sm-12 col-md-12'>",
     "      <div id='legend"+i+"' style='margin:auto; height:auto; max-width:100%; font-size: small;'></div>",
     "    </div>",
     "  </div>",
     "<br>",
     "</div>",
     "</div>"].join('\n');
    
    out.push("position")
    if (display === "phi") {
      out[i] += ",phi"
    } else {
      for (index in checkedSamples) {
        out[i] += ","+samples[index]
      }
    }
    out[i] += "\n"
    for (var j = 0; j < contigs[i].data.length; j++) {
      var line = "";
      line += contigs[i].data[j].position[0];
      if (display === "snps") {
        for (index in checkedSamples) {
          line += ","+contigs[i].data[j].SNPs[index]
        }
      } else if (display === "depth") {
        for (index in checkedSamples) {
          line += ","+contigs[i].data[j].depth[index]
        }
      } else if (display === "phi") {
        line += ","+Number.parseFloat(contigs[i].data[j].phi)
      }
      out[i] += line + "\n"
    }
  }
  // fill each graph with saved data
  var yLabel;
  if (display === "snps") {
    yLabel = "Number of SNPs"
  } else if (display === "depth") {
    yLabel = "Depth of SNPs"
  } else if (display === "phi") {
    yLabel = "PHI P-value"
  }
  for (var i = currentPage*maxContig; i < max; i++) {
    var g = new Dygraph(
      document.getElementById("graph"+i),
      out[i],
      {
        ylabel: yLabel,
        labelsShowZeroValues: true,
        labelsSeparateLines: true,
        labelsDiv: 'legend'+i,
        legend: 'never',
        valueRange: [null, null],
        showRangeSelector: false,
        rangeSelectorHeight: 25,
        maxNumberWidth: 4,
        highlightSeriesOpts: {
          strokeWidth: 3,
          strokeBorderWidth: 1,
          highlightCircleSize: 5
        },
        /*interactionModel : {
            'mousedown' : downV3,
            'mousemove' : moveV3,
            'mouseup' : upV3,
            'click' : clickV3,
            'dblclick' : dblClickV3,
            'mousewheel' : scrollV3
        }*/
        
        //highlightCallback: function(e, x, pts, row) {
        //  document.getElementById('status').innerHTML = pts[0].xval + ', ' + pts[0].name +', '+ pts[0].yval;
        //}
      }
    );
    // reverse y-axis for phi statistic
    if (display === "phi") {
      g.updateOptions({
        valueRange: [1,0]
      });
    }
  }
  buildNavigation();
  buildDisplay();
  buildSort();
}

// contig display buttons
function CONTIGdisplay(option) {
  display = option;
  CONTIGnavigation(0);
}
function buildDisplay() {
  document.getElementById("CONTIG-navigation").innerHTML += [
   "<div class='dropdown' style='margin-right:20px;'>",
   "  <button class='btn btn-secondary btn-sm dropdown-toggle' type='button' id='dropdownMenuButton' data-toggle='dropdown' aria-haspopup='true' aria-expanded='false'>Display Contigs</button>",
   "  <div class='dropdown-menu' aria-labelledby='dropdownMenuButton'>",
   "    <span class='dropdown-item' onclick='CONTIGdisplay(\"snps\")'>SNPs</span>",
   "    <span class='dropdown-item' onclick='CONTIGdisplay(\"depth\")'>Depth</span>",
   "    <span class='dropdown-item' onclick='CONTIGdisplay(\"phi\")'>PHI</span>",
   "  </div>",
   "</div>"].join('\n');
}

// contig sort buttons
function CONTIGsort(option) {
  if (option === "name") {
    contigs.sort(function(a,b){
      if(a.name < b.name)
        return -1;
      if(a.name > b.name)
        return 1;
      return 0;
    });
  } else if (option === "length") {
    contigs.sort(function(a,b){
      if(a.data.length > b.data.length)
        return -1;
      if(a.data.length < b.data.length)
        return 1;
      return 0;
    });
  } else if (option === "aggregate") {
    contigs.sort(function(a,b){
      var aTotal = 0;
      var bTotal = 0;
      for (var i = 0; i < a.data.length; i++) {
        aTotal += a.data[i].aggregate;
      }
      for (var i = 0; i < b.data.length; i++) {
        bTotal += b.data[i].aggregate;
      }
      if(aTotal/a.data.length > bTotal/b.data.length)
        return -1;
      if(aTotal/a.data.length < bTotal/b.data.length)
        return 1;
      return 0;
    });
  } else if (option === "depth") {
    contigs.sort(function(a,b){
      var aTotal = 0;
      var bTotal = 0;
      for (var i = 0; i < a.data.length; i++) {
        var aTemp;
        for (var j = 0; j < a.data[i].depth.length; j++) {
          aTemp += a.data[i].depth[j];
        }
        aTotal += aTemp/a.data[i].depth.length;
      }
      for (var i = 0; i < b.data.length; i++) {
        var bTemp;
        for (var j = 0; j < b.data[i].depth.length; j++) {
          bTemp += b.data[i].depth[j];
        }
        bTotal += bTemp/b.data[i].depth.length;
      }
      if(aTotal/a.data.length > bTotal/b.data.length)
        return -1;
      if(aTotal/a.data.length < bTotal/b.data.length)
        return 1;
      return 0;
    });
  } else if (option === "phi") {
    contigs.sort(function(a,b){
      var aTotal = 0;
      var bTotal = 0;
      for (var i = 0; i < a.data.length; i++) {
        aTotal += Number.parseFloat(a.data[i].phi);
      }
      for (var i = 0; i < b.data.length; i++) {
        bTotal += Number.parseFloat(b.data[i].phi);
      }
      if(aTotal/a.data.length > bTotal/b.data.length)
        return 1;
      if(aTotal/a.data.length < bTotal/b.data.length)
        return -1;
      return 0;
    });
  }
  CONTIGnavigation(0);
}
function buildSort() {
  document.getElementById("CONTIG-navigation").innerHTML += [
   "<div class='dropdown'>",
   "  <button class='btn btn-secondary btn-sm dropdown-toggle' type='button' id='dropdownMenuButton' data-toggle='dropdown' aria-haspopup='true' aria-expanded='false'>Sort Contigs</button>",
   "  <div class='dropdown-menu' aria-labelledby='dropdownMenuButton'>",
   "    <span class='dropdown-item' onclick='CONTIGsort(\"name\")'>Name</span>",
   "    <span class='dropdown-item' onclick='CONTIGsort(\"length\")'>Length</span>",
   "    <span class='dropdown-item' onclick='CONTIGsort(\"aggregate\")'>SNPs</span>",
   "    <span class='dropdown-item' onclick='CONTIGsort(\"depth\")'>Depth</span>",
   "    <span class='dropdown-item' onclick='CONTIGsort(\"phi\")'>PHI</span>",
   "  </div>",
   "</div>"].join('\n');
}

// contig navigation buttons
function CONTIGnavigation(page) {
  $('.loading').css('display', 'block');
  if (page === "-") {
    currentPage = parseInt(currentPage - 1);
  }
  else if(page === "+") {
    currentPage = parseInt(currentPage + 1);
  }
  else {
    currentPage = parseInt(page);
  }
  var timer = setTimeout(function() {
    drawContigs();
    $('.loading').css('display', 'none');
  }, 0);
}
function buildNavigation() {
  if (currentPage === 0) {
    document.getElementById("CONTIG-navigation").innerHTML = '<li class="page-item disabled"><span class="page-link" onclick="CONTIGnavigation(\'-\')">Prev</span></li>'
  } else {
    document.getElementById("CONTIG-navigation").innerHTML = '<li class="page-item"><span class="page-link" onclick="CONTIGnavigation(\'-\')">Prev</span></li>'
  }
  for (var i = 0; i < Math.ceil(contigs.length/maxContig); i++) {
    if (i === currentPage) {
      document.getElementById("CONTIG-navigation").innerHTML += '<li class="page-item disabled"><span class="page-link" onclick="CONTIGnavigation(\''+i+'\')">'+parseInt(i+1)+'</span></li>'
    } else {
      document.getElementById("CONTIG-navigation").innerHTML += '<li class="page-item"><span class="page-link" onclick="CONTIGnavigation(\''+i+'\')">'+parseInt(i+1)+'</span></li>'
    }
  }
  if (currentPage+1 === Math.ceil(contigs.length/maxContig)) {
    document.getElementById("CONTIG-navigation").innerHTML += '<li class="page-item disabled"><span class="page-link" onclick="CONTIGnavigation(\'+\')">Next</span></li>'

  } else {
    document.getElementById("CONTIG-navigation").innerHTML += '<li class="page-item"><span class="page-link" style="margin-right:20px;" onclick="CONTIGnavigation(\'+\')">Next</span></li>'
  }
}

// checkbox buttons
function CHECKBOXselectAll() {
  for (var i = 0; i < samples.length; i++) {
    document.getElementById("CHECKBOX-samples"+samples[i]).checked = true;
  }
}
function CHECKBOXselectNone() {
  for (var i = 0; i < samples.length; i++) {
    document.getElementById("CHECKBOX-samples"+samples[i]).checked = false;
  }
}
function CHECKBOXreturnChecked() {
  checkedSamples = [];
  for (var i = 0; i < samples.length; i++) {
    if (document.getElementById("CHECKBOX-samples"+samples[i]).checked === true) {
      checkedSamples.push(i);
    }
  }
  return checkedSamples;
}

// read in jsonp file and pass the buildGraphs
function readData() {
  var file = document.createElement("script");
  file.src = "snpDensityMatrix.jsonp?callback=buildGraphs";
  document.body.insertBefore(file, document.body.firstChild);
}

// save graph data to innerHTML
function buildGraphs(data) {
  samples = data.samples;
  contigs = data.contigs;
  for (var i = 0; i < data.samples.length; i++) {
    document.getElementById("SETTINGS-samples").innerHTML += '<input type="checkbox" id="CHECKBOX-samples'+data.samples[i]+'" checked> '+data.samples[i]+'<br>'
  }
}

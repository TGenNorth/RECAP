/*
TODO:
  1. make custom loading image
  3. add search option for samples & contigs
  4. fix contigs with single position (cannot see on graph)
*/

// ================================================== ( Global Variables )
// amount of contigs per page, and current page of contigs
// after sorting or changing graph type for contigs, current page defaults to 0
var maxContig = 20;
var currentPage = 0;

// variables to store imported JSON data
var samples;
var contigs;

// for current graph type on display
var display = "snps";

// ================================================== ( Page Display Function )
// handle multiple pages and navigation between them
jQuery(document).ready(function($) {
  // on window load (only runs once)
  window.history.replaceState(null, null, '#1') // push to default window state (home page)
  readData(); // read JSONP file and execute function
  // on window state change (runs every time a navigation bar link is clicked)
  if (window.history && window.history.replaceState) {
    $(window).on('popstate', function() {
      $('.loading').css('display', 'block');
      var url = window.location.hash;
      $(".active").removeClass("active");
      $('.page.current').css('opacity', 0);
      $('.page.current').removeClass('current');
      if(url == '#1') {
        var i=1;
        $(".button1").addClass("active");
      } else if(url == '#2') {
        var i=2;
        $(".button2").addClass("active");
      } /*else if(url == '#3') {
        var i=3;
        $(".button3").addClass("active");
      } else if(url == '#4') {
        var i=4;
        $(".button4").addClass("active");
      } else if(url == '#5') {
        var i=5;
        $(".button5").addClass("active");
      }*/ else if(url == "") {
        var i=1;
        $(".button1").addClass("active");
      } else {
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

// ================================================== ( Draw Contigs ) 
// create and build all graph objects
function drawContigs() {
  // only add data for samples that are checked off
  checkedSamples = CHECKBOXreturnChecked();
  var out = [];
  var max;
  if ((currentPage*maxContig)+maxContig > contigs.length) {
    max = contigs.length;
  } else {
    max = (currentPage*maxContig)+maxContig;
  }
  // empty and repopulate all graph objects
  document.getElementById("graphs").innerHTML = ''
  for (var i = currentPage*maxContig; i < max; i++) {
    document.getElementById("graphs").innerHTML +=[
     "  <div class='graph'>",
     "  <div class='container-fluid'>",
     "  <div class='row'>",
     "    <div class='col-xs-12 col-sm-12 col-md-12'>",
     "      <h2><font size='5'>Contig "+parseInt(i+1)+":</font><font size='4'> "+contigs[i].name+"</font></h2>",
     "    </div>",
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
    // get data for each graph object
    out[i] = ("position");
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
      if (contigs[i].data.length === 1) {
        out[i] += contigs[i].data[j].position[1]+line.substr(1, line.length) + "\n"
      }
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
  // graph plugin that resets graph zoom properly on dblclick,
  // resets to [0,1] range for phi, and [null,null] for other
  function getStuff(g) {
    var lowestPhi = 0.009;
    for (var i = 0; i < g.numRows(); i++) {
      if (lowestPhi > g.getValue(i,1)) {
        lowestPhi = g.getValue(i,1);
      }
    }
    return lowestPhi;
  }
  const zoomReset = {
    activate: function(g) {
      return {
        dblclick: e => {
          if (display === "phi") {
            e.dygraph.updateOptions({
              dateWindow: null,
              valueRange: [1,getStuff(g)/2],
              logscale: true,
              axes: {
                y: {
                  ticker: function(min, max, pixels) {
                    tickList = [
                      { v: 1.0 },
                      { v: 1.0, label: "1" },
                      { v: 0.05 },
                      { v: 0.05, label: "0.05" },
                    ]
                    // add custom ticks for phi graph formatting
                    if (getStuff(g) < 0.01) {
                      tickList.push({ v: 0.01 },{ v: 0.01, label: "0.01" })
                    }
                    if (getStuff(g) < 0.001) {
                      tickList.push({ v: 0.001 },{ v: 0.001, label: "0.001" })
                    }
                    if (getStuff(g) < 0.0001) {
                      tickList.push({ v: 0.0001 },{ v: 0.0001, label: "< 0.0001" })
                    } else if (getStuff(g) < 0.1) {
                      tickList.push({ v: 0.1 },{ v: 0.1, label: "0.1" })
                    }
                    return tickList;
                  }
                }
              }
            });
          } else {
            e.dygraph.updateOptions({
              dateWindow: null,
              valueRange: [null,null]
            });
          }
          e.preventDefault();
        }
      }
    }
  }
  // create each graph object
  for (var i = currentPage*maxContig; i < max; i++) {
    var g = new Dygraph(
      document.getElementById("graph"+i),
      out[i],
      {
        plugins: [
          zoomReset
        ],
        axes: {
          y: {
            axisLabelWidth: 80,
          }
        },
        ylabel: yLabel,
        labelsShowZeroValues: true,
        labelsSeparateLines: true,
        labelsDiv: 'legend'+i,
        legend: 'never',
        valueRange: [null, null],
        yRangePad: 1,
        digitsAfterDecimal: 3,
        animatedZooms: true,
        showRangeSelector: false,
        rangeSelectorHeight: 25,
        highlightSeriesOpts: {
          strokeWidth: 2,
          strokeBorderWidth: 1,
          highlightCircleSize: 4
        },
        
      }
    );
    // reverse y-axis for phi statistic
    if (display === "phi") {
      g.updateOptions({
        valueRange: [1,getStuff(g)/2],
        logscale: true,
        axes: {
          y: {
            ticker: function(min, max, pixels) {
              tickList = [
                { v: 1.0 },
                { v: 1.0, label: "1" },
                { v: 0.05 },
                { v: 0.05, label: "0.05" },
              ]
              // add custom ticks for phi graph formatting
              if (getStuff(g) < 0.01) {
                tickList.push({ v: 0.01 },{ v: 0.01, label: "0.01" })
              }
              if (getStuff(g) < 0.001) {
                tickList.push({ v: 0.001 },{ v: 0.001, label: "0.001" })
              }
              if (getStuff(g) < 0.0001) {
                tickList.push({ v: 0.0001 },{ v: 0.0001, label: "< 0.0001" })
              } else if (getStuff(g) < 0.1) {
                tickList.push({ v: 0.1 },{ v: 0.1, label: "0.1" })
              }
              return tickList;
            }
          }
        },
        // highlight portions of graph where phi statistic is under 0.05 "statistically significant"
        underlayCallback: function(canvas, area, g) {
          canvas.fillStyle = "rgba(255, 216, 0, 0.5)";
          function highlight(start, end) {
            var left = g.toDomXCoord(start);
            var right = g.toDomXCoord(end);
            var width = right - left;
            canvas.fillRect(left, area.y, width, area.h);
          }
          var minX = g.getValue(0,0);
          var maxX = g.getValue(g.numRows()-1,0);
          // check slope crossing 0.05 threshold and find X given Y=0.05
          for (var i = 0; i < g.numRows()-1; i++) {
            if (g.getValue(i,1) >= 0.05 && g.getValue(i+1,1) < 0.05) {
              var start = g.getValue(i,0);
              //console.log("slopeU:",Math.log10(g.getValue(i+1,1))-Math.log10(g.getValue(i,1)))
              //console.log("y    U:",g.getValue(i+1,1)*Math.pow((10),Math.log10(g.getValue(i,1))-Math.log10(g.getValue(i+1,1))))
              //console.log("phi  U:",g.getValue(i,1))
              //console.log("X1   U:",(Math.log10(0.05/g.getValue(i+1,1)))/(Math.log10(g.getValue(i+1,1))-Math.log10(g.getValue(i,1))))
              //console.log("X2   U:",(Math.log10(g.getValue(i+1,1)/g.getValue(i+1,1)))/(Math.log10(g.getValue(i+1,1))-Math.log10(g.getValue(i,1))))
              highlight(g.getValue(i+1,0)+((Math.log10(0.05/g.getValue(i+1,1)))/(Math.log10(g.getValue(i+1,1))-Math.log10(g.getValue(i,1)))*(g.getValue(i+1,0)-g.getValue(i,0))),g.getValue(i+1,0))
            } else if (g.getValue(i,1) < 0.05 && g.getValue(i+1,1) >= 0.05) {
              var start = g.getValue(i,0);
              //console.log("slopeD:",Math.log10(g.getValue(i,1))-Math.log10(g.getValue(i+1,1)))
              //console.log("y    D:",g.getValue(i,1)*Math.pow((10),Math.log10(g.getValue(i+1,1))-Math.log10(g.getValue(i,1))))
              //console.log("phi  D:",g.getValue(i+1,1))
              highlight(g.getValue(i+1,0)+((Math.log10(0.05/g.getValue(i+1,1)))/(Math.log10(g.getValue(i,1))-Math.log10(g.getValue(i+1,1)))*(g.getValue(i,0)-g.getValue(i+1,0))),g.getValue(i,0))
            } else if (g.getValue(i,1) < 0.05 && g.getValue(i+1,1) < 0.05) {
              highlight(g.getValue(i,0),g.getValue(i+1,0))
            }
          }
        }
      });
    }
  }
  // build all interactive buttons
  if (contigs.length > 5) {
    buildQuantity();
  }
  buildDisplay();
  buildSort();
  buildNavigation();
  buildBottom();
}

// ================================================== ( Build Buttons ) 
// drop down buton to change amount of contigs per page
function buildQuantity() {
  var buttonData = [
   "<div class='dropdown' style='margin-left:10px;margin-right:10px;'>",
   "  <button class='btn btn-secondary btn-sm dropdown-toggle' type='button' id='dropdownMenuButton' data-toggle='dropdown' aria-haspopup='true' aria-expanded='false'>Contigs Per Page</button>",
   "  <div class='dropdown-menu' aria-labelledby='dropdownMenuButton'>"]
  if (contigs.length > 5) {
    buttonData.push("<span class='dropdown-item' onclick='CONTIGquantity(5)'>5</span>");
  }
  if (contigs.length > 10) {
    buttonData.push("<span class='dropdown-item' onclick='CONTIGquantity(10)'>10</span>");
  }
  if (contigs.length > 20) {
    buttonData.push("<span class='dropdown-item' onclick='CONTIGquantity(20)'>20</span>");
  }
  if (contigs.length > 50) {
    buttonData.push("<span class='dropdown-item' onclick='CONTIGquantity(50)'>50</span>");
  }
  if (contigs.length > 100) {
    buttonData.push("<span class='dropdown-item' onclick='CONTIGquantity(100)'>100</span>");
  }
  if (contigs.length > 250) {
    buttonData.push("<span class='dropdown-item' onclick='CONTIGquantity(250)'>250</span>");
  }
  buttonData.push(
   "    <span class='dropdown-item' onclick='CONTIGquantity(\"all\")'>All</span>",
   "  </div>",
   "</div>" 
  );
  document.getElementById("CONTIG-options").innerHTML = buttonData.join('\n');
}

// drop down button to change graph type
function buildDisplay() {
  document.getElementById("CONTIG-options").innerHTML += [
   "<div class='dropdown' style='margin-left:10px;margin-right:10px;'>",
   "  <button class='btn btn-secondary btn-sm dropdown-toggle' type='button' id='dropdownMenuButton' data-toggle='dropdown' aria-haspopup='true' aria-expanded='false'>Graph Type</button>",
   "  <div class='dropdown-menu' aria-labelledby='dropdownMenuButton'>",
   "    <span class='dropdown-item' onclick='CONTIGdisplay(\"snps\")'>SNPs</span>",
   "    <span class='dropdown-item' onclick='CONTIGdisplay(\"depth\")'>Depth</span>",
   "    <span class='dropdown-item' onclick='CONTIGdisplay(\"phi\")'>PHI</span>",
   "  </div>",
   "</div>"].join('\n');
}

// drop down button to change contig sort order
function buildSort() {
  document.getElementById("CONTIG-options").innerHTML += [
   "<div class='dropdown' style='margin-left:10px;margin-right:10px;'>",
   "  <button class='btn btn-secondary btn-sm dropdown-toggle' type='button' id='dropdownMenuButton' data-toggle='dropdown' aria-haspopup='true' aria-expanded='false'>Sort</button>",
   "  <div class='dropdown-menu' aria-labelledby='dropdownMenuButton'>",
   "    <span class='dropdown-item' onclick='CONTIGsort(\"name\")'>Name</span>",
   "    <span class='dropdown-item' onclick='CONTIGsort(\"length\")'>Length</span>",
   "    <span class='dropdown-item' onclick='CONTIGsort(\"aggregate\")'>SNPs</span>",
   "    <span class='dropdown-item' onclick='CONTIGsort(\"depth\")'>Depth</span>",
   "    <span class='dropdown-item' onclick='CONTIGsort(\"phi\")'>PHI</span>",
   "  </div>",
   "</div>"].join('\n');
}

// buttons used for contig page navigation
function buildNavigation() {
  if (currentPage === 0) {
    document.getElementById("CONTIG-navigation").innerHTML = '<li class="page-item disabled"><span class="page-link" onclick="CONTIGnavigation(\'---\')">\<\<</span></li><li class="page-item disabled"><span class="page-link" onclick="CONTIGnavigation(\'-\')">Prev</span></li>'
  } else {
    document.getElementById("CONTIG-navigation").innerHTML = '<li class="page-item"><span class="page-link" onclick="CONTIGnavigation(\'---\')">\<\<</span></li><li class="page-item"><span class="page-link" onclick="CONTIGnavigation(\'-\')">Prev</span></li>'
  }
  var navStart = 0;
  if (currentPage-4 > 0) {
    navStart = currentPage - 4;
    document.getElementById("CONTIG-navigation").innerHTML += '<li class="page-item"><span class="page-link" onclick="CONTIGnavigation(\'--\')">...</span></li>'
  }
  var navEnd = Math.ceil(contigs.length/maxContig);
  if (currentPage+4 < Math.ceil(contigs.length/maxContig)) {
    navEnd = currentPage + 4;
  }
  for (var i = navStart; i < navEnd; i++) {
    var lastContig = contigs.length;
    if (lastContig > parseInt(maxContig*(i+1))) {
      lastContig = parseInt(maxContig*(i+1));
    }
    if (i === currentPage) {
      document.getElementById("CONTIG-navigation").innerHTML += '<li class="page-item disabled"><span class="page-link" onclick="CONTIGnavigation(\''+i+'\')">'+parseInt((maxContig*i)+1)+'-'+lastContig+'</span></li>'
    } else {
      document.getElementById("CONTIG-navigation").innerHTML += '<li class="page-item"><span class="page-link" onclick="CONTIGnavigation(\''+i+'\')">'+parseInt((maxContig*i)+1)+'-'+lastContig+'</span></li>'
    }
  }
  if (currentPage+4 < Math.ceil(contigs.length/maxContig)) {
    document.getElementById("CONTIG-navigation").innerHTML += '<li class="page-item"><span class="page-link" onclick="CONTIGnavigation(\'++\')">...</span></li>'
  }
  if (currentPage+1 === Math.ceil(contigs.length/maxContig)) {
    document.getElementById("CONTIG-navigation").innerHTML += '<li class="page-item disabled"><span class="page-link" onclick="CONTIGnavigation(\'+\')">Next</span></li><li class="page-item disabled"><span class="page-link" onclick="CONTIGnavigation(\'+++\')">\>\></span></li>'

  } else {
    document.getElementById("CONTIG-navigation").innerHTML += '<li class="page-item"><span class="page-link" onclick="CONTIGnavigation(\'+\')">Next</span></li><li class="page-item"><span class="page-link" onclick="CONTIGnavigation(\'+++\')">\>\></span></li>'
  }
}

function buildBottom() {
  document.getElementById("CONTIG-navigationBottom").innerHTML = document.getElementById("CONTIG-navigation").innerHTML
}

// ================================================== ( Button Functions )
// change page and redraw contig graphs
// used to selectively draw only a select few contigs
function CONTIGnavigation(page) {
  // show loading icon
  $('.loading').css('display', 'block');
  // prev, next, and manual selection for contig navigation
  if (page === "-") {
    currentPage = parseInt(currentPage - 1);
  } else if(page === "--") {
    currentPage = parseInt(currentPage - 4);
  } else if(page === "---") {
    currentPage = 0;
  } else if(page === "+") {
    currentPage = parseInt(currentPage + 1);
  } else if(page === "++") {
    currentPage = parseInt(currentPage + 4);
  } else if(page === "+++") {
    currentPage = Math.floor((contigs.length-1)/maxContig);
  } else {
    currentPage = parseInt(page);
  }
  // remove loading icon
  var timer = setTimeout(function() {
    drawContigs();
    $('.loading').css('display', 'none');
  }, 0);
}

// set graph type and reset pages
function CONTIGdisplay(option) {
  display = option;
  // reset pages
  CONTIGnavigation(0);
}

// set amount of contigs per page
function CONTIGquantity(option) {
  if (option === "all") {
    maxContig = contigs.length;
  } else {
  maxContig = option;
  }
  CONTIGnavigation(0);
}

// sort contig data and reset pages
function CONTIGsort(option) {
  // sort contigs by name (default option)
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
  // sort contigs by aggregate (SNP density)
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
  // sort contigs by depth
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
  // sort contigs by PHI p-value (lowest first)
  } else if (option === "phi") {
    contigs.sort(function(a,b){
      var aBestPhi = 1;
      var bBestPhi = 1;
      for (var i = 0; i < a.data.length; i++) {
        if (Number.parseFloat(a.data[i].phi) < aBestPhi) {
          aBestPhi = Number.parseFloat(a.data[i].phi);
        }
      }
      for (var i = 0; i < b.data.length; i++) {
        if (Number.parseFloat(b.data[i].phi) < bBestPhi) {
          bBestPhi = Number.parseFloat(b.data[i].phi);
        }
      }
      if(aBestPhi > bBestPhi)
        return 1;
      if(aBestPhi < bBestPhi)
        return -1;
      return 0;
    });
  }
  // reset pages
  CONTIGnavigation(0);
}

// ================================================== ( Sample Checkbox Buttons )
// select all samples
function CHECKBOXselectAll() {
  for (var i = 0; i < samples.length; i++) {
    document.getElementById("CHECKBOX-samples"+samples[i]).checked = true;
  }
}

// select no samples
function CHECKBOXselectNone() {
  for (var i = 0; i < samples.length; i++) {
    document.getElementById("CHECKBOX-samples"+samples[i]).checked = false;
  }
}

// return list of only selected samples
// used for drawing contigs with only relevant samples
function CHECKBOXreturnChecked() {
  checkedSamples = [];
  for (var i = 0; i < samples.length; i++) {
    if (document.getElementById("CHECKBOX-samples"+samples[i]).checked === true) {
      checkedSamples.push(i);
    }
  }
  return checkedSamples;
}

// ================================================== ( download file )
function thing() {
  console.log("hello");
}

// ================================================== ( Process JSON )
// callback JSON file and initiate 
function readData() {
  var file = document.createElement("script");
  file.src = "snpDensityMatrix.jsonp?callback=buildGraphs";
  document.body.insertBefore(file, document.body.firstChild);
}

// import sample and contig data from JSON
function buildGraphs(data) {
  samples = data.samples;
  // pre-sort the contigs by length, this is the best default order
  data.contigs.sort(function(a,b){
    if(a.data.length > b.data.length)
      return -1;
    if(a.data.length < b.data.length)
      return 1;
    return 0;
  });
  contigs = data.contigs;
  for (var i = 0; i < data.samples.length; i++) {
    document.getElementById("SETTINGS-samples").innerHTML += '<input type="checkbox" id="CHECKBOX-samples'+data.samples[i]+'" checked> '+data.samples[i]+'<br>'
  }
}

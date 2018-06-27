/*
TODO:
  1. make custom loading image
  3. rebuild multi-fasta and SNP matrix files
*/

// ================================================== ( Global Variables )
// amount of contigs per page, and current page of contigs
// after sorting or changing graph type for contigs, current page defaults to 0
var maxContig = 20;
var currentPage = 0;
var sideWidth = 0;

// variables to store imported JSON data
var samples;
var contigs;
var dygraphs = {};

// for current graph type on display
var display = "snps";
var out = [];

// pulls graph page aside to reveal samples menu
function toggleMenu() {
  $("#sidebar-wrapper").css('width',sideWidth);
  $("#graph-wrapper").toggleClass("toggle");
  $("#graph-wrapper").css("left", sideWidth);
  $("#graph-wrapper.toggle").css("left", 0);
  $("#menu-btn").toggleClass("toggle");
  $("#menu-btn").css("left", sideWidth+10);
  $("#menu-btn.toggle").css("left", 10);
  setTimeout(function(){ window.dispatchEvent(new Event('resize')); }, 250);
}

// ================================================== ( Update Graphs )
// update all graph data based on checkboxes clicked
function updateOut() {
  checkedSamples = CHECKBOXreturnChecked();
  var max;
  if ((currentPage*maxContig)+maxContig > contigs.length) {
    max = contigs.length;
  } else {
    max = (currentPage*maxContig)+maxContig;
  }
  for (var i = currentPage*maxContig; i < max; i++) {
    // get data for each graph object
    out[i] = ("position");
    if (display === "phi") {
      out[i] += ",phi"
    } else {
      for (var k = 0; k < checkedSamples.length; k++) { 
      //for (index in checkedSamples) {
        out[i] += ","+samples[checkedSamples[k]]
      }
    }
    out[i] += "\n"
    for (var j = 0; j < contigs[i].data.length; j++) {
      var line = "";
      line += contigs[i].data[j].position[0];
      if (display === "snps") {
        for (var k = 0; k < checkedSamples.length; k++) { 
        //for (index in checkedSamples) {
          line += ","+contigs[i].data[j].SNPs[checkedSamples[k]]
        }
      } else if (display === "depth") {
        for (var k = 0; k < checkedSamples.length; k++) { 
        //for (index in checkedSamples) {
          line += ","+contigs[i].data[j].depth[checkedSamples[k]]
        }
      } else if (display === "phi") {
        line += ","+Number.parseFloat(contigs[i].data[j].phi)
      }
      out[i] += line + "\n"
      if (contigs[i].data.length === 1) {
        out[i] += contigs[i].data[j].position[1]+line.substr(1, line.length) + "\n"
      }
    }
    dygraphs[i].updateOptions(
      {
        'file':out[i]
      }
    );
  }
}

// ================================================== ( Data Export )
// create download element
function download(filename, text) {
  var element = document.createElement('a');
  element.setAttribute('href', 'data:text/plain;charset=utf-8,' + encodeURIComponent(text));
  element.setAttribute('download', filename);
  element.style.display = 'none';
  document.body.appendChild(element);
  element.click();
  document.body.removeChild(element);
}

// push data in graph range to download element
function exportContig(i){
  range = dygraphs[i].xAxisRange()
  left = Math.floor(range[0])
  right = Math.ceil(range[1])
  download("range.txt", contigs[i].name+"\t"+String(left)+"\t"+String(right))
}

// ================================================== ( Draw Contigs ) 
// create and build all graph objects
function drawContigs() {
  // only add data for samples that are checked off
  checkedSamples = CHECKBOXreturnChecked();
  
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
     "<div class='graph'>",
     "  <div class='col-12'>",
     "    <button type='button' class='btn btn-outline-secondary btn-sm' style='top:5px;right:15px;position:absolute;' onclick='exportContig("+i+")'>Export</button>",
     "    <font size='5'>Contig "+parseInt(i+1)+":</font><font size='4'> <span id='link' onclick='console.log("+i+")'>"+contigs[i].name+"</span></font>",
     "  </div>",
     "  <div class='col-12'>",
     "    <div id='legend"+i+"' style='margin:auto;width:100%; font-size: small;'></div>",
     "    <div id='graph"+i+"' style='margin:auto;width:100%; height:200px;'></div>",
     "  </div>",
     "</div>"].join('\n');
    // get data for each graph object
    out[i] = ("position");
    if (display === "phi") {
      out[i] += ",phi"
    } else {
      for (var k = 0; k < checkedSamples.length; k++) { 
      //for (index in checkedSamples) {
        out[i] += ","+samples[checkedSamples[k]]
      }
    }
    out[i] += "\n"
    for (var j = 0; j < contigs[i].data.length; j++) {
      var line = "";
      line += contigs[i].data[j].position[0];
      if (display === "snps") {
        for (var k = 0; k < checkedSamples.length; k++) { 
        //for (index in checkedSamples) {
          line += ","+contigs[i].data[j].SNPs[checkedSamples[k]]
        }
      } else if (display === "depth") {
        for (var k = 0; k < checkedSamples.length; k++) { 
        //for (index in checkedSamples) {
          line += ","+contigs[i].data[j].depth[checkedSamples[k]]
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
  function getLowestPhi(g) {
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
              valueRange: [1,getLowestPhi(g)/2],
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
                    if (getLowestPhi(g) < 0.01) {
                      tickList.push({ v: 0.01 },{ v: 0.01, label: "0.01" })
                    }
                    if (getLowestPhi(g) < 0.001) {
                      tickList.push({ v: 0.001 },{ v: 0.001, label: "0.001" })
                    }
                    if (getLowestPhi(g) < 0.0001) {
                      tickList.push({ v: 0.0001 },{ v: 0.0001, label: "< 0.0001" })
                    } else if (getLowestPhi(g) < 0.1) {
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

  function legendFormatter(data) {
    var legend = '<b>Position: </b>';
    if (data.x >= 0) {
      legend += data.x;
    } else {
      legend += "-"
    }
    var highlightTrigger = 0;
    for (index in data.series) {
      if (data.series[index].labelHTML != "phi") {
        if (data.series[index].yHTML >= 0) {
          document.getElementById("legend-"+data.series[index].labelHTML).innerHTML = data.series[index].yHTML;
        } else {
          document.getElementById("legend-"+data.series[index].labelHTML).innerHTML = 0;
        }
      }
      if (data.series[index].isHighlighted) {
        highlightTrigger = 1;
        legend += '<br><b>' + data.series[index].labelHTML + ': </b>' + data.series[index].yHTML;
      }
    }
    if (highlightTrigger === 0) {
      legend += '<br><b>Sample: </b> -';
    }
    return legend
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
        highlightCallback: function(e, x, pts, row) {
          $('#mouseFollow').css('display','block')
        },
        unhighlightCallback: function(e) {
          $('#mouseFollow').css('display','none')
        },
        ylabel: yLabel,
        labelsShowZeroValues: true,
        labelsSeparateLines: true,
        labelsDiv: 'legend'+i,
        legend: 'always',
        legendFormatter: legendFormatter,
        valueRange: [null, null],
        yRangePad: 5,
        xRangePad: 5,
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
        valueRange: [1,getLowestPhi(g)/2],
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
              if (getLowestPhi(g) < 0.01) {
                tickList.push({ v: 0.01 },{ v: 0.01, label: "0.01" })
              }
              if (getLowestPhi(g) < 0.001) {
                tickList.push({ v: 0.001 },{ v: 0.001, label: "0.001" })
              }
              if (getLowestPhi(g) < 0.0001) {
                tickList.push({ v: 0.0001 },{ v: 0.0001, label: "< 0.0001" })
              } else if (getLowestPhi(g) < 0.1) {
                tickList.push({ v: 0.1 },{ v: 0.1, label: "0.1" })
              }
              return tickList;
            }
          }
        },
        // highlight portions of graph where phi statistic is under 0.05 "statistically significant"
        underlayCallback: function(canvas, area, g) {
          canvas.fillStyle = "rgba(0, 255, 216, 0.3)";
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
    dygraphs[i] = g
  }
  // build all interactive buttons
  
  buildQuantity();
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
  updateOut();
}

// select no samples
function CHECKBOXselectNone() {
  for (var i = 0; i < samples.length; i++) {
    document.getElementById("CHECKBOX-samples"+samples[i]).checked = false;
  }
  updateOut();
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
    document.getElementById("SETTINGS-samples").innerHTML += '<input type="checkbox" id="CHECKBOX-samples'+data.samples[i]+'" onclick="updateOut()" checked> '+data.samples[i]+': <div style="display:inline-block;" id="legend-'+data.samples[i]+'">0</div><br>';
  }
  sideWidth = $("#sidebar-wrapper").width() + 30;
  drawContigs()
}

// main function runs on page load
jQuery(document).ready(function($) {
  readData();
});










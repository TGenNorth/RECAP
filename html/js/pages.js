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
var loadContig = Math.min(10,maxContig-1);
var unloadContig = Math.max(loadContig-10,0);
var sideWidth = 0;
var topHeight = 0;
var currentContig = null;

// variables to store imported JSON data
var samples;
var contigs;
var dygraphs = [];
var gs = [];
var sync;

// for current graph type on display
var display = "SNPs";
var sort = "length";

$(window).resize(function () {
  resizeContig();
  loadHandler();
});

function closeContig() {
  currentContig = null;
  gs = [];
  $("#bottombar-wrapper").addClass("toggle");
  $("#dark-overlay").addClass("toggle");
  
  resizeContig();
}

function resizeContig() {
  topHeight = $("#graph-wrapper").height() + 30;
  $("#bottombar-wrapper").css("height", "auto");
  if ($("#bottombar-wrapper").height() > $("#main-wrapper").height()) {
    $("#bottombar-wrapper").css("height", $("#main-wrapper").height());
  }
  $("#bottombar-wrapper").css("bottom", 0);
  $("#bottombar-wrapper.toggle").css("bottom", -topHeight);
  $("#bottombar-wrapper").css("z-index", 200);
  $("#bottombar-wrapper.toggle").css("z-index", 99);
  $("#dark-overlay").css("opacity", 0.5);
  $("#dark-overlay.toggle").css("opacity", 0.0);
  $("#dark-overlay").css("z-index", 199);
  $("#dark-overlay.toggle").css("z-index", 99);
}

// pulls graph page aside to reveal samples menu
function openContig(i) {
  currentContig = i;
  document.getElementById("CONTIG-single").innerHTML = [
     "<div class='graph'>",
     "  <div class='col-12'>",
     "    <font size='6'>Contig "+parseInt(i+1)+": </font><font size='5'>"+contigs[i].name+"</font>",
     "    <div style='top:10px;right:10px;position:absolute;'>",
     "      <button type='button' id='exp-btn' class='btn btn-info btn-sm' onclick='exciseData()'>Excise</button>",
     "      <button type='button' id='exp-btn' class='btn btn-info btn-sm' onclick='exportData()'>Export</button>",
     "      <button type='button' id='close-btn' class='btn btn-danger btn-sm' onclick='closeContig()'>✕</button>",
     "    </div>",
     "  </div>",
     "  <div class='col-12'>",
     "  <div style='right:20px;position:absolute;'><font size='5'>SNPs</font></div>",
     "    <div id='snplegend"+i+"' style='width:100%; font-size: small;'></div>",
     "    <div id='snpgraph"+i+"' style='width:100%; height:200px;'></div>",
     "  </div>",
     "</div>",
     "<div class='graph'>",
     "  <div class='col-12'>",
     "  <div style='right:20px;position:absolute;padding-top:10px;'><font size='5'>Depth</font></div>",
     "    <div id='depthlegend"+i+"' style='padding-top:10px;width:100%; font-size: small;'></div>",
     "    <div id='depthgraph"+i+"' style='width:100%; height:200px;'></div>",
     "  </div>",
     "</div>",
     "<div class='graph'>",
     "  <div class='col-12'>",
     "  <div style='right:20px;position:absolute;padding-top:10px;'><font size='5'>Phi Statistic</font></div>",
     "    <div id='philegend"+i+"' style='padding-top:10px;width:100%; font-size: small;'></div>",
     "    <div id='phigraph"+i+"' style='width:100%; height:200px;'></div>",
     "  </div>",
     "</div>"].join('\n');
  if (gs.length == 0) {
    gs.push(fillGraphs("snp", i, "SNPs"));
    gs.push(fillGraphs("depth", i, "depth"));
    gs.push(fillGraphs("phi", i, "phi"));
  }
  $("#bottombar-wrapper").removeClass("toggle");
  $("#dark-overlay").removeClass("toggle");
  resizeContig();
  sync = Dygraph.synchronize(gs, {range: false, selection: false})
}

// pulls graph page aside to reveal samples menu
function toggleSamples() {
  $("#sidebar-wrapper").css("width",sideWidth);
  $("#main-wrapper").toggleClass("toggle");
  $("#main-wrapper").css("left", sideWidth);
  $("#main-wrapper.toggle").css("left", 0);
  setTimeout(function(){ window.dispatchEvent(new Event('resize')); }, 300);
}

// ================================================== ( Update Graphs )
// update all graph data based on checkboxes clicked

// get the highest graph value on current page
function getMax() {
  if ((currentPage*maxContig)+maxContig > contigs.length) {
    max = contigs.length;
  } else {
    max = (currentPage*maxContig)+maxContig;
  }
  return max
}

// update data in all loaded graphs
// data will be recalculated before calling this function
function updateOut() {
  for (var i = currentPage*maxContig; i < getMax(); i++) {
    if (dygraphs[i]) {
    dygraphs[i].updateOptions(
      {
        'file':getOut(i, display)
      }
    );
  }}
  if (currentContig != null) {
    gs[2].updateOptions(
      {
        'file':getOut(currentContig, "phi")
      }
    );
    gs[1].updateOptions(
      {
        'file':getOut(currentContig, "depth")
      }
    );
    gs[0].updateOptions(
      {
        'file':getOut(currentContig, "SNPs")
      }
    );
  }
}

// load and unload data based on scroll distance from the top
// dramatically increases speed when many contigs are present
function progressiveLoad() {
  while ((currentPage*maxContig)+loadContig < getMax()-1 && document.getElementById("graph-wrapper").clientHeight > $("#graph"+parseInt((currentPage*maxContig)+loadContig)).offset().top-$("#graph"+parseInt((currentPage*maxContig)+loadContig)).height()*2 && document.getElementById("graph-wrapper").clientHeight >= $("#graph"+parseInt((currentPage*maxContig)+loadContig)).offset().top-$("#graph"+(currentPage*maxContig)+parseInt(loadContig)).height()*3) {
    loadContig = Math.min(loadContig+1, maxContig)
    unloadContig = Math.max(unloadContig, loadContig-10)
  }
  while ((currentPage*maxContig)+unloadContig > (currentPage*maxContig) && document.getElementById("graph-wrapper").clientHeight < $("#graph"+parseInt((currentPage*maxContig)+loadContig)).offset().top-$("#graph"+(currentPage*maxContig)+parseInt(loadContig)).height()*3 && document.getElementById("graph-wrapper").clientHeight <= $("#graph"+parseInt((currentPage*maxContig)+loadContig)).offset().top-$("#graph"+parseInt((currentPage*maxContig)+loadContig)).height()*2) {
    unloadContig = Math.max(unloadContig-1, 0)
    loadContig = Math.min(loadContig, unloadContig+10)
  }
  loadHandler()
}

// load all graphs between "unloadContig" and "loadContig" variables
// unload all else if graph is loaded
function loadHandler() {
  for (var i = (currentPage*maxContig); i < getMax(); i++) {
    if (i >= (currentPage*maxContig)+unloadContig && i <= (currentPage*maxContig)+loadContig) {
      $("#legend"+i).css("width",parseInt($("#graph-wrapper").width()-30))
      $("#graph"+i).css("width",parseInt($("#graph-wrapper").width()-30))
      dygraphs[i] = fillGraphs("",i, display);
    }
    else if (dygraphs[i]) {
      document.getElementById("graph"+i).innerHTML = "";
      dygraphs[i] = null
    }
  }
}

// ================================================== ( Data Export )
// create download element
function download(filename, text) {
  var element = document.createElement('a');
  document.body.appendChild(element);
  element.style.display = 'none';

  var blob = new Blob([text], {encoding:"UTF-8",type:"text/plain;charset=UTF-8"}),
      url = window.URL.createObjectURL(blob);
  element.href = url;
  element.download = filename;
  element.click();
  window.URL.revokeObjectURL(url);

  document.body.removeChild(element);
}

function exciseData() {
  var file = document.createElement("script");
  file.src = "excise.jsonp?callback=exciseContig";
  document.body.insertBefore(file, document.body.firstChild);
}

function exportData() {
  var file = document.createElement("script");
  file.src = "export.jsonp?callback=exportContig";
  document.body.insertBefore(file, document.body.firstChild);
}

function exportContig(data) {
  if (currentContig != null) {
    var checkedSamples = CHECKBOXreturnChecked();
    for (var i = 0; i < checkedSamples.length; i++) {
      checkedSamples[i] = checkedSamples[i] + 1;
    }
    checkedSamples.unshift(0)
    range = gs[2].xAxisRange();
    left = Math.max(Math.floor(range[0]),0);
    right = Math.min(Math.ceil(range[1]),contigs[currentContig].data.slice(-1)[0].position[1]);
    var out = [];
    var keys = Object.keys(data[contigs[currentContig].name]);
    for (var i = 0; i < checkedSamples.length; i++) {
      out.push(">" + keys[checkedSamples[i]] + "\n");
      var tempLst = []
      var tempStr = data[contigs[currentContig].name][keys[i]].substring(left,right)
      for (var j = 0; j < right-left; j += 80) {
        tempLst.push(tempStr.substring(j,j+80))
      }
      out[i] += tempLst.join("\n")
    }
    download("(export)_"+contigs[currentContig].name+"_("+left+"-"+right+")"+".txt", out.join("\n\n"));
  }
}

function exciseContig(data) {
  if (currentContig != null) {
    range = gs[2].xAxisRange();
    left = Math.max(Math.floor(range[0]),0);
    right = Math.min(Math.ceil(range[1]),contigs[currentContig].data.slice(-1)[0].position[1]);

    var keys = Object.keys(data);
    out = [data[keys[0]]]
    for (var i = 1; i < keys.length; i++) {
      var positions = Object.keys(data[contigs[i-1].name])
      for (var j = 0; j < positions.length; j++) {
        if (currentContig != i-1) {
          out.push(data[keys[i]][positions[j]])
        } else if (positions[j] < left || positions[j] > right) {
          out.push(data[keys[i]][positions[j]])
        }
      }
    }
    download("(excise)_"+contigs[currentContig].name+"_(0-"+left+")("+right+"-"+contigs[currentContig].data.slice(-1)[0].position[1]+")"+".txt", out.join('\n'));
  }
}

// ================================================== ( Draw Contigs )

function getOut(i, display) {
  var out;
  var checkedSamples = CHECKBOXreturnChecked();
    out = ("position");
    if (display === "phi") {
      out += ",phi"
    } else {
      for (var k = 0; k < checkedSamples.length; k++) { 
      //for (index in checkedSamples) {
        out += ","+samples[checkedSamples[k]]
      }
    }
    out += "\n"
    for (var j = 0; j < contigs[i].data.length; j++) {
      var line = "";
      line += contigs[i].data[j].position[0];
      if (display === "SNPs") {
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
      out += line + "\n"
      if (contigs[i].data.length === 1) {
        out += contigs[i].data[j].position[1]+line.substr(1, line.length) + "\n"
      }
    }

  return out;
}

// create and build all graph objects
function drawGraphs() {
  // only add data for samples that are checked off
  var checkedSamples = CHECKBOXreturnChecked();
  
  // empty and repopulate all graph objects
  document.getElementById("graphs").innerHTML = ''

  for (var i = currentPage*maxContig; i < getMax(); i++) {
    document.getElementById("graphs").innerHTML +=[
     "<div class='graph'>",
     "  <div class='col-12'>",
     "    <font size='5'>Contig "+parseInt(i+1)+":</font><font size='4'> <span id='link' onclick='openContig("+i+")'>"+contigs[i].name+"</span></font>",
     "  </div>",
     "  <div class='col-12'>",
     "    <div id='legend"+i+"' style='margin:auto;width:"+parseInt($("#graph-wrapper").width()-30)+"px; font-size: small;'></div>",
     "    <div id='graph"+i+"' style='margin:auto;width:"+parseInt($("#graph-wrapper").width()-30)+"px; height:200px;'></div>",
     "  </div>",
     "</div>"].join('\n');
  }

  buildQuantity();
  buildDisplay();
  buildSort();
  buildNavigation();
  buildBottom();

  loadContig = Math.min(10,getMax()-currentPage*maxContig-1);
  unloadContig = Math.max(loadContig-10,0);
  progressiveLoad();

  document.getElementById("graphSort").innerHTML = "Showing " + maxContig + " contigs: " + display + ", sorted by " + sort;
  $('.loading').css('display', 'none');
}

function fillGraphs(type, i, display) {
  // fill each graph with saved data
  var yLabel;
  if (display === "SNPs") {
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
    for (var j = 0; j < g.numRows(); j++) {
      if (lowestPhi > g.getValue(j,1)) {
        lowestPhi = g.getValue(j,1);
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


    var g = new Dygraph(
      document.getElementById(type+"graph"+i),
      getOut(i, display),
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
        labelsDiv: type+'legend'+i,
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
          for (var j = 0; j < g.numRows()-1; j++) {
            if (g.getValue(j,1) >= 0.05 && g.getValue(j+1,1) < 0.05) {
              var start = g.getValue(j,0);
              //console.log("slopeU:",Math.log10(g.getValue(i+1,1))-Math.log10(g.getValue(i,1)))
              //console.log("y    U:",g.getValue(i+1,1)*Math.pow((10),Math.log10(g.getValue(i,1))-Math.log10(g.getValue(i+1,1))))
              //console.log("phi  U:",g.getValue(i,1))
              //console.log("X1   U:",(Math.log10(0.05/g.getValue(i+1,1)))/(Math.log10(g.getValue(i+1,1))-Math.log10(g.getValue(i,1))))
              //console.log("X2   U:",(Math.log10(g.getValue(i+1,1)/g.getValue(i+1,1)))/(Math.log10(g.getValue(i+1,1))-Math.log10(g.getValue(i,1))))
              highlight(g.getValue(j+1,0)+((Math.log10(0.05/g.getValue(j+1,1)))/(Math.log10(g.getValue(j+1,1))-Math.log10(g.getValue(j,1)))*(g.getValue(j+1,0)-g.getValue(j,0))),g.getValue(j+1,0))
            } else if (g.getValue(j,1) < 0.05 && g.getValue(j+1,1) >= 0.05) {
              var start = g.getValue(j,0);
              //console.log("slopeD:",Math.log10(g.getValue(i,1))-Math.log10(g.getValue(i+1,1)))
              //console.log("y    D:",g.getValue(i,1)*Math.pow((10),Math.log10(g.getValue(i+1,1))-Math.log10(g.getValue(i,1))))
              //console.log("phi  D:",g.getValue(i+1,1))
              highlight(g.getValue(j+1,0)+((Math.log10(0.05/g.getValue(j+1,1)))/(Math.log10(g.getValue(j,1))-Math.log10(g.getValue(j+1,1)))*(g.getValue(j,0)-g.getValue(j+1,0))),g.getValue(j,0))
            } else if (g.getValue(j,1) < 0.05 && g.getValue(j+1,1) < 0.05) {
              highlight(g.getValue(j,0),g.getValue(j+1,0))
            }
          }
        }
      });
    }
    return g;
  // build all interactive buttons
  
}

// ================================================== ( Build Buttons ) 
// drop down buton to change amount of contigs per page
function buildQuantity() {
  var buttonData = [
   "<div class='dropdown' style='margin-left:10px;margin-right:10px;'>",
   "  <button class='btn btn-info btn-sm dropdown-toggle' type='button' id='dropdownMenuButton' data-toggle='dropdown' aria-haspopup='true' aria-expanded='false'>Contigs Per Page</button>",
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
   "  <button class='btn btn-info btn-sm dropdown-toggle' type='button' id='dropdownMenuButton' data-toggle='dropdown' aria-haspopup='true' aria-expanded='false'>Graph Type</button>",
   "  <div class='dropdown-menu' aria-labelledby='dropdownMenuButton'>",
   "    <span class='dropdown-item' onclick='CONTIGdisplay(\"SNPs\")'>SNPs</span>",
   "    <span class='dropdown-item' onclick='CONTIGdisplay(\"depth\")'>Depth</span>",
   "    <span class='dropdown-item' onclick='CONTIGdisplay(\"phi\")'>PHI</span>",
   "  </div>",
   "</div>"].join('\n');
}

// drop down button to change contig sort order
function buildSort() {
  document.getElementById("CONTIG-options").innerHTML += [
   "<div class='dropdown' style='margin-left:10px;margin-right:10px;'>",
   "  <button class='btn btn-info btn-sm dropdown-toggle' type='button' id='dropdownMenuButton' data-toggle='dropdown' aria-haspopup='true' aria-expanded='false'>Sort</button>",
   "  <div class='dropdown-menu' aria-labelledby='dropdownMenuButton'>",
   "    <span class='dropdown-item' onclick='CONTIGsort(\"name\")'>Name</span>",
   "    <span class='dropdown-item' onclick='CONTIGsort(\"length\")'>Length</span>",
   "    <span class='dropdown-item' onclick='CONTIGsort(\"SNPs\")'>SNPs</span>",
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

function updateStatus() {
  $('.loading').css('display', 'inline-block');
  setTimeout(function wait() {
    drawGraphs();
  }, 100)
}

function CONTIGnavigation(page) {
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
  updateStatus()
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
  sort = option;
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
  } else if (option === "SNPs") {
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
  var checkedSamples = [];
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
  if (maxContig > data.contigs.length) {
    maxContig = data.contigs.length
  }
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
  updateStatus();
  sideWidth = $("#sidebar-wrapper").width() + 30;
  topHeight = $("#graph-wrapper").height() + 30;
}

// main function runs on page load
jQuery(document).ready(function($) {
  closeContig();
  readData();
});










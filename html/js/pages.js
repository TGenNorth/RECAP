jQuery(document).ready(function($) {
  window.history.replaceState(null, null, '#1')
  var speed = 150;
  $('.page.current').css('opacity',0).css('marginTop', 100).animate({opacity: 1, marginTop: 30}, speed);
  readData();
  if (window.history && window.history.replaceState) {
    $(window).on('popstate', function() {
      var url = window.location.hash;

      $(".collapsed").click(function() {
      });
      $(".active").removeClass("active");
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
        console.log(window.location)
        var i=1;
        $(".button1").addClass("active");
        window.history.replaceState(null, null, '#1')
      }
      if($(".navbar-toggle").attr("class") == "navbar-toggle"){
        $(".navbar-toggle").addClass("collapsed");
        $(".navbar-toggle").attr("aria-expanded", "false");
        $("#navbar").removeClass("in");
      }
      $('.page.current').animate({opacity: 0, marginTop: 40}, 0, function(){
        $(this).removeClass('current');
        $('.page').eq(i).css('marginTop', 40).animate({opacity: 1, marginTop: 30}, speed).addClass('current');
      });
    });
  }
});

function buildGraphs(data) {
  for (var i = 1; i < data.length; i++) {
    // console.log(data[i].split(","))
  }
}

function readData() {
  var file = document.createElement("script");
  file.src = "snpDensityMatrix.jsonp?callback=buildGraphs";
  document.body.appendChild(file);
}



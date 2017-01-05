//var opts = {
//  lines: 9, // The number of lines to draw
//  length: 9, // The length of each line
//  width: 5, // The line thickness
//  radius: 14, // The radius of the inner circle
//  color: '#EEEEEE', // #rgb or #rrggbb or array of colors
//  speed: 1.9, // Rounds per second
//  trail: 40, // Afterglow percentage
//  className: 'spinner', // The CSS class to assign to the spinner
//};
//var spinner = new Spinner(opts).spin(document.getElementById('#dbs_match_sunburst'));

var width = 250,
    height = 250,
    radius = (Math.min(width, height) / 2) - 10;

var formatNumber = d3.format(",d");

var x2 = d3.scale.linear()
    .range([0, 2 * Math.PI]);

var y2 = d3.scale.sqrt()
    .range([0, radius]);

var color3 = {
  'undefined':"#999966",
  'Barcode found':"#879942",
  'Barcode match pattern':"#669900",
  'Barcode not found':"#EB4714",
  'Barcode don\'t match pattern':"#FF3300"
}

var partition2 = d3.layout.partition()
    .value(function(d) { return d.size; });
    
var tmpText2 = d3.select("#dbs_match_info").text();

var arc2 = d3.svg.arc()
    .startAngle(function(d) { return Math.max(0, Math.min(2 * Math.PI, x2(d.x))); })
    .endAngle(function(d) { return Math.max(0, Math.min(2 * Math.PI, x2(d.x + d.dx))); })
    .innerRadius(function(d) { return Math.max(0, y2(d.y)); })
    .outerRadius(function(d) { return Math.max(0, y2(d.y + d.dy)); });

var svg2 = d3.select("#dbs_match_sunburst").append("svg")
    .attr("width", width)
    .attr("height", height)
  .append("g")
    .attr("transform", "translate(" + width / 2 + "," + (height / 2) + ")");

d3.json("dbs_match.json", function(error, root) {
  
  //spinner.stop();
  
  if (error) throw error;

  svg2.selectAll("path")
      .data(partition2.nodes(root))
    .enter().append("path")
      .attr("d", arc2)
      //.style("fill", function(d) { return color((d.children ? d : d.parent).name); })
      .style("fill", function(d) { var res = d.name.split("% "); return color3[res[1]]; })
      .on("click", click2)
      .on('mouseover', function changeLable(d) { d3.select("#dbs_match_info").text( d.name + " (" + formatNumber(d.value) + ")" );})
      .on("mouseout", function () {d3.select("#dbs_match_info").text(tmpText2);})
    .append("title")
      .text(function(d) { return d.name + " (" + formatNumber(d.value) + ")"; });
});

function click2(d) {
  svg2.transition()
      .duration(750)
      .tween("scale", function() {
        var xd = d3.interpolate(x2.domain(), [d.x, d.x + d.dx]),
            yd = d3.interpolate(y2.domain(), [d.y, 1]),
            yr = d3.interpolate(y2.range(), [d.y ? 20 : 0, radius]);
        return function(t) { x2.domain(xd(t)); y2.domain(yd(t)).range(yr(t)); };
      })
    .selectAll("path")
      .attrTween("d", function(d) { return function() { return arc2(d); }; });
}

d3.select(self.frameElement).style("height", height + "px");
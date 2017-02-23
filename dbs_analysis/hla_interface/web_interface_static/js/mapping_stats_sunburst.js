var opts = {
  lines: 9, // The number of lines to draw
  length: 9, // The length of each line
  width: 5, // The line thickness
  radius: 14, // The radius of the inner circle
  color: '#999999', // #rgb or #rrggbb or array of colors
  speed: 1.9, // Rounds per second
  trail: 40, // Afterglow percentage
  className: 'spinner', // The CSS class to assign to the spinner
};
var target = document.getElementById('mapping_stats_sunburst')
var mapping_stats_spinner = new Spinner(opts).spin(target);

var width = 250,
    height = 250,
    radius = (Math.min(width, height) / 2) - 10;

var formatNumber = d3.format(",d");

var mapping_stats_x = d3.scale.linear()
    .range([0, 2 * Math.PI]);

var mapping_stats_y = d3.scale.sqrt()
    .range([0, radius]);

var mapping_stats_color = {
  'undefined':"#EEEEEE",
  //'in bam':"#879942",
  //'mapped concordantly':"#669900",
  //'not in bam':"#EB4714",
  //'unMapped':"#FF3300"
 
'not in bam':"#D65C29",
'in bam':"#879942",
'unMapped':"#EB4714",
'Mapped':"#789924",
'Mapped as Proper Pair':"#009900",
'multi mapped pair':"#669900",
'mapped discordantly':"#669900",
'mapped single end':"#669900",
'multi mapped single end':"#669900" 
  
}

var mapping_stats_partition = d3.layout.partition()
    .value(function(d) { return d.size; });
    
var mapping_stats_tmpText = d3.select("#mapping_info").text();

var mapping_stats_arc = d3.svg.arc()
    .startAngle(function(d) { return Math.max(0, Math.min(2 * Math.PI, mapping_stats_x(d.x))); })
    .endAngle(function(d) { return Math.max(0, Math.min(2 * Math.PI, mapping_stats_x(d.x + d.dx))); })
    .innerRadius(function(d) { return Math.max(0, mapping_stats_y(d.y)); })
    .outerRadius(function(d) { return Math.max(0, mapping_stats_y(d.y + d.dy)); });

var mapping_stats_svg = d3.select("#mapping_stats_sunburst").append("svg")
    .attr("width", width)
    .attr("height", height)
  .append("g")
    .attr("transform", "translate(" + width / 2 + "," + (height / 2) + ")");

d3.json("mapping_stats.json", function(error, root) {
  
  mapping_stats_spinner.stop();
  
  if (error) throw error;

  mapping_stats_svg.selectAll("path")
      .data(mapping_stats_partition.nodes(root))
    .enter().append("path")
      .attr("d", mapping_stats_arc)
      //.style("fill", function(d) { return color((d.children ? d : d.parent).name); })
      .style("fill", function(d) { var res = d.name.split("% "); return mapping_stats_color[res[1]]; })
      .on("click", mapping_stats_click)
      .on('mouseover', function changeLable(d) { d3.select("#mapping_info").text( d.name + " (" + formatNumber(d.value) + ")" );})
      .on("mouseout", function () {d3.select("#mapping_info").text(mapping_stats_tmpText);})
    .append("title")
      .text(function(d) { return d.name + " (" + formatNumber(d.value) + ")"; });
});

function mapping_stats_click(d) {
  mapping_stats_svg.transition()
      .duration(750)
      .tween("scale", function() {
        var xd = d3.interpolate(mapping_stats_x.domain(), [d.x, d.x + d.dx]),
            yd = d3.interpolate(mapping_stats_y.domain(), [d.y, 1]),
            yr = d3.interpolate(mapping_stats_y.range(), [d.y ? 20 : 0, radius]);
        return function(t) { mapping_stats_x.domain(xd(t)); y2.domain(yd(t)).range(yr(t)); };
      })
    .selectAll("path")
      .attrTween("d", function(d) { return function() { return mapping_stats_arc(d); }; });
}

d3.select(self.frameElement).style("height", height + "px");
d3.csv("readspercluster.csv", function (data) {
    
    var width = 600,
    height = 300,
    padding = 50,
    move_left = 60,
    move_down = 5;

    // when the input range changes update the angle 
    d3.select("#nAngle").on("input", function() {
      update(+this.value);
    });
    
    var map = data.map( function (i) { return parseInt(i.reads); } )
        //.filter(function(i) { return i > nAngle })

    var bin_width = d3.max(map)/100
    var bins = []
    for (var i = 0; i < d3.max(map)+1; i=i+bin_width) {
        bins.push(i);
    }
    
    var histogram = d3.layout.histogram()
        .bins(bins)
        (map)
        
    var y = d3.scale.linear()
        .domain([0, d3.max(histogram.map( function (i) { return i.length; } ) )])
        .range([height-move_down, 0]);

    var x = d3.scale.linear()
        .domain([0, d3.max(map)])
        .range([0, width-move_left]);
    d3.select("#nAngle").property("max", d3.max(map));
    d3.select("#nAngle").property("style", "width:" + ( width - move_left ) + "px");
    d3.select("#nReadsMin").property("style","display: inline-block; width: " + ( move_left ) + "px; text-align: right");
        
    var xAxis = d3.svg.axis()
        .scale(x)
        .orient("bottom");

    var yAxis = d3.svg.axis()
        .scale( y )
        .orient("left");

    var tip = d3.tip()
      .attr('class', 'd3-tip')
      .offset([-10, 0])
      .html(function(d) {
        return d.y + " clusters has " + Math.ceil(Math.round(d.x*10)/10) + " to " + Math.floor(Math.round((d.x + d.dx)*10)/10) + " reads";
      })

    var canvas = d3.select("#image").append("svg")
        .attr("width",width)
        .attr("height",height + padding )
        .append("g")
            .attr( "transform", "translate( " + move_left + " ,  " + move_down + "  )" ) ;
    
    canvas.call(tip)
        
    var group = canvas.append("g")
        .attr("transform","translate(0," + (height-move_down) + ")" )
        .style({ 'stroke': 'Black', 'fill': 'none', 'stroke-width': '1px'})
        .call(xAxis)
    .selectAll("text")
        .attr("y", 0)
        .attr("x", 9)
        .attr("dy", ".35em")
        .attr("transform", "rotate(90)")
        .style("text-anchor", "start");

    var group2 = canvas.append("g")
        .style({ 'stroke': 'Black', 'fill': 'none', 'stroke-width': '1px'})
        .call(yAxis);
        
    var bars = canvas.selectAll("rect")
        .data(histogram)
        .enter()
        .append("svg:rect")
    
    //bars.append("rect")
        .attr("x", function (d) { return x(d.x) } )
        .attr("y", function (d) { return y(d.y) } )
        .attr("width", function (d) { return x(d.dx) } )
        .attr("height", function (d) { return y(0) - y(d.y) } )
        //.attr("fill", "steelblue")
        .on('mouseover', tip.show)
        .on('mouseout', tip.hide);
        
    // Initial starting angle of the text 
    update(0);
    
    // update the element
    function update(nAngle) {
    
      // adjust the text on the range slider
      d3.select("#nReadsMin-value").text(nAngle);
      d3.select("#nAngle").property("value", nAngle);
    
      // rotate the text
      //holder.select("text") 
      //  .attr("transform", "translate(300,150) rotate("+nAngle+")");
      redraw(nAngle)
    }
    
    function redraw(nAngle) {
        // Update
        
        var map = data.map( function (i) { return parseInt(i.reads); } )
            .filter(function(i) { return i >= nAngle })

        var histogram = d3.layout.histogram()
            .bins(bins)
            (map)
            
        y.domain([0, d3.max(histogram.map( function (i) { return i.length; } ) )])

        x.domain([0, d3.max(map)])
            
        console.log(map)
        
        group2.transition()
            .duration(1000)
            .call(yAxis);
        
        canvas.selectAll("rect")
                .data(histogram)
            .transition()
                .duration(1000)
                .attr("x", function (d) { return x(d.x) } )
                .attr("y", function (d) { return y(d.y) } )
                .attr("width", function (d) { return x(d.dx) } )
                .attr("height", function (d) { return y(0) - y(d.y) } )
                //.attr("fill", "steelblue");
    }
} )

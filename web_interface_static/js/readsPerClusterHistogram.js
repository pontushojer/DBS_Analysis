// makes the reads per cluster graph

d3.csv("readspercluster.csv", function (data) {
    
    // initiate variables defining the size of the image
    var width = 600,
    height = 300,
    padding = 50,
    move_left = 60,
    move_down = 5;

    // when the cutoff used change update the image 
    d3.select("#nRPCcutoff").on("input", function() {
      update(+this.value);
    });
    
    // get all the numbers as integers
    var map = data.map( function (i) { return parseInt(i.reads); } )

    // based on the values found define the bins so that there always are 100 bins
    var bin_width = d3.max(map)/100
    var bins = []
    for (var i = 0; i < d3.max(map)+1; i=i+bin_width) {
        bins.push(i);
    }
    
    // get the histogram
    var histogram = d3.layout.histogram()
        .bins(bins)
        (map)
    
    // create the scales for the x and y axes
    var y = d3.scale.linear()
        .domain([0, d3.max(histogram.map( function (i) { return i.length; } ) )])
        .range([height-move_down, 0]);

    var x = d3.scale.linear()
        .domain([0, d3.max(map)])
        .range([0, width-move_left]);
    d3.select("#nRPCcutoff").property("max", d3.max(map));
    d3.select("#nRPCcutoff").property("style", "width:" + ( width - move_left ) + "px");
    d3.select("#nReadsMin").property("style","display: inline-block; width: " + ( move_left ) + "px; text-align: right");
        
    var xAxis = d3.svg.axis()
        .scale(x)
        .orient("bottom");

    var yAxis = d3.svg.axis()
        .scale( y )
        .orient("left");

    // define how the mouse over tip should look/appear
    var tip = d3.tip()
      .attr('class', 'd3-tip')
      .offset([-10, 0])
      .html(function(d) {
        return d.y + " clusters has " + Math.ceil(Math.round(d.x*10)/10) + " to " + Math.floor(Math.round((d.x + d.dx)*10)/10) + " reads";
      })

    // define the image canvas to paint on
    var canvas = d3.select("#image").append("svg")
        .attr("width",width)
        .attr("height",height + padding )
        .append("g")
            .attr( "transform", "translate( " + move_left + " ,  " + move_down + "  )" ) ;
    
    // add the tip "layout"
    canvas.call(tip)
        
    // add x axis text
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

    // add y axis text
    var group2 = canvas.append("g")
        .style({ 'stroke': 'Black', 'fill': 'none', 'stroke-width': '1px'})
        .call(yAxis);
        
    // add the bars for the bins in the histogram
    var bars = canvas.selectAll("rect")
        .data(histogram)
        .enter()
        .append("svg:rect")
        .attr("x", function (d) { return x(d.x) } )
        .attr("y", function (d) { return y(d.y) } )
        .attr("width", function (d) { return x(d.dx) } )
        .attr("height", function (d) { return y(0) - y(d.y) } )
        .on('mouseover', tip.show)
        .on('mouseout', tip.hide);
        
    // set the initial cutof value to zero reads per cluster
    update(0);
    
    // the update the funtion to call when cutof changes
    function update(nRPCcutoff) {
    
        // adjust the text on the range slider
        d3.select("#nReadsMin-value").text(nRPCcutoff);
        d3.select("#nnRPCcutoff").property("value", nRPCcutoff);
      
        // redraw the image with the new cutof
        redraw(nRPCcutoff)
    }
    
    // the redraw function that actually redraws the image
    function redraw(nRPCcutoff) {
        
        // filter all the integers based on the defined cutof
        var map = data.map( function (i) { return parseInt(i.reads); } )
            .filter(function(i) { return i >= nRPCcutoff })

        // get the histogram data
        var histogram = d3.layout.histogram()
            .bins(bins)
            (map)
            
        // update the axes max values
        y.domain([0, d3.max(histogram.map( function (i) { return i.length; } ) )])
        x.domain([0, d3.max(map)]) // shouldn't change currently only have a lower cutoff

        // make the y axis transition
        group2.transition()
            .duration(1000)
            .call(yAxis);
        
        // repaint the bars with the new histogram data
        canvas.selectAll("rect")
                .data(histogram)
            .transition()
                .duration(1000)
                .attr("x", function (d) { return x(d.x) } )
                .attr("y", function (d) { return y(d.y) } )
                .attr("width", function (d) { return x(d.dx) } )
                .attr("height", function (d) { return y(0) - y(d.y) } )
    }
} )

// makes the reads per cluster graph

d3.csv("first1k.readspercluster.csv", function (data) {
    
    data.forEach(function(d) {
        d.reads = parseInt(d.reads);
        d.id = parseInt(d.id);
    });
    
    // initiate variables defining the size of the image
    var width = 600,
    height = 300,
    padding = 50,
    move_left = 60,
    move_down = 5;
    
    // get all the numbers as integers
    var reads_map = data.map( function (i) { return i.reads; } )
    var id_map = data.map( function (i) { return i.id; } )

    // create the scales for the x and y axes
    var y = d3.scale.linear()
        .domain([0, d3.max(reads_map)])
        .range([height-move_down, 0]);

    var x = d3.scale.linear()
        .domain([0, d3.max(id_map)])
        .range([0, width-move_left]);
        
    var xAxis = d3.svg.axis()
        .scale(x)
        .orient("bottom");

    var yAxis = d3.svg.axis()
        .scale( y )
        .orient("left");
        
    // define the image canvas to paint on
    var canvas = d3.select("#image").append("svg")
        .attr("width",width)
        .attr("height",height + padding )
        .append("g")
            .attr( "transform", "translate( " + move_left + " ,  " + move_down + "  )" ) ;
            
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

    // Define the line
    var valueline = d3.svg.line()
        .x(function(d) { return x(d.id); })
        .y(function(d) { return y(d.reads); });
    
    // Add the valueline path.
    canvas.append("path")
        .attr("id","cluster_line")
        .attr("class", "line")
        .attr("d", valueline(data));
    
} )

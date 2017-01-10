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
var target = document.getElementById('graph_container_1')
var spinner_1 = new Spinner(opts).spin(target);

function add(a, b) {
  return a + b;
}

function cumulative_array( array ) {
  return array.map( function(currentValue, index, arr) {
    return arr.slice(index).reduce(add, 0);
  })
}

function brushFilter(in_data,extent) {
    if ( extent[0] !=  extent[1]  ) {
      filtered_data = in_data.filter( function(i) { return (i.total || 0) >= extent[0] && (i.total || 0) <= extent[1]+1})
      return filtered_data
      }
    return in_data
}

d3.json("barcode_clusters.json", function(error, data) {
  
  spinner_1.stop();
  
  if (error) throw error;
  
  
  function initiate_plots( data ) {

    // initiate variables defining the size of the images
    var width = 750,
    height = 300,
    padding = 50,
    move_left = 80,
    move_down = 5;  
  
// HERE

    // get all the numbers as integers
    var hist_map = data.map( function (i) { return i.total; } )

    // based on the values found define the bins so that there always are 100 bins
    var bin_width = d3.max(hist_map)/100
    var bins = []
    for (var i = 0; i < d3.max(hist_map)+1; i=i+bin_width) {
        bins.push(i);
    }
    
    // get the histogram
    var histogram = d3.layout.histogram()
        .bins(bins)
        (hist_map)
    
    // create the scales for the x and y axes
    var hist_y = d3.scale.linear()
        .domain([0, d3.max(histogram.map( function (i) { return i.length; } ) )])
        .range([height-move_down, 0]);

    var hist_x = d3.scale.linear()
        .domain([0, d3.max(hist_map)])
        .range([0, width-move_left]);
        
    var hist_xAxis = d3.svg.axis()
        .scale(hist_x)
        .orient("bottom")
        .ticks(10, ",.0s")
        .tickSize(6, 0);

    var hist_yAxis = d3.svg.axis()
        .scale( hist_y )
        .orient("left")
        .ticks(10, ",.0s")
        .tickSize(6, 0);

    //// define how the mouse over tip should look/appear
    //var tip = d3.tip()
    //  .attr('class', 'd3-tip')
    //  .offset([-10, 0])
    //  .html(function(d) {
    //    return d.y + " clusters has " + Math.ceil(Math.round(d.x*10)/10) + " to " + Math.floor(Math.round((d.x + d.dx)*10)/10) + " reads";
    //  })

    // define the image canvas to paint on
    var hist_canvas = d3.select("#graph_container_2").append("svg")
        .attr("width",width)
        .attr("height",height + padding )
        .append("g")
            .attr( "transform", "translate( " + move_left + " ,  " + move_down + "  )" ) ;
    
    // add the tip "layout"
    //hist_canvas.call(tip)
        
    // add x axis text
    var hist_group = hist_canvas.append("g")
        .attr("transform","translate(0," + (height-move_down) + ")" )
        .attr('stroke', 'Black')
        //.style({ 'stroke': 'Black', 'fill': 'none', 'stroke-width': '1px'})
        .call(hist_xAxis)
    
    //hist_group.selectAll("text")
    //    .attr("y", 0)
    //    .attr("x", 9)
    //    .attr("dy", ".35em")
    //    .attr("transform", "rotate(90)")
    //    .style("text-anchor", "start");

    // add y axis text
    var hist_group2 = hist_canvas.append("g")
        .attr('stroke', 'Black')
        //.style({ 'stroke': 'Black', 'fill': 'none', 'stroke-width': '1px'})
        .call(hist_yAxis);
        
    // add the bars for the bins in the histogram
    var bars = hist_canvas.selectAll("rect")
        .data(histogram)
        .enter()
        .append("svg:rect")
        .attr("x", function (d) { return hist_x(d.x) } )
        .attr("y", function (d) { return hist_y(d.y) } )
        .attr("width", function (d) { return hist_x(d.dx) } )
        .attr("height", function (d) { return hist_y(0) - hist_y(d.y) } )
        .attr("fill", "#EAEAEA")
        .attr("stroke", "#000000");
        //.on('mouseover', tip.show)
        //.on('mouseout', tip.hide);

    function redraw(extent) {

        filtered_data = brushFilter(data,extent)

        // filter all the integers based on the defined cutof
        var hist_map_filt = filtered_data.map( function (i) { return i.total; } )

        //var bin_width_filt = (d3.max(hist_map_filt)-d3.min(hist_map_filt)) / 100
        //var bins_filt = []
        //for (var i = d3.min(hist_map_filt); i < d3.max(hist_map_filt)+1; i=i+bin_width_filt) {
        //    bins_filt.push(i);
        //}
        
        // get the histogram data
        var histogram_filt = d3.layout.histogram()
            .bins(bins)
            (hist_map_filt)
            
        // update the axes max values
        hist_y.domain([0, d3.max(histogram_filt.map( function (i) { return i.length; } ) )])
        hist_x.domain([0, d3.max(hist_map)]) // shouldn't change currently only have a lower cutoff
        //hist_x.domain([d3.min(hist_map_filt), d3.max(hist_map_filt)]) // shouldn't change currently only have a lower cutoff

        //hist_group.transition()
        //    .duration(1000)
        //    .call(hist_xAxis);

        // make the y axis transition
        hist_group2.transition()
            .duration(1000)
            .call(hist_yAxis);
        
        // repaint the bars with the new histogram data
        hist_canvas.selectAll("rect")
                .data(histogram_filt)
            .transition()
                .duration(1000)
                .attr("x", function (d) { return hist_x(d.x) } )
                .attr("y", function (d) { return hist_y(d.y) } )
                .attr("width", function (d) { return hist_x(d.dx) } )
                .attr("height", function (d) { return hist_y(0) - hist_y(d.y) } )
    }
//TO HERE
  
    var total_read_count = data.map( function (i) { return i.total; } ).reduce(add, 0)
    d3.select("#total_read_count").text( total_read_count);
    var total_cluster_count = data.length
    var read_counts_per_cluster = data.map( function (i) { return i.total; } );
    var y_data_1 = cumulative_array(read_counts_per_cluster).reverse()
    var y_data_2 = cumulative_array(read_counts_per_cluster.reverse())
    var x_data_1 = read_counts_per_cluster
    d3.select("#p_low").text( 2);
    d3.select("#p_high").text(d3.max(x_data_1));
    d3.select("#p_of_reads").text(100);
    d3.select("#p_num_clusters").text(data.length);
    d3.select("#p_of_clusters").text(100);
    d3.select("#p_average").text(Math.round(total_read_count/total_cluster_count));
    d3.select("#p_median").text(data[parseInt(Math.round(data.length/2))].total);

    //console.log( y_data_1 )
    //console.log( y_data_2 )
    //console.log( x_data_1 )

    // create the scales for the x and y axes
    var lines_y = d3.scale.linear()
        .domain([0, d3.max(y_data_1)])
        .range([height-move_down, 0]);

    var lines_x = d3.scale.log()
        .base(10)
        .domain([d3.min(x_data_1), d3.max(x_data_1)*2])
        .range([0, width-move_left]);

    var lines_xAxis = d3.svg.axis()
        .scale(lines_x)
        .orient("bottom")
        .ticks(10, ",.1s")
        .tickSize(6, 0);

    var lines_yAxis = d3.svg.axis()
        .scale( lines_y )
        .orient("left")
        .ticks(10, ",.2s")
        .tickSize(6, 0);

    // define the image canvas to paint on
    var lines_canvas = d3.select("#graph_container_1").append("svg")
        .attr("width",width)
        .attr("height",height + padding )
        .append("g")
            .attr( "transform", "translate( " + move_left + " ,  " + move_down + "  )" ) ;
            
    // add x axis text
    var lines_x_group = lines_canvas.append("g")
        .attr("transform","translate(0," + (height-move_down) + ")" )
        .attr('stroke', 'Black')
        //.style({ 'stroke': 'Black', 'fill': 'none', 'stroke-width': '1px'})
        .call(lines_xAxis)
    .selectAll("text")
        .attr("y", 10)
        //.attr("x", 0)
        //.attr("dy", ".35em")
        //.attr("transform", "rotate(90)")
        .style("text-anchor", "middle");

    lines_canvas.append("text")
        .attr("class", "h4")
        .attr("text-anchor", "middle")
        .attr("x", width/2)
        .attr("y", height + 40)
        .text("Read Pairs per Barcode Cluster");
    
    lines_canvas.append("text")
        .attr("class", "h4")
        .attr("text-anchor", "middle")
        .attr("x", -height/2)
        .attr("y", -move_left)
        .attr("dy", ".75em")
        .attr("transform", "rotate(-90)")
        .text("Read Pairs");
    
    // add y axis text
    var group2 = lines_canvas.append("g")
        //.style({ 'stroke': 'Black', 'fill': 'none', 'stroke-width': '1px'})
        .attr('stroke', 'Black')
        .call(lines_yAxis);

    // Define the line
    var valueline = d3.svg.line()
        .x(function(d) { return lines_x(d[0]); })
        .y(function(d) { return lines_y(d[1]); });
    
    // Add the valueline path.
    lines_canvas.append("path")
        .attr("id","cluster_line_1")
        .attr("class", "line")
        .attr("d", valueline( x_data_1.map(function (e, i) { return [e, y_data_1[i]]; }) ));

    lines_canvas.append("path")
        .attr("id","cluster_line_2")
        .attr("class", "line")
        .attr("d", valueline( x_data_1.map(function (e, i) { return [e, y_data_2[i]]; }) ));

    var brush = d3.svg.brush()
        .x(lines_x)
        .on("brushend", brushed)

    var gBrush = lines_canvas.append("g")
        .call(brush)
        .selectAll("rect")
            .attr("y", height-300)
            .attr("height", 300)
            .attr("fill", "#EAEAEA");

    function tabulate(data, columns, header_text) {
      d3.select('#graph_container_3').html('')
      var table = d3.select('#graph_container_3').append('table')
      table.attr("class","table")
      var thead = table.append('thead')
      var	tbody = table.append('tbody');

      // append the header row
      thead.append('tr')
        .selectAll('th')
        .data(columns).enter()
        .append('th')
          .text(function (column) { return header_text[column]; });

      // create a row for each object in the data
      var rows = tbody.selectAll('tr')
        .data(data)
        .enter()
        .append('tr');

      // create a cell in each row for each column
      var cells = rows.selectAll('td')
        .data(function (row) {
          return columns.map(function (column) {
            return {column: column, value: row[column]};
          });
        })
        .enter()
        .append('td')
          .text(function (d) { return d.value; });

      return table;
    }
            
    function brushed() {
    //redraw(age_chart_info.brush.extent(),nSexCutoff)
    //update()
      console.log(brush.extent().map( function (i) { return parseInt(i)}))
      filtered_data = brushFilter(data,brush.extent())
      tmp_percentage = filtered_data.map( function (i) { return i.total; } ).reduce(add, 0) / total_read_count
      tmp_percentage2 = filtered_data.length / total_cluster_count
      console.log( tmp_percentage  )
      //d3.select("#info_lable_1").html(Math.round(tmp_percentage*10000)/100 + "% of reads<br>"+Math.round(tmp_percentage2*10000)/100 + "% of clusters");
      //d3.select("#info_lable_2").text(Math.round(tmp_percentage2*10000)/100 + "%");
      if ( brush.extent()[0] !=  brush.extent()[1]  ) {
        d3.select("#p_low").text( Math.round(brush.extent()[0]));
        d3.select("#p_high").text(Math.round(brush.extent()[1]));
      }else{
        d3.select("#p_low").text( 2);
        d3.select("#p_high").text(d3.max(x_data_1));
      }
      d3.select("#p_of_reads").text(   Math.round(tmp_percentage*10000)/100);
      d3.select("#p_num_clusters").text(filtered_data.length);
      d3.select("#p_of_clusters").text(Math.round(tmp_percentage2*10000)/100);
      d3.select("#p_average").text(Math.round(filtered_data.map( function (i) { return i.total; } ).reduce(add, 0)/filtered_data.length));
      if ( filtered_data.length > 2  ) {
        d3.select("#p_median").text(filtered_data[parseInt(Math.round(filtered_data.length/2))].total);
      }else{
        d3.select("#p_median").text("NaN");
      }
      redraw(brush.extent())

	// render the table(s)
	tabulate(filtered_data, ['id', 'total','inbam','mapped'], {'id':'Cluster Id', 'total':'Read Pairs in cluster','inbam':'SE Read in bam','mapped':'SE Reads Mapped'})
      
    };
 
  }
  
  initiate_plots( data )
    
  //var sum = [1, 2, 3].reduce(add, 0);
  //array.forEach(function(currentValue, index, arr))
  //array.reverse()
  
});
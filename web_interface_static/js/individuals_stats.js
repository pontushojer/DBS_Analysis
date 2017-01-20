
function tabulate(container_id, data, columns, header_bool,header_text) {
  d3.select(container_id).html('')
  var table = d3.select(container_id).append('table')
  table.attr("class","table")
  var thead = table.append('thead')
  var	tbody = table.append('tbody');

  // append the header row
  if (header_bool) {
    thead.append('tr')
      .selectAll('th')
      .data(columns).enter()
      .append('th')
        .text(function (column) { return header_text[column]; });
  }

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

function appendTable(container_id, data, columns, header_bool,header_text) {
  d3.select(container_id).html(d3.select(container_id).html()+"<br>")
  var table = d3.select(container_id).append('table')
  table.attr("class","table")
  var thead = table.append('thead')
  var	tbody = table.append('tbody');

  // append the header row
  if (header_bool) {
    thead.append('tr')
      .selectAll('th')
      .data(columns).enter()
      .append('th')
        .text(function (column) { return header_text[column]; });
  }

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

d3.json("individuals.json", function(error, data) {
    tabulate( '#ind_stats_table', data,
             ['id', 'a1','a2'],
             true,
             {
                'id':'individual Id',
                'a1':'Allele 1, % of all reads, % of reads in ind.',
                'a2':'Allele 2, % of all reads, % of reads in ind.',
                }
            )
    }
)

d3.json("ind_details.json", function(error, data) {

    d3.select('#ind_details').html("")
    for (i in data) {
        
        individual = data[i]
        d3.select('#ind_details').html(d3.select('#ind_details').html()+"<h3>Individual "+individual.ind_id+":</h3>")
        appendTable( '#ind_details', individual.alleles,
                 ['id', 'barcode_cluster_count','barcode_cluster_perc','read_count',"read_perc_tot","read_perc_ind"],
                 true,
                 {
                        'id':'Id',
                        "barcode_cluster_count": "Barcode Clusters",
                        "barcode_cluster_perc": "% of Barcode Clusters",
                        "read_count":"Reads",
                        "read_perc_tot":"% of Total Reads",
                        "read_perc_ind":"% of Reads in ind."+individual.ind_id
                    }
                )
        }
    }
    

)

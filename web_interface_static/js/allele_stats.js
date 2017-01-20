console.log('hej')

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

d3.json("allele_matches.json", function(error, data) {
    console.log('before')
    console.log(data)
    console.log('after')

    tabulate( '#allele_stats_table', data,
             ['id', 'barcode_clusters','reads','in_sample','match'],
             true,
             {
                'id':'Allele Id',
                'barcode_clusters':'Barcode Clusters',
                'reads':'% of Total Reads*',
                'in_sample':'Found in sample ID',
                'match':'IPD-IMGT/HLA Correspondence'
                }
            )
    }
)

d3.json("cluster_trash.json", function(error, data) {
    console.log('before2')
    console.log(data)
    console.log('after2')

    tabulate( '#cluster_trash_table', data,
             ['id', 'barcode_clusters','reads'],
             true,
             {
                'id':'Reason',
                'barcode_clusters':'Barcode Clusters',
                'reads':'% of Total Reads*',
                }
            )
    }
)

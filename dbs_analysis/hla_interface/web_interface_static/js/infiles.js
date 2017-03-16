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

d3.json("infiles.json", function(error, data) {

    tabulate( '#infiles_table', data,
             ['file_pair_id', 'read_count','r1_path','r2_path'],
             true,
             {
                'file_pair_id':'Id',
                'read_count':'Read Pairs',
                'r1_path':'Read 1 Filename',
                'r2_path':'Read 2 Filename'                }
            )
    }
)


#! /usr/bin/env python

#
# Imports
#
from flask import Flask
from flask import render_template
from flask import json
from flask import Response
from dbs_analysis import metadata
import sys
import time
import os
from flask import request

#
# initiate the flask application
#
static_folder_path = '/'.join(os.path.realpath(__file__).split('/')[:-1])+'/web_interface_static/'
template_folder_path = '/'.join(os.path.realpath(__file__).split('/')[:-1])+'/templates/'
app = Flask(
    __name__,
    static_folder   = static_folder_path,
    static_url_path = '/static',
    template_folder=template_folder_path,
    )

# #
# # check input and get commandline args
# #
# try:
#     analysisfolder = metadata.AnalysisFolder(sys.argv[1])
#     #analysisfolder.readindexTsv()
# except IndexError: sys.stderr.write('please supply a commandline on format:\n'+os.path.basename(__file__)+' <analysis-output-folder>\n');sys.exit()
# app.analysisfolder = analysisfolder
# app.analysisfolder.settings.mapqCutOff = int(app.analysisfolder.settings.mapqCutOff)
# app.analysisfolder.dataPath=os.path.abspath(app.analysisfolder.dataPath) 
# 
# #
# # check analysis folder
# #
# if not analysisfolder.checkIntegrity() == 'PASS': print analysisfolder.checkIntegrity()+'\nERROR: Now exiting'
# 
# #
# # create a logfile
# #
# logfile = open(analysisfolder.logpath+'/dbs_hla_server.log.txt','a',1)
# logfile.write(time.strftime("%y%m%d-%H:%M:%S",time.localtime())+'\n')
# logfile.write('cmd: '+' '.join(sys.argv)+'\n')
# analysisfolder.logfile = logfile

@app.route("/")
def start():
    """ initial page does nothing but redirecting the user to the index.html route """
    
    page =  """<!DOCTYPE HTML>
        <html lang="en-US">
        <head>
            <meta charset="UTF-8">
            <meta http-equiv="refresh" content="1;url=index.html">
            <script type="text/javascript">
              window.location.href = "index.html"
            </script>
            <title>Page Redirection</title>
        </head>
        <body>
            If you are not redirected automatically, follow <a href='index.html'>this link</a>
        </body>
      </html>"""
    
    return page

@app.route("/index.html")
def index():
    """ function that builds the index page of the webinterface """
    return render_template('index.html',analysispath=app.analysisfolder.path)

@app.route("/read_pairs")
def read_pairs():
    from dbs_analysis.misc import thousandString
    from dbs_analysis.seqdata import revcomp
    from dbs_analysis.sequences import HLA_DBS as DBS
    from dbs_analysis.sequences import HLA_H1 as H1
    from dbs_analysis.sequences import HLA_H2 as H2
    from dbs_analysis.sequences import HLA_H3 as H3
    # from seqdata import revcomp
    # from misc import HtmlColors
    # layout_html = HtmlColors.BlueIntense+H1+HtmlColors.Color_Off+'-'+HtmlColors.CyanIntense+DBS+HtmlColors.Color_Off+'-'+HtmlColors.PurpleIntense+revcomp(H2)+HtmlColors.Color_Off+HtmlColors.BlackIntense+'- DNA Sequence Of Targeted Amplicon -'+HtmlColors.Color_Off+HtmlColors.YellowIntense+revcomp(H3)+HtmlColors.Color_Off
    return render_template('read_pairs.html',total_read_pairs=thousandString(app.analysisfolder.results.totalReadCount),h1=H1,h2=revcomp(H2),h3=revcomp(H3),dbs=DBS,bt2=str(app.analysisfolder.results.bt2AlignmentRate),bt2q20=str(app.analysisfolder.results.alignmentRateQ20),analysispath=app.analysisfolder.path)

@app.route("/barcode_clusters")
def barcode_clusters():
    from dbs_analysis import misc
    return render_template('barcode_clusters.html', non_singleton_clusters=misc.thousandString(app.analysisfolder.results.barcodeClusterCount-app.analysisfolder.results.singeltonBarcodeClusters),analysispath=app.analysisfolder.path)

@app.route("/alleles")
def alleles():
    import json
    from dbs_analysis import misc
    import operator
    with open( app.analysisfolder.dataPath+'/find_alleles_info.dict' ) as infile: info_vars = eval(infile.read())
    with open( app.analysisfolder.dataPath+'/allele_matches.dict' ) as infile: allele_matches = eval(infile.read())
    with open( app.analysisfolder.dataPath+'/allele_types_dict' ) as infile: allele_types = eval(infile.read())
    print info_vars.keys()
    print allele_matches
    print allele_types.keys()
    return render_template(
        'alleles.html',
        reference_id='HLA-A',
        allele_ids=[allele_id for readcount, allele_id in sorted([(sum([count for barcode_cluster_id, count in list(info_vars['allele_cluster_info'][allele_id]['barcode_clusters_in_this_allele_cluster'].iteritems())+list(info_vars['allele_cluster_info'][allele_id]['supporting_barcode_clusters'].iteritems())]), allele_id) for allele_id in info_vars['allele_cluster_info'].keys()],key=operator.itemgetter(0), reverse=True)],#info_vars['allele_cluster_info'].keys(),
        allele_count=len(info_vars['allele_cluster_info'].keys()),
        barcode_cluster_count=misc.thousandString(info_vars['barcode_clusters_loaded_for_find_alleles_step']),
        read_pair_count=misc.thousandString(info_vars['read_pairs_loaded_for_find_alleles_step']),
        read_pair_percentage=misc.percentage(info_vars['read_pairs_loaded_for_find_alleles_step'],app.analysisfolder.results.totalReadCount),
        min_reads = info_vars['read_cutoff_for_variable_position_info'],
        analysispath=app.analysisfolder.path
        )

@app.route("/individuals")
def individuals():
    with open( app.analysisfolder.dataPath+'/allele_types_dict' ) as infile: allele_types = eval(infile.read())
    return render_template('individuals.html',ind_count=sum([1 for ind_id in allele_types if ind_id not in ['No info','Non clonal','unknown id','nosiy signal']]),analysispath=app.analysisfolder.path)

@app.route("/ind_details.json")
def ind_details_json():
    import json
    from dbs_analysis import misc
    import operator
    import re
    
# trc = total read count
# bc = barcode cluster
# bci = barcodes cluster in individual
# rci = read counts in individual
# 
# Individual id "No info" has 976 barcode clusters (22.14% ot total read count):
# 	#bc:	%trc:	%bci:	%rci:	what:
    with open( app.analysisfolder.dataPath+'/find_alleles_info.dict' ) as infile: info_vars = eval(infile.read())
    with open( app.analysisfolder.dataPath+'/allele_types_dict' ) as infile: allele_types = eval(infile.read())

    out_list = []
    in_sample = []
    
    for ind_id, total_reads_in_this_ind in sorted([ (ind_id, sum([sum(read_counts_list) for read_counts_list in allele_types_in_this_ind.values()])) for ind_id, allele_types_in_this_ind in allele_types.iteritems()], key=operator.itemgetter(1),reverse=True):
        allele_types_in_this_ind = allele_types[ind_id]
        barcode_clusters = sum([len(read_counts_list) for read_counts_list in allele_types_in_this_ind.values()])
        
        alleles = []
        for allele_id, read_count in  sorted([ ( allele_id, sum(read_counts_list) ) for allele_id, read_counts_list in allele_types_in_this_ind.iteritems() ], key=operator.itemgetter(1),reverse=True):
            read_counts_list = allele_types_in_this_ind[allele_id]

            alleles.append(
                    {
                        'id':allele_id,
                        "barcode_cluster_count": str(len(read_counts_list)),
                        "barcode_cluster_perc": str(misc.percentage(len(read_counts_list), barcode_clusters))+'%',
                        "read_count":misc.thousandString(read_count),
                        "read_perc_tot":str(misc.percentage(read_count,app.analysisfolder.results.totalReadCount))+'%',
                        "read_perc_ind":str(misc.percentage(read_count,total_reads_in_this_ind))+'%'
                    }
                )
        alleles.append(
            {
                'id':"Total",
                "barcode_cluster_count": str(barcode_clusters),
                "barcode_cluster_perc": str(misc.percentage(barcode_clusters, barcode_clusters))+'%',
                "read_count":misc.thousandString(total_reads_in_this_ind),
                "read_perc_tot":str(misc.percentage(total_reads_in_this_ind,app.analysisfolder.results.totalReadCount))+'%',
                "read_perc_ind":str(misc.percentage(total_reads_in_this_ind,total_reads_in_this_ind))+'%'
            }
        )

        out_list.append( {'ind_id': ind_id, 'alleles':alleles} )

    return json.dumps(out_list)

@app.route("/individuals.json")
def individuals_json():
    import json
    from dbs_analysis import misc
    import operator
    import re
    
    with open( app.analysisfolder.dataPath+'/find_alleles_info.dict' ) as infile: info_vars = eval(infile.read())
    with open( app.analysisfolder.dataPath+'/allele_types_dict' ) as infile: allele_types = eval(infile.read())

    out_list = []
    in_sample = []
    
    for ind_id, total_reads_in_this_ind in sorted([ (ind_id, sum([sum(read_counts_list) for read_counts_list in allele_types_in_this_ind.values()])) for ind_id, allele_types_in_this_ind in allele_types.iteritems()], key=operator.itemgetter(1),reverse=True):
        allele_types_in_this_ind = allele_types[ind_id]

        #total_reads_in_this_ind = sum([sum(read_counts_list) for read_counts_list in allele_types_in_this_ind.values()])
        
        alleles = []
        for allele_id, read_count in  sorted([ ( allele_id, sum(read_counts_list) ) for allele_id, read_counts_list in allele_types_in_this_ind.iteritems() ], key=operator.itemgetter(1),reverse=True):
            if ind_id not in ['No info', 'nosiy signal', 'unknown id','Non clonal'] and not re.search('trashed',allele_id):
                alleles.append({'id':allele_id.replace('allele ',''), 'perc_of_total':str(misc.percentage(read_count,app.analysisfolder.results.totalReadCount))+'%', 'perc_of_self':str(misc.percentage(read_count,total_reads_in_this_ind))+'%'})
                #str(allele_id).replace('allele ','')+', '+str(misc.percentage(read_count,app.analysisfolder.results.totalReadCount))+'%, '+str(misc.percentage(read_count,total_reads_in_this_ind))+'%')

        if len(alleles) :
            out_list.append(
                    {
                        "id":ind_id,
                        "a1":alleles[0]['id'],
                        "a1_perc_tot":alleles[0]['perc_of_total'],
                        "a1_perc_self":alleles[0]['perc_of_self'],
                        "a2":alleles[1]['id'] if len(alleles) >=2 else '-',
                        "a2_perc_tot":alleles[1]['perc_of_total'] if len(alleles) >=2 else '-',
                        "a2_perc_self":alleles[1]['perc_of_self'] if len(alleles) >=2 else '-',
                    }
                )

    return json.dumps(out_list)

@app.route("/cluster_trash.json")
def cluster_trash_json():
    import json
    from dbs_analysis import misc
    
    with open( app.analysisfolder.dataPath+'/find_alleles_info.dict' ) as infile: info_vars = eval(infile.read())
    with open( app.analysisfolder.dataPath+'/allele_types_dict' ) as infile: allele_types = eval(infile.read())

    # info_vars['trash_bin'] = {
    #     'No_similarity_to_allele_cluster': { barcode_cluster.id:barcode_cluster.readPairCount }, 
    #     'non clonal': {
    #         'hetro_bases':  { barcode_cluster.id:barcode_cluster.readPairCount }, 
    #         'individual_id':{ barcode_cluster.id:barcode_cluster.readPairCount }
    #         }, 
    #     'matching_more_than_1_allele_cluster': { barcode_cluster.id:barcode_cluster.readPairCount  }, 
    #     'Low_read_count': { barcode_cluster.id:barcode_cluster.readPairCount  }
    # }

    out_list = []
    in_sample = []
    total_barcode_cluster_count = 0
    total_read_count = 0
    table_name = {'No_similarity_to_allele_cluster':'Match no allele', 'matching_more_than_1_allele_cluster':'Match several alleles', 'Low_read_count':'Below read pair count cut off','hetro_bases':'Non clonal, has hetrozygous base overlapping potential variant','individual_id':'Non clonal, has several individual ids'}
    for key in ['No_similarity_to_allele_cluster', 'matching_more_than_1_allele_cluster', 'Low_read_count']:
        readcount = sum(info_vars['trash_bin'][key].values())
        reads = misc.percentage(readcount, app.analysisfolder.results.totalReadCount)
        total_read_count += readcount
        total_barcode_cluster_count += len(info_vars['trash_bin'][key])
        out_list.append(
                {
                    "id":table_name[key],
                    "barcode_clusters":str(len(info_vars['trash_bin'][key])) + ' ('+str(misc.percentage(len(info_vars['trash_bin'][key]),info_vars['barcode_clusters_loaded_for_find_alleles_step']))+'%)',
                    "reads":str(reads)+'%',
                    "in_sample":" ,".join(in_sample),
                    "match":''
                }
            )
    for key in ['hetro_bases','individual_id']:
        readcount = sum(info_vars['trash_bin']['non clonal'][key].values())
        reads = misc.percentage(readcount, app.analysisfolder.results.totalReadCount)
        total_read_count += readcount
        total_barcode_cluster_count += len(info_vars['trash_bin']['non clonal'][key])
        out_list.append(
                {
                    "id":table_name[key],
                    "barcode_clusters":str(len(info_vars['trash_bin']['non clonal'][key])) + ' ('+str(misc.percentage(len(info_vars['trash_bin']['non clonal'][key]),info_vars['barcode_clusters_loaded_for_find_alleles_step']))+'%)',
                    "reads":str(reads)+'%',
                    "in_sample":" ,".join(in_sample),
                    "match":''
                }
            )
    out_list.append(
            {
                "id":'In Total',
                "barcode_clusters":str(total_barcode_cluster_count) + ' ('+str(misc.percentage(total_barcode_cluster_count,info_vars['barcode_clusters_loaded_for_find_alleles_step']))+'%)',
                "reads":str(misc.percentage(total_read_count, app.analysisfolder.results.totalReadCount))+'%',
                "in_sample":"",
                "match":''
            }
        )    
    return json.dumps(out_list)

@app.route("/allele_matches.json")
def allele_matches_json():
    
    from dbs_analysis import misc
    import operator
    
    with open( app.analysisfolder.dataPath+'/find_alleles_info.dict' ) as infile: info_vars = eval(infile.read())
    with open( app.analysisfolder.dataPath+'/allele_matches.dict' ) as infile: allele_matches = eval(infile.read())
    with open( app.analysisfolder.dataPath+'/allele_types_dict' ) as infile: allele_types = eval(infile.read())
    
    total_read_count = 0
    total_barcode_cluster_count = 0
    
    out_list = []
    #for allele_id in sorted(info_vars['allele_cluster_info'].keys()):
    for readcount, allele_id in sorted([(sum([count for barcode_cluster_id, count in list(info_vars['allele_cluster_info'][allele_id]['barcode_clusters_in_this_allele_cluster'].iteritems())+list(info_vars['allele_cluster_info'][allele_id]['supporting_barcode_clusters'].iteritems())]), allele_id) for allele_id in info_vars['allele_cluster_info'].keys()],key=operator.itemgetter(0), reverse=True):
        
        matches = [tmp_string.replace('A*','').split(' ')[1] for tmp_string in allele_matches[allele_id].keys()]
        barcode_clusters = list(info_vars['allele_cluster_info'][allele_id]['barcode_clusters_in_this_allele_cluster'].iteritems())+list(info_vars['allele_cluster_info'][allele_id]['supporting_barcode_clusters'].iteritems())
        #readcount = sum([count for barcode_cluster_id, count in barcode_clusters])
        reads = misc.percentage(readcount, app.analysisfolder.results.totalReadCount)
        total_read_count += readcount
        total_barcode_cluster_count += len(barcode_clusters)
        
        in_sample = []
        for ind_id in allele_types:
            allele_types_in_this_ind = allele_types[ind_id]
            if 'allele '+str(allele_id) in allele_types_in_this_ind.keys():
                total_reads_in_this_ind = sum([sum(read_counts_list) for read_counts_list in allele_types_in_this_ind.values()])
                if misc.percentage(sum(allele_types_in_this_ind['allele '+str(allele_id)]),app.analysisfolder.results.totalReadCount) > 1.0 and ind_id != 'No info': in_sample.append(ind_id)
        
        out_list.append(
                {
                    "id":allele_id,
                    "barcode_clusters":str(len(barcode_clusters)) + ' ('+str(misc.percentage(len(barcode_clusters),info_vars['barcode_clusters_loaded_for_find_alleles_step']))+'%)',
                    "reads":str(reads)+'%',
                    "in_sample":" ,".join(in_sample),
                    "match":', '.join(matches)
                }
            )

    out_list.append(
            {
                "id":'In Total',
                "barcode_clusters":str(total_barcode_cluster_count) + ' ('+str(misc.percentage(total_barcode_cluster_count,info_vars['barcode_clusters_loaded_for_find_alleles_step']))+'%)',
                "reads":str(misc.percentage(total_read_count, app.analysisfolder.results.totalReadCount))+'%',
                "in_sample":"",
                "match":''
            }
        )    
    import json
    return json.dumps(out_list)#'[{"id":1, "barcode_clusters":1, "reads":10, "in_sample":1,"match":1},{"id":2, "barcode_clusters":7489,"reads":10,"in_sample":"1 , 8 , 9 , 0","match":"matchar pelle och kurt"}]'

@app.route("/barcode_clusters.json")
def barcode_clusters_json():
    with open(app.analysisfolder.dataPath+'/barcode_cluster.json') as infile: return infile.read()

@app.route("/reference.fasta.fai")
def reference_fasta_fai():
   with open(app.analysisfolder.settings.bowtie2Reference+'.fai') as infile: return infile.read()

@app.route("/reference.fasta")
def reference_fasta():
#   with open(app.analysisfolder.settings.bowtie2Reference) as infile: return infile.read()
    
    import os
    import sys
    import re
    from flask import request, send_from_directory

    full_path = os.path.abspath(os.path.join(app.analysisfolder.settings.bowtie2Reference, app.analysisfolder.settings.bowtie2Reference))
    print full_path

    # security check - only files under READ_VIZ_DIR should be accsessible
    if not full_path.startswith(app.analysisfolder.settings.bowtie2Reference):
        return "Invalid path: %s" % app.analysisfolder.settings.bowtie2Reference

    # handle igv.js Range header which it uses to request a subset of a .bam
    range_header = request.headers.get('Range', None)
    if not range_header:
        return send_from_directory(app.analysisfolder.settings.bowtie2Reference, app.analysisfolder.settings.bowtie2Reference)

    m = re.search('(\d+)-(\d*)', range_header)
    if not m:
        error_msg = "ERROR: unexpected range header syntax: %s" % range_header
        sys.stderr.write(error_msg+'\n')
        return error_msg

    size = os.path.getsize(full_path)
    offset = int(m.group(1))
    length = int(m.group(2) or size) - offset

    data = None
    with open(full_path, 'rb') as f:
        f.seek(offset)
        data = f.read(length)

    rv = Response(data, 206, mimetype="application/octet-stream", direct_passthrough=True)
    rv.headers.add('Content-Range', 'bytes {0}-{1}/{2}'.format(offset, offset + length - 1, size))

    sys.stderr.write("readviz: range request: %s-%s %s" % (m.group(1), m.group(2), full_path)+'\n')
    return rv

@app.route("/target_regions.bed")
def target_regions_bed():
    with open(app.analysisfolder.settings.targetRegionBed) as infile: return infile.read()

@app.route('/read_viz/<path:path>')
def read_viz_files(path):
    
    import os
    import sys
    import re
    from flask import request, send_from_directory

    full_path = os.path.abspath(os.path.join(app.analysisfolder.dataPath, path))

    # security check - only files under READ_VIZ_DIR should be accsessible
    if not full_path.startswith(app.analysisfolder.dataPath): return "Invalid path: %s" % path

    # handle igv.js Range header which it uses to request a subset of a .bam
    range_header = request.headers.get('Range', None)
    if not range_header: return send_from_directory(app.analysisfolder.dataPath, path)

    m = re.search('(\d+)-(\d*)', range_header)
#    print 'range_header=',range_header
    if not m:
        error_msg = "ERROR: unexpected range header syntax: %s" % range_header
        sys.stderr.write(error_msg+'\n')
        return error_msg

    size = os.path.getsize(full_path)
    offset = int(m.group(1))
    length = int(m.group(2) or size) - offset
#    print 'offset=',offset
#   print 'length=',length
#    if length == size: print 'whole file requested!'

    data = None
    with open(full_path, 'rb') as f:
        f.seek(offset)
        data = f.read(length)

    rv = Response(data, 206, mimetype="application/octet-stream", direct_passthrough=True)
    rv.headers.add('Content-Range', 'bytes {0}-{1}/{2}'.format(offset, offset + length - 1, size))
    rv.headers.add('Content-Transfer-Encoding','binary')
    rv.headers.add('Accept-Ranges','bytes')
    
    sys.stderr.write("readviz: range request: %s-%s %s" % (m.group(1), m.group(2), full_path)+'\n')

    print "\n\n-request-"
    print request.headers
    print "-response-"
    print rv.headers
    return rv

@app.route("/mapping_stats.json")
def make_mapping_stats_json():
    #pattern="""(?P<totalReads>\d+) reads; of these:\n\s+(?P<pairedReads>\d+) \(\d+.\d+\%\) were paired; of these:\n\s+(?P<notPropMapedPair>\d+) \(\d+.\d+\%\) aligned concordantly 0 times\n\s+(?P<properPairs>\d+) \(\d+.\d+\%\) aligned concordantly exactly 1 time\n\s+(?P<properPairsMultiMap>\d+) \(\d+.\d+\%\) aligned concordantly >1 times\n\s+----\n\s+(?P<notPropMapedPair2>\d+) pairs aligned concordantly 0 times; of these:\n\s+(?P<discordantPairs>\d+) \(\d+.\d+\%\) aligned discordantly 1 time\n\s+----\n\s+(?P<unMappedPair>\d+) pairs aligned 0 times concordantly or discordantly; of these:\n\s+(?P<possibleSingletons>\d+) mates make up the pairs; of these:\n\s+(?P<unMappedReads>\d+) \(\d+.\d+\%\) aligned 0 times\n\s+(?P<singleSingleMap>\d+) \(\d+.\d+\%\) aligned exactly 1 time\n\s+(?P<singleMultiMap>\d+) \(\d+.\d+\%\) aligned >1 times\n(?P<overallAlignmentRate>\d+.\d+)\% overall alignment rate"""
    pattern="""(?P<totalReads>\d+) reads; of these:
\s+(?P<pairedReads>\d+) \(\d+.\d+\%\) were paired; of these:
\s+(?P<notPropMapedPair>\d+) \(\d+.\d+\%\) aligned concordantly 0 times
\s+(?P<properPairs>\d+) \(\d+.\d+\%\) aligned concordantly exactly 1 time
\s+(?P<properPairsMultiMap>\d+) \(\d+.\d+\%\) aligned concordantly >1 times
\s+----
\s+(?P<notPropMapedPair2>\d+) pairs aligned concordantly 0 times; of these:
\s+(?P<discordantPairs>\d+) \(\d+.\d+\%\) aligned discordantly 1 time
\s+----
\s+(?P<unMappedPair>\d+) pairs aligned 0 times concordantly or discordantly; of these:
\s+(?P<possibleSingletons>\d+) mates make up the pairs; of these:
\s+(?P<unMappedReads>\d+) \(\d+.\d+\%\) aligned 0 times
\s+(?P<singleSingleMap>\d+) \(\d+.\d+\%\) aligned exactly 1 time
\s+(?P<singleMultiMap>\d+) \(\d+.\d+\%\) aligned >1 times
(?P<overallAlignmentRate>\d+.\d+)\% overall alignment rate"""
    import re
    from dbs_analysis.misc import percentage
    with open(app.analysisfolder.dataPath+'/bt2.stat.txt') as infile:
        data = infile.read()
        p = re.compile(pattern)
        m = p.search(data)
        m.groupdict()['overallAlignmentRate'] #m.groupdict()['']
        
        total = int(app.analysisfolder.results.totalReadCount)*2
        not_in_bam = int(app.analysisfolder.results.totalReadCount)*2 - int(m.groupdict()['pairedReads'])*2
        
        # mapped
        proper_pair = int(m.groupdict()['properPairs'])*2
        proper_pair_mulit = int(m.groupdict()['properPairsMultiMap'])*2
        not_pp_aligned = int(m.groupdict()['discordantPairs'])*2
        mapped_SE_reads = int(m.groupdict()['singleSingleMap'])
        multim_SE_reads = int(m.groupdict()['singleMultiMap'])
        mapped_total = sum([proper_pair,proper_pair_mulit,not_pp_aligned,mapped_SE_reads,multim_SE_reads])
        
        # not mapped
        unmapped_SE_reads = int(m.groupdict()['unMappedReads'])
        
        in_bam = sum([proper_pair,proper_pair_mulit,not_pp_aligned,mapped_SE_reads,multim_SE_reads,unmapped_SE_reads])
    
    json  ='     {"name": "Total SE reads ", "children": ['
    json +='        {"name": "'+str(percentage(not_in_bam,total))+'%'+' not in bam", "size": '+str(not_in_bam)+'},'
    json +='        {"name": "'+str(percentage(in_bam,total))+'%'+' in bam", "children": ['
    json +='            {"name": "'+str(percentage(unmapped_SE_reads,total))+'%'+' unMapped", "size": '+str(unmapped_SE_reads)+'},'
#    json +='            {"name": "'+str(percentage(mapped_total,total))+'%'+' mapped concordantly", "size": '+str(mapped_total)+'}'
    json +='            {"name": "'+str(percentage(mapped_total,total))+'%'+' Mapped", "children": ['
    json +='                {"name": "'+str(percentage(proper_pair,total))+'%'+' Mapped as Proper Pair", "size": '+str(proper_pair)+'},'
    json +='                {"name": "'+str(percentage(proper_pair_mulit,total))+'%'+' multi mapped pair", "size": '+str(proper_pair_mulit)+'},'
    json +='                {"name": "'+str(percentage(not_pp_aligned,total))+'%'+' mapped discordantly", "size": '+str(not_pp_aligned)+'},'
    json +='                {"name": "'+str(percentage(mapped_SE_reads,total))+'%'+' mapped single end", "size": '+str(mapped_SE_reads)+'},'
    json +='                {"name": "'+str(percentage(multim_SE_reads,total))+'%'+' multi mapped single end", "size": '+str(multim_SE_reads)+'}'
    json +='            ] }'
    json +='        ] }'
    json +='     ] }'
    return json

@app.route("/dbs_match.json")
def make_dbs_match_json():
    import time
    from dbs_analysis.misc import percentage
    
    if app.analysisfolder.results.readsWithDbsPatternMatch:
        tmp = eval(app.analysisfolder.results.readsWithDbsPatternMatch)
        if False not in tmp: tmp[False] = 0
        total = app.analysisfolder.results.totalReadCount
    else: tmp = {False:0, True:0, None:0}

    json  ='     {"name": "Total read pairs", "children": ['
    json +='        {"name": "'+str(percentage(tmp[None],total))+'%'+' Barcode not found", "size": '+str(tmp[None])+'},'
    json +='        {"name": "'+str(percentage(tmp[True]+tmp[False],total))+'%'+' Barcode found", "children": ['
    json +='            {"name": "'+str(percentage(tmp[True],total))+'%'+' Barcode match pattern", "size": '+str(tmp[True])+'},'
    json +='            {"name": "'+str(percentage(tmp[False],total))+'%'+' Barcode don\'t match pattern", "size": '+str(tmp[False])+'}'
    json +='        ] }'
    json +='     ] }'
    return json

@app.route("/handles.json")
def makehandlejson():
    """ gets the handle content information and dums it to json """
    
    import operator
    from dbs_analysis.misc import percentage
    from dbs_analysis.misc import thousandString
    constructOK = 0
    missing_h3 = 0
    missing_h2 = 0
    missing_h3_h2 = 0
    missing_h1_h2_h3 = 0
    missing_h1_h2 = 0
    missing_h1 = 0
    missing_h1_h3 = 0
    if app.analysisfolder.results.constructTypes:
        tmpDict = eval(app.analysisfolder.results.constructTypes)
        for what, count in sorted(eval(app.analysisfolder.results.constructTypes).iteritems(), key=operator.itemgetter(1))[::-1]:
            what = what.split(' ')
#            print ' '+str(what)+' => '+str(percentage(count,app.analysisfolder.results.totalReadCount))+'%, ('+thousandString(count)+')'
            if what == ['constructOK']: constructOK = count
            elif 'h1' in what:
                if 'h2' in what:
                    if 'h3' in what: missing_h1_h2_h3 = count
                    else: missing_h1_h2 = count
                elif 'h3' in what: missing_h1_h3 = count
                else: missing_h1 = count
            elif 'h2' in what:
                if 'h3' in what: missing_h3_h2 = count
                else: missing_h2 = count
            elif 'h3' in what: missing_h3 = count
            else: print 'unexpected', what
        total = constructOK +missing_h3 +missing_h2 +missing_h3_h2 +missing_h1_h2_h3 +missing_h1_h2 +missing_h1 +missing_h1_h3
    else: total = 0

    json  ='     {"name": "Total read pairs", "children": ['
    json +='        {"name": "'+str(percentage(constructOK+missing_h3+missing_h2+missing_h3_h2,total))+'%'+' Has h1", "children": ['
    json +='            {"name": "'+str(percentage(constructOK+missing_h3,total))+'%'+' Has h1 and h2", "children": ['
    json +='                {"name": "'+str(percentage(constructOK,total))+'%'+' Has h1, h2 and h3", "size": '+str(constructOK)+'},'
    json +='                {"name": "'+str(percentage(missing_h3,total))+'%'+' Has h1, h2 and miss h3", "size": '+str(missing_h3)+'}   '
    json +='            ] },'
    json +='            {"name": "'+str(percentage(missing_h2+missing_h3_h2,total))+'%'+' Has h1 and miss h2", "children": ['
    json +='                {"name": "'+str(percentage(missing_h2,total))+'%'+' Has h1, miss h2 and has h3", "size": '+str(missing_h2)+'},  '
    json +='                {"name": "'+str(percentage(missing_h3_h2,total))+'%'+' Has h1, miss h2 and miss h3", "size": '+str(missing_h3_h2)+'}'    
    json +='            ] }'
    json +='        ] },'
    json +='        {"name": "'+str(percentage(missing_h1_h2_h3+missing_h1_h2+missing_h1+missing_h1_h3,total))+'%'+' Miss h1", "children": ['
    json +='            {"name": "'+str(percentage(missing_h1_h2_h3+missing_h1_h2,total))+'%'+' Miss h1 and h2", "children": ['
    json +='                {"name": "'+str(percentage(missing_h1_h2_h3,total))+'%'+' Miss h1, h2 and h3", "size": '+str(missing_h1_h2_h3)+'},         '
    json +='                {"name": "'+str(percentage(missing_h1_h2,total))+'%'+' Miss h1, miss h2 and has h3", "size": '+str(missing_h1_h2)+'}   '
    json +='            ] },'
    json +='            {"name": "'+str(percentage(missing_h1+missing_h1_h3,total))+'%'+' Miss h1 and has h2", "children": ['
    json +='                {"name": "'+str(percentage(missing_h1,total))+'%'+' Miss h1, has h2 and h3", "size": '+str(missing_h1)+'},      '
    json +='                {"name": "'+str(percentage(missing_h1_h3,total))+'%'+' Miss h1, has h2 and miss h3", "size": '+str(missing_h1_h3)+'}     '
    json +='            ] }'
    json +='        ] }'
    json +='     ] }'

    return json

def get_bootstrap(static_folder_path):
    
    import sys
    import urllib2
    hdr = {'User-Agent': 'Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.11 (KHTML, like Gecko) Chrome/23.0.1271.64 Safari/537.11',
       'Accept': 'text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8',
       'Accept-Charset': 'ISO-8859-1,utf-8;q=0.7,*;q=0.3',
       'Accept-Encoding': 'none',
       'Accept-Language': 'en-US,en;q=0.8',
       'Connection': 'keep-alive'}
    
    sys.stderr.write('Looking for bootstrap files ... \n')
    if not os.path.isfile(static_folder_path+'js/bootstrap.js'):
        sys.stderr.write('bootstrap not found downloading ... \n')
        url = 'https://github.com/twbs/bootstrap/releases/download/v3.3.7/bootstrap-3.3.7-dist.zip'
        req = urllib2.Request(url, headers=hdr)
        response = urllib2.urlopen(req)
        data = response.read()
        zip_file_name = static_folder_path+'bootstrap-3.3.7-dist.zip'
        with open(zip_file_name,'w') as outfile: outfile.write(data)
        sys.stderr.write('download complete, unzipping ... \n')
        import zipfile
        zip_ref = zipfile.ZipFile(zip_file_name, 'r')
        zip_ref.extractall(static_folder_path)
        zip_ref.close()
        sys.stderr.write('extraction complete, installing ... \n')
        for folder in ['css','fonts','js']:
            if not os.path.isdir(static_folder_path+folder): os.mkdir(static_folder_path+folder)
        for tmp in [
                'css/bootstrap-theme.css',
                'css/bootstrap-theme.css.map',
                'css/bootstrap-theme.min.css',
                'css/bootstrap-theme.min.css.map',
                'css/bootstrap.css',
                'css/bootstrap.css.map',
                'css/bootstrap.min.css',
                'css/bootstrap.min.css.map',
                'fonts/glyphicons-halflings-regular.eot',
                'fonts/glyphicons-halflings-regular.svg',
                'fonts/glyphicons-halflings-regular.ttf',
                'fonts/glyphicons-halflings-regular.woff',
                'fonts/glyphicons-halflings-regular.woff2',
                'js/bootstrap.js',
                'js/bootstrap.min.js',
                'js/npm.js']:
            src = static_folder_path+'bootstrap-3.3.7-dist/'+tmp
            dst = static_folder_path+tmp
            os.rename(src,dst)
        sys.stderr.write('installation complete, removing temporary files ... \n')
        for folder in ['css','fonts','js']: os.rmdir(static_folder_path+'bootstrap-3.3.7-dist/'+folder)
        os.rmdir(static_folder_path+'bootstrap-3.3.7-dist/')
        os.remove(zip_file_name)
        sys.stderr.write('bootstrap now present in '+static_folder_path+'.\n')
    else:
        sys.stderr.write('bootstrap files found in '+static_folder_path+'.\n')

    sys.stderr.write('Looking for jquery files ... \n')
    if not os.path.isfile(static_folder_path+'js/jquery-3.1.1.js'):
        sys.stderr.write('jquery not found downloading ... \n')
        url = 'https://code.jquery.com/jquery-3.1.1.js'
        req = urllib2.Request(url, headers=hdr)
        response = urllib2.urlopen(req)
        data = response.read()
        zip_file_name = static_folder_path+'js/jquery-3.1.1.js'
        with open(zip_file_name,'w') as outfile: outfile.write(data)
        sys.stderr.write('download complete.\n')
    else:
        sys.stderr.write('jquery files found in '+static_folder_path+'js/jquery-3.1.1.js .\n')
        
    sys.stderr.write('Looking for d3 files ... \n')
    if not os.path.isfile(static_folder_path+'js/d3.v3.js'):
        sys.stderr.write('d3 not found downloading ... \n')
        url = 'http://d3js.org/d3.v3.js'
        # https://github.com/d3/d3/releases/download/v4.4.0/d3.zip
        req = urllib2.Request(url, headers=hdr)
        response = urllib2.urlopen(req)
        data = response.read()
        zip_file_name = static_folder_path+'js/d3.v3.js'
        with open(zip_file_name,'w') as outfile: outfile.write(data)
        sys.stderr.write('download complete.\n')
    else:
        sys.stderr.write('d3 files found in '+static_folder_path+'js/d3.v3.js .\n')

    sys.stderr.write('Looking for spin.js files ... \n')
    if not os.path.isfile(static_folder_path+'js/spin.js'):
        sys.stderr.write('spin.js not found downloading ... \n')
        url = 'http://spin.js.org/spin.js'
        # https://github.com/d3/d3/releases/download/v4.4.0/d3.zip
        req = urllib2.Request(url, headers=hdr)
        response = urllib2.urlopen(req)
        data = response.read()
        zip_file_name = static_folder_path+'js/spin.js'
        with open(zip_file_name,'w') as outfile: outfile.write(data)
        sys.stderr.write('download complete.\n')
    else:
        sys.stderr.write('spin.js files found in '+static_folder_path+'js/spin.js .\n')

    # https://igv.org/web/release/1.0.6/
    # [TXT]	igv-1.0.6.css	22-Dec-2016 13:37	25K	 
    # [TXT]	igv-1.0.6.js	22-Dec-2016 13:37	813K	 
    # [TXT]	igv-1.0.6.min.js	22-Dec-2016 13:37	381K	 
    # [   ]	igv-1.0.6.min.map	22-Dec-2016 13:37	414K	 
    # [DIR]	img/	22-Dec-2016 13:37	-
    #      [IMG]	cursor_logo.png	22-Dec-2016 13:37	5.3K	 
    #      [IMG]	cursor_logo.svg	22-Dec-2016 13:37	3.7K	 
    #      [IMG]	igv_logo_letters_paths.svg	22-Dec-2016 13:37	2.8K	  

    sys.stderr.write('Looking for igv.js and dependancies ... \n')
    for filename, web_folder,local_folder in [
        ('igv-1.0.6.js','http://igv.org/web/release/1.0.6/','js/'),
        ('igv-1.0.6.css','http://igv.org/web/release/1.0.6/','css/'),
        ('igv-1.0.6.min.js','http://igv.org/web/release/1.0.6/','js/'),
        ('igv-1.0.6.min.map','http://igv.org/web/release/1.0.6/','js/'),
        ('cursor_logo.png','http://igv.org/web/release/1.0.6/img/','img/'),
        ('cursor_logo.svg','http://igv.org/web/release/1.0.6/img/','img/'),
        ('igv_logo_letters_paths.svg','http://igv.org/web/release/1.0.6/img/','img/'),
        ('jquery-ui.css','https://ajax.googleapis.com/ajax/libs/jqueryui/1.11.2/themes/smoothness/','css/'),
        ('font-awesome.min.css','https://maxcdn.bootstrapcdn.com/font-awesome/4.2.0/css/','css/'),
        ('jquery-ui.min.js','https://ajax.googleapis.com/ajax/libs/jqueryui/1.11.2/','js/'),
        ('fontawesome-webfont.ttf','https://maxcdn.bootstrapcdn.com/font-awesome/4.2.0/fonts/','fonts/'),
        ('fontawesome-webfont.woff','https://maxcdn.bootstrapcdn.com/font-awesome/4.2.0/fonts/','fonts/'),
        ('jquery.min.js','https://ajax.googleapis.com/ajax/libs/jquery/1.11.1/','js/')]:
        if not os.path.isfile(static_folder_path+local_folder+filename):
            sys.stderr.write(filename +' not found downloading ... \n')
            url = web_folder+filename
            req = urllib2.Request(url, headers=hdr)
            response = urllib2.urlopen(req)
            data = response.read()
            outfile_name = static_folder_path+local_folder+filename
            with open(outfile_name,'w') as outfile: outfile.write(data)
            sys.stderr.write(url+' => '+outfile_name+' ... downloaded.\n')
    sys.stderr.write('all files present.\n')

def make_cluster_stats_json(analysisfolder, update=True):
    
    import operator
    import json
    import os
    
    if not os.path.isfile(analysisfolder.dataPath+'/barcode_cluster.json') or update:
    
        analysisfolder.database.getConnection()
        data = analysisfolder.database.c.execute('SELECT clusterId, clusterTotalReadCount, readPairsInBamFile, mappedSEReads, targetInfo FROM barcodeClusters WHERE clusterTotalReadCount > 1').fetchall()
        data = analysisfolder.database.c.execute('SELECT clusterId, clusterTotalReadCount, readPairsInBamFile, mappedSEReads FROM barcodeClusters WHERE clusterTotalReadCount > 1').fetchall()
        #tmp_list = [{'id':item[0],'total':item[1],'inbam':item[2]*2,'mapped':item[3],'ti':[{'name':entry['entry_name'],'ard':entry['averageReadDepth']} for entry in eval(item[4])] } for item in data ]
        tmp_list = [{'id':item[0],'total':item[1],'inbam':item[2],'mapped':item[3] } for item in data ]
        # analysisfolder.results.totalReadCount
        # analysisfolder.results.singeltonBarcodeClusters
        # analysisfolder.results.barcodeClusterCount
    
        barcode_cluster_jsons = []
        sort_dict = {}    
        for item in tmp_list:
            try:
                sort_dict[ item['total'] ].append( item )
            except KeyError:
                sort_dict[ item['total'] ] = [ item ]
        for count, sort_list in sorted(sort_dict.iteritems(), key=operator.itemgetter(0), reverse=True):
            for item in sort_list: barcode_cluster_jsons.append( item )
    
        #print '\n'.join( [ str(item['total'])+':'+str(item) for item in barcode_cluster_jsons ] )
        #print json.dumps(barcode_cluster_jsons)
    
        with open(analysisfolder.dataPath+'/barcode_cluster.json','w') as outfile: json.dump(barcode_cluster_jsons, outfile)

if __name__ == "__main__":

    import random
    import platform
    import subprocess

    get_bootstrap(static_folder_path)
    make_cluster_stats_json(app.analysisfolder)

    if app.analysisfolder.settings.port == 'random':
        app.run(host='0.0.0.0',port=random.randint(1000,9999),debug=True)
    else:
        app.run(host='0.0.0.0',port=int(app.analysisfolder.settings.port),debug=True)
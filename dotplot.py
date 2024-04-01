from lib.KJH_SVG.KJH_SVG import element
from optparse import OptionParser
import sys

#option parser
parser = OptionParser(usage="""Run annotation.py \n Usage: %prog [options]""")
parser.add_option("-i","--input",action = 'store',type = 'string',dest = 'INPUT',help = "")
parser.add_option("-m","--identity",action = 'store',type = 'float',dest = 'IDENTITY',help = "")
parser.add_option("-l","--length",action = 'store',type = 'int',dest = 'LENGTH',help = "")
(opt, args) = parser.parse_args()
if opt.INPUT == None or opt.IDENTITY == None or opt.LENGTH == None:
    print('     python dotplot.py --input test.blastn_m6 --identity 99 --length 10000')
    sys.exit()


def get_color(identity):
    if 99 < identity:
        return 'black'
    elif 90 < identity:
        return 'red'
    else:
        return 'green'





infile = opt.INPUT
threshold_identity = opt.IDENTITY 
threshold_length = opt.LENGTH

imageLength = 800
margin = 80



length_DICT = {}
fin = open(infile)
for line in fin:
    mylist = line.rstrip('\n').split('\t')
    if len(mylist) != 12: break
    query, sbjct, identity, alignment_length, mismatches, gap_openings, query_start, query_end, sbjct_start, sbjct_end, e_value, bit_score = mylist

    identity = float(identity)
    if identity < threshold_identity: continue

    alignment_length = int(alignment_length)
    if alignment_length < threshold_length: continue

    query_start = int(query_start)
    query_end = int(query_end)
    sbjct_start = int(sbjct_start)
    sbjct_end = int(sbjct_end)
    
    if not query in length_DICT: length_DICT[query] = [min(query_start, query_end), max(query_start, query_end)]
    if not sbjct in length_DICT: length_DICT[sbjct] = [min(sbjct_start, sbjct_end), max(sbjct_start, sbjct_end)]

    length_DICT[query][0] = min([length_DICT[query][0], query_start, query_end])
    length_DICT[query][1] = max([length_DICT[query][1], query_start, query_end])

    length_DICT[sbjct][0] = min([length_DICT[sbjct][0], sbjct_start, sbjct_end])
    length_DICT[sbjct][1] = max([length_DICT[sbjct][1], sbjct_start, sbjct_end])
fin.close()

fin = open(infile)
image_DICT = {}
ratio_DICT = {}
for line in fin:
    mylist = line.rstrip('\n').split('\t')
    if len(mylist) != 12: break
    query, sbjct, identity, alignment_length, mismatches, gap_openings, query_start, query_end, sbjct_start, sbjct_end, e_value, bit_score = mylist

    identity = float(identity)
    if identity < threshold_identity: continue

    alignment_length = int(alignment_length)
    if alignment_length < threshold_length: continue

    query_start = int(query_start)
    query_end = int(query_end)
    sbjct_start = int(sbjct_start)
    sbjct_end = int(sbjct_end)

    key = (query, sbjct)

    if not key in image_DICT: 
        svg = element('svg', None)

        svg.attr('viewBox', '0 0 ' + str(imageLength) + ' ' + str(imageLength))
        svg.attr('height', imageLength)
        svg.attr('width', imageLength)
        svg.style('background', 'white')


        



        image_DICT[key] = svg
        ratio_DICT[key] = {}

        queryMinPos, queryMaxPos = length_DICT[query]
        sbjctMinPos, sbjctMaxPos = length_DICT[sbjct]

        queryLength = queryMaxPos - queryMinPos + 1
        sbjctLength = sbjctMaxPos - sbjctMinPos + 1

        ratio_DICT[key][query] = float(imageLength - margin * 2) / max([queryLength, sbjctLength])
        ratio_DICT[key][sbjct] = float(imageLength - margin * 2) / max([queryLength, sbjctLength])
        

        #query
        x1 = margin
        x2 = margin + ratio_DICT[key][query] * queryLength
        y1 = imageLength - margin
        y2 = imageLength - margin
        line = element('line', svg)
        line.attr('x1', x1)
        line.attr('y1', y1)
        line.attr('x2', x2)
        line.attr('y2', y2)
        line.attr('stroke', 'black')
        line.attr('stroke-width', '1')

        x = margin
        y = imageLength - margin + 15
        text = element('text', svg)
        text.attr('x', x)
        text.attr('y', y)
        text.attr('font-size', '10')
        text.attr('fill', 'red')
        text.add(queryMinPos)

        x = margin + ratio_DICT[key][query] * queryLength  - 15
        y = imageLength - margin + 15
        text = element('text', svg)
        text.attr('x', x)
        text.attr('y', y)
        text.attr('font-size', '10')
        text.attr('fill', 'red')
        text.add(queryMaxPos)

        x = (margin * 2 + ratio_DICT[key][query] * queryLength  - 15)/2
        y = imageLength - margin + 30
        text = element('text', svg)
        text.attr('x', x)
        text.attr('y', y)
        text.attr('font-size', '10')
        text.attr('fill', 'red')
        text.add(query)

        #sbjct
        y1 = imageLength - margin
        y2 = imageLength - margin - ratio_DICT[key][sbjct] * sbjctLength
        x1 = margin
        x2 = margin
        line = element('line', svg)
        line.attr('x1', x1)
        line.attr('y1', y1)
        line.attr('x2', x2)
        line.attr('y2', y2)
        line.attr('stroke', 'black')
        line.attr('stroke-width', '1')

        x = margin - 50
        y = imageLength - margin
        text = element('text', svg)
        text.attr('x', x)
        text.attr('y', y)
        text.attr('font-size', '10')
        text.attr('fill', 'red')
        text.add(sbjctMinPos)

        x = margin - 50
        y = imageLength - margin - ratio_DICT[key][sbjct] * sbjctLength
        text = element('text', svg)
        text.attr('x', x)
        text.attr('y', y)
        text.attr('font-size', '10')
        text.attr('fill', 'red')
        text.add(sbjctMaxPos)

        x = margin - 60
        y = (imageLength * 2 - margin * 2 - ratio_DICT[key][sbjct] * sbjctLength) /2
        text = element('text', svg)
        text.attr('x', x)
        text.attr('y', y)
        text.attr('font-size', '10')
        text.attr('fill', 'red')
        text.add(sbjct)
        

        



    
    x1 = margin + ratio_DICT[key][query] * (query_start - length_DICT[query][0])
    y1 = imageLength - margin - ratio_DICT[key][sbjct] * (sbjct_start - length_DICT[sbjct][0])

    x2 = margin + ratio_DICT[key][query] * (query_end - length_DICT[query][0])
    y2 = imageLength - margin - ratio_DICT[key][sbjct] * (sbjct_end - length_DICT[sbjct][0])

    
    line = element('line', image_DICT[key])
    line.attr('x1', x1)
    line.attr('y1', y1)
    line.attr('x2', x2)
    line.attr('y2', y2)
    line.attr('stroke', get_color(identity))
    line.attr('stroke-width', '1')

    line.attr('onmouseover', 'console.log([' + str(query_start) + ',' + str(query_end) + ',' + str(sbjct_start) + ',' + str(sbjct_end) + ']' + ')')


for key, image in image_DICT.items():
    fout = open('-'.join(key) + ".html", 'w')
    fout.write(str(image))
    fout.close()



html = element('html', None)
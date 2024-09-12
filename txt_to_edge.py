#this program reads a txt-file from stdin and outputs edgecode
# The txt should be in the following format: 
# every vertex and every edge is represented by a unique number
# for every vertex, write a line starting with its number, then a :, and then the numbers of its incident edges in clockwise order, separated by commas (no spaces)

# An example: the octahedron
#1:4,3,2,1
#2:3,7,6,5
#3:2,5,8,11
#4:6,9,10,8
#5:7,4,12,9
#6:10,12,1,11

'''
Created on Nov 3, 2022
@authors: hvdncamp, based on txt_to_planar.py by nvcleemp
'''

import sys

def writeGraph(lines, count, zeroBased):
    nr_of_bytes = len(lines) + count - 1
    sys.stdout.flush()
    sys.stdout.buffer.write(nr_of_bytes.to_bytes(1,'big'))
    for index, line in enumerate(lines):
        _,line = line.split(':')
        for edge in [int(n) + (1 if zeroBased else 0) for n in line.split(',')]:
            sys.stdout.write('{:c}'.format(edge))
        sys.stdout.flush()
        if(index != len(lines) - 1):
            sys.stdout.buffer.write(b"\xFF")

#parse command-line arguments
zeroBased = '-0' in sys.argv

sys.stdout.write('>>edge_code<<')

lines = []
count = 0;
for line in sys.stdin:
    line = line.strip()
    if line:
        lines.append(line)
        count += line.count(",") + 1
    else:
        writeGraph(lines, count, zeroBased)
        lines=[]
        count = 0

if lines:
    writeGraph(lines, count, zeroBased)

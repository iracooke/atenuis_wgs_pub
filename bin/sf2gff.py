#!/usr/bin/env python3

import argparse
import sys
import vcf
import pdb
import math

parser = argparse.ArgumentParser(description='Convert high LR regions in sf2 output to gff')

parser.add_argument('seqid',type=str,help='Scaffold')

parser.add_argument('insf2',metavar='FILE',nargs='?',type=argparse.FileType('r'),help='Input SF2 File',default=sys.stdin)
parser.add_argument('-t','--lr-threshold',type=int,help='Keep regions with LR above this value',default=20)


args  = parser.parse_args()

start_pos = None
end_pos = None
scores = []

for line in args.insf2:

	line_values = line.split()

	if ( line_values[0] != 'location' ):
		lr = float(line_values[1])
		if ( lr > args.lr_threshold ):
			if start_pos == None:
				start_pos = float(line_values[0])

			end_pos = float(line_values[0])
			scores.append(lr)

		else:
			if ( start_pos != None ):
				sys.stdout.write(args.seqid+"\tSF2\tsweep\t"+str(math.floor(start_pos))+"\t"+
											str(math.floor(end_pos))+"\t"+str(round(max(scores),1))+
											"\t+\t.\tNote=lr_threshold"+str(args.lr_threshold)+"\n")
			start_pos = None
			scores=[]

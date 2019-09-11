#!/usr/bin/env python
import sys
from rblib import cnvproc
def parsecnvfile(cnvfile,hchr,ploidy):
	"""
	#Sample(tissue/cell)	Patient	        chrom   start           end             probes  log2ratio 
	TCGA-OR-A5LA-01A        TCGA-OR-A5LA    chr1    3301765         5128894         1270    -0.2177
	TCGA-OR-A5LA-01A        TCGA-OR-A5LA    chr1    5132643         47878718        22491   -0.179
	TCGA-OR-A5LA-01A        TCGA-OR-A5LA    chr1    161959375       247650984       54615   -0.2341
	"""
	h = {}
	f = open(cnvfile,"r")
	for line in f:
		if line.startswith("#"):continue
		sn,pn,chrom,start,end,nums,logratio = line.rstrip("\n").split("\t")
		assert chrom in hchr
		if sn not in h: h[sn] = []
		start = int(start); end = int(end)
		ratio  = 2**float(logratio) * ploidy
		fploidy = 2**float(logratio) * ploidy
		h[sn].append([chrom,start,end,ratio,ratio,fploidy])
	f.close()
	return h
def parsechrfile(chrfile):
	h = {}
	f = open(chrfile,"r")
	totlen = 0.0
	for line in f:
		arr = line.rstrip("\n").split("\t")
		tmplen = float(arr[1])
		h[arr[0]] = tmplen
		totlen += tmplen
	f.close()
	return h,totlen
def runscript(cnvfile,chrfile,ploidy,upcut,lowcut):
	hchr,totlen = parsechrfile(chrfile)
	hsn  = parsecnvfile(cnvfile,hchr,ploidy)
	sys.stdout.write("#SN\twploidy\twgii\n")
	for sn in hsn:
		segments = hsn[sn]
		cin_ins = cnvproc.CNVproc(segments,totseglen=totlen)
		wploidy = cin_ins.wploidy()
		#print 2**upcut * ploidy,2**lowcut * ploidy
		wgii    = cin_ins.wgii(tgain=2**upcut * ploidy,tloss=2**lowcut * ploidy,chrlen=hchr,numchrs=len(hchr))
		sys.stdout.write("%s\t%.5f\t%.5f\n"%(sn,wploidy,wgii))
	return 0

from optparse import OptionParser,OptionGroup
import time
import os

def checkfile(fns):
	for fn in fns:
		if not os.path.isfile(fn): return 2
	return 0

def __main():
	usage = "usage: %prog CNVfile"
	description = "Contact: Rong Zhengqin <rongzhengqin@basepedia.com>"
	parser = OptionParser(usage,version="%prog 1.0.0",description = description)
	Required_group = OptionGroup(parser,'Required Options')
	Required_group.add_option('-r',dest='chromfile',help="chromosome size file",metavar='FILE',type='string',default=None)
	Required_group.add_option('-p',dest='ploidy',help="Species ploidy [2]",metavar='INT',type='int',default=2.0)

	Other_group    = OptionGroup(parser,'Threshold Options')
	Other_group.add_option('-l',dest='low',help="Low cutoff for loss fragments",metavar='FLOAT',type='float',default = -0.4)
	Other_group.add_option('-u',dest='up',help="Up cutoff for gain fragments",metavar='FLOAT',type='float',default = 0.4)

	parser.add_option_group(Required_group)
	parser.add_option_group(Other_group)

	(options, args) = parser.parse_args()
	if len(args) < 1:
		parser.print_help()
		return -1

	chrfile = str(options.chromfile)
	ploidy  = float(options.ploidy)
	lowcut  = float(options.low)
	upcut   = float(options.up)
	cnvfile = args[0]
	for fn in [cnvfile,chrfile]:
		if not os.path.isfile(fn):
			sys.stderr.write("file '%s' not found!"%fn)
			return 2
	ret = runscript(cnvfile,chrfile,ploidy,upcut,lowcut)
	return ret


if __name__ == "__main__":
	start_time = time.time()
	ret        = __main()
	cost_time  = time.time()-start_time
	if ret: sys.stderr.write("[ERROR] Task interrupt, Code: %d\n"%ret)
	else: sys.stderr.write("[INFO] Time consumed: %.2fs, Code: %d\n"%(cost_time,ret))
	exit(ret)


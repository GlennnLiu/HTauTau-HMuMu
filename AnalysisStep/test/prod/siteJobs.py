import os,glob

log="9120328"
chunks=glob.glob("*Chunk*")

def getFile(run):
	for l in open(run).readlines():
		if "fileNames = cms.untracked.vstring(" in l:
			return l.split("'")[1]

remaining=[]
for c in chunks:
	if os.path.exists("{0}/log/{1}.log".format(c,log)):
		os.system('grep "terminated" {0}/log/{1}.log > tmp.txt'.format(c,log))
		if len(open('tmp.txt').readlines())==0:
			continue
	#os.system('rm {0}/log/{1}*'.format(c,log))
	f=getFile("{0}/run_cfg.py".format(c))
	if "/store/test/xrootd/" in f:
		continue
	os.system('dasgoclient -query="site file={0}" > tmp.txt'.format(f))
	site=""
	for l in open('tmp.txt').readlines():
		if "Tape" in l:
			continue
		site=l.split("_Disk")[0]
		break
	if site=="":
		print "NO SITE!! {0}".format(c)
	else:
		site=site.rstrip("\n")
		print site
		os.system('sed -i "s#/store#/store/test/xrootd/{0}/store#g" {1}/run_cfg.py'.format(site,c))
		remaining.append(c)

print " ".join(remaining)

import os,glob

chunks=glob.glob("*Chunk*")

for c in chunks:
	print c
	cfg=open("{0}/run_cfg.py".format(c))
	for l in cfg.readlines():
		if "fileNames = cms.untracked.vstring(" in l:
			allFile=l
			break
	allFile=allFile.split("(")[-1].split(")")[0]
	files=allFile.split(',')
	cfg.close()
	for i,f in enumerate(files):
		os.system("cp -r {0} {0}_{1}".format(c,i))
		os.system('sed -i "s#{0}#{1}#" {2}_{3}/run_cfg.py'.format(allFile,f,c,i))
	os.system("rm -r {0}".format(c))
	os.system("rm -r /eos/home-g/geliu/LepUni/BigTrees/UL/2016/{0}".format(c))
	

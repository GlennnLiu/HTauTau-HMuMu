import os,glob

chunkList=glob.glob("*Chunk*")
eosPath="/eos/home-g/geliu/LepUni/BigTrees/UL/2016/"

for chunk in chunkList:
	remove=False
	if os.path.exists("{0}{1}/log.txt".format(eosPath,chunk)):
		os.system('grep "T---Report end!" {0}{1}/log.txt > tmp.txt'.format(eosPath,chunk))
	elif os.path.exists("{0}{1}/log.txt.gz".format(eosPath,chunk)):
                os.system('zgrep "T---Report end!" {0}{1}/log.txt.gz > tmp.txt'.format(eosPath,chunk))
	else:
		remove=True
	if len(open("tmp.txt").readlines())==0:
		remove=True
	if remove:
		os.system('rm {0}{1}/*root'.format(eosPath,chunk))
		print 'rm {0}{1}/*root'.format(eosPath,chunk)	

import os,glob

i=1
log="9468244"
f=open("leftJobs.txt","w")
for l in glob.glob("*Chunk*"):
	os.system("grep 'terminated' {0}/log/{1}.log > tmp.txt".format(l,log))
	nter=len(open("tmp.txt").readlines())
	
	nerr=0
	err=glob.glob("{0}/log/{1}*err".format(l,log))
	if len(err)>0:
		os.system("grep 'cannot create' {0}/log/{1}*err > tmp.txt".format(l,log))
		nerr=len(open("tmp.txt").readlines())
		#if nerr>0:
		#	continue
	os.system("{1} 'Successful' /eos/home-g/geliu/LepUni/BigTrees/UL/2016/{0}/log* > tmp.txt".format(l,"zgrep" if not os.path.exists("/eos/home-g/geliu/LepUni/BigTrees/UL/2016/{0}/log.txt".format(l)) else "grep"))
	nsucc=len(open("tmp.txt").readlines())
	os.system("{1} 'T---Report end!' /eos/home-g/geliu/LepUni/BigTrees/UL/2016/{0}/log* > tmp.txt".format(l,"zgrep" if not os.path.exists("/eos/home-g/geliu/LepUni/BigTrees/UL/2016/{0}/log.txt".format(l)) else "grep"))
        nfinish=len(open("tmp.txt").readlines())
	os.system("{1} 'FileReadError' /eos/home-g/geliu/LepUni/BigTrees/UL/2016/{0}/log* > tmp.txt".format(l,"zgrep" if not os.path.exists("/eos/home-g/geliu/LepUni/BigTrees/UL/2016/{0}/log.txt".format(l)) else "grep"))
        os.system("{1} 'FileOpenError' /eos/home-g/geliu/LepUni/BigTrees/UL/2016/{0}/log* >> tmp.txt".format(l,"zgrep" if not os.path.exists("/eos/home-g/geliu/LepUni/BigTrees/UL/2016/{0}/log.txt".format(l)) else "grep"))
	os.system("{1} 'A fatal system signal has occurred' /eos/home-g/geliu/LepUni/BigTrees/UL/2016/{0}/log* >> tmp.txt".format(l,"zgrep" if not os.path.exists("/eos/home-g/geliu/LepUni/BigTrees/UL/2016/{0}/log.txt".format(l)) else "grep"))
	nerror=len(open("tmp.txt").readlines())
	#print(l,nter,nerror,os.path.exists("/eos/home-g/geliu/LepUni/BigTrees/UL/2016/{0}/log.txt".format(l)))
	if nter>0 and nsucc>0 and nfinish>0 and nerr==0 and nerror==0:
		print("{0} {1} finished".format(i,l))
		i=i+1
		os.system("mv {0} /eos/home-g/geliu/LepUni/BigTrees/UL/2016/CondorMC/".format(l))
	print(l,nter,nsucc,nerr if len(err)>0 else 0,nerror)
	if nter>0 and nerror>0:
		f.write("{0} has file open/reading errors.\n".format(l))
	if nter==0 or nsucc==0:
		f.write("{0} is not finished.\n".format(l))


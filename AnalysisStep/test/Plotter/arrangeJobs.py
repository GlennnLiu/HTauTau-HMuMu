import os,glob

regions=["SRTree","CRQCDTree","CRWJTree","CRQCDvSRTree","CRQCDvSROSTree","CRQCDvSRSSTree","CRAPPOSTree","CRAPPSSTree"]
os.system("mkdir Condor")

for re in regions:
	os.system("mkdir Condor/{0}".format(re))
	os.system("cp DistributionPlots.sh Condor/{0}".format(re))
	os.system("mkdir Condor/{0}/log".format(re))
	os.system('sed -i "s#DistributionPlots --region#DistributionPlots --region {0} --addData#" Condor/{0}/DistributionPlots.sh'.format(re))
	os.system('sed -i "s#Plots/Plots#Plots/{0}.Plots#" Condor/{0}/DistributionPlots.sh'.format(re))

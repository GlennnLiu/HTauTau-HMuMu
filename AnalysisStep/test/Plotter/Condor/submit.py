import os,glob

f=glob.glob("*Tree")

print "condor_submit condor.sub -queue directory in "+" ".join(f)

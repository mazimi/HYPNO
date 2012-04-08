#!/opt/python/bin/python2.7

import subprocess, os

class GenTree:

    # __init__: Constructor method
    # input: self
    # output: none
    def __init__( self ):
        pass

    def makeTree( self, MSA , name):
    	print 'got into MakeTree'
    	flag = 0
    	if flag == 1:
    		os.system("FastTree -gtr -nt " + MSA + " > " + name)
    	else: 
    		subprocess.call(["FastTree -gtr -nt " + MSA + " > " + name], shell=True)

if __name__ == '__main__':
	sys.exit(1)
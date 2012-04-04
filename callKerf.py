#!/opt/python/bin/python2.7

from kerf import Kerf
import sys

# main: Assign protein MSA (a2m), DNA sequences (FASTA) and aligned DNA output file (FASTA)
# input: none
# output: none
def main():
	fileMSA = "MSA.a2m"
	fileTree = "tree.ml"
	threshold = 80.0
	mainTree = Kerf()
	mainTree.kerfRun(fileTree, fileMSA, threshold)

if __name__ == '__main__':
  main()
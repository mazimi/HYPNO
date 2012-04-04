#!/usr/bin/env python

from Bio import Phylo # from Biopython
from Bio import SeqIO # from Biopython
import sys



class Kerf:

    # __init__: Constructor method
    # input: self
    # output: none
    def __init__( self ):
        pass


    def kerfRun( self, tree2kerf, MSA2kerf, kerfThresh, outputDir ):

        global thresh  # some variable declarations
        global path
        global done
        global subfamily
        global csvBuffer
        global MSA
        global tree
        global node_list
        ###########################################################################
        #
        # calc_minPID - This function calculates the minimum pairwise percent 
        #           identity for all leaves in the subtree pointed to by the parameter 
        #           clade.    
        #          
        ###########################################################################

        def calc_minPID( clade ):
            mini = 100
            terminal_list = clade.get_terminals()
            for iseq1 in range(len(terminal_list)-1):
                for iseq2 in range(iseq1+1, len(terminal_list)):
                    # do a pairwise comparison bewteen all leaves
                    seq1 = getSeq(terminal_list[iseq1].name)
                    seq2 = getSeq(terminal_list[iseq2].name)
                    # calculate percent identity
                    PID = calc_PID(seq1,seq2)
                    # update mini if necessary
                    if (PID < mini):
                        mini = PID
            return mini


        ###########################################################################
        #
        # calc_PID - This function calculates the percent identity as I understand
        #            it to be.  Specifically, it searches for matches between the
        #            seq1 and seq2 at each position along their lengths.  If the
        #            matches aren't the '-' character, it scores an identity.  The
        #            percent ID is calculated and returned after iteration through
        #            the entire string.
        #
        ###########################################################################


        def calc_PID( seq1,seq2 ):
            identity = 0  # number of identities
            num_align = 0 # number of single alignments
            if (len(seq1) == len(seq2)):  
                for amino_acid in range(len(seq1)):
                    if (seq1[amino_acid] == seq2[amino_acid]):
                        if (seq1[amino_acid] != '-'):
                        # if both sequences match, but aren't both -
                            identity+=1
                            num_align+=1
                        else:
                            num_align+=1
                return (100.0*identity)/num_align # return % ID
            else:
                return -1 # sequences don't match in length


        ###########################################################################
        #
        # getSeq - This function takes a string as an input parameter and attempts to 
        #           return the corresponding sequence whose name matches the string. 
        #           This is used to tie the tree nodes to their corresponding MSA seq.    
        #          
        ###########################################################################

        def getSeq( rec ):
            global MSA
            for record in MSA:
                if rec in record.name:
                    return record.seq

        ###########################################################################
        #
        # getParent - This function returns the parent node of the parameter clade. 
        #         If no parent exists, -1 is returned.   
        #               
        ###########################################################################

        def getParent( clade ):
            global node_list
            foundParent = False
            for node in node_list:
                if clade in node:
                    foundParent = True
                    return node
            if not foundParent:
                return -1


        ###########################################################################
        #
        # traverse - This function recursively traverses through the tree, as
        #       as described in the problem statement    
        #
        ###########################################################################

        def traverse( clade ):
            global thresh

            if (getParent( clade ) != -1):
                if (calc_minPID(getParent(clade)) < thresh):
                    # this clade is the root of a subfamily, write to file
                    writeClade(clade)
                else:
                    traverse( getParent(clade) )
            else:
                # clade is root, all nodes are in the same subfamily
                writeClade(clade)
            return 0


        ###########################################################################
        #
        # writeClade - This function calculates the minimum pairwise percent 
        #           identity for all leaves in a subtree pointed to by the parameter 
        #           clade.    
        #          
        #            
        #            
        #
        ###########################################################################

        def writeClade( clade ):

            # initialize
            global done
            global subfamily
            global MSA
            global path
            global csvBuffer

            terminal_list = clade.get_terminals()
            writeBuffer = []
            subfamily+=1  # got a new subfamily
            for leaf in terminal_list:
                for record in MSA:
                    if leaf.name in record.name:
                        # we found the record, mark it done and update the 2 write 
                        done.append(leaf.name)  # buffers
                        writeBuffer.append(record)
                        csvBuffer+='%d, %s\n' % (subfamily, leaf.name)
            # now write the subfamily to file
            filename = path.split('.')[0]+'sf%d.' % subfamily + path.split('.')[-1]
            f = open(outputDir + '/' + filename, 'w')
            SeqIO.write(writeBuffer, f, 'fasta')
            f.close()
            return 0

        ###########################################################################
        #
        # isDone - This function returns a Boolean corresponding to if node is in
        #           the global done list.       
        #
        ###########################################################################

        def isDone( node ):
            global done
            for nodes in done:
                if (node.name == nodes):
                    return True
            return False


        ###############################################################
        #
        #           Inputs
        #
        ###############################################################

        #treepath = raw_input('path to tree --> ')
        #MSApath = raw_input('path to MSA --> ')
        #thresh = float(raw_input('threshold --> '))
        treepath = tree2kerf
        MSApath = MSA2kerf
        thresh = kerfThresh

        ###############################################################
        #
        #           Variable initializations
        #
        ###############################################################

        path = MSApath
        done=[]
        subfamily=0
        csvBuffer = ''

        # open and read MSA file
        f = open(outputDir + '/' + MSApath)
        MSA = list(SeqIO.parse(f,'fasta'))
        f.close()

        # open and read tree
        tree = Phylo.read(outputDir + '/' + treepath,'newick')

        leaves = tree.get_terminals()
        node_list = tree.get_terminals() + tree.get_nonterminals()

        # loop  while there are still leaves that aren't done
        #print "%d\n" % len(leaves)

        while (len(done) < len(leaves)):  #loop termination
            for leaf in leaves:
                if not isDone(leaf):    # if a leaf isn't done, traverse it 
                    traverse(leaf)

        # Now write the CSV file

        fn = path.split('.')[0] + '.csv'
        f = open(outputDir + '/' + fn, 'w')
        f.write(csvBuffer)
        f.close()

        sys.exit(0)

if __name__ == '__main__':
    sys.exit(1)
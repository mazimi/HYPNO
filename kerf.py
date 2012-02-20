#!/opt/python/bin/python2.7
"""Module to cut a tree into subtrees to define subfamilies"""

from __future__ import division
import argparse
import csv
import itertools
import logging
from Bio import AlignIO, Phylo
from Bio.Align import MultipleSeqAlignment


class MSATree:
    """An MSA with a Phylo Tree"""
    def __init__(self, msa, tree):
        self.msa = msa
        self.tree = tree

        def get_key(align):
            """Gets a key for an alignment

            Given the Key: tr|B3SB47|B3...
            Returns: B3SB47

            The join logic is definitely not universal,
            but works with bpg087857

            It looks like ETE2 forces the sequence IDs to be equal
            to the tree node names"""
            # makes a predictable lookup key for our msa

            return align.name.split('|')[1]

        self.msa_by_name = dict((get_key(alignment), alignment)
                                for alignment in self.msa)

        #place to store the phylo sub trees for output
        self.sub_trees = []

        # stolen from http://biopython.org/wiki/Phylo_cookbook
        self.parents = {}
        for clade in self.tree.find_clades(order='level'):
            for child in clade:
                self.parents[child] = clade

        self._prune_tree()

        # Flag for if we add the root to the output trees
        self.added_root = False

        # Cache for sequence distances
        # memoizing with a decorator might be more pythonic, but this
        # makes it easy to use the names as the key instead of hashing
        self.distance_cache = {}

    def dump(self, output_dir):
        """Dumps the split contents into files

        returns nothing"""

        # Files:
        #   outputdir/msa_map.csv: same size as the input msa,
        #     but has the output msa number
        #   outputdir/n.a2m where n is the msa number - msa files
        #   outputdir/n.nf where n is the msa number - output trees

        def get_alignment_from(tree):
            """Fetches the MSA for a Phylo Tree"""
            msa = []
            for node in tree.get_terminals():
                alignment = self.msa_by_name[node.name.split('|')[1]]
                if msa:
                    msa.append(alignment)
                else:
                    msa = MultipleSeqAlignment([alignment])

            return msa

        msa_map = csv.writer(open(output_dir + '/msa_map.csv', 'wb'))

        for index, tree in enumerate(self.sub_trees):
            Phylo.write(tree, output_dir + "/" + str(index) + ".nh", "newick")
            AlignIO.write(
                get_alignment_from(tree),
                output_dir + "/" + str(index) + ".a2m",
                "fasta")
            for node in tree.get_terminals():
                msa_map.writerow([index, node.name])

    def _prune_tree(self):
        """Prunes nodes that don't have a sequence in the MSA
        Returns Nothing"""

        for leaf in self.tree.get_terminals():
            if not self.msa_by_name.get(leaf.name.split('|')[1]):
                logging.warning("Pruning: " + leaf.name)
                self.tree.prune(leaf)

    def _node_distance(self, first, second):
        """Takes in two clade nodes

        Returns the pairwise identity percentage

        Nodes are remembered by name, and the results cached
        """

        name_1 = first.name.split('|')[1]
        name_2 = second.name.split('|')[1]

        if self.distance_cache.get(name_1 + name_2):
            return self.distance_cache[name_1 + name_2]
	
        seq1 = self.msa_by_name[name_1]
        seq2 = self.msa_by_name[name_2]

        distance = self._seq_distance(seq1, seq2)
        self.distance_cache[name_1 + name_2] = distance
        return distance

    # This is split out to make testing sane
    # might still be an annoying performance hit
    def _seq_distance(self, seq1, seq2):
        """Takes two sequence strings

        Returns the pairwise identity percentage

        From: http://openwetware.org/wiki/Wikiomics:Percentage_identity
        'there are many different ways of calculating percentage identity'

        From the Spec:
        Note that percent identity is defined as: #exact matches of amino acids
        / #positions where at least one sequence aligns. (i.e., do not consider
        positions where both sequences have gaps).

        Sample MSA from The Spec:
        MSTPP----W
        -TTPPPPP-W

        Based on the spec's definition, there are 9 spots where one of the
        sequences align, and 4 identities, so the example seems to be based
        on the number of positions being the longest sequence string after the
        gaps have been removed

        This method implements the written description from the spec, and not
        the example.
        """

        num_positions = 0  # the number of positions with at least one amino
        num_matches = 0    # the number of matches seen

        for pair in zip(seq1, seq2):
            if not (pair[0] == "-" and pair[1] == "-"):
                num_positions += 1

                if (pair[0] == pair[1]):
                    num_matches += 1

        return num_matches / num_positions * 100

    def subtree_distances(self, root):
        """Takes a node

        returns a list of the percent identities between all leaves

        This isn't needed for the program to run, but can help with debugging
        The other method is faster, as it only computes distances until one
        exeeds the threshold
        """

        nodes = root.get_terminals()
        nodes.reverse()
        node_pairs = itertools.ifilter(
            lambda (a1, a2): a1.name < a2.name,
            itertools.product(nodes, nodes))

        distances = [self._node_distance(pair[0], pair[1])
            for pair in node_pairs]

        return distances

    def _subtree_above_threshold(self, root, threshold=0):
        """Takes in a clade node and a threshold

        Returns True if all children have a pairwise sequence
        alignment below the threshold"""

        # BFS will grab the most distant leaves last
        # reverse should guarantee that we get the smallest subtrees
        nodes = root.get_terminals()
        nodes.reverse()

        # always making n1 < n2 so that we can memoize
        # this all could be done in parallel
        node_pairs = itertools.ifilter(
            lambda (a1, a2): a1.name < a2.name,
            itertools.product(nodes, nodes))

        for pair in node_pairs:
            distance = self._node_distance(pair[0], pair[1])
            if distance and distance < threshold:
                return False

        return True

    def _subtree_from(self, clade):
        """Creates a subtree for output and removes it from the main tree

        returns nothing"""
        parent = self.parents.get(clade)

        for child in clade.clades:
            # Seems odd to waste the copy, but saves mucking about in its
            # internals.
            self.sub_trees.append(Phylo.BaseTree.Tree.from_clade(child))

        if parent:
            parent.clades.remove(clade)

        if clade == self.tree.root:
            self.added_root = True

    def create_subtrees(self, threshold):
        """Breaks the main tree into sub trees based on the threshold

        returns self"""

        nodes = self.tree.get_nonterminals(order='level')
        nodes.reverse()

        for node in nodes:
            if not self._subtree_above_threshold(node, threshold):
                self._subtree_from(node)

        #If we have to split on the root node, we do it.
        if not self.added_root:
            if not self._subtree_above_threshold(self.tree.root, threshold):
                self._subtree_from(self.tree.root)
            else:
                self.sub_trees.append(self.tree.root)

        return self

if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description='Method to cut a tree into subtrees to define subfamilies')
    parser.add_argument('--msa-file', metavar='file', help='msa file',
        type=argparse.FileType('r'), required=True)
    parser.add_argument('--tree-file', metavar='file', help='tree file',
        type=argparse.FileType('r'), required=True)
    parser.add_argument('--output-directory', metavar='dir',
        help='output directory', required=True)
    parser.add_argument('--threshold', metavar='float', help='threshold',
        type=float, required=True)

    args = parser.parse_args()

    msa_tree = MSATree(
        AlignIO.read(args.msa_file, "fasta"),
        Phylo.read(args.tree_file, "newick"))

    msa_tree.create_subtrees(args.threshold)
    msa_tree.dump(args.output_directory)
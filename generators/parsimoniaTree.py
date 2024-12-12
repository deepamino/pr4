from Bio.Phylo.TreeConstruction import ParsimonyScorer, NNITreeSearcher, ParsimonyTreeConstructor
from Bio import Phylo

class ParsimoniaTree:

    @staticmethod
    def generate_tree(sequences):
        scorer = ParsimonyScorer()
        searcher = NNITreeSearcher(scorer)
        constructor = ParsimonyTreeConstructor(searcher)
        return constructor.build_tree(sequences)

    @staticmethod
    def show_tree(tree):
        return Phylo.draw(tree)


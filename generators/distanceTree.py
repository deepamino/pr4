from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio import Phylo

class DistanceTree:

    @staticmethod
    def generate_tree_upgma(sequences, method):
        calculator = DistanceCalculator(method)
        dm = calculator.get_distance(sequences)
        constructor = DistanceTreeConstructor()
        return constructor.upgma(dm)

    @staticmethod
    def generate_tree_nj(sequences, method):
        calculator = DistanceCalculator(method)
        dm = calculator.get_distance(sequences)
        constructor = DistanceTreeConstructor()
        return constructor.nj(dm)
    
    @staticmethod
    def show_tree(tree):
        return Phylo.draw(tree)
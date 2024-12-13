from .parsimoniaTree import ParsimoniaTree
from .distanceTree import DistanceTree

class GeneratorTreeFactory:

    __generators = {
        'ParsimoniaTree': ParsimoniaTree(),
        'DistanceTree': DistanceTree
    }

    @staticmethod
    def initialize_generator(key):
        identifier_class = GeneratorTreeFactory.__generators.get(key)
        if identifier_class is not None:
            return identifier_class
        return None
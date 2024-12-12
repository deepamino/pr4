from .parsimoniaTree import ParsimoniaTree

class GeneratorTreeFactory:

    __generators = {
        'ParsimoniaTree': ParsimoniaTree()
    }

    @staticmethod
    def initialize_generator(key):
        identifier_class = GeneratorTreeFactory.__generators.get(key)
        if identifier_class is not None:
            return identifier_class
        return None
#!/usr/bin/env python

import networkx

class braidedTree:

    def __init__(self, n_levels, alphabet):

        self.alphabet = set(alphabet)
        self.n_levels = n_levels
        self.levels = []

        for level in range(n_levels):
            self.add_level(self)

    def add_level(self):
        self.levels.append([])
        for node in self.levels[-2]:
            self.levels.extend(node.generate_children())
        

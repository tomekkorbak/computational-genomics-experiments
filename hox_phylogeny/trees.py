class Node(object):
    
    def __init__(self, children=[], species=None):
        self.children = children
        self.species = species
        if self.species is not None:
            self.species = set([self.species])
        self.parent = None
        for child in children:
            child.parent = self
            
    def __repr__(self):
        return f'Node(species={self.species})'
    
    def get_species(self):
        if self.is_terminal:
            return set([self.species])
        else:
            _species = set()
            for child in self.children:
                 _species.update(child.get_species())
            return _species
        
    def get_path_to_root(self):
        if self.is_root:
            return [self]
        else:
            path = [self]
            while True:
                path.append(path[-1].parent)
                if path[-1].is_root:
                    return path
        
    def get_path_length(self, other):
        if self == other:
            return 0
        self_to_root = self.get_path_to_root()
        other_to_root = other.get_path_to_root()
        if self.is_root:
            return len(other_to_root)
        if other.is_root:
            return len(self_to_root)
        for i in range(1, max(len(self_to_root), len(other_to_root))):
            if self_to_root[-i] != other_to_root[-i]:
                return len(self_to_root) + len(other_to_root) - 2*i
        return 0

    def get_terminals(self):
        if self.is_terminal:
            return [self]
        else:
            leaves = []
            for child in self.children:
                leaves += child.get_terminals()
            return leaves
        
    def get_nonterminals(self):
        nonterminals = [self]
        for child in self.children:
            if not child.is_terminal:
                nonterminals += child.get_nonterminals()
        return nonterminals
    
    @property
    def is_root(self):
        return self.parent is None
    
    @property
    def is_terminal(self):
        return self.children == []


def find_least_common_ancestor(species_tree, gene_tree, source_gene_node):
    paths = []
    for leaf in species_tree.get_terminals():
        path_to_species_root = leaf.get_path_to_root()
        for i, node in enumerate(path_to_species_root):
            if source_gene_node.species.issubset(node.species):
                paths.append(path_to_species_root[:i+1])
    return sorted(paths, key=lambda path: len(path), reverse=True)[-1][-1]


def reconcile(species_tree, gene_tree):
    for node in gene_tree.get_nonterminals():
        lca = find_least_common_ancestor(species_tree, gene_tree, node)
        node.mapped_species_node = lca
        print(node, '-->', lca)
    species_leaves = {list(node.species)[0]: node 
                      for node in species_tree.get_terminals()}
    for node in gene_tree.get_terminals():
        node.mapped_species_node = species_leaves[list(node.species)[0]]


def compute_deep_coalescence(species_tree, gene_tree):
    nodes = set(gene_tree.get_terminals() + gene_tree.get_nonterminals())
    score = 0
    for node in nodes:
        if node != gene_tree.root:
            species_node_a = node.mapped_species_node
            species_node_b = node.parent.mapped_species_node
            score += species_node_a.get_path_length(species_node_b)
    return score


def label_terminal_node_with_species(node):
    # Monkey-patching a species and class labels
    node.species, _, _, node.family = node.name.split('|')
    node.species = set([node.species])

    
def label_nonterminal_clade_with_species(node):
    # Monkey-patching a set of species of node's leaves
    species = set()
    for leaf in node.get_terminals():
        species |= leaf.species
    node.species = species


def label_with_species(tree):
    for leaf in tree.root.get_terminals():
        label_terminal_node_with_species(leaf)
    
    for clade in tree.root.get_nonterminals():
        label_nonterminal_clade_with_species(clade)

   
def get_parent(node, tree):
    # A hacky way to get parent of a node in a Bio.Phylo tree object
    return tree.root.trace(node, tree.root)[0]


def label_with_parents(tree):
    for node in tree.get_terminals() + tree.get_nonterminals():
        if node != tree.root:
            node.parent = get_parent(node, tree)


def find_duplications(gene_tree):
    duplications = set()
    for node in gene_tree.get_nonterminals():
        if node == gene_tree.root:
            continue
        if node.mapped_species_node == node.parent.mapped_species_node:
            duplications.add(node.parent)
    return duplications


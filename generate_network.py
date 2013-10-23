#!/usr/bin/env python

from __future__ import division

from argparse import ArgumentParser
from os import makedirs
from os.path import isdir, join
from shutil import rmtree
from collections import defaultdict

parser = ArgumentParser()
parser.add_argument('-r', '--reactions', help='The reactions file from '
    'which the output edges and nodes files will be generated.',
    required=True)
parser.add_argument('-c', '--cofactors', help='The cofactors file from which '
    'the list of cofacors will be read. One cofactor per line.')
parser.add_argument('-o', '--output', help='The path to the output '
    'directory. The directory must not already exist.', required=True)
parser.add_argument('-f', '--force', help='Force replacement of output '
    'directory, even if it already exists.', action="store_true")
parser.add_argument('-s', '--strip_localizations', help='Strips off cellular '
    'localizations [c], [p], and [e].', action="store_true")

# FLAGS
REACTANT = 1
PRODUCT = 2
SAME_SIDE = 1
DIFFERENT_SIDE = 2

NODE_COLOR_DICT = {
    REACTANT: 'cyan',
    PRODUCT: 'yellow',
    REACTANT | PRODUCT: 'green'}

EDGE_WEIGHT_DICT = {
    SAME_SIDE: 1,
    DIFFERENT_SIDE: 10,
    SAME_SIDE | DIFFERENT_SIDE: 1}

EDGE_COLOR_DICT = {
    # side flag, reversible
    (SAME_SIDE, True): 'black',
    (DIFFERENT_SIDE, True): 'blue',
    (SAME_SIDE | DIFFERENT_SIDE, True): 'blue',

    # side flag, not reversible
    (SAME_SIDE, False): 'black',
    (DIFFERENT_SIDE, False): 'red',
    (SAME_SIDE | DIFFERENT_SIDE, False): 'red'}

def remove_localizations(input_string):
    """Removes cellular localizations from a string

    Replaces [c], [p], and [e] with the empty string
    """
    return input_string.replace('[c]','').replace('[e]','').replace('[p]','')

def parse_cofactors(cofactors_fp, strip_localizations):
    """Parses a list of cofactors from a file

    Returns a set of cofactors
    """
    cofactors = []
    for line in open(cofactors_fp, 'U'):
        if strip_localizations:
            line = remove_localizations(line)
        line = line.strip()
        cofactors.append(line)

    return set(cofactors)

def iter_reactions(reactions_fp, cofactors, strip_localizations):
    """Iterates over a reaction file, returning a tuples for each reaction

    Each returned tuple contains three items, a set of reactants on the left,
    a set of reactants on the right, and a boolean "reversible"

    cofactors will be removed from all reactions. Reactions will NOT be
    returned if either the right side or the left side is empty after removal
    of cofactors.
    """
    for line in open(reactions_fp, 'U'):
        if strip_localizations:
            line = remove_localizations(line)
        line = line.strip()
        reversible = '<=>' in line
        divider = ' <=> ' if reversible else ' -> '
        left, right = line.split(divider)
        left_reactants = left.split(' + ')
        right_reactants = right.split(' + ')

        # filter out cofactors from the reaction
        filter_fn = lambda x: x not in cofactors
        left_reactants = filter(filter_fn, left_reactants)
        right_reactants = filter(filter_fn, right_reactants)

        # if there are reactants remaining on both sides of the equation, yield
        # the reaction
        if left_reactants and right_reactants:
            yield (set(left_reactants), set(right_reactants), reversible)

            # also yield the reverse reaction if this reaction is reversible
            if reversible:
                yield (set(right_reactants), set(left_reactants), reversible)

def generate_network(reactions_fp, cofactors_fp, output_dir,
    strip_localizations):
    """Generates edge and node attribute dicts from an input list of reactions.

    Given reactions in the format:

    etoh[c] + nad[c] <=> acald[c] + h[c] + nadh[c]

    Returns 5 dicts:
        edge_colors             Keys are frozensets of pairs of metabolites,
                                and values represent colors that can be
                                looked up in the EDGE_COLOR_DICT.
        
        edge_weights            Keys are as above, and values are represent
                                weights that can be looked up in the
                                EDGE_WEIGHT_DICT

        edge_reversibilities    Keys are as above, and values are booleans,
                                True if the pair are ever observed on opposite
                                sides of a reversible reaction, False if not.

        node_colors             Keys are metabolite names, values represent
                                colors that can be looked up in the
                                NODE_COLOR_DICT

        node_sizes              Keys are metabolites, values record how many
                                different reactions the metabolite is observed
                                in.
    """
    edge_colors = defaultdict(int)
    edge_weights = defaultdict(int)
    edge_reversibilities = defaultdict(bool)
    node_colors = defaultdict(int)
    node_sizes = defaultdict(int)

    cofactors = parse_cofactors(cofactors_fp, strip_localizations)

    for left, right, reversible in iter_reactions(reactions_fp, cofactors,
            strip_localizations):

        # update node color and size trackers
        for metabolite in left:
            node_colors[metabolite] |= REACTANT
            node_sizes[metabolite] += 1
        for metabolite in right:
            node_colors[metabolite] |= PRODUCT
            node_sizes[metabolite] += 1

        # convert to lists to do index-based metabolite pairing
        left, right = list(left), list(right)

        # update edge color and size trackers
        # same side, on the left
        for i, first in enumerate(left[:len(left)-1]):
            for second in left[i+1:]:
                key = frozenset([first, second])
                edge_colors[key] |= SAME_SIDE
                edge_weights[key] |= SAME_SIDE

        # same side, on the right
        for i, first in enumerate(right[:len(right)-1]):
            for second in right[i+1:]:
                key = frozenset([first, second])
                edge_colors[key] |= SAME_SIDE
                edge_weights[key] |= SAME_SIDE

        # different sides
        for l in left:
            for r in right:
                key = frozenset([l, r])
                edge_colors[key] |= DIFFERENT_SIDE
                edge_weights[key] |= DIFFERENT_SIDE
                edge_reversibilities[key] |= reversible

    return (edge_colors, edge_weights, edge_reversibilities, node_colors, 
        node_sizes)

def write_node_table(node_colors, node_sizes, output_fp):
    """Writes the node table

    The node table has three columns: node_name, node_size, and node_color.
    The node_name is the name of the metabolite, and the size and color are
    self-explanatory fields.
    """
    node_table = open(output_fp, 'w')
    node_cols = ['node_name', 'node_size', 'node_color']
    node_table.write('\t'.join(node_cols) + '\n')

    for metabolite, size in node_sizes.iteritems():
        color = NODE_COLOR_DICT[node_colors[metabolite]]
        size = str(size)
        node_table.write('\t'.join([metabolite, size, color]) + '\n')
    node_table.close()

def write_edge_table(edge_colors, edge_weights, edge_reversibilities,
        output_fp):
    """
    """
    # write edge table
    edge_table = open(output_fp, 'w')
    edge_cols = ['reactant', 'product', 'distance', 'color']
    edge_table.write('\t'.join(edge_cols) + '\n')

    for reactant_product, weight in edge_weights.iteritems():
        color = EDGE_COLOR_DICT[(edge_colors[reactant_product],
            edge_reversibilities[reactant_product])]
        weight = str(EDGE_WEIGHT_DICT[weight])
        reactant, product = reactant_product
        edge_table.write('\t'.join([reactant, product, weight, color]) + '\n')
    edge_table.close()

if __name__ == '__main__':
    args = parser.parse_args()
    reactions_fp = args.reactions
    cofactors_fp = args.cofactors
    output_dir = args.output
    force = args.force
    strip_localizations = args.strip_localizations

    if isdir(output_dir):
        if force:
            rmtree(output_dir)
        else:
            raise IOError, ("Output directory already exists. Please remove "
                "it before continuing. You may pass -f or --force to omit "
                "this check.")

    makedirs(output_dir)

    (edge_colors,
    edge_weights,
    edge_reversibilities,
    node_colors, 
    node_sizes) = generate_network(reactions_fp, cofactors_fp,
        output_dir, strip_localizations)

    write_node_table(node_colors, node_sizes,
        join(output_dir, 'node_table.tsv'))
    write_edge_table(edge_colors, edge_weights, edge_reversibilities,
        join(output_dir, 'edge_table.tsv'))

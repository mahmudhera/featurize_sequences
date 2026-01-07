import os
import sys

sys.path.append('/home/jl2791/scripts')
import featurizer
from featurizer_util import *
from featurizer_paths import *

example_dir = os.path.dirname(os.path.abspath("__file__"))

Features = ['deepbind', '5mer', 'fimo_summary', 'polyA_polyT_GC', 'dna_shape']


# epigenetic
# 'encode_matrix'

featurizer.generate_features(Features, example_dir + '/', max_num_threads=6)
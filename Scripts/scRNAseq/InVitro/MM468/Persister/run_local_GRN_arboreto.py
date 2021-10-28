import os
import pandas as pd

from arboreto.algo import grnboost2, genie3
from arboreto.utils import load_tf_names

wd = os.getcwd() + '/'

net1_ex_path = wd + 'int/1.1_exprMatrix_filtered_t.txt'
net1_tf_path = wd + 'int/1.1_inputTFs.txt'

ex_matrix = ex_matrix = pd.read_csv(net1_ex_path, sep='\t', index_col=0, header=None).T


ex_matrix.shape

tf_names = load_tf_names(net1_tf_path)

network = grnboost2(expression_data=ex_matrix,
                    tf_names=tf_names)

network.head()

network.to_csv('ex_01_network.tsv', sep='\t', header=False, index=False)


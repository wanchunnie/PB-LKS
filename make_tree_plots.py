# draw trees
from sklearn import tree
from itertools import product
import pickle
import matplotlib.pyplot as plt

RF_PATH = './PBLKS/PBLKS_model.pkl'

clf = pickle.load(open(RF_PATH, 'rb'))
fn = [''.join(kmer) for kmer in product(['A', 'T', 'G', 'C'], repeat=4)]

cn = ['0(not host)', '1(is host)']
for i in range(len(clf.estimators_)):
    fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(10, 10), dpi=1000)
    tree.plot_tree(clf.estimators_[i],
                   feature_names=fn,
                   class_names=cn,
                   filled=True)
    fig.savefig('./PBLKS/Trees/decision_tree_' +
                str(i) + '_visualization.png')

from pytsne import tsne
import numpy as np

a = np.random.random((1000,5))
tsne.run(a, 2, 5, 0, 0, 0)

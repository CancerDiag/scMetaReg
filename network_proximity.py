import numpy as np
import pandas as pd
import numpy.ma as ma
import networkx as nx
from time import sleep
from tqdm import tqdm
import seaborn as sns
import matplotlib.pyplot as plt

# ### Network Proximity
import sys
import random
import numpy as np
import numpy.ma as ma
import networkx as nx
import pandas as pd


class Interactome(object):
    #def __init__(self, pathG="HumanInteractome.tsv", pathSD="a.npy.npz", binSize=300):
    def __init__(self, pathG=pcnet , SD= pcnet_SD, binSize=300):
        #self.G = nx.read_edgelist(pathG, delimiter="\t", data=[("src", str), ("typ", str)])
        #self.G.remove_edges_from(nx.selfloop_edges(self.G))
        self.G = pathG
        self.SD = SD
        self.nodes = sorted(self.G.nodes())
        self.i2n = {index: node for index, node in enumerate(self.nodes)}
        self.n2i = {node: index for index, node in enumerate(self.nodes)}
        self.i2d = {}
        self.d2i = {}
        for node, degree in self.G.degree():
            index = self.n2i[node]
            if degree in self.d2i:
                self.d2i[degree].append(index)
            else:
                self.d2i[degree] = [index]
            self.i2d[index] = degree
        self.dmin = min(self.d2i.keys())
        self.dmax = max(self.d2i.keys())
        # ------------------------------
        self.binSize = binSize
        self.d2b = {}
        self.b2i = {}
        self.d2b = {}
        self.b2i = {0: []}
        degrees = sorted(self.d2i.keys())
        b = 0
        for curr, till in zip(degrees, degrees[1:] + [degrees[-1] + 1]):
            for d in range(curr, till):
                self.d2b[d] = b
            self.b2i[b].extend(self.d2i[curr])
            if curr != degrees[-1] and len(self.b2i[b]) >= binSize:
                b += 1
                self.b2i[b] = []
        if len(self.b2i[b]) < binSize and b > 0:
            for d in range(degrees[-1], -1, -1):
                if self.d2b[d] != b:
                    break
                self.d2b[d] = b - 1
            self.b2i[b - 1].extend(self.b2i[b])
            del self.b2i[b]

    def Name2Index(self, names, skipUnknown=True):
        if skipUnknown:
            return [self.n2i[n] for n in names if n in self.n2i]
        else:
            return [self.n2i[n] for n in names]

    def DegreePreserveSampling(self, indexes):
        b2n = {}
        for index in indexes:
            b = self.d2b[max(self.dmin, min(self.dmax, self.i2d[index]))]
            if b not in b2n:
                b2n[b] = 0
            b2n[b] += 1
        for b in b2n:
            if b2n[b] > len(self.b2i[b]):
                raise ValueError("Number to sample > size of bin. Try to increase binSize")
        while True:
            yield sum([random.sample(self.b2i[b], b2n[b]) for b in sorted(b2n)], [])

    def Proximity(self, mod1, mod2):
        SD = self.SD[np.ix_(mod1, mod2)]
        closest1 = SD.min(0)
        closest2 = SD.min(1)
        return (closest1.sum() + closest2.sum()) / (closest1.count() + closest2.count())

    def ProximityRandom(self, mod1, mod2, repeat):
        result = np.zeros(repeat)
        index = 0
        for mod1r, mod2r in zip(self.DegreePreserveSampling(mod1), self.DegreePreserveSampling(mod2)):
            v = self.Proximity(mod1r, mod2r)
            if not ma.is_masked(v):
                result[index] = v
                index += 1
                if index == repeat:
                    break
        return result

    def ProximityZ(self, mod1, mod2, repeat):
        d = self.Proximity(mod1, mod2)
        b = self.ProximityRandom(mod1, mod2, repeat=repeat)
        z, p = Z_Score(d, b)
        return d, z, p


def Z_Score(real, background):
    m = background.mean()
    s = background.std(ddof=1)
    z = (real - m) / s
    p = np.mean(background < real)
    return z, p


def network_proximity(input1, input2, repeat = 1000, random_seed = 11096):
    random.seed(int(11096))
        
    interactome = Interactome()
    
    genes1 = interactome.Name2Index(input1)
    genes2 = interactome.Name2Index(input2)
    
    d, z, p = interactome.ProximityZ(genes1, genes2, repeat=repeat)
    print(d, z, p)
    return [d, z, p]


# Shortest Distance
def compute_shortest_distance(interactome, mod1, mod2):
    SD = np.zeros((len(mod1), len(mod2)))
    for row, source in enumerate(mod1):
        for col, target in tqdm(enumerate(mod2)):
            
            try:
                SD[row, col] = nx.shortest_path_length(interactome, source, target)
            except nx.NetworkXNoPath:
                SD[row, col] = 9999
        #sleep(0.0001)
    SD = ma.masked_values(SD, 9999)
    return SD

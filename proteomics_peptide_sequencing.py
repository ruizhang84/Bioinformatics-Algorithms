"""
    Solve the Peptide Sequencing Problem.
    """
# general imports
import time
import random

# list of amino acid molecular weights (one-letter)
MASS_TOY = {4:'X', 5:'Z'}
MASS_PROTEIN = {57:'G', 71:'A', 87:'S', 97:'P', 99:'V', 101:'T',
                103:'C', 113:'IL', 114:'N', 115:'D', 128:'KQ',
                129:'E', 131:'M', 137:'H', 147:'F', 156:'R',
                163:'Y', 186:'W'}

#Sample Input:
#0 0 0 4 -2 -3 -1 -7 6 5 3 2 1 9 3 -8 0 3 1 2 1 8
#Sample Output:
#XZZXX
def path_to_peptide(path):
    """
        take a path (0, 131, ..485)
        return the corresponding peptide
        """
    peptide = ''   #peptide
    for id in range(len(path)-1):
        #key->path_j-path_i, equal to an amino acid mass
        key = abs(path[id+1]-path[id])
        if key == 113 or key == 128:
            peptide += random.choice(MASS_PROTEIN[key])
        else:
            peptide += MASS_PROTEIN[key]
    return peptide

class TxtExtract:
    """
        Given: A space-delimited spectral vector 'Spectrum'.
        Return: An amino acid string with maximum score against 'Spectrum'. For masses
        with more than one amino acid, any choice may be used.
        """
    def __init__(self, txt_file):
        f = open(txt_file, 'r')
        self._vector = [0]+[int(x) for x in f.read().split(" ")]
        
        self.graph = {}
        self.edge = {}
        self.in_graph = {}
        self.weight = self._vector[:]
        f.close()
    
    def spectrum_graph(self):
        """
            connecting nodes si and sj by a directed edge labeled 
            by an amino acid a if sj-si is equal to the mass of a.
            """
        for i in range(len(self._vector)):
            for j in range(i+1, len(self._vector)):
                key = abs(j - i)
                if key in MASS_PROTEIN:
                    self.graph[i] = self.graph.get(i,[]) + [j]
                    self.edge[(i, j)] = 1
    
    def sequence(self):
        """
            Any path connecting source to sink corresponds to an amino acid string Peptide
            assign weight si to node i and assign weight zero to node 0
            return a maximum-weight path from source to sink in a node-weighted DAG.
            """
        self.spectrum_graph()
        score = [ [float('-inf') for x in range(len(self.weight))] for x in range(len(self.weight)-1)]
        
        for j in self.graph[0]:                     # initilize
            score[0][j] = self.weight[j]
        for i in  range(1, len(self.weight)-1):     # dynamic programming, vector source -> sink
            score[i] = score[i-1][:]
            if i in self.graph:
                for j in self.graph[i]:
                    score[i][j] = max(self.weight[j]+score[i-1][i], score[i-1][j])

        path = [len(self.weight)-1]                 # path (0, ... sink)
        score_idx = len(self.weight)-2              # score index
        vector_idx = len(self.weight)-1             # vector index
        
        while score_idx > 0:                        # tace back
            if score[score_idx][vector_idx] == score[score_idx-1][vector_idx]:
                score_idx -= 1
            else:
                current = vector_idx
                while vector_idx > 0:
                    vector_idx -= 1
                    if score[score_idx][current] == score[score_idx-1][vector_idx] + self.weight[current]:
                        if (vector_idx, current) in self.edge:
                            score_idx = vector_idx
                            path.insert(0, vector_idx)
                            break
        path.insert(0, 0)
        print path_to_peptide(path)

if __name__=="__main__":
    test=TxtExtract('dataset_11813_10.txt')
    #start_time = time.time()
    test.sequence()
    #print("---running time is %s seconds ---" % (time.time() - start_time))



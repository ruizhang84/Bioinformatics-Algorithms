"""
    Solve the Peptide Search Problem.
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

TOY_MASS = {'X':4, 'Z':5}
PROTEIN_MASS = {'G':57, 'A':71, 'S':87, 'P':97, 'V':99, 'T':101,
                'C':103, 'I':113, 'L':113, 'N':114, 'D': 115, 'K':128,
                'Q':128, 'E':129, 'M':131, 'H':137, 'F':147, 'R':156,
                'Y':163, 'W':186}

#Sample Input:
#0 0 0 4 -2 -3 -1 -7 6 5 3 2 1 9 3 -8 0 3 1 2 1 8
#XZZXZXXXZXZZXZXXZ
#Sample Output:
#ZXZXX

class TxtExtract:
    """
        Given: A set of space-delimited spectral vectors SpectralVectors, 
        an amino acid string Proteome, and an integer threshold.
        Return: The set PSM_threshold (Proteome, SpectralVectors).
        """
    def __init__(self, txt_file):
        f = open(txt_file, 'r')
        
        self._vectors = []
        self.weights = []
        proteome_read = 0
        for line in f:
            if len(line.strip().split(" ")) > 1:
                line = [0] + [int(x) for x in line.strip().split(" ")]
                self.weights.append(line)
                self._vectors.append(line)
            elif proteome_read == 0:
                self.proteome = line.strip()
                proteome_read = 1
            elif proteome_read == 1:
                self.threshold = int(line.strip())
        self.graph = {}
        f.close()
    
    def spectrum_graph(self, graph_idx):
        """
            connecting nodes si and sj by a directed edge labeled 
            by an amino acid a if sj-si is equal to the mass of a.
            """
        self.graph = {}         # initialize
        for i in range(len(self._vectors[graph_idx])):
            for j in range(i+1, len(self._vectors[graph_idx])):
                key = abs(j - i)
                if key in MASS_PROTEIN:
                    self.graph[i] = self.graph.get(i,[]) + [j]
    
    def sequence(self):
        """
            take substring from proteome, match with spectrum, and score
            return the maximum substring
            """
        peptide_set = set([])
        for graph_idx in range(len(self._vectors)):
            score_max = float('-inf')
            peptide_max = ''
            self.spectrum_graph(graph_idx)       # construct the spectrum graph
            for i in range(len(self.proteome)):
                peptide, score= self.peptide_to_path(i, graph_idx)
                if score > score_max:
                    score_max = score
                    peptide_max = peptide
        
            if score_max >= self.threshold:
                peptide_set.add(peptide_max)
                    
        for element in peptide_set:
            print element

    def peptide_to_path(self, idx, graph_idx):
        """
            A match function
            take a peptide substring starting at idx in proteome,
            return the corresponding substring and score
            """
        score = 0
        spectrum_idx = 0
        spectrum_idx_path = [0]
        spectrum_peak = PROTEIN_MASS[self.proteome[idx]]
        spectrum_length = len(self.weights[graph_idx])-1
        peptide = ''
        peptide_temp = self.proteome[idx]

        while (spectrum_peak in self.graph[spectrum_idx]):
            if spectrum_peak == spectrum_length:
                peptide = peptide_temp
                spectrum_idx_path.append(spectrum_peak)
                break
            elif spectrum_peak > spectrum_length:
                break
            idx += 1
            if idx == len(self.proteome):
                break
            spectrum_idx = spectrum_peak
            if spectrum_idx not in self.graph:
                break
            spectrum_idx_path.append(spectrum_idx)
            spectrum_peak += PROTEIN_MASS[self.proteome[idx]]
            peptide_temp += self.proteome[idx]

        if len(peptide) == 0:
            return ('',float('-inf'))
        
        for i in spectrum_idx_path:
            score += self.weights[graph_idx][i]
        return (peptide, score)


if __name__=="__main__":
    test=TxtExtract('dataset_11866_5.txt')
    #start_time = time.time()
    test.sequence()
    #print("---running time is %s seconds ---" % (time.time() - start_time))



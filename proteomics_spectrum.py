"""
    Solve the Peptide Sequencing.
    """
# general imports
import time
import random

# list of amino acid molecular weights (one-letter)
MASS_PROTEIN = {57:'G', 71:'A', 87:'S', 97:'P', 99:'V', 101:'T',
                103:'C', 113:'IL', 114:'N', 115:'D', 128:'KQ',
                129:'E', 131:'M', 137:'H', 147:'F', 156:'R',
                163:'Y', 186:'W'}

#Sample Input:
#57 71 154 185 301 332 415 429 486
#Sample Output:
#0->57:G
#0->71:A
#57->154:P
#57->185:K
#71->185:N
#154->301:F
#185->332:F
#301->415:N
#301->429:K
#332->429:P
#415->486:A
#429->486:G

class TxtExtract:
    """
        Given: A space-delimited list of integers Spectrum.
        Return: Graph(Spectrum).
        """
    def __init__(self, txt_file):
        f = open(txt_file, 'r')
        self._mass = [0] + [int(x) for x in f.read().split(" ")]
        self.graph = {} #  0: [57, 71]..
        self.edge = {}  # (0, 57): G, (0, 71): A..
        f.close()

    def spectrum_graph(self):
        """
            connecting nodes si and sj by a directed edge labeled 
            by an amino acid a if sj-si is equal to the mass of a.
            """
        for x in range(len(self._mass)):
            for y in range(x+1, len(self._mass)):
                key = self._mass[y] - self._mass[x]
                if key in MASS_PROTEIN:
                    self.graph[self._mass[x]] = self.graph.get(self._mass[x],[]) + [self._mass[y]]
                    if key == 113 or key == 128:
                        self.edge[(self._mass[x], self._mass[y])] = random.choice(MASS_PROTEIN[key])
                    else:
                        self.edge[(self._mass[x], self._mass[y])] = MASS_PROTEIN[key]
    
        for key in sorted(self.edge):  #print out graph
            print "%d->%d:%s" %(key[0], key[1], self.edge[key])


if __name__=="__main__":
    test=TxtExtract('dataset_11813_2.txt')
    #start_time = time.time()
    test.spectrum_graph()
    #print("---running time is %s seconds ---" % (time.time() - start_time))

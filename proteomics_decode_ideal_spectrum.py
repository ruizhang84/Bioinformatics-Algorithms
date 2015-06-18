"""
    Decoding an Ideal Spectrum.
    """
# general imports
import time
import random

# list of amino acid molecular weights (one-letter)
MASS_PROTEIN = {57:'G', 71:'A', 87:'S', 97:'P', 99:'V', 101:'T',
                103:'C', 113:'IL', 114:'N', 115:'D', 128:'KQ',
                129:'E', 131:'M', 137:'H', 147:'F', 156:'R',
                163:'Y', 186:'W'}

PROTEIN_MASS = {'G':57, 'A':71, 'S':87, 'P':97, 'V':99, 'T':101,
                'C':103, 'I':113, 'L':113, 'N':114, 'D': 115, 'K':128,
                'Q':128, 'E':129, 'M':131, 'H':137, 'F':147, 'R':156,
                'Y':163, 'W':186}

#Sample Input:
#57 71 154 185 301 332 415 429 486
#Sample Output:
#GPFNA

def ideal_spectrum(peptide):
    """
        A helper function: take the peptide sequence,
        return the ideal spectrum
        """
    spectrum = set([])
    prefix_mass, surffix_mass = 0, 0
    for amino_id in range(len(peptide)):
        #peptide->prefix+suffix, ABCD-> A,D; (A)B, (C)D...
        prefix_peptide = peptide[amino_id]
        surffix_peptide = peptide[-amino_id-1]
        #iterations to add mass of the new amino acid
        prefix_mass += PROTEIN_MASS[prefix_peptide]
        surffix_mass += PROTEIN_MASS[surffix_peptide]
        spectrum.update([prefix_mass, surffix_mass])
    return spectrum

def path_to_peptide(path):
    """
        take a path (0, 131, ..485)
        return the corresponding peptide
        """
    peptide = ''   #peptide
    for id in range(len(path)-1):
        #key->path_j-path_i, equal to an amino acid mass
        key = path[id+1]-path[id]
        if key == 113 or key == 128:
            peptide += random.choice(MASS_PROTEIN[key])
        else:
            peptide += MASS_PROTEIN[key]
    return peptide

def bread_first_search(graph, start, end):
    """
        A bread first algorithm,
        return all the path that correspoding to graph
        """
    all_path = []
    # maintain a queue of paths
    queue = []
    # push the first path into the queue
    queue.append([start])
    while queue:
        # get the first path from the queue
        path = queue.pop(0)
        # get the last node from the path
        node = path[-1]
        # path found
        if node == end:
            all_path.append(path)
        # enumerate all adjacent nodes, construct a new path and push it into the queue
        for adjacent in graph.get(node, []):
            new_path = list(path)
            new_path.append(adjacent)
            queue.append(new_path)
    return all_path


class TxtExtract:
    """
        Given: A space-delimited list of integers Spectrum.
        Return: An amino acid string that explains Spectrum.
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

    def decoding_ideal_spectrum(self):
        """
            for each path Path from source to sink in Graph(Spectrum)
                Peptide <- the amino acid string spelled by the edge labels of Path
                if IdealSpectrum(Peptide) = Spectrum then return Peptide
            """
        self.spectrum_graph()
        spectrum = set(self._mass[1:])          # spectrum
        last = self._mass[-1]
        for path in bread_first_search(self.graph, 0, last):
            peptide = path_to_peptide(path)
            if ideal_spectrum(peptide) == spectrum:
                return path_to_peptide(path)


if __name__=="__main__":
    test=TxtExtract('dataset_11813_4.txt')
    #start_time = time.time()
    print test.decoding_ideal_spectrum()
    #print("---running time is %s seconds ---" % (time.time() - start_time))



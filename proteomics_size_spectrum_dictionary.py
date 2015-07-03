"""
    Solve the Size of Spectral Dictionary Problem.
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
MASS_PROTEIN_LIST = [57, 71, 87, 97, 99, 101, 103, 113, 113, 114, 115, 128, 128,
                129, 131, 137, 147, 156, 163, 186]
#Sample Input:
#4 -3 -2 3 3 -4 5 -3 -1 -1 3 4 1 3
#1
#8
#Sample Output:
#3

class TxtExtract:
    """
        Given: A spectral vector Spectrum', an integer threshold, and an integer max_score.
        Return: The size of the dictionary Dictionarythreshold(Spectrum').
        """
    def __init__(self, txt_file):
        f = open(txt_file, 'r')
        self.vector = [0] + [int (x) for x in f.readline().strip().split(" ")]
        self.threshold = int(f.readline().strip())
        self.max_score = int(f.readline().strip())  # provided max_score for the height of the table.
        f.close()

    def dictionary(self):
        """
            introduce a variable Size(i, t) as the number of peptides Peptide of mass i 
            such that Score(Peptide, Spectrum'i) is equal to t.
            return the size of a spectral dictionary
            """
        size_dp = [[0 for t in range(self.max_score+1)] for a in range(len(self.vector))]
        size_dp[0][0] = 1
        
        for i in range(1, len(self.vector)):
            for t in range(self.max_score+1):
                for x in MASS_PROTEIN_LIST:
                    if i-x < 0:
                        continue
                    elif 0 <= t-self.vector[i] <= self.max_score:
                        size_dp[i][t] += size_dp[i-x][t-self.vector[i]]
        
        print sum(size_dp[len(self.vector)-1][self.threshold:])



if __name__=="__main__":
    test=TxtExtract('dataset_11866_8.txt')
    #start_time = time.time()
    test.dictionary()
    #print("---running time is %s seconds ---" % (time.time() - start_time))



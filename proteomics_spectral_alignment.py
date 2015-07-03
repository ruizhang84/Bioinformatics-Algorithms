"""
    Solve the Spectral Alignment Problem.
    """
# general imports
import time
import random

# list of amino acid molecular weights (one-letter)
PROTEIN_MASS = {'G':57, 'A':71, 'S':87, 'P':97, 'V':99, 'T':101,
                'C':103, 'I':113, 'L':113, 'N':114, 'D': 115, 'K':128,
                'Q':128, 'E':129, 'M':131, 'H':137, 'F':147, 'R':156,
                'Y':163, 'W':186}

#Sample Input:
#XXZ
#4 -3 -2 3 3 -4 5 -3 -1 -1 3 4 1 -1
#2
#Sample Output:
#XX(-1)Z(+2)

class TxtExtract:
    """
        Given: A peptide Peptide, a spectral vector Spectrum', and an integer k.
        Return: A peptide Peptide' related to Peptide by up to k modifications with
        maximal score against Spectrum' out of all possibilities.
        """
    def __init__(self, txt_file):
        f = open(txt_file, 'r')
        self.peptide = f.readline().strip()
        self.vector = [0] + [int (x) for x in f.readline().strip().split(" ")]
        self.post_modif = int(f.readline().strip())
        self.diff = {0 : 0}     # difference array implemented with dict{}
        f.close()

    def spectrum_alignment(self):
        """
            dynamic programming for score(i,j,t)
            score(i,j,t) = sj + max(score(i-diff(i), j',t), score(i-diff(i),j-diff(j),t))
            returnt the post-translation modification sequence
            """
        self.diff_PROTEIN()
        
        score = []  #node->>([t][i][j])
        for t in range(self.post_modif+1):
            pos = 0             # position of peptide for converting mass
            score_ij = {0: [ float('-inf') for t in range(len(self.vector))]}
            for amino in self.peptide:
                score_j = [ float('-inf') for t in range(len(self.vector))]
                pos += PROTEIN_MASS[amino]
                score_ij[pos] = score_j
            score.append(score_ij)
        
        score[0][0][0] = 0
        # score for node(i,j,t)
        for t in range(self.post_modif+1):
            for i in sorted(score[t]):
                if i > 0:  # i-self.diff[i]
                    for j in range(len(self.vector)):
                        temp_max = float('-inf')
                        if j >= self.diff[i]:
                            temp_max = score[t][i-self.diff[i]][j-self.diff[i]]
                        if t > 0:
                            for j_p in range(j):
                                if temp_max < score[t-1][i-self.diff[i]][j_p]:
                                    temp_max = score[t-1][i-self.diff[i]][j_p]
                
                        score[t][i][j] = self.vector[j] + temp_max
    
        # trace back --> the longest path
        max_score = float('-inf')
        layer = 0                         # modify
        row = pos                           # mass
        column = len(self.vector)-1       # vector
        modify = []
        for t in range(self.post_modif+1):
            if max_score < score[t][pos][-1] :
                max_score = score[t][pos][-1]
                layer = t
        
        while layer > 0:
            score_temp = score[layer][row][column] - self.vector[column]
            if score_temp == score[layer][row-self.diff[row]][column-self.diff[row]]:
                column -= self.diff[row]
                row -= self.diff[row]
            else:
                for j_p in range(column-1):
                    if score_temp == score[layer-1][row-self.diff[row]][j_p]:
                        modify.append((row, column-row))
                        row -= self.diff[row]
                        column = j_p
                        layer -= 1
                        break
    

        # print out the sequence
        modify.sort()
        sequence = ""
        pos = 0
        i = 0
        mass = 0
        for amino in self.peptide:
            pos += PROTEIN_MASS[amino]
            sequence += str(amino)
            if pos == modify[i][0]:
                if i == 0:
                    mass = modify[i][1]
                else:
                    mass = modify[i][1]-modify[i-1][1]
                
                if mass > 0:
                    sequence += "(+"+str(mass)+")"
                else:
                    sequence += "("+str(mass)+")"
                i += 1
        
        print sequence
        

    def diff_PROTEIN(self):
        """
            A helper function to calculate difference array
            diff(2) = 2 - 0, diff(5) = 5 - 2, etc
            """
        pos = 0
        for i in range(len(self.peptide)):
            pos += PROTEIN_MASS[self.peptide[i]]
            self.diff[pos] = PROTEIN_MASS[self.peptide[i]]



if __name__=="__main__":
    test=TxtExtract('dataset_11866_14.txt')
    #start_time = time.time()
    test.spectrum_alignment()
    #print("---running time is %s seconds ---" % (time.time() - start_time))



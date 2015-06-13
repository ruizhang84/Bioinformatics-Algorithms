"""
    Solve the Soft Decoding Problem
    """

import math
import re
class TxtExtract:
    """
        Input: A string x, followed by the alphabet from which x was constructed,
        followed by the states States, transition matrix Transition, and emission matrix
        Emission of an HMM (alphabet, States, Transition, Emission).
        Output: An |x| x |States| matrix whose (i, k)-th element holds the conditional probability.
        """
    def __init__(self,txt_file):
        f = open(txt_file,'r')
        
        change_line = 0  #count "----"
        self.transit_matrix = []
        self.emission_matrix = []
        for line in f:
            line = line.strip()
            if line == '--------':
                change_line += 1
            elif change_line == 0:
                self.string = line
            elif change_line == 1:
                self.alphabet = line.split(" ")
            elif change_line == 2:
                self.state = line.split(" ")
            elif change_line == 3:
                self.transit_matrix.append( line.split("\t") )
            else:
                self.emission_matrix.append( line.split("\t") )
        
        self.transit_matrix.pop(0) #pop out lable
        self.emission_matrix.pop(0)
        for i, row in enumerate(self.transit_matrix):
            row.pop(0)
            self.transit_matrix[i] = [float(x) for x in row ]
        for i, row in enumerate(self.emission_matrix):
            row.pop(0)
            self.emission_matrix[i] = [float(x) for x in row ]

        self.alphabet_index = {} #index {'x'->0,'y'->1,..}
        for i, alphabet in enumerate(self.alphabet):
            self.alphabet_index[alphabet] = i

        self.sink = 0
        self.forward = [ [1.0 for x in range(len(self.state))] for x in range(len(self.string))]
        self.backward = [ [1.0 for x in range(len(self.state))] for x in range(len(self.string))]
        self.condit_prob = [ [1.0 for x in range(len(self.state))] for x in range(len(self.string))]
        f.close()

    def soft_decoding(self):
        """
            return the conditional probability that the HMM was in state k
            at time i given that it emitted string x.
            """
        self.forward_prob()
        self.backward_prob()
        for i in range(len(self.string)):
            for j in range(len(self.state)):
                self.condit_prob[i][j] = self.forward[i][j]*self.backward[i][j]/self.sink

        for i in self.state:
            print i,
        print
        for prob in self.condit_prob:
            for j in prob:
                print round(j,4),
            print


    def forward_prob(self):
        """
            calculate the forward probability P(forward)_i
            """
        for i in range(len(self.state)):
            self.forward[0][i] = self.emission_matrix[i][self.alphabet_index[self.string[0]]]/(len(self.state))
        for x in range(1, len(self.string)):
            for i in range(len(self.state)):
                emit  = self.emission_matrix[i][self.alphabet_index[self.string[x]]]
                transit = 0
                for j in range(len(self.state)):
                    transit += self.transit_matrix[j][i]*self.forward[x-1][j]
                self.forward[x][i] = transit*emit
        for i in range(len(self.state)):
            self.sink += self.forward[-1][i]


    def backward_prob(self):
        """
            calculate the backward probability P(backward)_i
            """
        for i in range(len(self.state)):
            self.backward[-1][i] = 1.0
        for i in range(1, len(self.string)):
            for k in range(len(self.state)):
                backward = 0
                for l in range(len(self.state)):
                    weight  = self.emission_matrix[l][self.alphabet_index[self.string[-i]]]*self.transit_matrix[k][l]
                    backward += self.backward[-i][l]*weight
                self.backward[-i-1][k] = backward


if __name__=="__main__":
    test=TxtExtract('dataset_11632_12.txt')
    test.soft_decoding()
    #test=TXT_Extract('dataset_10928_3.txt')



"""
    Implement Viterbi learning for estimating the parameters of an HMM.
    """

import math
import re
class TxtExtract:
    """
        Input: A number of iterations j, followed by a string x of symbols emitted by an HMM,
        followed by the HMM's alphabet, followed by the HMM's states, followed by initial transition
        and emission matrices for the HMM.
        Output: Emission and transition matrices resulting from applying Viterbi learning for j iterations.
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
                self.num_iter = int(line)
            elif change_line == 1:
                self.string = line
            elif change_line == 2:
                self.alphabet = line.split(" ")
            elif change_line == 3:
                self.state = line.split(" ")
            elif change_line == 4:
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

        self.state_index = {} #index{'A'->0, 'B'->1, ..}
        self.alphabet_index = {}
        for i, state in enumerate(self.state):
            self.state_index[state] = i
        for i, alphabet in enumerate(self.alphabet):
            self.alphabet_index[alphabet] = i
        self.path = ''
        f.close()

    def iter_learning(self):
        """
            iteration of viterbi learning
            """
        for i in range(self.num_iter):
            self.prob_path()
            self.parameter_estimate()
    
        print "",
        for symbol in self.state:
            print symbol,
        print
        for i in range(len(self.transit_matrix)):
            print self.state[i],
            for fraction in self.transit_matrix[i]:
                print round(fraction,3),
            print
        
        print "--------"
     
        print "",
        for symbol in self.alphabet:
            print symbol,
        print
        for i in range(len(self.emission_matrix)):
            print self.state[i],
            for fraction in self.emission_matrix[i]:
                print round(fraction,3),
            print

    def parameter_estimate(self):
        """
            return the transition matrix followed 
            by the emission matrix of HMM.
            """
        self.transit_matrix_build()
        self.emission_matrix_build()
            
    def emission_matrix_build(self):
        """
            return the emission matrix.
            """
        self.emission_matrix = []
        # count the number of xyz
        for abc in self.state:
            pos_index = []
            count_xyz = [0 for x in range(len(self.alphabet))]
            pattern = re.compile(abc)  #pattern search
            for m in pattern.finditer(self.path):
                pos_index.append(m.start())
        
            for pos in pos_index:
                for xyz in range(len(self.alphabet)):
                    if self.alphabet[xyz] == self.string[pos]:
                        count_xyz[xyz] += 1
            self.emission_matrix.append(count_xyz)
        #convert to fraction
        matrix_len = len(self.emission_matrix)
        sum_temp = 0
        for i in range(matrix_len):
            sum_temp = sum(self.emission_matrix[i])
            if sum_temp == 0:
                for j in range(len(self.emission_matrix[0])):
                    self.emission_matrix[i][j] = 1.0/len(self.emission_matrix[0])
            else:
                for j in range(len(self.emission_matrix[0])):
                    self.emission_matrix[i][j] = self.emission_matrix[i][j]/float(sum_temp)

    
    def transit_matrix_build(self):
        """
            return the transition matrix.
            """
        self.transit_matrix = []
        for abc in self.state:
            temp_transit = []
            for pair in self.state:
                pattern = str(abc) + str(pair)
                num = countoverlapping(pattern, self.path)
                temp_transit.append(num)
            self.transit_matrix.append(temp_transit)

        #convert to fraction
        matrix_len = len(self.transit_matrix)
        sum_temp = 0
        for i in range(matrix_len):
            sum_temp = sum(self.transit_matrix[i])
            if sum_temp == 0:
                for j in range(matrix_len):
                    self.transit_matrix[i][j] = 1.0/matrix_len
            else:
                for j in range(matrix_len):
                    self.transit_matrix[i][j] = self.transit_matrix[i][j]/float(sum_temp)

    def prob_path(self):
        """
            a dynamic programming alogrithm to solve the Decoding Problem.
            
            maximizes Pr(outcome,state) over all possible paths
            """
        # dynamic programming array [A,B,...] x string
        dp = [ [0.0 for n in range(len(self.state))] for x in range(len(self.string)) ]
        for i in range(len(self.string)):
            dp_temp = dp[i-1][:]
            for j in range(len(self.state)):
                dp[i][j] = self.max_state(self.state[j], self.string[i], dp_temp)
        
        dp_max = max(dp[len(self.string)-1]) # traceback the path
        for i, val in enumerate(dp[len(self.string)-1]):
            if dp_max == val:
                hidden_path = [self.state[i]]
                break

        for i in range(len(self.string)-1):
            temp_dp = []                      # backtracking the dp
            for st in self.state:
                temp_dp.append(dp_max-self.transition(st, hidden_path[0])-self.emission(hidden_path[0], self.string[-i-1]))
               
            for j, dp_i in enumerate(temp_dp):
                if compare_float(dp_i, dp[-i-2]):
                    hidden_path = [self.state[j]]+hidden_path
                    dp_max = dp_i
                    break
                        
        self.path = "".join(x for x in hidden_path)
        

    def max_state(self, state, outcome, dp):
        """
            return the probability Pr(outcome, state) at the current state
            """
        temp_dp = []
        for i, st in enumerate(self.state):
            temp_dp.append( dp[i]+self.transition(st, state)+self.emission(state, outcome) )
        return max(temp_dp)

    def transition(self, prev_state, state):
        """
            A helper function return the transition
            probability (logP) given states.
            """
        transit_temp = self.transit_matrix[self.state_index[prev_state]][self.state_index[state]]
        if transit_temp == 0:
            return float('-inf')
        return math.log(transit_temp)
        
    def emission(self, state, outcome):
        """
            A helper function return a emission 
            probability (logP)given states and outcome.
            """
        emitt_temp = self.emission_matrix[self.state_index[state]][self.alphabet_index[outcome]]
        if emitt_temp == 0:
            return float('-inf')
        return math.log(self.emission_matrix[self.state_index[state]][self.alphabet_index[outcome]])

def compare_float(num_1, num_list):
    """ 
        return true if number is in the list
        """
    for number in num_list:
        if round(num_1, 4) == round(number, 4):
            return True

def countoverlapping(pattern, string):
    """
        A helper function
        counting overlapping matches
        """
    total = 0
    start = 0
    pat = re.compile(pattern)
    while True:
        pos = pat.search(string, start)
        if pos is None:
            return total
        total += 1
        start = 1 + pos.start()

if __name__=="__main__":
    test=TxtExtract('dataset_11632_10.txt')
    test.iter_learning()
    #test=TXT_Extract('dataset_10928_3.txt')



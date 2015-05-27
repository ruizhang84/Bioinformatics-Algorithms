import math
class TXT_Extract:
    def __init__(self,dna_file):
        f = open(dna_file,'r')
        self.lists = []
        for line in f:
            line = line.strip()
            line = line.split("\t")
            self.lists.append(line)
        
        self.outcome = tuple(self.lists.pop(0)[0]) #xyzxyzzz

        self.probAA = math.log(float(self.lists[6][1]))   #transition probability Pr(A->A)
        self.probAB = math.log(float(self.lists[6][2]))
        self.probAC = math.log(float(self.lists[6][3]))
        self.probBA = math.log(float(self.lists[7][1]))
        self.probBB = math.log(float(self.lists[7][2]))
        self.probBC = math.log(float(self.lists[7][3]))
        self.probCA = math.log(float(self.lists[8][1]))
        self.probCB = math.log(float(self.lists[8][2]))
        self.probCC = math.log(float(self.lists[8][3]))
        
        self.probAX = math.log(float(self.lists[11][1]))  #emission probability Pr(x|A)
        self.probAY = math.log(float(self.lists[11][2]))
        self.probAZ = math.log(float(self.lists[11][3]))
        self.probBX = math.log(float(self.lists[12][1]))
        self.probBY = math.log(float(self.lists[12][2]))
        self.probBZ = math.log(float(self.lists[12][3]))
        self.probCX = math.log(float(self.lists[13][1]))
        self.probCY = math.log(float(self.lists[13][2]))
        self.probCZ = math.log(float(self.lists[13][3]))
        f.close()

    def prob_path(self):
        """
            a dynamic programming alogrithm to solve the Decoding Problem.
            
            maximizes Pr(outcome,state) over all possible paths
            """
        dp = [ [0.0, 0.0, 0.0] for x in range(len(self.outcome)) ] # dynamic programming array [A,B]
        for i in range(len(self.outcome)):
            dp_temp = dp[i-1][:]
            dp[i][0] = self.max_state('A', self.outcome[i], dp_temp)
            dp[i][1] = self.max_state('B', self.outcome[i], dp_temp)
            dp[i][2] = self.max_state('C', self.outcome[i], dp_temp)
        
        dp_max = max(dp[len(self.outcome)-1]) # traceback the path
        if dp_max == dp[len(self.outcome)-1][0]:
            hidden_path = ['A']
        elif dp_max == dp[len(self.outcome)-1][1]:
            hidden_path = ['B']
        else:
            hidden_path = ['C']
        for i in range(len(self.outcome)-1):
            dA = dp_max-self.transition('A', hidden_path[0])-self.emission(hidden_path[0], self.outcome[-i-1])
            dB = dp_max-self.transition('B', hidden_path[0])-self.emission(hidden_path[0], self.outcome[-i-1])
            dC = dp_max-self.transition('C', hidden_path[0])-self.emission(hidden_path[0], self.outcome[-i-1])
            if compare_float(dA, dp[-i-2]):  ##set the path
                hidden_path = ['A']+hidden_path
                dp_max = dA
            elif compare_float(dB, dp[-i-2]):
                hidden_path = ['B']+hidden_path
                dp_max = dB
            else:
                hidden_path = ['C']+hidden_path
                dp_max = dC
        
        print "".join(x for x in hidden_path)

    def max_state(self,state, outcome, dp):
        """
            return the probability Pr(outcome, state) at the current state
            """
        dA = dp[0]+self.transition('A', state)+self.emission(state, outcome)
        dB = dp[1]+self.transition('B', state)+self.emission(state, outcome)
        dC = dp[2]+self.transition('C', state)+self.emission(state, outcome)
        return max(dA, dB, dC)

    def transition(self, prev_state, state):
        """
            return a transitoin probability given a state and
            previous state
            if transition at initial state, return 0.5
                                                        
            """
        if prev_state == 'A' and state == 'A':
            return self.probAA
        elif prev_state == 'A' and state == 'B':
            return self.probAB
        elif prev_state == 'A' and state == 'C':
            return self.probAC
        elif prev_state == 'B' and state == 'A':
            return self.probBA
        elif prev_state == 'B' and state == 'B':
            return self.probBB
        elif prev_state == 'B' and state == 'C':
            return self.probBC
        elif prev_state == 'C' and state == 'A':
            return self.probCA
        elif prev_state == 'C' and state == 'B':
            return self.probCB
        elif prev_state == 'C' and state == 'C':
            return self.probCC
        
    def emission(self, state, outcome):
        """
            return a emission probability given a state and
            an observed outcome
        
            """
        if state == 'A' and outcome == 'x':
            return self.probAX
        elif state == 'A' and outcome == 'y':
            return self.probAY
        elif state == 'A' and outcome == 'z':
            return self.probAZ
        elif state == 'B' and outcome == 'x':
            return self.probBX
        elif state == 'B' and outcome == 'y':
            return self.probBY
        elif state == 'B' and outcome == 'z':
            return self.probBZ
        elif state == 'C' and outcome == 'x':
            return self.probCX
        elif state == 'C' and outcome == 'y':
            return self.probCY
        elif state == 'C' and outcome == 'z':
            return self.probCZ

def compare_float(num_1, num_list):
    if round(num_1, 4) == round(num_list[0], 4):
        return True
    elif round(num_1, 4) == round(num_list[1], 4):
        return True
    elif round(num_1, 4) == round(num_list[2], 4):
        return True
    return False

if __name__=="__main__":
    test=TXT_Extract('dataset_11594_6.txt')
    test.prob_path()
    #test=TXT_Extract('dataset_10928_3.txt')



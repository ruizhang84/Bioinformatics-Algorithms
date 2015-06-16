"""
    Solve the Sequence Alignment with Profile HMM Problem.
    """
import math

class TxtExtract:
    """
        Input: A string x followed by a threshold theta and a pseudocount sigma, followed by an
        alphabet, followed by a multiple alignment Alignment whose strings are formed from alphabet.
        return: An optimal hidden path emitting x in HMM(Alignment, theta, sigma).
        """
    def __init__(self,txt_file):
        f = open(txt_file,'r')

        change_line = 0  #count "----"
        self.alphabet = []
        self.alignment = []
        for line in f:
            line = line.strip()
            if line == '--------':
                change_line += 1
            elif change_line == 0:
                self.string = line
            elif change_line == 1:
                line = line.split(" ")
                self.theta = float(line[0])
                self.sigma = float(line[1])
            elif change_line == 2:
                self.alphabet = line.split(" ")
            else:
                self.alignment.append(list(line))
    
        self.hidden_path = self.alignment[:]
        self.hidden_path_insert = []
        self.transit_matrix = []
        self.emission_matrix = []
        self.dp = []
        self.optimal_path = []
        self.state_num = 0 # state number
        f.close()
    
    def sequence_alignment(self):
        """
            return an optimal hidden path emitting a string x
            in HMM(Alignment, threshold, pseudocount).
            """
        self.profile_HMM()
        self.dynamic_max()
        self.trace_back()
    
    def trace_back(self):
        """
            trace back the dynamic programming
            return a state path
            """
        trace_index = {0: 'I', 1: 'M', 2:'D'}
        state_path = []                    #store path
        temp_I = self.dp[-1][0][-1] + self.transit_matrix[-2][-1]
        temp_M = self.dp[-1][1][-1] + self.transit_matrix[-4][-1]
        temp_D = self.dp[-1][2][-1] + self.transit_matrix[-3][-1]
        
        my_list = [temp_I, temp_M, temp_D] #traceback the path
        dp_max = max(my_list)
        max_index = my_list.index(dp_max)
        dp_max = self.dp[-1][max_index][-1]
        
        #test row, column
        #for row in self.dp:
        #    for col in row[0]:
        #        print math.exp(col),
        #    print
        #    for col in row[1]:
        #        print math.exp(col),
        #    print
        #    for col in row[2]:
        #        print math.exp(col),
        #    print
        #    print
        
        self.state_num = 0                                 #number of M/D state
        for i in range(1, len(self.hidden_path[0])-1):
            if type(self.hidden_path[0][i]) is str:
                self.state_num += 1
        #for i in range(state_num):
        #    Insert_symbol = 'I'+ str(i+1)
        #    Match_symbol = 'M'+str(i+1)
        #    Dele_symbol ='D'+str(i+1)
        #    state_symbol += [Match_symbol]+[Dele_symbol]+[Insert_symbol]
        #state_symbol = ['S','I0']+state_symbol+['E']

        ##print out the transtion matrix and emission matrix.
        #print "\t",
        #for symbol in state_symbol:
        #   print symbol,"\t",
        #print
        #print trace_index[max_index], state_num
        
        if max_index == 0:
            string = 'I' + str(self.state_num)
            self.optimal_path.insert(0, string)
            self.trace_I(dp_max, 0, 0)
        elif max_index == 1:
            string = 'M' + str(self.state_num)
            self.optimal_path.insert(0, string)
            self.trace_M(dp_max, 0, 0)
        else:
            string = 'D' + str(self.state_num)
            self.optimal_path.insert(0, string)
            self.trace_D(dp_max, 0, 0)

        for string in self.optimal_path:
            print string,

    def trace_I(self, dp_max, i, j):
        """
            trace back the M_i state, 
            with reversed state i, string j
            return max value and emitting string.
            """
        k = {}                      #{'A':0, 'B':1, 'C':2, ...}
        for l,s in enumerate(self.alphabet):
            k[s] = l                # alphabtic order index
        symbol = self.string[-j-1]
        if i*3+4 > len(self.transit_matrix):
            #print 'END'
            return 0
        else:
            temp_I = dp_max - self.transit_matrix[-i*3-2][-i*3-2] - self.emission_matrix[-i*3-2][ k[symbol] ]
            temp_M = dp_max - self.transit_matrix[-i*3-4][-i*3-2] - self.emission_matrix[-i*3-2][ k[symbol] ]
            temp_D = dp_max - self.transit_matrix[-i*3-3][-i*3-2] - self.emission_matrix[-i*3-2][ k[symbol] ]

            if compare_float(temp_I, self.dp[-i-1][0][-j-2]):
                #print 'I', len(self.string)-j-1
                string  = 'I' + str(self.state_num)
                self.optimal_path.insert(0, string)
                return self.trace_I(temp_I, i, j+1)
            elif compare_float(temp_M, self.dp[-i-1][1][-j-2]):
                #print 'M', len(self.string)-j-1
                string  = 'M' + str(self.state_num)
                self.optimal_path.insert(0, string)
                return self.trace_M(temp_M, i, j+1)
            elif compare_float(temp_D, self.dp[-i-1][2][-j-2]):
                #print 'D', len(self.string)-j-1
                string  = 'D' + str(self.state_num)
                self.optimal_path.insert(0, string)
                return self.trace_D(temp_D, i, j+1)


    def trace_M(self, dp_max, i, j):
        """
            trace back the I_i state
            return max value and emitting string.
            """
        k = {}                      #{'A':0, 'B':1, 'C':2, ...}
        for l,s in enumerate(self.alphabet):
            k[s] = l                # alphabtic order index
        symbol = self.string[-j-1]
        
        if i*3+7 > len(self.transit_matrix):
            #print 'END'
            return 0
        else:
            temp_I = dp_max - self.transit_matrix[-i*3-4-1][-i*3-4] - self.emission_matrix[-i*3-4][ k[symbol] ]
            temp_M = dp_max - self.transit_matrix[-i*3-4-3][-i*3-4] - self.emission_matrix[-i*3-4][ k[symbol] ]
            temp_D = dp_max - self.transit_matrix[-i*3-4-2][-i*3-4] - self.emission_matrix[-i*3-4][ k[symbol] ]
            self.state_num -= 1  #state_number -1
            if compare_float(temp_I, self.dp[-i-2][0][-j-2]):
                #print 'I', len(self.string)-j-1
                string  = 'I' + str(self.state_num)
                self.optimal_path.insert(0, string)
                return self.trace_I(temp_I, i+1, j+1)
            elif compare_float(temp_M, self.dp[-i-2][1][-j-2]):
                #print 'M', len(self.string)-j-1
                string  = 'M' + str(self.state_num)
                self.optimal_path.insert(0, string)
                return self.trace_M(temp_M, i+1, j+1)
            elif compare_float(temp_D, self.dp[-i-2][2][-j-2]):
                #print 'D', len(self.string)-j-1
                string  = 'D' + str(self.state_num)
                self.optimal_path.insert(0, string)
                return self.trace_D(temp_D, i+1, j+1)
                
    def trace_D(self, dp_max, i, j):
        """
            trace back the I_i state
            return max value and emitting string.
            """
        k = {}                      #{'A':0, 'B':1, 'C':2, ...}
        for l,s in enumerate(self.alphabet):
            k[s] = l                # alphabtic order index
        symbol = self.string[-j-1]
        if i*3+7 > len(self.transit_matrix):
            #print 'END'
            return 0
        else:
            temp_I = dp_max - self.transit_matrix[-i*3-3-2][-i*3-3] - self.emission_matrix[-i*3-3][ k[symbol] ]
            temp_M = dp_max - self.transit_matrix[-i*3-3-4][-i*3-3] - self.emission_matrix[-i*3-3][ k[symbol] ]
            temp_D = dp_max - self.transit_matrix[-i*3-3-3][-i*3-3] - self.emission_matrix[-i*3-3][ k[symbol] ]
            self.state_num -= 1
            if compare_float(temp_I, self.dp[-i-2][0][-j-1]):
                string  = 'I' + str(self.state_num)
                self.optimal_path.insert(0, string)
                #print 'I', len(self.string)-j-1
                return self.trace_I(temp_I, i+1, j)
            elif compare_float(temp_M, self.dp[-i-2][1][-j-1]):
                string  = 'M' + str(self.state_num)
                self.optimal_path.insert(0, string)
                #print 'M', len(self.string)-j-1
                return self.trace_M(temp_M, i+1, j)
            elif compare_float(temp_D, self.dp[-i-2][2][-j-1]):
                string  = 'D' + str(self.state_num)
                self.optimal_path.insert(0, string)
                #print 'D', len(self.string)-j-1
                return self.trace_D(temp_D, i+1, j)


    def dynamic_max(self):
        """
            dynamic programming
                            |- S_I(j-1),(i-1) * weight(I(j-1), M(j),i-1) -> j string, i state
            S_M(j),i = max -|- S_D(j-1),(i-1) * weight(D(j-1), M(j),i-1)
                            |- S_M(j-1),(i-1) * weight(M(j-1), M(j),i-1)
                            
            return a state prob value
            """
        # dynamic programming array [ [sequence_align -> silent + A, B, ...],.. x state]
        self.dp = [ [[0.0 for x in range(len(self.string)+1) ] for x in range(3)] for x in range(len(self.transit_matrix)/3) ]
        # alphabtic order index
        k = {}    #{'A':0, 'B':1, 'C':2, ...}
        for i,s in enumerate(self.alphabet):
            k[s] = i
        
        #initilize j->string, 0/1/2-> I/M/D
        symbol = self.string[0]
        temp_S = self.transit_matrix[0][1] + self.emission_matrix[1][ k[symbol] ]
        self.dp[0][0][1] = temp_S
        
        for j in range(1, len(self.string)):
            symbol = self.string[j]
            temp_I = self.dp[0][0][j] + self.transit_matrix[1][1] + self.emission_matrix[1][ k[symbol] ]
            self.dp[0][0][j+1] = temp_I
        
        #Adding 0-th column that contains only silent states
        temp_S =  self.transit_matrix[0][3]
        self.dp[1][2][0] = temp_S
        
        for j in range(len(self.string)):
            symbol = self.string[j]
            start = self.string[0]
            # S_M(j),i
            if j == 0:
                temp_I = self.dp[0][0][j] + self.transit_matrix[1][2] + self.emission_matrix[2][ k[symbol] ]
                temp_S = self.transit_matrix[0][2] + self.emission_matrix[2][ k[start] ]
                self.dp[1][1][j+1] = max(temp_I, temp_S)
            else:
                temp_I = self.dp[0][0][j] + self.transit_matrix[1][2] + self.emission_matrix[2][ k[symbol] ]
                self.dp[1][1][j+1] = temp_I
            # S_I(j),i
            if j == 0:
                temp_D = self.dp[1][2][0] + self.transit_matrix[3][4] + self.emission_matrix[4][ k[symbol] ]
                self.dp[1][0][1] = temp_D
            elif j == 1:
                temp_M = self.dp[1][1][j] + self.transit_matrix[2][4] + self.emission_matrix[4][ k[symbol] ]
                temp_D = self.dp[1][2][j] + self.transit_matrix[3][4] + self.emission_matrix[4][ k[symbol] ]
                self.dp[1][0][j+1]= max(temp_M, temp_D)
            else:
                temp_I = self.dp[1][0][j] + self.transit_matrix[4][4] + self.emission_matrix[4][ k[symbol] ]
                temp_M = self.dp[1][1][j] + self.transit_matrix[2][4] + self.emission_matrix[4][ k[symbol] ]
                temp_D = self.dp[1][2][j] + self.transit_matrix[3][4] + self.emission_matrix[4][ k[symbol] ]
                self.dp[1][0][j+1]= max(temp_M, temp_D, temp_I)
            
            # S_D(j),i
            if j == 0:
                temp_D = self.dp[1][2][0] + self.transit_matrix[1][3] + self.emission_matrix[1][ k[symbol] ]
                self.dp[1][2][1] = temp_D
            else:
                temp_D = self.dp[1][2][j] + self.transit_matrix[1][3] + self.emission_matrix[1][ k[symbol] ]
                self.dp[1][2][j+1] = temp_D

        # 0-th column, S_D(j),i 0/1/2-> I/M/D
        for i in range(2,len(self.transit_matrix)/3):
            temp_D = self.dp[i-1][2][0] + self.transit_matrix[i*3-3][i*3]
            self.dp[i][2][0]= temp_D
        
        # i-> state I_i-1/M_i/D_i, j->string, 0/1/2-> I/M/D
        
        for i in range(2,len(self.transit_matrix)/3):
            for j in range(len(self.string)+1):
                if j > 0:
                    self.cal_M_i(i,j) # update i, 1, j+1 from i-1, 0/1/2, j
                self.cal_D_i(i,j)     # update i,2, j from i-1, 0/1/2, j
                if j > 0:
                    self.cal_I_i(i,j) # update i, 0, j+1 from i, 0/1/2, j
            
    def cal_M_i(self, i, j):
        """
            calculate the prob of  M_i state,
            i state, j string, return max value.
            """
        # alphabtic order index
        k = {}    #{'A':0, 'B':1, 'C':2, ...}
        for l,s in enumerate(self.alphabet):
            k[s] = l
        string = ' ' + self.string
        symbol = string[j]
        if symbol in k:
            emitt = self.emission_matrix[i*3-1][ k[symbol] ]
        # S_M(j),i
        temp_I = float('-inf')
        temp_M = float('-inf')
        if j > 1:
            temp_I = self.dp[i-1][0][j-1] + self.transit_matrix[i*3-2][i*3-1] + emitt #I(i-1) -> M(i)
            temp_M = self.dp[i-1][1][j-1] + self.transit_matrix[i*3-4][i*3-1] + emitt
        temp_D = self.dp[i-1][2][j-1] + self.transit_matrix[i*3-3][i*3-1] + emitt
        self.dp[i][1][j] = max(temp_I, temp_M, temp_D)

    def cal_D_i(self, i, j):
        """
            calculate the prob of  D_i state,
            i state, j string, return max value.
            """
        # alphabtic order index
        k = {}    #{'A':0, 'B':1, 'C':2, ...}
        for l,s in enumerate(self.alphabet):
            k[s] = l
        # S_D(j),i
        string = ' ' + self.string
        symbol = string[j]
        temp_I = float('-inf')
        temp_M = float('-inf')
        if j > 0:
            temp_I = self.dp[i-1][0][j] + self.transit_matrix[i*3-2][i*3]
            temp_M = self.dp[i-1][1][j] + self.transit_matrix[i*3-4][i*3]
        temp_D = self.dp[i-1][2][j] + self.transit_matrix[i*3-3][i*3]
        self.dp[i][2][j]= max(temp_I, temp_M, temp_D)

    def cal_I_i(self, i, j):
        """
            calculate the prob of  I_i state,
            i state, j string, return max value.
            """
        # alphabtic order index
        k = {}    #{'A':0, 'B':1, 'C':2, ...}
        for l,s in enumerate(self.alphabet):
            k[s] = l
        string = ' ' + self.string
        symbol = string[j]
        if symbol in k:
            emitt = self.emission_matrix[i*3+1][ k[symbol] ]
        # S_I(j),i
        temp_I = float('-inf')
        temp_M = float('-inf')
        if j > 1:
            temp_I = self.dp[i][0][j-1] + self.transit_matrix[i*3+1][i*3+1] + emitt
            temp_M = self.dp[i][1][j-1] + self.transit_matrix[i*3-1][i*3+1] + emitt
        temp_D = self.dp[i][2][j-1] + self.transit_matrix[i*3][i*3+1] + emitt
        self.dp[i][0][j] = max(temp_I, temp_M, temp_D)
            
    def profile_HMM(self):
        """
            return the transition matrix followed 
            by the emission matrix of HMM.
            """
        self.remove_col()
        self.transit_matrix_build()
        self.emission_matrix_build()
        self.pseudocount_build()
        
        state_num = 0                                   #number of M/D state
        state_symbol = []
        for i in range(1, len(self.hidden_path[0])-1):
            if type(self.hidden_path[0][i]) is str:
                state_num += 1
        for i in range(state_num):
            Insert_symbol = 'I'+ str(i+1)
            Match_symbol = 'M'+str(i+1)
            Dele_symbol ='D'+str(i+1)
            state_symbol += [Match_symbol]+[Dele_symbol]+[Insert_symbol]
        state_symbol = ['S','I0']+state_symbol+['E']
    
    def pseudocount_build(self):
        """
            Introduce pseudocounts sigma to transition matrix
            and emission matrix in logrithmic value.
            """
        state_num = len(self.transit_matrix)/3-1 #number of M/D/I triple state.
        #update S, I_0
        for i in range(2):
            temp_total = []
            for j in [1,2,3]:
                temp_total += [self.transit_matrix[i][j]+self.sigma]
            for j in [1,2,3]:
                self.transit_matrix[i][j] = math.log(temp_total[j-1]/sum(temp_total))
        #update E, M_n, D_n, I_n
        for i in range(3):
            temp_total = []
            for j in [1,2]:
                temp_total += [self.transit_matrix[-i-2][-j]+self.sigma]
            for j in [1,2]:
                self.transit_matrix[-i-2][-j] = math.log(temp_total[j-1]/sum(temp_total))
        #update M_i, D_i, I_i
        for i in range(state_num-1):
            temp_total = []
            for j in [0, 1, 2]:
                temp_total += [self.transit_matrix[3*i+2][3*i+4+j]+self.sigma]
            for j in [0, 1, 2]:
                self.transit_matrix[3*i+2][3*i+4+j] = math.log(temp_total[j]/sum(temp_total))
            temp_total = []
            for j in [0, 1, 2]:
                temp_total += [self.transit_matrix[3*i+3][3*i+4+j]+self.sigma]
            for j in [0, 1, 2]:
                self.transit_matrix[3*i+3][3*i+4+j] = math.log(temp_total[j]/sum(temp_total))
            temp_total = []
            for j in [0, 1, 2]:
                temp_total += [self.transit_matrix[3*i+4][3*i+4+j]+self.sigma]
            for j in [0, 1, 2]:
                self.transit_matrix[3*i+4][3*i+4+j] = math.log(temp_total[j]/sum(temp_total))
        
        #update emission matrix
        for i in range(state_num):
            temp_total = []
            for j in range(len(self.alphabet)):
                temp_total += [ self.emission_matrix[3*i+1][j]+self.sigma ]
            for j in range(len(self.alphabet)):
                self.emission_matrix[3*i+1][j] = math.log(temp_total[j]/sum(temp_total))
            temp_total = []
            for j in range(len(self.alphabet)):
                temp_total += [ self.emission_matrix[3*i+2][j]+self.sigma ]
            for j in range(len(self.alphabet)):
                self.emission_matrix[3*i+2][j] = math.log(temp_total[j]/sum(temp_total))
        #update I_n
        temp_total = []
        for j in range(len(self.alphabet)):
            temp_total += [ self.emission_matrix[-2][j]+self.sigma ]
        for j in range(len(self.alphabet)):
            self.emission_matrix[-2][j] = math.log(temp_total[j]/sum(temp_total))

    def emission_matrix_build(self):
        """
            return the emission matrix.
            """
        align_num = len(self.alignment)                 #number of alignment
        path_len  = len(self.hidden_path_insert[0])     #number of hidden path states
        state_num = 0                                   #number of M/D state
        for i in range(1, len(self.hidden_path[0])-1):
            if type(self.hidden_path[0][i]) is str:
                state_num += 1
        
        emission_matrix_row = [0 for x in range(len(self.alphabet))]
        emission_matrix = [ emission_matrix_row[:] for x in range(state_num*3+3) ]

        for j in range(path_len):                         #j current state number
            for i in range(align_num):                    #i current alignment number
                element = self.hidden_path_insert[i][j]   #current element
                cur_point = j//2                          #current state number 'i'-> M_i, D_i, I_i
                if j == 1:
                    if type(element) is list:
                        if len(element) > 0:
                            for symbol in element:
                                for k in range(len(self.alphabet)):
                                    if symbol == self.alphabet[k]:
                                        emission_matrix[1][k] += 1
                    else:
                        for k in range(len(self.alphabet)):
                            if element == self.alphabet[k]:
                                emission_matrix[1][k] += 1
            
                elif type(element) is list:
                    if len(element) > 0:
                        for symbol in element:
                            for k in range(len(self.alphabet)):
                                if symbol == self.alphabet[k]:
                                    emission_matrix[(cur_point-1)*3+2+2][k] += 1
                else:
                    for k in range(len(self.alphabet)):
                        if element == self.alphabet[k]:
                            emission_matrix[(cur_point-1)*3+2][k] += 1

        #convert the number to the fraction
        for emission_row in emission_matrix:
            row_temp = []
            for state in emission_row:
                if sum(emission_row) > 0 and state != 0:
                    row_temp.append(float(state)/sum(emission_row))
                else:
                    row_temp.append(state)
            self.emission_matrix.append(row_temp)
    
    def transit_matrix_build(self):
        """
            return the transition matrix.
            """
        align_num = len(self.alignment)         #number of alignment
        path_len  = len(self.hidden_path[0])    #number of hidden path states
        state_num = 0                           #number of M/D state
        for i in range(1, len(self.hidden_path[0])-1):
            if type(self.hidden_path[0][i]) is str:
                state_num += 1

        # transition matrix
        # [S -> I_0 -> (M_i, D_i, I_i) -> E ](#column) x number of state (#row)
        transit_matrix_row = [0,0] + [0 for x in range(state_num*3)] +[0]
        transit_matrix = [ transit_matrix_row[:] for x in range(state_num*3+3) ]
        
        # handle I_i case, insert []
        hidden_path_insert = [ [] for x in range(align_num) ]
        for i in range(align_num):
            for j in range(path_len-1):
                element = self.hidden_path[i][j]
                next_element = self.hidden_path[i][j+1]
                if (type(element) is str) and (type(next_element) is str):
                    hidden_path_insert[i].append(element)
                    hidden_path_insert[i].append([])
                else:
                    hidden_path_insert[i].append(element)
            hidden_path_insert[i].append('END')
        self.hidden_path_insert = hidden_path_insert[:][:]   #update hidden path with inerstion
        path_len  = len(hidden_path_insert[0])               #update, number of hidden path states
        
        # count the number of transition from hidden path
        # the transiton defined case by case, from state to states
        for j in range(path_len-1):                         #j current state number
            for i in range(align_num):                      #i current alignment number
                element = hidden_path_insert[i][j]             #current element
                next_element = hidden_path_insert[i][j+1]      #next element
            
                cur_point = j//2                              #current state number 'i'-> M_i, D_i, I_i
                # transition from the start state
                if j == 0:
                    if type(next_element) is list and len(next_element) > 0:
                        num_d = 0
                        for symbol in next_element:
                            if symbol == '-' :
                                num_d += 1
                        if num_d < len(next_element):
                            transit_matrix[0][1] += 1
                        elif hidden_path_insert[i][2] == '-':
                            transit_matrix[0][3] += 1
                        elif hidden_path_insert[i][2] in self.alphabet:
                            transit_matrix[0][2] += 1
                        elif hidden_path_insert[i][2] == 'END':
                            transit_matrix[0][-1] += 1
                    else:
                        next_element = hidden_path_insert[i][2]
                        if next_element == '-' :
                            transit_matrix[0][3] += 1
                        elif next_element == 'END':
                            transit_matrix[0][-1] += 1
                        elif next_element in self.alphabet:
                            transit_matrix[0][2] += 1
            
                # transition I_0 case
                elif j == 1 and len(element) > 0:
                    num_d = 0    # number of '-'
                    for symbol in element:
                        if symbol == '-' :
                            num_d += 1
                    if num_d < len(element):
                        # I_i -> I_i
                        transit_matrix[1][1] += len(element)-num_d-1
                        # I_i -> M_i/D_i
                        if hidden_path_insert[i][j+1] == '-':
                            transit_matrix[1][3] += 1
                        else:
                            transit_matrix[1][2] += 1
            
            
                #transition from a deletion state
                elif element == '-':
                    if type(next_element) is list and len(next_element) > 0:
                        num_d = 0
                        for symbol in next_element:
                            if symbol == '-' :
                                num_d += 1
                        if num_d < len(next_element):
                            transit_matrix[3*(cur_point-1)+2+1][3*(cur_point-1)+2+2] += 1
                        elif hidden_path_insert[i][j+2] == '-':
                            transit_matrix[3*(cur_point-1)+2+1][3*(cur_point)+2+1] += 1
                        elif hidden_path_insert[i][j+2] in self.alphabet:
                            transit_matrix[3*(cur_point-1)+2+1][3*(cur_point)+2] += 1
                        elif hidden_path_insert[i][j+2] == 'END':
                            transit_matrix[3*(cur_point-1)+2+1][-1] += 1
                    else:
                        next_element = hidden_path_insert[i][j+2]
                        if next_element == '-' :
                            transit_matrix[3*(cur_point-1)+2+1][3*(cur_point)+2+1] += 1
                        elif next_element == 'END':
                            transit_matrix[3*(cur_point-1)+2+1][-1] += 1
                        elif next_element in self.alphabet:
                            transit_matrix[3*(cur_point-1)+2+1][3*(cur_point)+2] += 1
        
                #transition from a insertion state
                elif type(element) is list and len(element) > 0:
                    num_d = 0    # number of '-'
                    for symbol in element:
                        if symbol == '-' :
                            num_d += 1
                    if num_d < len(element):
                        # I_i -> I_i, len(element) - num -1, more then one element
                        transit_matrix[3*(cur_point-1)+2+2][3*(cur_point-1)+2+2] += len(element)-num_d-1
                        # I_i -> M_i/D_i/'END'
                        if  hidden_path_insert[i][j+1] == 'END':
                            transit_matrix[3*(cur_point-1)+2+2][-1] += 1
                        elif hidden_path_insert[i][j+1] == '-':
                            transit_matrix[3*(cur_point-1)+2+2][3*(cur_point)+2+1] += 1
                        elif hidden_path_insert[i][j+1] in self.alphabet:
                            transit_matrix[3*(cur_point-1)+2+2][3*(cur_point)+2] += 1

                #transition from a M state xxx-> I/M/D/E
                elif element in self.alphabet:
                    if type(next_element) is list and len(next_element) > 0:
                        num_d = 0
                        for symbol in next_element:
                            if symbol == '-' :
                                num_d += 1
                        if num_d < len(next_element):
                            transit_matrix[3*(cur_point-1)+2][3*(cur_point-1)+2+2] += 1
                        elif hidden_path_insert[i][j+2] == '-':
                            transit_matrix[3*(cur_point-1)+2][3*(cur_point)+2+1] += 1
                        elif hidden_path_insert[i][j+2] in self.alphabet:
                            transit_matrix[3*(cur_point-1)+2][3*(cur_point)+2] += 1
                        elif hidden_path_insert[i][j+2] == 'END':
                            transit_matrix[3*(cur_point-1)+2][-1] += 1
                    else:
                        next_element = hidden_path_insert[i][j+2]
                        if next_element == '-' :
                            transit_matrix[3*(cur_point-1)+2][3*(cur_point)+2+1] += 1
                        elif next_element == 'END':
                            transit_matrix[3*(cur_point-1)+2][-1] += 1
                        elif next_element in self.alphabet:
                            transit_matrix[3*(cur_point-1)+2][3*(cur_point)+2] += 1

        #convert the number to the fraction
        for transit_row in transit_matrix:
            row_temp = []
            for state in transit_row:
                if sum(transit_row) > 0 and state != 0:
                    row_temp.append(float(state)/sum(transit_row))
                else:
                    row_temp.append(state)
            self.transit_matrix.append(row_temp)

    def remove_col(self):
        """
            Remove columns if the fraction of "-" exceed theta
            construct a hidden_path.
            """
        string_len = len(self.hidden_path[0])
        align_num = len(self.alignment)
        tol_col = len(self.hidden_path)*self.theta #exceed limit
        pos_d = [0 for x in range(string_len)]    #i postion and its number of '-'
        
        for i in range(string_len):
            for j in range(align_num):
                if self.alignment[j][i] == '-':
                    pos_d[i] += 1
        
        self.hidden_path = [ [] for x in range(align_num) ]
        for pos in range(len(pos_d)):
            if pos_d[pos] >= tol_col:
                for i in range(align_num):
                    if len(self.hidden_path[i]) > 0 and type(self.hidden_path[i][-1]) is list:
                        self.hidden_path[i][-1].append(self.alignment[i][pos])
                    else:
                        self.hidden_path[i].append([ self.alignment[i][pos] ])
            else:
                for i in range(align_num):
                    self.hidden_path[i].append(self.alignment[i][pos])
        # add the start and end state
        for i in range(align_num):
            self.hidden_path[i] = ['START'] + self.hidden_path[i] + ['END']

def compare_float(num_1, num_2):
    """
        helper function, compare two number
        return True if equal
        """
    if round(num_1, 3) == round(num_2, 3):
        return True
    return False


if __name__=="__main__":
    test = TxtExtract('dataset_11632_6.txt')
    test.sequence_alignment()
    #test.profile_HMM()
    #test=TXT_Extract('dataset_10928_3.txt')



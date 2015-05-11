import math
import numpy as np

class TXT_Extract:
    def __init__(self,dna_file):
        f = open(dna_file,'r')
        
        self.lists = []
        for line in f:
            line = line.strip()
            line = line.split(" ")
            for n,x in enumerate(line):
                line[n] = float(x)
            self.lists.append(line)
    
        self.num_center, self.dimensions = [int(x) for x in self.lists.pop(0)]
        self.beta = self.lists.pop(0)[0]

        self.matrix = []
        self.center  = []
        
        for i in range(self.num_center):
            self.center += [self.lists[i]]  # cooridnate center

        f.close()

    def search_cluster(self):
        self.matrix = partition(self.lists, self.center, self.beta)
    
    def search_center(self):
        return Mstep(self.lists, self.matrix, self.dimensions)

    def iteration(self):
        center = []
        prev_center = self.center
        
        while compare_float(prev_center, center):
            prev_center = self.center[:]
            self.search_cluster()
            center = self.search_center()
            self.center = center[:]

    def print_center(self):
        self.iteration()
        for center in self.center:
            for j in center:
                print "%.3f" %( j ),
            print ""

def compare_float(list1, list2):
    if sorted(round_up(list1)) != sorted(round_up(list2)):
        return True
    return False

def round_up(list):
    round_list = []
    for i in list:
        temp = []
        for j in i:
            temp += [round(j, 4)]
        round_list += [temp]
    return round_list

def partition(datapoint, center, stiff):
    matrix = []  ##HiddenMatrix_ij -> matrix[i][j]
    matrix_j = []

    for j in datapoint:
        matrix_jk = []
        for k in center:
            matrix_jk += [math.exp(-stiff*euclidean(k,j))]
        matrix_j += [sum(matrix_jk)]
    

    for i in center:
        matrix_i = []
        for n,j in enumerate(datapoint):
            matrix_i += [ math.exp(-stiff*euclidean(i,j))/matrix_j[n]]
        matrix += [matrix_i]

    return matrix

def Mstep(datapoint, matrix, dimension):
    center = []
    
    for m in matrix:
        center_i = []
        for i in range(dimension):
            weight = []
            for n,j in enumerate(datapoint):
                weight += [ j[i] * m[n] ]
            center_i += [sum(weight)/sum(m)]
        center += [center_i]
    
    return center

def euclidean(c1, c2):
    distance =0
    for i,j in zip(c1,c2):
        distance += (i-j)**2
    return math.sqrt(distance)

if __name__=="__main__":
    #test=TXT_Extract('test.txt')
    test=TXT_Extract('dataset_10933_7.txt')
    test.print_center()


import math

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
    
        self.num_center = int(self.lists.pop(0)[0])
        
        self.cluster = {}
        self.center  = []
        
        for i in range(self.num_center):
            self.center += [self.lists[i]]  # cooridnate center
        
        f.close()

    def search_cluster(self):
        for i in self.lists:
            center = min_center(i, self.center)
            self.cluster[ tuple(center) ] = self.cluster.get(tuple(center), []) + [i]
    
    def search_center(self):
        center = []
        for i in self.center:
            center += [ gravity(self.cluster[tuple(i)]) ]
        return center

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

def gravity(datapoint):
    center = []
    num = len(datapoint[0])
    for i in range(num):
        distance = []
        for j in datapoint:
            distance += [ j[i] ]
        center += [ sum(distance)/float(len(distance))]
    return center

def min_center(datapoint, center):
    """calculate the d(DataPoint,Center),
        return the minimum distance center"""
    distance = []
    for i in center:
        distance += [ euclidean(datapoint,i) ]
    index = distance.index(min(distance))
    return center[index]


def euclidean(c1, c2):
    distance =0
    for i,j in zip(c1,c2):
        distance += (i-j)**2
    return math.sqrt(distance)

if __name__=="__main__":
    #test=TXT_Extract('test.txt')
    test=TXT_Extract('dataset_10928_3.txt')
    test.print_center()


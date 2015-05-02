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
        
        self.center = [self.lists[1]]  # cooridnate center
        self.num_center = int(self.lists.pop(0)[0])

        f.close()

    def search_center(self):
        distance = []
 
        for i in self.lists:
            if i not in self.center:
                distance += [ self.min_center(i,self.center) ]
            elif i in self.center:
                distance += [0.0]
    
        return search_max(distance)
                
    def print_center(self):
        while len(self.center) < self.num_center:
            center_lists = self.search_center()
            for i in center_lists:
                self.center += [ self.lists[i] ]
    
        for i,j in enumerate(self.center):
            if i< self.num_center:
                for x in j:
                    print x,
                print ""

    def min_center(self, datapoint, center):
        """calculate the d(DataPoint,Center),
            return the minimum distance"""
        distance = []
        for i in center:
            distance += [ euclidean(datapoint,i) ]

        return min(distance)


def euclidean(c1, c2):
    distance =0
    for i,j in zip(c1,c2):
        distance += (i-j)**2
    
    return math.sqrt(distance)

def search_max(lists):
    d_max= max(lists)
    
    return [i for i,j in enumerate(lists) if j == d_max]


if __name__=="__main__":
    #test=TXT_Extract('test.txt')
    test=TXT_Extract('dataset_10926_14.txt')
    test.print_center()



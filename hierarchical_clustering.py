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
    
        self.index = {}
        self.cluster = {}
        self.dim = int (self.lists.pop(0)[0] )
        f.close()

    def print_cluster(self):
        self.index = {x:x+1 for x in range(self.dim)}
        self.cluster = {x:set([x]) for x in range(self.dim+1) }
        
        for n in range(self.dim-1):
            search_i, search_j = search_min(self.lists)
            interpret_i = self.index[search_i]
            interpret_j = self.index[search_j]
            
            self.cluster[interpret_i].update(self.cluster[interpret_j])
            self.cluster[interpret_j].update(self.cluster[interpret_i])
            C1 = float(len(self.cluster[interpret_i]))
            C2 = float(len(self.cluster[interpret_i]))
            
            self.dim = len(self.lists)
            new_lists = [row[:] for row in self.lists]
            self.lists = []
            for i in range(self.dim):           #update
                temp_list = []
                for j in range(self.dim):
                    if i == search_i and j != search_j:
                        temp_list += [ (new_lists[search_i][j]*C1+new_lists[search_j][j]*C2)/(C1+C2) ]
                    elif i != search_j and j == search_i:
                        temp_list += [ (new_lists[i][search_i]*C1+new_lists[i][search_j]*C2)/(C1+C2) ]
                    elif i != search_j and j != search_j:
                        temp_list += [ new_lists[i][j] ]
                if len(temp_list) != 0:
                    self.lists.append(temp_list)
                        
            temp_dic = self.index.copy()
            for m in range(search_j,self.dim-1):
                self.index[m] = temp_dic[m+1]

            for i in self.cluster[interpret_i]:
                print i,
            print ""



def search_min(lists):
    min = float("inf")
    search_i = 0
    search_j = 1
    lists_range = len(lists)
    
    for i in range(lists_range):
        for j in range(i+1, lists_range):
            if lists[i][j] < min:
                search_i = i
                search_j = j
                min = lists[i][j]
    if search_j < search_i:
        return (serch_j, search_i)
    return (search_i, search_j)


if __name__=="__main__":
    test=TXT_Extract('dataset_10934_7.txt')
    test.print_cluster()



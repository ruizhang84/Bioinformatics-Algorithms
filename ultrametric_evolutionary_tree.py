class TXT_Extract:
    def __init__(self,dna_file):
        f = open(dna_file,'r')
        
        self.lists = []
        for line in f:
            line = line.strip()
            line = line.split("\t")
            for n,x in enumerate(line):
                line[n] = float(x)
            self.lists.append(line)
        
        self.age = {}  # node => age
        self.root = {} # node1, node2 => node
        self.count_node = {} # leaf by node
        
        self.num_node=int (self.lists.pop(0)[0] )
        
        f.close()

    def build_tree(self):
        self.age = { x: 0 for x in range(self.num_node) }
        self.count_node = { x: 1 for x in range(self.num_node) }
        interpret = { x: x for x in range(self.num_node) }
        
        for n in range(self.num_node-1):
            search_i=0; search_j=1
            small = self.lists[0][1]
            lists_range = len(self.lists)
            
            for i in range(lists_range):    #search cluster
                for j in range(i+1, lists_range):
                    if self.lists[i][j] < small:
                        search_i = i
                        search_j = j
                        small = self.lists[i][j]


            interpret_i = interpret[search_i]  #covert i to origninal order
            interpret_j = interpret[search_j]

            self.age[self.num_node+n] = small/2.0    #build dependence in trees
            self.root[interpret_i] = self.num_node+n
            self.root[interpret_j] = self.num_node+n
            self.count_node[self.num_node+n] = self.count_node[interpret_i]+self.count_node[interpret_j]
            
            new_lists = [row[:] for row in self.lists]
            for i in range(lists_range):           #update value
                for j in range(lists_range):
                    if i == search_i or i ==  search_j:
                        if j == search_i or j == search_j:
                            new_lists[i][j] = 0.0
                        else:
                            C1 = float(self.count_node[ interpret_i ])
                            C2 = float(self.count_node[ interpret_j ])
                            new_lists[i][j] = (self.lists[search_i][j]*C1+self.lists[search_j][j]*C2)/(C1+C2)
                    else:
                        if j == search_i or j == search_j:
                            C1 = float(self.count_node[ interpret_i ])
                            C2 = float(self.count_node[ interpret_j ])
                            new_lists[i][j] = (self.lists[i][search_i]*C1+self.lists[i][search_j]*C2)/(C1+C2)
                            


            self.lists=[]
            for i in range(lists_range):           #update tree
                temp_list = []
                for j in range(lists_range):
                    if i != search_j and j != search_j:
                        temp_list += [ new_lists[i][j] ]
                if len(temp_list) != 0:
                    self.lists.append(temp_list)
            
            interpret[search_i] = self.num_node+n   # update list order
            temp_dic = interpret.copy()
            for m in range(search_j,lists_range-1):
                interpret[m] = temp_dic[m+1]
                

    def print_tree(self):
        leaf_set = {}
        for i in self.root:
            root = self.root[i]
            leaf_set[i] = leaf_set.get(i,[])+[root]
            leaf_set[root] = leaf_set.get(root,[])+[i]
        
        for i in leaf_set:
            leaf_group = leaf_set[i]
            for leaf in leaf_group:
                edge = abs(self.age[leaf]-self.age[i])
                if edge != 0:
                    print "%d->%d:%.3f" %(i,leaf,edge)



if __name__=="__main__":
    test=TXT_Extract('dataset_10332_8.txt')
    test.build_tree()
    test.print_tree()



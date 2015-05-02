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
        
        self.tree = {}  # node1 <-> node2 => length
        self.root = {} # node <= node1, node2
        self.leaf = {}# node => node1, node2
        
        self.num_node = int(self.lists.pop(0)[0])
        
        f.close()

    def build_tree(self):
        interpret = { x: x for x in range(len(self.lists)) }
        search_i = 0
        search_j = 1
        
        for n in range(len(self.lists)-3) :
            (search_i, search_j) = search_min(self.lists)# neighbor search
            (limb_1,limb_2) = self.nj_tree(search_i, search_j)

            self.root[interpret[search_i]] = self.num_node+n
            self.root[interpret[search_j]] = self.num_node+n
            self.leaf[self.num_node+n] = (interpret[search_i],interpret[search_j])
            self.tree[interpret[search_i]] = limb_1
            self.tree[interpret[search_j]] = limb_2
        
            interpret[search_i] = self.num_node+n   # update list order
            temp_dic = interpret.copy()
            for m in range(search_j,len(self.lists)):
                interpret[m] = temp_dic[m+1]

        search_i,search_j = [x for x in range(len(self.lists)) if x != search_i] ##the other side node
        (limb_1,limb_2) = self.nj_tree(search_i, search_j)
        self.root[interpret[search_i]] = self.num_node*2-3
        self.root[interpret[search_j]] = self.num_node*2-3
        self.leaf[self.num_node*2-3] = (interpret[search_i],interpret[search_j])
        self.tree[interpret[search_i]] = limb_1
        self.tree[interpret[search_j]] = limb_2
        
        self.root[self.num_node*2-3] = self.num_node*2-4 #node1 <==> node2
        self.root[self.num_node*2-4] = self.num_node*2-3
        self.tree[self.num_node*2-3] = self.lists[0][1]
        self.tree[self.num_node*2-4] = self.lists[0][1]
                
    
    def print_tree(self):
        
        for i in range(self.num_node*2-2):
            if i in self.leaf:
                print "%d->%d:%.3f" %(i,self.leaf[i][0],self.tree[self.leaf[i][0]])
                print "%d->%d:%.3f" %(i,self.leaf[i][1],self.tree[self.leaf[i][1]])
            if i in self.root:
                print "%d->%d:%.3f" %(i, self.root[i],self.tree[i])
    
    

    def nj_tree(self,search_i,search_j):
        """D_k,m = (1/2)(D_k,i + D_k,j - Di,j)"""
        new_lists = []
        lists_range = len(self.lists)
        for i in range(lists_range):
            row_lists = []
            for j in range(lists_range):
                if i == search_i :
                    if j == search_i:
                        row_lists.append(0)
                    elif j == search_j:
                        continue
                    else:
                        d_km = (self.lists[search_i][j] + self.lists[search_j][j] - self.lists[search_i][search_j])/2
                        row_lists.append(d_km)
                elif i == search_j:
                    continue
                elif j == search_i:
                    d_km = (self.lists[i][search_i] + self.lists[i][search_j] - self.lists[search_i][search_j])/2
                    row_lists.append(d_km)
                elif j == search_j:
                    continue
                else:
                    row_lists.append(self.lists[i][j])
            if len(row_lists) != 0:
                new_lists.append(row_lists)
        
        total_dist = total_distance(self.lists)
        d_ij = (total_dist[search_i] - total_dist[search_j])/(float(lists_range) - 2.0)
        limb_1 = (self.lists[search_i][search_j] + d_ij)/2.0
        limb_2 = (self.lists[search_i][search_j] - d_ij)/2.0
        
        self.lists = [ row[:] for row in new_lists ]
        return limb_1, limb_2


def total_distance(lists):
    total_dist = [ float(sum(lists[x])) for x in range(len(lists)) ]
    
    return total_dist

def search_min(lists):
    search_i = 0
    search_j = 1
    
    lists_t = nj_matrix(lists) ##neighbor-joining matrix
    
    small = lists_t[0][1]
    for i in range(len(lists_t)):
        for j in range(i+1, len(lists_t)):
            if lists_t[i][j] < small:
                search_i = i
                search_j = j
                small = lists_t[i][j]

    return (search_i,search_j)

def nj_matrix(lists):
    """d_ij =(n-2)*d_ij-total_i-total_j"""
    nj_len = len(lists)
    nj_lists = [ row[:] for row in lists ]

    total_dist = total_distance(lists)
    for i in range(nj_len):
        for j in range(nj_len):
            if i != j:
                nj_lists[i][j] = (nj_len-2)*lists[i][j]-total_dist[i]-total_dist[j]

    return nj_lists


if __name__=="__main__":
    test=TXT_Extract('dataset_10333_6.txt')
    test.build_tree()
    test.print_tree()



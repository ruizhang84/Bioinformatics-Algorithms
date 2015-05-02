import heapq
import itertools

class TXT_Extract:
    def __init__(self,dna_file):
        f = open(dna_file,'r')
        self.lists = []
        for line in f:
            line = line.strip()
            line = line.split(" ")
            for n,x in enumerate(line):
                line[n] = int(x)
            
            self.lists.append(line)
    
        self.edge = {}        #key (1,2) => length
        self.pair = {}        #connection 1=>2
        self.node = set()     #node

        self.limb_pair ={}    # n => i, j
        self.limb_group = {}  # distance for n => n, i,j pairs
        
        self.num_node=int(self.lists.pop(0)[0])
        f.close()

    def search(self,node_1,node_2):
        """given node_1,node_2,dynamically search for the path"""
        index = {}
        search_path = []
        node_next = node_1
        
        for node in self.node:
            index[node] = 0
        
        while node_next!=node_2:  #depth first search
            if index[node_next] < len(self.pair[node_next]):
                index_n = index[node_next]
                index[node_next] += 1
                if self.pair[node_next][index_n] not in search_path:
                    search_path += [node_next]
                    node_next = self.pair[node_next][index_n]
            else:
                node_next = search_path.pop()
        
        search_path += [node_2]
        return search_path

    def cal_length(self,cal_node,node_list):
        """given a node,calculate limb length(node) and return i,j"""
        dis = []
        n = cal_node
        limb_legnth = 0
        
        for (i,j) in itertools.combinations(node_list,2):
            distance_ni = self.lists[n][i]
            distance_nj = self.lists[n][j]
            distance_ij = self.lists[i][j]
                    
            limb = (distance_ni+distance_nj-distance_ij)/2
            heapq.heappush(dis,(limb,n,i,j))

        (limb_length,n,i,j) = heapq.heappop(dis)

        self.limb_pair[n] = (i,j)
        distance_bald_i = self.lists[n][i]-limb_length
        distance_bald_j = self.lists[n][j]-limb_length
        
        return (i,j,limb_length, distance_bald_i,distance_bald_j)

    def sub_tree(self):
        """decrease the distance matrix to 2-D
            only keep 1st 2nd column/row    """
        node_list=[ x for x in range(self.num_node) ]

        while len(node_list) >2 :
            node = node_list.pop(1)
            
            (i,j,limb, distance_bald_i,distance_bald_j) = self.cal_length(node,node_list)
            self.limb_group[(node,i)] = distance_bald_i
            self.limb_group[(node,j)] = distance_bald_j
            self.limb_group[(node,node)] = limb

        node_distance = self.lists[ node_list[0] ][ node_list[1]  ]
        return node_list, node_distance

    def attach_limb(self):
        """attach node back to the tree"""
        (node_list,distance) = self.sub_tree()  #initilize
        self.node.update(node_list)
 
        node_1,node_2=node_list                 #1st edge
        self.pair[node_1] = [node_2]
        self.pair[node_2] = [node_1]
        self.edge[(node_1,node_2)] = distance
        self.edge[(node_2,node_1)] = self.edge[(node_1,node_2)]

        for n in range(self.num_node-2):
            node_in  = self.num_node-2-n             #exist node inserted
            node_new = self.num_node+n               #new   node inserted
            
            i,j = self.limb_pair[node_in]
            path = self.search(i,j)                 #i -> j
            distance_bald_i = self.limb_group[(node_in,i)]
            distance_bald_j = self.limb_group[(node_in,j)]

            distance_sum=0
            for k in range(len(path)):
                node_pre = path[k]
                node_pass = path[k+1]

                distance_sum += self.edge[(node_pre,node_pass)]
                if distance_bald_i < distance_sum:
                    self.pair[node_new] = [node_pre]
                    self.pair[node_new] += [node_pass]
                    self.pair[node_pre] += [node_new]
                    self.pair[node_pre].remove(node_pass)
                    self.pair[node_pass] += [node_new]
                    self.pair[node_pass].remove(node_pre)
                    self.edge[(node_pass,node_new)] = distance_sum-distance_bald_i
                    self.edge[(node_new,node_pass)] = self.edge[(node_pass,node_new)]
                    self.edge[(node_pre,node_new)] = self.edge[(node_pre,node_pass)]-self.edge[(node_new,node_pass)]
                    self.edge[(node_new,node_pre)] = self.edge[(node_pre,node_new)]
                    del self.edge[(node_pre,node_pass)]
                    del self.edge[(node_pass,node_pre)]
                    break
    
            self.pair[node_in]  = [node_new]
            self.pair[node_new] += [node_in]
            self.edge[(node_in,node_new)] = self.limb_group[(node_in,node_in)]
            self.edge[(node_new,node_in)] = self.edge[(node_in,node_new)]
            self.node.update([node_in,node_new])

    def print_tree(self):
        for i in self.node:
            for j in self.pair[i]:
                print "%d->%d:%d" %(i,j,self.edge[(i,j)])

if __name__=="__main__":
    test=TXT_Extract('dataset_10330_6.txt')
    test.attach_limb()
    test.print_tree()



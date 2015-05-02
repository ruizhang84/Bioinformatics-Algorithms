class TXT_Extract:
    def __init__(self,dna_file):
        import re
        regx=re.compile(r"(\d+)->(\d+):(\d+)")
        self.node=[]  #node [1,2,3..]
        self.inner_node=[]
        self.vertice=[]
        self.edge={}  #key (1,2) => length
        self.pair={}  #connection 1=>2

        f=open(dna_file,'r')
        for line in open(dna_file):
            line=f.readline()
            line=line.strip()
            
            r=regx.search(line)
            if r:
                node_1=int(r.group(1))
                node_2=int(r.group(2))
                length=int(r.group(3))
                self.edge[(node_1,node_2)]=length
                self.node+=[node_1,node_2]
                
                if node_1 not in self.pair:
                    self.pair[node_1]=[node_2]
                else:
                    self.pair[node_1]+=[node_2]
                    self.inner_node+=[node_1]
        
            else:
             self.num_vertice=int(line)
        
        self.node=set(self.node)
        self.inner_node=set(self.inner_node)
        f.close()

    def search(self,node_1,node_2):
        index={}
        search_path=[]
        node_next=node_1
        
        for node in self.node:
            index[node]=0
        
        while node_next!=node_2:  #depth first search
            if index[node_next]<len(self.pair[node_next]):
                index_n=index[node_next]
                index[node_next]+=1
                if self.pair[node_next][index_n] not in search_path:
                    search_path+=[node_next]
                    node_next=self.pair[node_next][index_n]
            else:
                node_next=search_path.pop()
        
        search_path+=[node_2]

        return search_path

    def cal_distance(self,node_1,node_2):
        distance=0
        
        search_path=self.search(node_1,node_2)
        for i in range(len(search_path)-1):
            node=search_path[i]
            node_next=search_path[i+1]
            distance+=self.edge[(node,node_next)]

        return distance
    
    
    def build_matrix(self):
        for node_1 in range(self.num_vertice):
            for node_2 in range(self.num_vertice):
                print self.cal_distance(node_1, node_2),
            print ''

if __name__=="__main__":
    test=TXT_Extract('dataset_10328_11.txt')
    test.build_matrix()

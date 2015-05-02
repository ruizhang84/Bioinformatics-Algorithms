class TXT_Extract:
    def __init__(self,dna_file):
        f=open(dna_file,'r')
        self.lists=[]
        for line in open(dna_file):
            line=f.readline()
            line=line.strip()
            line=line.split(" ")
            for n,x in enumerate(line):
                line[n]=int(x)
            
            self.lists.append(line)
        
        self.num_node=int(self.lists.pop(0)[0])
        self.cal_node=int(self.lists.pop(0)[0])
        
        f.close()
    def cal_length(self):
        import heapq
        dis=[]
        limb_legnth=0
        for i in range(self.num_node):
            for j in range(self.num_node):
                if i != self.cal_node and j !=self.cal_node and i !=j:
                    distance_ni=self.lists[self.cal_node][i]
                    distance_nj=self.lists[self.cal_node][j]
                    distance_ij=self.lists[i][j]
                    
                    limb=(distance_ni+distance_nj-distance_ij)/2
                    heapq.heappush(dis,limb)

        limb_legnth=heapq.heappop(dis)
        print limb_legnth
        return limb_legnth


if __name__=="__main__":
    test=TXT_Extract('dataset_10329_11.txt')
    test.cal_length()

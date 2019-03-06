import numpy as np
def process_file(filename):
    """
    data format:
    4
    0   13  21  22
    13  0   12  13
    21  12  0   13
    22  13  13  0
    """
    with open(filename) as f:
        for i, line in enumerate(f):
            line = line.strip()
            if i == 0:
                dim = int(line)
                matrix = np.zeros((dim, dim), dtype=np.int32)
            else:
                rowData = np.array(line.split()).astype(np.int32)
                matrix[i-1] += rowData
    return (dim, matrix)

def limbLength(matrix, j):
    """
    Compute the Limb Length by
    limb = min_ik(d_ij + d_jk - d_ik)/2
    """
    limb = np.inf
    mask = [i for i in range(matrix.shape[0]) if i != j]
    for k in range(matrix.shape[0]):
        if k == j:
            continue
        limb = min(limb, min(matrix[mask, j]-matrix[mask,k]+matrix[j,k]))
    return limb//2

def find(matrix):
    """
    find i, k s.t. Di,k = Di,n + Dn,k
    """
    for k in range(matrix.shape[0]-1):
        arr = matrix[k] - matrix[-1]
        index = np.where(arr == matrix[k, -1])
        if len(index[0]) > 0:
            return (index[0][0], k)
    return None

def nearest(edge, weight, x, i, k):
    """
    find the nearest two nodes on path i->k
    to insert new node, BFS
    """
    queue = [[i]]
    visited = set([i])
    findPath = []
    while len(queue) > 0:
        path = queue.pop()
        node = path[-1]
        visited.add(node)
        if node == k:
            findPath = path
            break
        for next_node in edge[node]:
            if next_node not in visited:
                queue.append(path+[next_node])
    
    # distance
    dist = 0
    for k in range(len(findPath)-1):
        i, j = findPath[k], findPath[k+1]
        if dist+weight[(i, j)] > x:
            return (i, j, x-dist, dist+weight[(i, j)]-x)
        dist += weight[(i, j)]
    
    

def additivePhylogeny(matrix, n, inner_n):
    """
    finds the simple tree fitting an n x n 
    additive distance matrix D. 
    Tree: Edge(u, v) as [... u-th[v, ...]..]
    Weight(u, v) as {(u, v): w}
    inner_n, count the index of inner node
    """
    if n == 2:
        edge = {}
        edge[0] = [1]
        edge[1] = [0]
        weight = {}
        weight[(0, 1)] = matrix[0, 1]
        weight[(1, 0)] = matrix[0, 1]
        return (edge, weight, inner_n)

    limb = limbLength(matrix, n-1)
    matrix[:-1,-1] -= limb
    matrix[-1,:-1] -= limb
    #i, n, k three leaves such that Di,k = Di,n + Dn,k
    i, k = find(matrix)
    x = matrix[i, -1]
    # remove row n and column n from D
    edge, weight, inner_n = additivePhylogeny(matrix[:-1, :-1], n-1, inner_n)
    # the (potentially new) node in T at distance x from i on the path between i and k
    
     # find the insert node
    i_near, k_near, i_x, n_x = nearest(edge, weight, x, i, k)
    new_node = i_near

    # check if need to creat new node
    if i_x != 0:
        new_node = inner_n
        inner_n += 1
        # insert between i, k
        edge[i_near].remove(k_near)
        edge[k_near].remove(i_near)
        edge[i_near].append(new_node)
        edge[k_near].append(new_node)
        edge[new_node] = [i_near, k_near]

        weight[(new_node, i_near)] = i_x
        weight[(i_near, new_node)] = i_x
        weight[(new_node, k_near)] = n_x
        weight[(k_near, new_node)] = n_x
        del weight[(i_near, k_near)]
        del weight[(k_near, i_near)]
    # add leaf n back to T by creating a limb (v, n) of length limbLength
    edge[new_node].append(n-1)
    edge[n-1] = [new_node]
    weight[(n-1, new_node)] = limb
    weight[(new_node, n-1)] = limb
    return (edge, weight, inner_n)


def printOut(edge, weight):
    for i in sorted(edge):
        for j in sorted(edge[i]):
            print (str(i)+"->"+str(j)+":"+str(weight[(i, j)]))


n, matrix = process_file("rosalind_ba7c.txt")
edge, weight, _ = additivePhylogeny(matrix, n, n)
printOut(edge, weight)
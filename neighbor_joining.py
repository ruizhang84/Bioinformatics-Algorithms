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
                matrix = {}
            else:
                row = line.split()
                matrix[i-1] = {}
                for j in range(len(row)):
                    matrix[i-1][j] = float(row[j])
    return (dim, matrix)

def neighborJoining(matrix, n):
    if n == 2:
        edge = {}
        weight = {}
        for i in matrix:
            for j in matrix[i]:
                if  j == i:
                    continue
                edge[i] = [j]
                edge[j] = [i]
                weight[(i, j)] = matrix[i][j]
                weight[(j, i)] = matrix[i][j]
        return (edge, weight)
    # Construct neighbor-joining matrix D* from D.
    totalDistance = {}
    for i in matrix:
        totalDistance[i] = 0
        for j in matrix[i]:
            totalDistance[i] += matrix[i][j]
    MatrixD = {}
    for i in matrix:
        MatrixD[i] = {}
        for j in matrix[i]:
            if i == j:
                MatrixD[i][j] = 0
            else:
                MatrixD[i][j] = (n-2)*matrix[i][j] - totalDistance[i] - totalDistance[j]
    # Find a minimum element D*i,j of D*.
    mini = float("inf")
    _i = _j = None
    for i in MatrixD:
        for j in MatrixD[i]:
            if MatrixD[i][j] < mini:
                mini = MatrixD[i][j]
                _i, _j = i, j
    # Compute delta(i,j) = (TotalDistanceD(i) –TotalDistanceD(j)) / (n – 2).
    delta = (totalDistance[_i] - totalDistance[_j])/(n-2.0)
    # Set LimbLength(i) equal to 1/2(delta_i,j + delta_i,j) and LimbLength(j) equal to 1/2(Di,j – delta_i,j).
    limbLength_i = (matrix[_i][_j] + delta)/2.0
    limbLength_j = (matrix[_i][_j] - delta)/2.0
    # Form a matrix D' by removing i-th and j-th row/column from D and adding an m-th row/column 
    # such that for any k, Dk,m = (Di,k + Dj,k – Di,j) / 2.
    rows = list(matrix)
    m = max(matrix) + 1
    matrix[m] = {}
    matrix[m][m] = 0.0
    for k in rows:
        if k != _i and k != _j:
            matrix[m][k] = (matrix[_i][k] + matrix[_j][k] - matrix[_i][_j])/2.0
            matrix[k][m] = matrix[m][k]
            del matrix[k][_i]
            del matrix[k][_j]
    del matrix[_i]
    del matrix[_j]
    # Apply NeighborJoining to D' to obtain Tree(D').
    edge, weight = neighborJoining(matrix, n-1)
    # Reattach limbs of i and j to obtain Tree(D).
    edge[_i] = [m]
    edge[m].append(_i)
    edge[_j] = [m]
    edge[m].append(_j)
    weight[(_i, m)] = limbLength_i
    weight[(m, _i)] = limbLength_i
    weight[(_j, m)] = limbLength_j
    weight[(m, _j)] = limbLength_j
    return edge, weight




n, matrix = process_file("rosalind_ba7e.txt")
edge, weiight = neighborJoining(matrix, n)
for i in sorted(edge):
    for j in edge[i]:
        w = "%.3f" %weiight[(i, j)]
        print (str(i)+"->"+str(j)+":"+w)
import copy
def process_file(filename):
    s = []
    with open(filename) as f:
        for i, line in enumerate(f):
            line = line.strip()
            s.append(line)
    return s


# s = process_file('2BreakDistance.txt')
# s = process_file('result.txt')
# s = process_file('rosalind_ba6c.txt')

def chromeToCycle(chrome):
    nodes = []
    temp = ''
    for s in chrome:
        if s == '(':
            temp = ''
        elif s == ')':
            nums = temp.split()
            for i in range(len(nums)):
                n = int(nums[i])
                if n > 0:
                    head = n*2-1
                    tail = n*2
                else:
                    head = -n*2
                    tail = -n*2-1
                nodes.extend([head, tail])
        else:
            temp += s
    return nodes

def cycletoChrome(nodes):
    chrome = []
    for i in range(len(nodes)//2):
        if nodes[2*i] < nodes[2*i+1]:
            chrome.append(nodes[2*i+1]//2)
        else:
            chrome.append(-nodes[2*i]//2)
    return chrome

def colorEdges(chrome_s):
    edges = {}
    temp = ''
    for s in chrome_s:
        if s == '(':
            temp = '('
        elif s == ')':
            temp += ')'
            nodes = chromeToCycle(temp)
            for i in range(len(nodes)//2):
                n1, n2 = nodes[2*i+1], nodes[(2*i+2)%len(nodes)]
                edges[n1] = n2
                edges[n2] = n1
        else:
            temp += s
    return edges

def blackEdges(chrome_s):
    edges = {}
    temp = ''
    for s in chrome_s:
        if s == '(':
            temp = '('
        elif s == ')':
            temp += ')'
            nodes = chromeToCycle(temp)
            for i in range(len(nodes)//2):
                n1, n2 = nodes[2*i], nodes[2*i+1]
                edges[n1] = n2
                edges[n2] = n1
        else:
            temp += s
    return edges

def graphEdges(chrome_s):
    black_edges = blackEdges(chrome_s)
    color_edges = colorEdges(chrome_s)
    return (black_edges, color_edges)

def findCycle(graph):
    black_edges, color_edges = graph
    black_turn, color_turn = 0, 1

    cycles = []
    visited = set()
    node_list = sorted(list(black_edges))

    while len(node_list) > 0:
        node = node_list.pop(0)
        if node in visited:
            continue
        path = []
        turn = black_turn
        while node not in visited:
            visited.add(node)
            path.append(node)
            if turn == black_turn:
                node = black_edges[node]
                turn = color_turn
            else:
                node = color_edges[node]
                turn = black_turn                   
        cycles.append(path)
    return cycles

def graphToGenome(graph):
    path = ''
    for cycle in findCycle(graph):
        nodes = cycle
        chrome = cycletoChrome(nodes)
        path += '('
        temp = []
        for n in chrome:
            if n > 0:
                temp.append('+'+str(n))
            else:
                temp.append(str(n))
        path += ' '.join(temp)
        path += ')'
    return path

def twoBreakOnGraph(graph, i1, i2, i3, i4):
    _, color_edges = graph
    # remove (i1, i2), (i3, i4)
    # add (i1, i3) (i2, i4)
    del color_edges[i1]
    del color_edges[i2]
    del color_edges[i3]
    del color_edges[i4]

    color_edges[i1] = i3
    color_edges[i3] = i1
    color_edges[i2] = i4
    color_edges[i4] = i2
    return

def twoBreakOnGenome(P, a, b, c, d):
    graph = graphEdges(P)
    twoBreakOnGraph(graph, a, b, c, d)
    return graphToGenome(graph)

# print (twoBreakOnGenome('(+1 -2 -3 +4)', 2, 4, 5, 7))

def findNonTrivialCycle(breakpointGraph):
    red_edges, blue_edges = breakpointGraph
   
    red_turn, blue_turn = 0, 1
    visited = set()

    node_list = sorted(list(blue_edges))

    while len(node_list) > 0:
        root = node_list.pop(0)
        node = root
        if node in visited:
            continue
        path = []
        turn =blue_turn
        while node not in visited:
            visited.add(node)
            path.append(node)
            if turn == blue_turn:
                node = blue_edges[node]
                turn = red_turn
            else:
                node = red_edges[node]
                turn = blue_turn                   
        if len(path) > 2:
            next_root = blue_edges[root]
            return [root, red_edges[root], next_root, red_edges[next_root]]
    return None   


def shortestRearrangementScenario(P, Q):
    print (P)
    red_edges = colorEdges(P)
    blue_edges = colorEdges(Q)
    breakpointGraph = (red_edges, blue_edges)
    while findNonTrivialCycle(breakpointGraph) is not None:
        i1, i2, i3, i4 = findNonTrivialCycle(breakpointGraph)
        P = twoBreakOnGenome(P, i1, i2, i3, i4)
        red_edges = colorEdges(P)
        breakpointGraph = (red_edges, blue_edges)
        print (P)
    return 

# print (shortestRearrangementScenario('(+1 -2 -3 +4)', '(+1 +2 -4 -3)'))
s = process_file("rosalind_ba6d.txt")
shortestRearrangementScenario(s[0], s[1])
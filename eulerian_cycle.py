import random
def process_file(filename):
    edge = {}
    with open(filename) as f:
        for _, line in enumerate(f):
            line = line.strip()
            pref, suf = line.split(' -> ')
            edge[pref] = suf.split(',')
    return edge

def oneCycle(node, graph):
    cycle = [node]
    while node in graph:
        node_list = graph[node]
        next_node = node_list.pop()
        if len(node_list) == 0:
            del graph[node]
        node = next_node
        cycle.append(node)
    return cycle

def nextNode(cycle, graph):
    for node in cycle:
        if node in graph:
            return node
    return None

def shift(path, node):
    for i in range(len(path)):
        if path[i] == node:
            return path[i:-1] + path[:i]

def eulerianCycle(graph):
    path = []
    node = random.choice(list(graph))
    while len(graph) > 0:
        cycle = oneCycle(node, graph)
        path.extend(cycle)
        node = nextNode(path, graph)
        if node is not None:
            path = shift(path, node)
    return ('->'.join(path))

# edge= process_file('results.txt')
edge = process_file('rosalind_ba3f.txt')
print (eulerianCycle(edge))
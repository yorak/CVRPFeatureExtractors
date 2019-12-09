from heapq import heappush, heappop

def dijkstra(graph, start):
    dist = {}
    queue = [(0, start)]
    while queue:
        l, v = heappop(queue)
        if v not in dist: #not visited
            dist[v] = l
            if v in graph: # nodes leaving from
                for w, c in graph[v].items():
                    if w not in dist:
                        heappush(queue, (l + c, w))
    return dist

def dijkstra_path(graph, start):
    prev = {}
    queue = [(0, None, start)]
    while queue:
        l, from_v, to_v = heappop(queue)
        if to_v not in prev: #not visited
            prev[to_v] = from_v
            if to_v in graph: # nodes leaving from
                for next_v, c in graph[to_v].items():
                    if next_v not in prev:
                        heappush(queue, (l + c, to_v, next_v))
    return prev
    
def test():
    
    G = {
        'a':{'b':7,'c':9,'f':14},
        'b':{'c':10,'d':15},
        'c':{'d':11,'f':2},
        'd':{'e':6},
        'e':{'f':9}
    }
    print dijkstra(G, 'a')
    
    
if __name__=="__main__":
    test()
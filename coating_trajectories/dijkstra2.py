import heapq

def dijkstra(adj, costs, s, t):
    ''' Return predecessors and min distance if there exists a shortest path 
        from s to t; Otherwise, return None '''
    Q = []     # priority queue of items; note item is mutable.
    d = {s: 0} # vertex -> minimal distance
    Qd = {}    # vertex -> [d[v], parent_v, v]
    p = {}     # predecessor
    visited_set = set([s])

    for v in adj.get(s, []):
        d[v] = costs[s, v]
        item = [d[v], s, v]
        heapq.heappush(Q, item)
        Qd[v] = item

    while Q:
        cost, parent, u = heapq.heappop(Q)
        if u not in visited_set:
            p[u]= parent
            visited_set.add(u)
            if u == t:
                return p, d[u]
            for v in adj.get(u, []):
                if d.get(v):
                    if d[v] > costs[u, v] + d[u]:
                        d[v] =  costs[u, v] + d[u]
                        Qd[v][0] = d[v]    # decrease key
                        Qd[v][1] = u       # update predecessor
                        heapq._siftdown(Q, 0, Q.index(Qd[v]))
                else:
                    d[v] = costs[u, v] + d[u]
                    item = [d[v], u, v]
                    heapq.heappush(Q, item)
                    Qd[v] = item

    return None

def make_undirected(cost):
    ucost = {}
    for k, w in cost.iteritems():
        ucost[k] = w
        ucost[(k[1],k[0])] = w
    return ucost

class dijkstra_cost:
    joint_distance = None
    vs = None
    vt = None
    stored = dict()
    limits = None
    joints = None
    def __init__(self,distance_cost,joints,vs,vt,limits = None):
        self.joint_distance = distance_cost
        self.vs = vs
        self.vt = vt
        self.limits = limits
        self.joints = joints
        return

    def __getitem__(self, item):
        from_node = item[0]
        to_node = item[1]

        if (self.vs in item) or (self.vt in item):
            return 0

        if not self.stored.has_key(item):
            mdistance = self.joint_distance(self.joints[from_node[0]][from_node[1]], self.joints[to_node[0]], self.limits)
            for target in range(0, len(mdistance)):
                self.stored[(item[0],(to_node[0],target))] = mdistance[target]


        return self.stored.get(item)
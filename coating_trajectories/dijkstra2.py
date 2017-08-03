import heapq
from numpy import finfo, array, inf, zeros,round

def dijkstra(adj, costs, s, t):
    ''' Return predecessors and min distance if there exists a shortest path 
        from s to t; Otherwise, return None '''
    Q = []     # priority queue of items; note item is mutable.
    d = {s: 0} # vertex -> minimal distance
    Qd = {}    # vertex -> [d[v], parent_v, v]
    p = {}     # predecessor
    visited_set = {s}

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
                    if d[v] > costs[u, v] + d[u]:#max(costs[u, v],d[u])
                        d[v] =   costs[u, v] + d[u] #max(costs[u, v],d[u])
                        Qd[v][0] = d[v]    # decrease key
                        Qd[v][1] = u       # update predecessor
                        heapq._siftdown(Q, 0, Q.index(Qd[v]))
                else:
                    d[v] = costs[u, v] + d[u] #max(costs[u, v],d[u])
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

class virtual_node:
    def __init__(self,point,velocity):
        self.content = (point, velocity)

    def __getitem__(self, item):
        return self.content[item]

    def __eq__(self, other):
        return (self.content[0]==other[0])

    def __hash__(self):
        return hash(self.content[0])


class dijkstra_adj:
    def __init__(self,joints,dtimes):
        self.adj = dict()
        self.joints = joints
        self.dtimes = dtimes
        return

    def add_link(self,from_point,to_point):
        adjacents = self.adj.get(from_point, [])
        adjacents.append(to_point)
        self.adj[from_point] = adjacents

    def get_linked(self,from_point):
        return self.adj.get(from_point,[])

    def __getitem__(self, item):
        return self.get(item)

    def get(self, item, std_ans = None):

        from_node = item[0]
        full_states = []

        for to_node in self.get_linked(from_node):
            velocity =  (self.joints[to_node[0]][to_node[1]]-self.joints[from_node[0]][from_node[1]])/self.dtimes[to_node[0]]
            full_states += [(to_node, tuple(round(velocity, 9)))]

        return full_states

class dijkstra_acc_cost:
    def __init__(self, acc_cost, vs, vt, dtimes,vel_limits,acc_limits):
        self.acc_cost = acc_cost
        self.vs = vs
        self.vt = vt
        self.vel_limits = vel_limits
        self.acc_limits = acc_limits
        self.dtimes = dtimes
        self.stored = dict()
        return

    def __getitem__(self, item):
        from_node, from_velocity = item[0]
        from_velocity = array(from_velocity)
        to_node, to_velocity = item[1]
        to_velocity = array(to_velocity)

        if (self.vs in item) or (self.vt in item):
            return finfo(float).eps

        if (from_node[0] == 0) or (to_node[0] == len(self.dtimes)-1):
            if not self.stored.has_key(item):
                self.stored[item] = self.acc_cost(from_velocity, to_velocity, (self.dtimes[to_node[0]]+self.dtimes[from_node[0]])/2., self.vel_limits, self.acc_limits*inf)

        if not self.stored.has_key(item):
            self.stored[item] = self.acc_cost(from_velocity, to_velocity, (self.dtimes[to_node[0]]+self.dtimes[from_node[0]])/2., self.vel_limits, self.acc_limits)

        return self.stored[item]

class dijkstra_vel_cost:
    def __init__(self,distance_cost,joints,vs,vt,dtimes,vel_limits):
        self.distance_cost = distance_cost
        self.dtimes = dtimes
        self.vs = vs
        self.vt = vt
        self.vel_limits = vel_limits
        self.joints = joints
        self.stored = dict()
        return

    def __getitem__(self, item):
        from_node = item[0]
        to_node = item[1]

        if (self.vs in item) or (self.vt in item):
            return finfo(float).eps

        if not self.stored.has_key(item):
            mdistance = self.distance_cost(self.joints[from_node[0]][from_node[1]], self.joints[to_node[0]],self.dtimes[to_node[0]], self.vel_limits)
            for target in range(len(mdistance)):
                self.stored[(item[0],(to_node[0],target))] = mdistance[target]

        return self.stored[item]
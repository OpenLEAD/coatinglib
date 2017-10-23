import heapq
from numpy import finfo, array, inf, zeros,round

def dijkstra(adj, costs, s, t):
    ''' Return predecessors and min distance if there exists a shortest path 
        from s to t; Otherwise, raise exception
    '''
    Q = []     # priority queue of items; note item is mutable.
    d = {s: 0} # vertex -> minimal distance
    p = {}     # predecessor
    visited_set = set()
    successors = {}

    heapq.heappush(Q, (0, None, s))
    while Q:
        u_cost, parent, u = heapq.heappop(Q)
        if u not in visited_set:
            p[u]= parent
            visited_set.add(u)
            if u == t:
                return p, d[u]

            u_geometric  = u[0]
            u_successors = successors.get(u_geometric)
            if not u_successors:
                successors[u_geometric] = u_successors = adj.get(u_geometric)

            for v in u_successors:
                v_cost_new = u_cost + costs[u, v]
                v_cost_old = d.get(v)
                if (not v_cost_old) or (v_cost_old > v_cost_new):
                    d[v] = v_cost_new
                    item = [v_cost_new, u, v]
                    heapq.heappush(Q, item)

    raise ValueError('No shortest path to target.')

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

    def get(self, from_node):
        full_states = []

        from_joints = self.joints[from_node[0]][from_node[1]]
        for to_node in self.get_linked(from_node):
            dtime = self.dtimes[to_node[0]]
            velocity =  (self.joints[to_node[0]][to_node[1]]-from_joints)/dtime
            velocity = round(velocity, 9)
            state = (to_node, tuple(velocity))
            full_states.append(state)

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
        to_node, to_velocity = item[1]

        if (self.vs in item) or (self.vt in item):
            return finfo(float).eps

        if not self.stored.has_key(item):
            if (from_node[0] == 0) or (to_node[0] == len(self.dtimes)-1):
                acc_limits = self.acc_limits*inf
            else:
                acc_limits = self.acc_limits
            from_velocity = array(from_velocity)
            to_velocity = array(to_velocity)
            self.stored[item] = self.acc_cost(from_velocity, to_velocity,
                    (self.dtimes[to_node[0]]+self.dtimes[from_node[0]])/2.,
                    self.vel_limits, acc_limits)

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

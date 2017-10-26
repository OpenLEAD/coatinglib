import heapq
from numpy import finfo, array, inf, zeros,round, amax
from numpy.core.multiarray import zeros, array
from numpy.core.umath import absolute


def dijkstra_vel(adj, costs, s, t):
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

            u_geometric  = u
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


class VirtualNode:
    def __init__(self,point,velocity):
        self.content = (point, velocity)

    def __getitem__(self, item):
        return self.content[item]

    def __eq__(self, other):
        return (self.content[0]==other[0])

    def __hash__(self):
        return hash(self.content[0])


class DijkstraAdj:
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


def make_dijkstra_vel(joints, dtimes, vel_limits, verbose = False):
    virtual_start = (-1,-1)
    virtual_end = (-2,-2)
    adj = dict()

    for jointsi in range(0,len(joints)-1):
         for u in range(0,len(joints[jointsi])):
             for v in range(0,len(joints[jointsi+1])):
                 l = adj.get((jointsi,u),[])
                 l.append((jointsi+1,v))
                 adj[(jointsi,u)] = l

    for joints0i in range(0,len(joints[0])):
        l = adj.get(virtual_start,[])
        l.append((0,joints0i))
        adj[virtual_start] = l

    for jointsi in range(0,len(joints[-1])):
        l = adj.get((len(joints)-1,jointsi),[])
        l.append(virtual_end)
        adj[(len(joints)-1,jointsi)] = l

    cost = dijkstra_vel.dijkstra_vel_cost(joint_distance_mh12, joints, virtual_start, virtual_end, dtimes, vel_limits)

    predecessors, min_cost = dijkstra_vel.dijkstra_vel(adj, cost, virtual_start, virtual_end)

    c = virtual_end
    path = [c]

    while predecessors.get(c):
        path.insert(0, predecessors[c])
        c = predecessors[c]

    joint_path = []
    for i in range(1,len(path)-1):
        joint_index = path[i][0]
        joint_configuration = path[i][1]
        joint_path.append(joints[joint_index][joint_configuration])

    if verbose:
        return joint_path, path, min_cost, adj, cost
    else:
        return joint_path


def compute_foward_cost(joints0, joints1, limits):
    cost = zeros((len(joints0),len(joints1)))
    for i in range(0,len(joints0)):
        cost[i] = joint_distance_mh12(joints0[i], joints1, limits)
    return cost


def joint_distance_mh12(joint1, joint2, dt, vel_limits):
    dif = abs(array(joint1) - array(joint2))
    if dt == 0:
        dif = max(dif,1)
        blown = (dif > 1e-5)
        dif[blown] = inf
        dif[~blown] = 0
        return dif
    dif /= dt
    percent_dif = dif/vel_limits
    return max(percent_dif,1)
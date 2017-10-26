import heapq
from numpy import finfo, array, inf, zeros,round, amax
from numpy import max as npmax
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


class DijkstraAdj:
    def __init__(self,joints):
        self.joints = joints
        return

    def get(self,from_point):
        return [(from_point[0]+1,joint) for joint in range(len(self.joints[from_point[0]+1]))]

    def __getitem__(self, item):
        return self.get(item)


class DijkstraVelCost:
    def __init__(self,distance_cost,joints,dtimes,vel_limits):
        self.distance_cost = distance_cost
        self.dtimes = dtimes
        self.vel_limits = vel_limits
        self.joints = joints
        self.stored = dict()
        return

    def __getitem__(self, item):
        from_node = item[0]
        to_node = item[1]

        if  ((0,0) in item) or ((len(self.joints)-1,0) in item):
            return finfo(float).eps

        if not self.stored.has_key(item):
            mdistance = self.distance_cost(self.joints[from_node[0]][from_node[1]], self.joints[to_node[0]],self.dtimes[to_node[0]], self.vel_limits)
            for target in range(len(mdistance)):
                self.stored[(item[0],(to_node[0],target))] = mdistance[target]

        return self.stored[item]


def make_dijkstra_vel(joints, dtimes, vel_limits, verbose = False):

    exd_joints = [[zeros(len(joints[0][0]))]] + list(joints) + [[zeros(len(joints[0][0]))]]
    dtimes = [0] + list(dtimes) + [0]

    target = (len(exd_joints)-1,0)

    adj = DijkstraAdj(exd_joints)

    cost = DijkstraVelCost(joint_distance_mh12, exd_joints, dtimes, vel_limits)

    predecessors, min_cost = dijkstra_vel(adj, cost, (0,0), target)

    path = [target]

    while predecessors.get(target):
        path.insert(0, predecessors[target])
        target = predecessors[target]

    joint_path = []
    for i in range(1,len(path)-1):
        joint_index = path[i][0]
        joint_configuration = path[i][1]
        joint_path.append(exd_joints[joint_index][joint_configuration])

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
        dif = npmax(dif,1)
        blown = (dif > 1e-5)
        dif[blown] = inf
        dif[~blown] = 0
        return dif
    dif /= dt
    percent_dif = dif/vel_limits
    return npmax(percent_dif,1)
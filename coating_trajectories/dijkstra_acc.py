import heapq
from numpy import finfo, array, inf, zeros,round
from numpy.core.multiarray import zeros, array
from numpy.core.umath import absolute


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


class VirtualNode:
    def __init__(self,point,velocity):
        self.content = (point, velocity)
        self.hash = hash(self.content[0])

    def __getitem__(self, item):
        return self.content[item]

    def __eq__(self, other):
        return self.content[0]==other[0]

    def __hash__(self):
        return self.hash


class DijkstraAdj:
    def __init__(self,joints,dtimes):
        self.joints = joints
        self.dtimes = dtimes
        return

    def get_linked(self,from_point):
        return [(from_point[0]+1,joint) for joint in range(len(self.joints[from_point[0]+1]))]

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


class DijkstraAccCost:
    def __init__(self, acc_cost, dtimes,vel_limits,acc_limits):
        self.acc_cost = acc_cost
        self.vel_limits = vel_limits
        self.acc_limits = acc_limits
        self.dtimes = dtimes
        self.stored = dict()
        return

    def __getitem__(self, item):
        from_node, from_velocity = item[0]
        to_node, to_velocity = item[1]

        if (from_node[0]==0) or (to_node[0]==(len(self.dtimes)-1)):
            return finfo(float).eps

        if not self.stored.has_key(item):
            if (from_node[0] == 1) or (to_node[0] == (len(self.dtimes)-2)):
                acc_limits = self.acc_limits*inf
            else:
                acc_limits = self.acc_limits
            from_velocity = array(from_velocity)
            to_velocity = array(to_velocity)
            self.stored[item] = self.acc_cost(from_velocity, to_velocity,
                    (self.dtimes[to_node[0]]+self.dtimes[from_node[0]])/2.,
                    self.vel_limits, acc_limits)

        return self.stored[item]


def make_dijkstra(joints, dtimes, vel_limits, acc_limits, verbose = False):

    zero_joints = zeros(len(joints[0][0]))

    exd_joints = [[zero_joints]] + list(joints) + [[zero_joints]]
    dtimes = [1] + list(dtimes) + [1]
    dtimes[1] = 1

    start = VirtualNode((0, 0), tuple(zero_joints))
    target = VirtualNode((len(exd_joints)-1,0), tuple(zero_joints))

    cost = DijkstraAccCost(single_vel_distance, dtimes, vel_limits, acc_limits)

    adj = DijkstraAdj(exd_joints, dtimes)
    predecessors, min_cost = dijkstra(adj, cost, start, target)

    c = next(y for y in predecessors.keys() if y == target)

    path = [c]
    while predecessors.get(c):
        path.insert(0, predecessors[c])
        c = predecessors[c]

    joint_path = []
    for i in range(1,len(path)-1):
        joint_index = path[i][0][0]
        joint_configuration = path[i][0][1]
        joint_path.append(exd_joints[joint_index][joint_configuration])

    if verbose:
        return joint_path, path, min_cost, adj, cost
    else:
        return joint_path


def single_vel_distance( vel1, vel2, dt, vel_limits, acc_limits):
    dif = abs(array(vel1)-array(vel2))/dt
    return sum((dif/acc_limits)**2) + sum((vel2/vel_limits)**2)
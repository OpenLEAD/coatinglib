from numpy import array, zeros, max
from numpy import abs, inf
from openravepy import IkFilterOptions, interfaces, databases, IkParameterization
from openravepy import CollisionOptions, RaveCreateCollisionChecker, CollisionReport
import dijkstra2





def single_vel_distance( vel1, vel2, dt, vel_limits, acc_limits):
    dif = abs(array(vel1)-array(vel2))/dt
    return sum((dif/acc_limits)**2) + sum((abs(vel2)/vel_limits)**2)
       
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

def compute_foward_cost(joints0, joints1, limits):
    cost = zeros((len(joints0),len(joints1)))
    for i in range(0,len(joints0)):
        cost[i] = joint_distance_mh12(joints0[i], joints1, limits)
    return cost


def make_dijkstra(joints, dtimes, vel_limits, acc_limits, verbose = False):
    virtual_start = (-1,-1)
    virtual_end = (-2,-1)
    adj = dijkstra2.dijkstra_adj(joints,dtimes)

    for jointsi in range(len(joints)-1):
        for u in range(len(joints[jointsi])):
            for v in range(len(joints[jointsi+1])):
                adj.add_link((jointsi,u),(jointsi+1,v))

    for joints0i in range(len(joints[0])):
        adj.add_link(virtual_start,(0,joints0i))

    for jointsi in range(len(joints[-1])):
        adj.add_link((len(joints)-1,jointsi),virtual_end)

    vs = dijkstra2.virtual_node((-1,-1),tuple(zeros(len(joints[0][0]))))
    ve = dijkstra2.virtual_node((-2,-1),tuple(zeros(len(joints[0][0]))))
    dtimes[0] = 1

    cost = dijkstra2.dijkstra_acc_cost(single_vel_distance,vs,ve,dtimes,vel_limits,acc_limits)

    predecessors, min_cost = dijkstra2.dijkstra(adj, cost, vs, ve)

    c = next(y for y in predecessors.keys() if y == ve)

    path = [c]
    while predecessors.get(c):
        path.insert(0, predecessors[c])
        c = predecessors[c]

    joint_path = []
    for i in range(1,len(path)-1):
        joint_index = path[i][0][0]
        joint_configuration = path[i][0][1]
        joint_path.append(joints[joint_index][joint_configuration])

    if verbose:
        return joint_path, path, min_cost, adj, cost
    else:
        return joint_path


def make_dijkstra_vel(joints, dtimes, vel_limits, acc_limits, verbose = False):
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

    cost = dijkstra2.dijkstra_vel_cost(joint_distance_mh12,joints,virtual_start,virtual_end,dtimes,vel_limits)

    predecessors, min_cost = dijkstra2.dijkstra(adj, cost, virtual_start, virtual_end)

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
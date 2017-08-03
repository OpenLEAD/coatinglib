import graph_tool as gt
import planning
import dijkstra2


class DjkVisitor(gt.search.DijkstraVisitor):
    def __init__(self,weight,graph,adj,cost,pv):
        self.weight = weight
        self.graph = graph
        self.adj = adj
        self.cost = cost
        self.pv = pv

    def discover_vertex(self,u):
        adjs = self.adj.get(self.pv[u])
        for i,vertex in enumerate(self.graph.add_vertex(len(adjs))):
            self.pv[vertex] = adjs[i]
            self.graph.add_edge(u,vertex)

    def examine_edge(self,e):
        self.weight[e] = self.cost[self.pv[e.source()],self.pv[e.target()]]


def make_dijkstra(joints, dtimes, vel_limits, acc_limits, verbose = False):
    djvisitor = DjkVisitor(joints,dtimes)
    virtual_start = (-1,-1) #falta velocidade inicial
    virtual_end = (-2,-1) #falta velocidade final
    adj = dijkstra2.dijkstra_adj(joints,dtimes)

    for jointsi in range(len(joints)-1):
         for u in range(len(joints[jointsi])):
             for v in range(len(joints[jointsi+1])):
                 adj.add_link((jointsi,u),(jointsi+1,v))

    for joints0i in range(len(joints[0])):
        adj.add_link(virtual_start,(0,joints0i))

    for jointsi in range(len(joints[-1])):
        adj.add_link((len(joints)-1,jointsi),virtual_end)

    vs = dijkstra2.virtual_node((-1,-1),tuple(zeros(len(joints[0][0])))) #falta velocidade inicial
    ve = dijkstra2.virtual_node((-2,-1),tuple(zeros(len(joints[0][0])))) #falta velocidade final
    dtimes[0] = 1

    cost = dijkstra2.dijkstra_acc_cost(single_vel_distance,vs,ve,dtimes,vel_limits,acc_limits)

    min_costs, predecessors = gt.search.dijkstra_search(g, weight, source=None, visitor= djvisitor)

    predecessors, min_cost = dijkstra2.dijkstra(adj, cost, vs, ve)

    c = [y for x, y in enumerate(predecessors.keys()) if y == ve][0]
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
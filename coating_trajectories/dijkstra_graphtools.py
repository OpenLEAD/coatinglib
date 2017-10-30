import graph_tool as gt
from graph_tool.search import StopSearch
from numpy import zeros
import dijkstra_acc


class AstarVisitor(gt.search.DijkstraVisitor):
    def __init__(self,weight,graph,adj,cost,dist,graph_cost,point_velocity,virtual_end):
        self.graph_cost = graph_cost
        self.dist = dist
        self.weight = weight
        self.graph = graph
        self.adj = adj
        self.cost = cost
        self.point_velocity = point_velocity
        self.virtual_end = virtual_end
        self.target = None

    def examine_vertex(self,u):
        adjs = self.adj.get(self.point_velocity[u])

        if len(adjs) == 1:
            vertex = self.graph.add_vertex()
            self.point_velocity[vertex] = adjs[0]
            self.graph.add_edge(u,vertex)
            self.dist[vertex] = self.graph_cost[vertex] = float('inf')
        else:
            for i,vertex in enumerate(self.graph.add_vertex(len(adjs))):
                self.point_velocity[vertex] = adjs[i]
                self.graph.add_edge(u,vertex)
                self.dist[vertex] = self.graph_cost[vertex] = float('inf')


    def examine_edge(self,e):
        self.weight[e] = self.cost[self.point_velocity[e.source()],
                                   self.point_velocity[e.target()]]

    def edge_relaxed(self, e):
        if self.point_velocity[e.target()] == self.virtual_end:
            self.target = e.target()
            raise StopSearch()


def make_dijkstra(joints, dtimes, vel_limits, acc_limits, verbose = False):
    virtual_start = (-1,-1) #falta velocidade inicial
    virtual_end = (-2,-1) #falta velocidade final
    adj = dijkstra_acc.DijkstraAdj(joints, dtimes)

    for jointsi in range(len(joints)-1):
         for u in range(len(joints[jointsi])):
             for v in range(len(joints[jointsi+1])):
                 adj.add_link((jointsi,u),(jointsi+1,v))

    for joints0i in range(len(joints[0])):
        adj.add_link(virtual_start,(0,joints0i))

    for jointsi in range(len(joints[-1])):
        adj.add_link((len(joints)-1,jointsi),virtual_end)

    vs = dijkstra_acc.VirtualNode((-1, -1), tuple(zeros(len(joints[0][0])))) #falta velocidade inicial
    ve = dijkstra_acc.VirtualNode((-2, -1), tuple(zeros(len(joints[0][0])))) #falta velocidade final
    dtimes[0] = 1

    cost = dijkstra_acc.DijkstraAccCost(dijkstra_acc.single_vel_distance, vs, ve, dtimes, vel_limits, acc_limits)

    graph = gt.Graph()
    weight = graph.new_ep("double")
    dist = graph.new_vp("double")
    graph_cost = graph.new_vp("double")
    point_velocity = graph.new_vp("object")
    graph.add_vertex(1)
    point_velocity[graph.vertex(0)] = vs

    astarvisitor = AstarVisitor(weight,graph,adj,cost,dist,graph_cost,point_velocity,ve)

    min_costs, predecessors = gt.search.astar_search(g=graph, weight=weight, dist_map=dist,
                                                     source=graph.vertex(0), visitor=astarvisitor,
                                                     implicit=True, cost_map = graph_cost)

    path = []
    v = graph.vertex(predecessors[astarvisitor.target])
    while v != graph.vertex(0):
        path = [point_velocity[v][0]] + path
        v = graph.vertex(predecessors[v])

    joint_path = []
    for i in range(len(path)):
        joint_index = path[i][0]
        joint_configuration = path[i][1]
        joint_path.append(joints[joint_index][joint_configuration])

    if verbose:
        return joint_path, path, min_costs[astarvisitor.target], adj, cost
    else:
        return joint_path
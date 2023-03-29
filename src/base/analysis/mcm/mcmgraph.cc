/*
 *  Eindhoven University of Technology
 *  Eindhoven, The Netherlands
 *  Dept. of Electrical Engineering
 *  Electronics Systems Group
 *  Model Based Design Lab (https://computationalmodeling.info/)
 *
 *  Name            :   mcmgraph.cc
 *
 *  Author          :   Sander Stuijk (sander@ics.ele.tue.nl)
 *
 *  Date            :   November 7, 2005
 *
 *  Function        :   Convert HSDF graph to weighted directed graph for MCM
 *                      calculation.
 *
 *  History         :
 *      07-11-05    :   Initial version.
 *
 *
 *  Copyright 2023 Eindhoven University of Technology
 *
 *  Permission is hereby granted, free of charge, to any person obtaining
 *  a copy of this software and associated documentation files (the “Software”),
 *  to deal in the Software without restriction, including without limitation
 *  the rights to use, copy, modify, merge, publish, distribute, sublicense,
 *  and/or sell copies of the Software, and to permit persons to whom the
 *  Software is furnished to do so, subject to the following conditions:
 *
 *  The above copyright notice and this permission notice shall be included
 *  in all copies or substantial portions of the Software.
 *
 *  THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 *  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 *  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 *  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 *  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 *  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 *  SOFTWARE.
 */

#include "base/analysis/mcm/mcmgraph.h"
#include "base/analysis/mcm/mcm.h"
#include "base/analysis/mcm/mcmyto.h"
#include "base/base.h"

namespace Graphs {
/**
 * ~MCMgraph
 * Constructor.
 */
MCMgraph::MCMgraph() {}

/**
 * ~MCMgraph
 * Destructor.
 */
MCMgraph::~MCMgraph() {
}

// Copy Constructor
MCMgraph::MCMgraph(const MCMgraph &g) {
    for (auto i = g.nodes.cbegin(); i != g.nodes.cend(); i++) {
        const MCMnode &n = **i;
        this->addNode(std::make_shared<MCMnode>(n.id, n.visible));
    }
    for (auto i = g.edges.cbegin(); i != g.edges.cend(); i++) {
        const MCMedge &e = **i;
        this->addEdge(e.id, this->getNode(e.src->id), this->getNode(e.dst->id), e.w, e.d);
    }
}

/**
 * MCMedge
 * Constructor.
 */
MCMedge::MCMedge(CId eId, bool eVisible) :
    id(eId), visible(eVisible), src(nullptr), dst(nullptr), w(0.0), d(0.0) {}

/**
 * MCMnode
 * Constructor.
 */
MCMnode::MCMnode(CId nId, bool nVisible) : id(nId), visible(nVisible) {}

/**
 * splitMCMedgeToSequence ()
 * The function converts an MCM edge with more then one delay
 * into a sequence of edges with one delay (uses recursive call
 * to itself).
 */
static void splitMCMedgeToSequence(std::shared_ptr<MCMgraph> g, std::shared_ptr<MCMedge> e) {
    // Create dummy node n;
    std::shared_ptr<MCMnode> n = std::make_shared<MCMnode>((CId)g->getNodes().size(), true);
    g->addNode(n);

    // Create a new edge between the src node of e and a new
    // dummy node.
    std::shared_ptr<MCMedge> eN = std::make_shared<MCMedge>((CId)g->getEdges().size(), true);
    eN->src = e->src;
    eN->dst = n;
    eN->w = 0;
    eN->d = e->d - 1;
    g->addEdge(eN);

    // Remove e from the set of edges its source node is connected to
    // and add eN to this list
    for (MCMedgesIter iter = e->src->out.begin(); iter != e->src->out.end(); iter++) {
        if ((*iter)->id == e->id) {
            e->src->out.erase(iter);
            break;
        }
    }
    e->src->out.push_back(eN);

    // Connect e to node n
    e->src = n;
    n->out.push_back(e);
    n->in.push_back(eN);

    // One delay left on e
    e->d = 1;

    // More then one delay on the new edge?
    if (eN->d > 1)
        splitMCMedgeToSequence(g, eN);
}

/**
 * addLongestDelayEdgeForNode ()
 * The function add an additional edge between the node n and any node
 * reachable over e which only uses edges with no initial delays (except
 * for the initial edge e). The weight of the new edge is the longest
 * path from the node n to the destination node.
 * The function uses an adapted version of Dijkstra's shortest path
 * algorithm.
 * Note: algorithm assumes that edge weights are integer values !
 */
static void addLongestDelayEdgeForNode(std::shared_ptr<MCMgraph> g, std::shared_ptr<MCMnode> n, std::shared_ptr<MCMedge> e, int **visit) {
    int dMax;
    int *d = new int[g->getNodes().size()];
    std::shared_ptr<MCMnode> *pi = new std::shared_ptr<MCMnode>[g->getNodes().size()];
    MCMnodes S;
    MCMnodes Q = g->getNodes();

    // Initialize single source
    for (uint v = 0; v < g->getNodes().size(); v++) {
        d[v] = -1;
        pi[v] = nullptr;
    }
    d[n->id] = 0;

    // Initialize the node connected via the edge e to n
    d[e->dst->id] = (int)e->w;

    // Remove node n from Q and add it to S (only when edge e is no
    // self-edge)
    if (e->dst->id != n->id) {
        for (MCMnodesIter iter = Q.begin(); iter != Q.end(); iter++) {
            if ((*iter)->id == n->id) {
                // Is edge e
                Q.erase(iter);
                break;
            }
        }
        S.push_back(n);
    }

    // Find all longest paths till Q is empty or all reachable paths seen
    while (!Q.empty()) {
        std::shared_ptr<MCMnode> u = nullptr;

        // Find node u in Q with largest distance
        dMax = -1;
        for (MCMnodesIter iter = Q.begin(); iter != Q.end(); iter++) {
            if (d[(*iter)->id] > dMax) {
                u = *iter;
                dMax = d[u->id];
            }
        }

        // Found no node with a non-negative distance ?
        if (dMax < 0)
            break;

        // Remove node u from Q and add it to S
        for (MCMnodesIter iter = Q.begin(); iter != Q.end(); iter++) {
            if ((*iter)->id == u->id) {
                Q.erase(iter);
                break;
            }
        }
        S.push_back(u);

#if __VISITED_NODES__
        // Visited node u before?
        if (visit[u->id] != nullptr) {
            for (MCMnodesIter iter = g->nodes.begin(); iter != g->nodes.end(); iter++) {
                MCMnode *v = *iter;

                // Node v reachable from u in the past?
                if (visit[u->id][v->id] > visit[u->id][u->id]) {
                    // Distance to v is distance to v in the past minus
                    // the distance to u in the past plus the distance
                    // to u in the current situation.
                    d[v->id] = visit[u->id][v->id] - visit[u->id][u->id] + d[u->id];
                    S.push_back(v);
                }
            }
        }
#endif

        // Relax all nodes v adjacent to u (connected via edges with no tokens)
        for (MCMedgesIter iter = u->out.begin(); iter != u->out.end(); iter++) {
            std::shared_ptr<MCMedge> e = *iter;
            std::shared_ptr<MCMnode> v = e->dst;

            if (e->d == 0 && d[v->id] < d[u->id] + e->w) {
                if (d[v->id] != -1) {
                    // Found a longer path to v, add it to list of
                    // nodes to be checked (again).
                    Q.push_back(v);
                }

                d[v->id] = d[u->id] + (int)e->w;
                pi[v->id] = u;
            }
        }
    }

#if __VISITED_NODES__
    // Store the distance for next runs of the longest path function
    // Skip first node on longest path (this node may have outgoing edges
    // not explored in this run - i.e. algo start with specific edge)
    MCMnode *v = S.back();
    while (v != nullptr && pi[v->id] != nullptr) {
        visit[v->id] = d;
        v = pi[v->id];
    }
#endif

    // Add an edge between the node n and any node m reachable from n
    // with a weight equal to the longest path from n to m
    for (MCMnodesIter iter = S.begin(); iter != S.end(); iter++) {
        std::shared_ptr<MCMnode> m = *iter;

        // Node m reachable from n and not connected directly to n via e?
        if (d[m->id] > 0 && e->dst->id != m->id) {
            // Create an edge between n and m
            std::shared_ptr<MCMedge> eN = std::make_shared<MCMedge>((CId)g->getEdges().size(), false);
            eN->src = n;
            eN->dst = m;
            eN->w = d[m->id];
            eN->d = e->d; // This should always be 1
            n->out.push_back(eN);
            m->in.push_back(eN);
            g->addEdge(eN);
        }
    }
}

/**
 * addLongestDelayEdgesToMCMgraph ()
 * The function adds additional edges to the graph which express the
 * longest path between two nodes crossing one edge with a delay. Edges
 * with no delay are removed and edges with more then one delay element
 * are converted into a sequence of edges with one delay element.
 */
void addLongestDelayEdgesToMCMgraph(std::shared_ptr<MCMgraph> g) {
    int **visit;

    // Find longest path between a node n and a node m
    // over all sequences of edges in which only the
    // first edge may contain a delay
    for (MCMedgesCIter iter = g->getEdges().begin(); iter != g->getEdges().end(); iter++) {
        std::shared_ptr<MCMedge> e = *iter;

        // More then one delay on the edge?
        if (e->d > 1)
            splitMCMedgeToSequence(g, e);
    }

    // Initialize visit array which contains distances of already
    // computed paths
    visit = new int *[g->getNodes().size()];
    for (uint i = 0; i < g->getNodes().size(); i++)
        visit[i] = nullptr;

    // Find longest path between a node n and a node m
    // over all sequences of edges in which only the
    // first edge may contain a delay
    for (MCMedgesCIter iter = g->getEdges().begin(); iter != g->getEdges().end(); iter++) {
        std::shared_ptr<MCMedge> e = *iter;

        // Initial tokens on edge?
        if (e->d != 0 && e->visible) {
            // Find the longest path from n to any node m
            // and add an edge with the path wait to the graph
            addLongestDelayEdgeForNode(g, e->src, e, visit);

            // Seen this edge (set it to invisible)
            e->visible = false;
        }
    }

    // Hide all edges which do not contain a delay
    for (MCMedgesCIter iter = g->getEdges().begin(); iter != g->getEdges().end(); iter++) {
        std::shared_ptr<MCMedge> e = *iter;

        // No initial tokens on edge?
        if (e->d == 0)
            e->visible = false;
        else
            e->visible = true;
    }
    delete[] visit;
}

/**
 * getAdjacentNodes ()
 * The function returns a list with nodes directly reachable from
 * node a (in case transpose if false). If transpose is 'true', the
 * graph is transposed and the function returns a list with nodes
 * which are directly reachable from a in the transposed graph.
 */
static MCMnodes getAdjacentNodes(std::shared_ptr<MCMnode> a, bool transpose) {
    MCMnodes nodes;

    if (!transpose) {
        for (MCMedgesIter iter = a->out.begin(); iter != a->out.end(); iter++) {
            std::shared_ptr<MCMedge> e = *iter;

            if (e->visible)
                nodes.push_back(e->dst);
        }
    } else {
        for (MCMedgesIter iter = a->in.begin(); iter != a->in.end(); iter++) {
            std::shared_ptr<MCMedge> e = *iter;

            if (e->visible)
                nodes.push_back(e->src);
        }
    }

    return nodes;
}

/**
 * getNextNode ()
 * The function returns the node from the list of nodes with the highest
 * order. Its order is set to -1. If all nodes have order -1, a nullptr pointer is
 * returned.
 */
static std::shared_ptr<MCMnode> getNextNode(const MCMnodes &nodes, v_int &order) {
    std::shared_ptr<MCMnode> a = nullptr;
    int orderA = -1;

    // Find actor with largest order
    for (MCMnodesCIter iter = nodes.begin(); iter != nodes.end(); iter++) {
        std::shared_ptr<MCMnode> b = *iter;

        if (orderA < order[b->id]) {
            a = b;
            orderA = order[b->id];
        }
    }

    // All actors have order -1?
    if (orderA == -1)
        return nullptr;

    return a;
}

/**
 * dfsVisit ()
 * The visitor function of the DFS algorithm.
 */
static void
dfsVisit(std::shared_ptr<MCMnode> u, int &time, v_int &color, v_int &d, v_int &f, std::shared_ptr<MCMnode> *pi, bool transpose) {
    // color[u] <- gray
    color[u->id] = 1;

    // d[u] <- time <- time + 1
    time++;
    d[u->id] = time;

    // for each v in Adj(e)
    MCMnodes adj = getAdjacentNodes(u, transpose);
    for (MCMnodesIter iter = adj.begin(); iter != adj.end(); iter++) {
        std::shared_ptr<MCMnode> v = *iter;

        // do if color[v] = white
        if (color[v->id] == 0) {
            pi[v->id] = u;
            dfsVisit(v, time, color, d, f, pi, transpose);
        }
    }

    // color[u] <- black
    color[u->id] = 2;

    // f[u] <- time <- time + 1
    time++;
    f[u->id] = time;
}

/**
 * dfsMCMgraph ()
 * The function performs a depth first search on the graph. The discover
 * time (d), finish time (f) and predecessor tree (pi) are passed back.
 * The graph is transposed if the argument 't' is 'true'. Vertices are
 * considered in decreasing order of f[u].
 *
 * The function derives an 'order' vector from the vector 'f'. Actors are
 * considered in decreasing order according to this vector. When an actor
 * is visited, its order is set to -1. At the end, all actors must have
 * order -1 (i.e. all actors are visited).
 */
void dfsMCMgraph(const MCMgraph& g, v_int &d, v_int &f, std::shared_ptr<MCMnode> *pi, bool transpose) {
    int time;

    // for each u in G do order[u] <- f[u]
    v_int order(f);

    // for each u in G do color[u] <- white
    v_int color(g.nrVisibleNodes(), 0);

    // time <- 0
    time = 0;

    // for each u in G do pi[u] <- NIL
    for (uint u = 0; u < g.nrVisibleNodes(); u++)
        pi[u] = nullptr;

    // for each u in G (visit in order given by order)
    for (std::shared_ptr<MCMnode> a = getNextNode(g.getNodes(), order); a != nullptr;
         a = getNextNode(g.getNodes(), order)) {
        // Mark node as visited
        order[a->id] = -1;

        // if color[u] = white
        if (color[a->id] == 0)
            dfsVisit(a, time, color, d, f, pi, transpose);
    }
}

/**
 * addNodeToComponent ()
 * The function adds a new copy of a node to the component and it
 * also creates copies for all edges between this node and all nodes
 * already in the component.
 */
static void addNodeToComponent(std::shared_ptr<MCMnode> n, std::shared_ptr<MCMgraph> comp) {
    std::shared_ptr<MCMnode> m;

    // Create a copy of n and add it to the component
    m = std::make_shared<MCMnode>(n->id, true);
    comp->addNode(m);

    // Check all edges of n for inclusion in the component, first the outgoing edges...
    for (MCMedgesIter iter = n->out.begin(); iter != n->out.end(); iter++) {
        std::shared_ptr<MCMedge> e = *iter;

        if (!e->visible)
            continue;

        // Is destination node in the component?
        for (MCMnodesCIter iterN = comp->getNodes().begin(); iterN != comp->getNodes().end();
             iterN++) {
            if (e->dst->id == (*iterN)->id) {
                // Add a copy of the edge to the component
                std::shared_ptr<MCMedge> eN = std::make_shared<MCMedge>(e->id, e->visible);
                eN->d = e->d;
                eN->w = e->w;
                eN->src = m;
                eN->dst = *iterN;
                comp->addEdge(eN);
                break;
            }
        }
    }
    // ... and now the incoming edges. We must now skip the self edges since they have been included
    // in the previous loop already
    for (MCMedgesIter iter = n->in.begin(); iter != n->in.end(); iter++) {
        std::shared_ptr<MCMedge> e = *iter;

        // Is source node in the component?
        for (MCMnodesCIter iterN = comp->getNodes().begin(); iterN != comp->getNodes().end();
             iterN++) {
            // if the source node is in the component and it is not a self-edge
            if ((e->src->id == (*iterN)->id) && (e->src->id != e->dst->id)) {
                // Add a copy of the edge to the component
                std::shared_ptr<MCMedge> eN = std::make_shared<MCMedge>(e->id, e->visible);
                eN->d = e->d;
                eN->w = e->w;
                eN->src = *iterN;
                eN->dst = m;
                comp->addEdge(eN);
                break;
            }
        }
    }
}

/**
 * treeVisitChildren ()
 * The function visits all children of the actor 'u'. The parent-child
 * relation is given via the vector 'pi'.
 */
static bool treeVisitChildren(const MCMgraph& g, std::shared_ptr<MCMnode> *pi, std::shared_ptr<MCMnode> u, std::shared_ptr<MCMgraph> comp) {
    bool children = false;

    for (uint i = 0; i < g.getNodes().size(); i++) {
        // Node v points to node u (i.e. v is a child of u)?
        if (pi[i] != nullptr && pi[i]->id == u->id) {
            std::shared_ptr<MCMnode> v = nullptr;

            for (MCMnodesCIter iter = g.getNodes().begin(); iter != g.getNodes().end(); iter++) {
                if ((*iter)->id == i) {
                    v = *iter;
                    break;
                }
            }
            ASSERT(v != nullptr, "The must always be a node v.");

            // Add node v to the component
            v->visible = true;
            addNodeToComponent(v, comp);
            children = true;

            // Find all children of v in the tree
            treeVisitChildren(g, pi, v, comp);
        }
    }

    return children;
}

/**
 * findComponentsInMCMgraph ()
 * The function determines the strongly connected components in the graph. To do
 * this, it performs depth-first walk on the forest given by 'pi'.
 */
static void findComponentsInMCMgraph(const MCMgraph& g,
                                     std::shared_ptr<MCMnode> *pi,
                                     MCMgraphs &components,
                                     bool includeComponentsWithoutEdges = false) {
    std::shared_ptr<MCMgraph> comp;

    // Set all node as invisible
    for (MCMnodesCIter iter = g.getNodes().begin(); iter != g.getNodes().end(); iter++) {
        std::shared_ptr<MCMnode> n = *iter;

        n->visible = false;
    }

    // Find the strongly connected component and make all of its nodes visible
    for (MCMnodesCIter iter = g.getNodes().begin(); iter != g.getNodes().end(); iter++) {
        std::shared_ptr<MCMnode> n = *iter;

        if (pi[n->id] == nullptr) {
            // Create a new graph for the component
            comp = std::make_shared<MCMgraph>();

            // Find all children of n in the tree
            if (treeVisitChildren(g, pi, n, comp)) {
                // Node n has children, so it belongs to the strongly
                // connected component
                n->visible = true;
                addNodeToComponent(n, comp);
            } else {
                if (!includeComponentsWithoutEdges) {
                    // Node n may have a self-loop making it a strongly
                    // connected component
                    for (MCMedgesIter iterE = n->in.begin(); iterE != n->in.end(); iterE++) {
                        // Is edge a self-loop?
                        if ((*iterE)->src->id == n->id) {
                            n->visible = true;
                            addNodeToComponent(n, comp);
                        }
                    }
                } else {
                    n->visible = true;
                    addNodeToComponent(n, comp);
                }
            }

            // Found a strongly connected component (at least one edge)?
            if (includeComponentsWithoutEdges || comp->nrVisibleEdges() > 0)
                components.push_back(comp);
        }
    }

    // Make all edges to invisible nodes also invisible
    for (MCMedgesCIter iter = g.getEdges().begin(); iter != g.getEdges().end(); iter++) {
        std::shared_ptr<MCMedge> e = *iter;

        if (e->visible)
            if (!e->src->visible || !e->dst->visible)
                e->visible = false;
    }
}

/**
 * Extract the strongly connected components from the graph. These components
 * are returned as a set of MCM graphs. All nodes which belong to at least
 * one of the strongly connected components are set to visible in the graph g,
 * all other nodes are made invisible. Also edges between two nodes in (possibly
 * different) strongly connected components are made visible and all others
 * invisible. The graph g consists in the end of only nodes which are part of
 * a strongly connected component and all the edges between these nodes. Some
 * MCM algorithms work also on this graph (which reduces the execution time
 * needed in some of the conversion algorithms).
 */
void stronglyConnectedMCMgraph(const MCMgraph& g,
                               MCMgraphs &components,
                               bool includeComponentsWithoutEdges) {
    // Initialize
    v_int d(g.nrVisibleNodes());
    v_int f(g.nrVisibleNodes(), 0);
    std::shared_ptr<MCMnode> *pi = (std::shared_ptr<MCMnode> *)malloc(sizeof(std::shared_ptr<MCMnode>) * g.nrVisibleNodes());
    ASSERT(pi != nullptr, "Failed malloc");

    // Call dfs(g) to compute f[u] for each u
    dfsMCMgraph(g, d, f, pi, false);

    // Call dfs(gT), consider the vertices in order of decreasing f[u]
    dfsMCMgraph(g, d, f, pi, true);

    // Component is given by the vertices of the tree in the depth-first forest
    findComponentsInMCMgraph(g, pi, components, includeComponentsWithoutEdges);

    //// Print the MCM graph in dot-format
    // cerr << "digraph g {" << endl;
    // cerr << "    size=\"7,10\";" << endl;
    // for (MCMedgesIter iter = g->edges.begin(); iter != g->edges.end(); iter++)
    //{
    //     MCMedge *e = *iter;

    //    if (e->visible)
    //    {
    //        cerr << "    " << e->src->id << " -> " << e->dst->id;
    //        cerr << " [ label=" << e->w << " ];" << endl;
    //    }
    //}
    // cerr << "}" << endl;

    // Cleanup
    free(pi);
}

/**
 * relabelMCMgraph ()
 * The function removes all hidden nodes and edges from the graph. All visible
 * edges are assigned a new id starting in the range [0,nrNodes()).
 */
void relabelMCMgraph(std::shared_ptr<MCMgraph> g) {
    uint nodeId = 0, edgeId = 0;

    // Relabel nodes
    for (MCMnodesCIter iter = g->getNodes().begin(); iter != g->getNodes().end(); iter++) {
        std::shared_ptr<MCMnode> n = *iter;

        if (n->visible) {
            n->id = nodeId;
            nodeId++;
        }
    }

    // Relabel edges
    for (MCMedgesCIter iter = g->getEdges().begin(); iter != g->getEdges().end(); iter++) {
        std::shared_ptr<MCMedge> e = *iter;

        if (e->visible) {
            e->id = edgeId;
            edgeId++;
        }
    }

    // Remove edges from nodes
    for (MCMnodesCIter iter = g->getNodes().begin(); iter != g->getNodes().end(); iter++) {
        std::shared_ptr<MCMnode> n = *iter;

        if (!n->visible)
            continue;

        for (MCMedgesIter iterE = n->in.begin(); iterE != n->in.end();) {
            MCMedgesIter iterEN = iterE;

            // Next iterator
            iterE++;

            // Erase current iterator?
            if (!(*iterEN)->visible)
                n->in.erase(iterEN);
        }

        for (MCMedgesIter iterE = n->out.begin(); iterE != n->out.end();) {
            MCMedgesIter iterEN = iterE;

            // Next iterator
            iterE++;

            // Erase current iterator?
            if (!(*iterEN)->visible)
                n->out.erase(iterEN);
        }
    }

    // Remove nodes from graph
    for (MCMnodesCIter iter = g->getNodes().begin(); iter != g->getNodes().end();) {
        MCMnodesCIter iterN = iter;

        // Next iterator
        iter++;

        if (!(*iterN)->visible)
            g->removeNode(*iterN);
    }

    // Remove edges from graph
    for (MCMedgesCIter iter = g->getEdges().begin(); iter != g->getEdges().end();) {
        MCMedgesCIter iterE = iter;

        // Next iterator
        iter++;

        if (!(*iterE)->visible)
            g->removeEdge(*iterE);
    }
}

// Prune the edges in the MCMgraph to maintain only Pareto maximal combinations
// Of the weight and number of tokens pair pair of nodes
// Generates a new graph
// Note this algorithm does currently not distinguish visible and invisible edges!
std::shared_ptr<MCMgraph> MCMgraph::pruneEdges() {

    class _local {
    public:
        map<std::shared_ptr<MCMnode>, MCMedges, MCMNodeLess> paretoEdges;
        void insert(std::shared_ptr<MCMedge> e) {
            MCMedges &edges = paretoEdges[e->dst];
            // Simple Cull
            MCMedges::iterator i = edges.begin();
            bool cont = i != edges.end();
            bool add_e = true;
            while (cont) {
                std::shared_ptr<MCMedge> f = *i;
                // if e is worse than f
                if ((e->d <= f->d) && (e->w >= f->w)) {
                    // remove f from edges;
                    i = edges.erase(i);
                } else {
                    // if f is worse than e
                    if ((f->d <= e->d) && (f->w >= e->w)) {
                        // forget e and stop iterating;
                        add_e = false;
                        cont = false;
                    } else {
                        // incomparable, go on
                        i++;
                    }
                }
                cont = cont && (i != edges.end());
            }
            if (add_e)
                edges.push_back(e);
        }
    };

    _local local;

    // create new graph
    std::shared_ptr<MCMgraph> result = std::make_shared<MCMgraph>();

    // create all nodes.
    map<std::shared_ptr<MCMnode>, std::shared_ptr<MCMnode>> newNodeMap;
    std::shared_ptr<MCMnode> u;
    for (MCMnodes::iterator iter = this->nodes.begin(); iter != this->nodes.end(); iter++) {
        u = *iter;
        std::shared_ptr<MCMnode> n = std::make_shared<MCMnode>(u->id, u->visible);
        newNodeMap[u] = n;
        result->addNode(n);
    }

    // for every node
    for (MCMnodes::iterator iter = this->nodes.begin(); iter != this->nodes.end(); iter++) {
        u = *iter;

        // for every outgoing edges to a simple cull Pareto filtering
        MCMedges::const_iterator i;
        for (i = u->out.begin(); i != u->out.end(); i++) {
            local.insert(*i);
        }

        // add Pareto Edges to new Graph.
        // for every dst node
        map<std::shared_ptr<MCMnode>, MCMedges>::const_iterator j;
        for (j = local.paretoEdges.begin(); j != local.paretoEdges.end(); j++) {
            const MCMedges &edges = (*j).second;
            // for every edge
            MCMedges::const_iterator k;
            for (k = edges.begin(); k != edges.end(); k++) {
                const std::shared_ptr<MCMedge>& e = *k;
                std::shared_ptr<MCMedge> ne = std::make_shared<MCMedge>(e->id, true);
                ne->d = e->d;
                ne->w = e->w;
                ne->src = newNodeMap[u];
                ne->dst = newNodeMap[e->dst];
                result->addEdge(ne);
            }
        }
    }

    return result;
}

CDouble MCMgraph::calculateMaximumCycleMeanKarp() const { return maximumCycleMeanKarp(*this); }

CDouble MCMgraph::calculateMaximumCycleMeanKarpDouble(const MCMnode *criticalNode) const {
    return maximumCycleMeanKarpDouble(*this, criticalNode);
}

CDouble MCMgraph::calculateMaximumCycleRatioAndCriticalCycleYoungTarjanOrlin(std::shared_ptr<MCMedge> **cycle,
                                                                             uint *len) {
    return maxCycleRatioAndCriticalCycleYoungTarjanOrlin(*this, cycle, len);
}

void MCMgraph::relabelNodeIds(std::map<int, int> &nodeIdMap) {
    int k = 0;
    for (auto i = this->nodes.begin(); i != this->getNodes().end(); i++, k++) {
        MCMnode &n = **i;
        nodeIdMap[k] = n.id;
        n.id = k;
    }
}

MCMgraph MCMgraph::normalize(CDouble mu) const {
    MCMgraph result(*this);
    for (auto i = result.getEdges().begin(); i != result.getEdges().end(); i++) {
        MCMedge &e = **i;
        e.w -= mu;
    }
    return result;
}

MCMgraph MCMgraph::normalize(const std::map<CId, CDouble> &mu) const {
    MCMgraph result(*this);
    for (auto i = result.getEdges().begin(); i != result.getEdges().end(); i++) {
        MCMedge &e = **i;
        CDouble nc = mu.at(e.src->id);
        if (nc != -DBL_MAX) {
            e.w -= nc;
        }
    }
    return result;
}

// longest path computation, may not be implemented optimally.
std::map<CId, CDouble> MCMgraph::longestPaths(const CId rootNodeId) const {
    std::map<CId, CDouble> result;
    for (auto i = this->nodes.begin(); i != this->nodes.end(); i++) {
        const MCMnode &n = **i;
        result[n.id] = (n.id == rootNodeId) ? 0.0 : -DBL_MAX;
    }
    bool changed = true;
    while (changed) {
        changed = false;
        for (auto i = this->edges.begin(); i != this->edges.end(); i++) {
            MCMedge &e = **i;
            if (result[e.src->id] != -DBL_MAX) {
                if (result[e.src->id] + e.w > result[e.dst->id]) {
                    changed = true;
                    result[e.dst->id] = result[e.src->id] + e.w;
                }
            }
        }
    }
    return result;
}

std::map<CId, CDouble> MCMgraph::normalizedLongestPaths(const CId rootNodeId,
                                                        const CDouble mu) const {
    MCMgraph normalizedGraph = this->normalize(mu);
    const std::map<CId, CDouble> result =
            normalizedGraph.longestPaths(rootNodeId); // const keeps gcc happy
    return result;
}

std::map<CId, CDouble> MCMgraph::normalizedLongestPaths(const CId rootNodeId,
                                                        const std::map<CId, CDouble> &mu) const {
    MCMgraph normalizedGraph = this->normalize(mu);
    const std::map<CId, CDouble> result =
            normalizedGraph.longestPaths(rootNodeId); // const keeps gcc happy
    return result;
}

} // namespace Graphs

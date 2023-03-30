/*
 *  Eindhoven University of Technology
 *  Eindhoven, The Netherlands
 *  Dept. of Electrical Engineering
 *  Electronics Systems Group
 *  Model Based Design Lab (https://computationalmodeling.info/)
 *
 *  Name            :   mcmgraph.h
 *
 *  Author          :   Sander Stuijk (sander@ics.ele.tue.nl)
 *
 *  Date            :   November 7, 2005
 *
 *  Function        :   Convert graph to weighted directed graph for MCM
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

#ifndef BASE_ANALYSIS_MCM_MCMGRAPH_H_INCLUDED
#define BASE_ANALYSIS_MCM_MCMGRAPH_H_INCLUDED

#include "base/basic_types.h"
#include <memory>
#include <utility>

namespace Graphs {
class MCMnode;

class MCMedge {
public:
    // Constructor
    MCMedge(CId eId, bool eVisible);
    CId id;
    bool visible;
    std::shared_ptr<MCMnode> src;
    std::shared_ptr<MCMnode> dst;
    CDouble w;
    CDouble d;
};

using MCMedges = list<std::shared_ptr<MCMedge>>;

class MCMnode {
public:
    // Constructor
    MCMnode(CId nId, bool nVisible);
    CId id;
    bool visible;
    MCMedges in;
    MCMedges out;
};

struct MCMNodeLess {
    bool operator()(std::shared_ptr<MCMnode> const &lhs, std::shared_ptr<MCMnode> const &rhs) const { return lhs->id < rhs->id; };
};

using MCMnodes = list<std::shared_ptr<MCMnode>>;

class MCMgraph {
public:
    // Constructor
    MCMgraph();

    // Destructor
    ~MCMgraph();

    // Copy Constructor
    MCMgraph(const MCMgraph &g);

    const MCMnodes &getNodes() const { return nodes; };

    uint nrVisibleNodes() const {
        uint nrNodes = 0;
        for (auto &node : nodes) {
            if (node->visible) {
                nrNodes++;
            }
        }
        return nrNodes;
    };
    std::shared_ptr<MCMnode> getNode(CId id) {
        for (auto &node : nodes) {
            if (node->id == id) {
                return node;
            }
        }
        return nullptr;
    };

    const MCMedges &getEdges() const { return edges; };

    std::shared_ptr<MCMedge> getEdge(CId id) {
        for (auto &edge : edges) {
            if (edge->id == id) {
                return edge;
            }
        }
        return nullptr;
    };

    std::shared_ptr<MCMedge> getEdge(CId srcId, CId dstId) {
        for (auto &edge : edges) {
            if (edge->src->id == srcId) {
                if (edge->dst->id == dstId) {
                    return edge;
                }
            }
        }
        return nullptr;
    };

    [[nodiscard]] uint nrVisibleEdges() const {
        uint nrEdges = 0;
        for (const auto &edge : edges) {
            if (edge->visible) {
                nrEdges++;
            }
        }
        return nrEdges;
    };

    // Construction

    // Add a node to the MCM graph
    void addNode(const std::shared_ptr<MCMnode> &n) {
        // Add the node to the MCM graph
        this->nodes.push_back(n);
    }

    // Remove a node from the MCMgraph.
    // Note: containers of nodes are lists, so remove is expensive!
    void removeNode(const std::shared_ptr<MCMnode> &n) {
        // remove any remaining edges
        while (!n->in.empty()) {
            this->removeEdge(*(n->in.begin()));
        }
        while (!n->out.empty()) {
            this->removeEdge(*(n->out.begin()));
        }

        this->nodes.remove(n);
    }

    // Add an edge to the MCMgraph.
    std::shared_ptr<MCMedge> addEdge(CId id,
                                     std::shared_ptr<MCMnode> src,
                                     std::shared_ptr<MCMnode> dst,
                                     CDouble w,
                                     CDouble d) {
        std::shared_ptr<MCMedge> e = std::make_shared<MCMedge>(id, true);
        e->src = std::move(src);
        e->dst = std::move(dst);
        e->w = w;
        e->d = d;
        this->addEdge(e);
        return e;
    }

    // Add an edge to the MCMgraph.
    void addEdge(std::shared_ptr<MCMedge> e) {
        this->edges.push_back(e);
        e->src->out.push_back(e);
        e->dst->in.push_back(e);
    }

    // Remove an edge from the MCMgraph.
    // Note: containers of edges are lists, so remove is expensive!
    void removeEdge(std::shared_ptr<MCMedge> e) {
        this->edges.remove(e);
        e->src->out.remove(e);
        e->dst->in.remove(e);
    }

    void relabelNodeIds(std::map<int, int> &nodeIdMap);

    // reduce the MCM graph by removing obviously redundant edges
    // in particular if there are multiple edges between the same pair
    // of nodes and for some edge (w1, d1) there exists a different edge
    // (w2, d2) such that d2<=d1 and w2>=w1, then (w2, d2) is removed
    // Note this algorithm does currently not distinguish visible and invisible edges!
    std::shared_ptr<MCMgraph> pruneEdges();

    [[nodiscard]] CDouble calculateMaximumCycleMeanKarp() const;
    [[nodiscard]] CDouble calculateMaximumCycleMeanKarpDouble(MCMnode** criticalNode = nullptr) const;

    [[nodiscard]] CDouble calculateMaximumCycleRatioAndCriticalCycleYoungTarjanOrlin(std::shared_ptr<MCMedge> **cycle = nullptr,
                                                                       uint *len = nullptr);

    [[nodiscard]] MCMgraph normalize(CDouble mu) const;
    [[nodiscard]] MCMgraph normalize(const std::map<CId, CDouble> &mu) const;
    [[nodiscard]] std::map<CId, CDouble> longestPaths( CId rootNodeId) const;
    [[nodiscard]] std::map<CId, CDouble> normalizedLongestPaths( CId rootNodeId,  CDouble mu) const;
    [[nodiscard]] std::map<CId, CDouble> normalizedLongestPaths( CId rootNodeId,
                                                  const std::map<CId, CDouble> &) const;

private:
    // Nodes
    MCMnodes nodes;

    // Edges
    MCMedges edges;
};

using MCMgraphs = list<std::shared_ptr<MCMgraph>>;
using MCMgraphsIter = MCMgraphs::iterator;

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
                               bool includeComponentsWithoutEdges = false);

/**
 * relabelMCMgraph ()
 * The function removes all hidden nodes and edges from the graph. All visible
 * edges are assigned a new id starting in the range [0,nrNodes()).
 */
void relabelMCMgraph(std::shared_ptr<MCMgraph> g);

/**
 * addLongestDelayEdgesToMCMgraph ()
 * The function adds additional edges to the graph which express the
 * longest path between two nodes crossing one edge with a delay. Edges
 * with no delay are removed and edges with more then one delay element
 * are converted into a sequence of edges with one delay element.
 */
void addLongestDelayEdgesToMCMgraph(std::shared_ptr<MCMgraph> g);

} // namespace Graphs
#endif

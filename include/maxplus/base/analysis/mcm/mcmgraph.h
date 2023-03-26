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

namespace Graphs
{
    class MCMnode;

    class MCMedge
    {
        public:
            // Constructor
            MCMedge(CId eId, bool eVisible);
            CId id;
            bool visible;
            MCMnode *src;
            MCMnode *dst;
            CDouble w;
            CDouble d;
    };


    typedef list<MCMedge *>             MCMedges;
    typedef MCMedges::iterator          MCMedgesIter;
    typedef MCMedges::const_iterator    MCMedgesCIter;

    class MCMnode
    {
        public:
            // Constructor
            MCMnode(CId nId, bool nVisible);
            CId id;
            bool visible;
            MCMedges in;
            MCMedges out;
    };

	struct MCMNodeLess {
		bool operator ()(MCMnode* const& lhs, MCMnode* const& rhs) const
		{
			return lhs->id < rhs->id;
		};
	};

    typedef list<MCMnode *>             MCMnodes;
    typedef MCMnodes::iterator          MCMnodesIter;
    typedef MCMnodes::const_iterator    MCMnodesCIter;

    class MCMgraph
    {
        public:
            // Constructor
            MCMgraph();

            // Destructor
            ~MCMgraph();

            // Copy Constructor
            MCMgraph(const MCMgraph& g);

            const MCMnodes &getNodes()
            {
                return nodes;
            };

            uint nrVisibleNodes()
            {
                uint nrNodes = 0;
                for (MCMnodesIter iter = nodes.begin(); iter != nodes.end(); iter++)
                    if ((*iter)->visible) nrNodes++;
                return nrNodes;
            };
            MCMnode *getNode(CId id)
            {
                for (MCMnodes::iterator i = nodes.begin(); i != nodes.end(); i++)
                    if ((*i)->id == id)
                        return (*i);
                return NULL;
            };

            const MCMedges &getEdges()
            {
                return edges;
            };

            MCMedge *getEdge(CId id)
            {
                for (MCMedges::iterator i = edges.begin(); i != edges.end(); i++)
                    if ((*i)->id == id)
                        return (*i);
                return NULL;
            };

			MCMedge *getEdge(CId srcId, CId dstId)
			{
				for (MCMedges::iterator i = edges.begin(); i != edges.end(); i++)
					if ((*i)->src->id == srcId) {
						if ((*i)->dst->id == dstId) {
							return (*i);
						}
					}
				return NULL;
			};

            uint nrVisibleEdges()
            {
                uint nrEdges = 0;
                for (MCMedgesIter iter = edges.begin(); iter != edges.end(); iter++)
                    if ((*iter)->visible) nrEdges++;
                return nrEdges;
            };

            // Construction

            // Add a node to the MCM graph
            void addNode(MCMnode *n)
            {
                // Add the node to the MCM graph
                this->nodes.push_back(n);
            }

            // Remove a node from the MCMgraph.
            // Note: containers of nodes are lists, so remove is expensive!
            void removeNode(MCMnode *n)
            {
                // remove any remaining edges
                while (! n->in.empty())
                {
                    this->removeEdge(*(n->in.begin()));
                }
                while (! n->out.empty())
                {
                    this->removeEdge(*(n->out.begin()));
                }

                this->nodes.remove(n);
            }

            // Add an edge to the MCMgraph.
            MCMedge *addEdge(CId id, MCMnode *src, MCMnode *dst, CDouble w, CDouble d)
            {
                MCMedge *e = new MCMedge(id, true);
                e->src = src;
                e->dst = dst;
                e->w = w;
                e->d = d;
                this->addEdge(e);
                return e;
            }


            // Add an edge to the MCMgraph.
            void addEdge(MCMedge *e)
            {
                this->edges.push_back(e);
                e->src->out.push_back(e);
                e->dst->in.push_back(e);
            }

            // Remove an edge from the MCMgraph.
            // Note: containers of edges are lists, so remove is expensive!
            void removeEdge(MCMedge *e)
            {
                this->edges.remove(e);
                e->src->out.remove(e);
                e->dst->in.remove(e);
            }

            void relabelNodeIds(std::map<int, int>& nodeIdMap);

            // reduce the MCM graph by removing obviously redundant edges
            // in particular if there are multiple edges between the same pair
            // of nodes and for some edge (w1, d1) there exists a different edge
            // (w2, d2) such that d2<=d1 and w2>=w1, then (w2, d2) is removed
            // Note this algorithm does currently not distinguish visible and invisible edges!
            MCMgraph *pruneEdges(void);

            CDouble calculateMaximumCycleMeanKarp();
            CDouble calculateMaximumCycleMeanKarpDouble(MCMnode** criticalNode = NULL);

            CDouble calculateMaximumCycleRatioAndCriticalCycleYoungTarjanOrlin(MCMedge*** cycle = NULL, uint* len = NULL);

            MCMgraph normalize(CDouble mu) const;
            MCMgraph normalize(const std::map<CId,CDouble>& mu) const;
            std::map<CId, CDouble> longestPaths(const CId rootNodeId) const;
            std::map<CId, CDouble> normalizedLongestPaths(const CId rootNodeId, const CDouble mu) const;
            std::map<CId, CDouble> normalizedLongestPaths(const CId rootNodeId, const std::map<CId, CDouble>&) const;

        private:
			// Nodes
            MCMnodes nodes;

            // Edges
            MCMedges edges;

    };

    typedef list<MCMgraph *>      MCMgraphs;
    typedef MCMgraphs::iterator  MCMgraphsIter;


    /**
     * Extract the strongly connected components from the graph. These components
     * are returned as a set of MCM graphs. All nodes which belong to at least
     * one of the strongly connected components are set to visible in the graph g,
     * all other nodes are made invisible. Also edges between two nodes in (possibly
     * different) strongly connected components are made visible and all others
     * invisible. The graph g consists in the end of only nodes which are part of
     * a strongly connnected component and all the edges between these nodes. Some
     * MCM algorithms work also on this graph (which reduces the execution time
     * needed in some of the conversion algorithms).
     */
    void stronglyConnectedMCMgraph(MCMgraph *g, MCMgraphs &components, bool includeComponentsWithoutEdges = false);

    /**
     * relabelMCMgraph ()
     * The function removes all hidden nodes and edges from the graph. All visible
     * edges are assigned a new id starting in the range [0,nrNodes()).
     */
    void relabelMCMgraph(MCMgraph *g);

    /**
    * addLongestDelayEdgesToMCMgraph ()
    * The function adds additional edges to the graph which express the
    * longest path between two nodes crossing one edge with a delay. Edges
    * with no delay are removed and edges with more then one delay element
    * are converted into a sequence of edges with one delay element.
    */
    void addLongestDelayEdgesToMCMgraph(MCMgraph *g);

}//namespace SDF
#endif

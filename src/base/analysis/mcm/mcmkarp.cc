/*
 *  Eindhoven University of Technology
 *  Eindhoven, The Netherlands
 *  Dept. of Electrical Engineering
 *  Electronics Systems Group
 *  Model Based Design Lab (https://computationalmodeling.info/)
 *
 *  Name            :   mcmkarp.cc
 *
 *  Author          :   Sander Stuijk (sander@ics.ele.tue.nl)
 *
 *  Date            :   November 8, 2005
 *
 *  Function        :   Compute the MCM for an HSDF graph using Karp's
 *                      algorithm.
 *
 *  History         :
 *      08-11-05    :   Initial version.
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
#include "base/math/cmath.h"
#include <climits>


namespace Graphs {
/**
 * mcmKarp ()
 * The function computes the maximum cycle mean of an MCMgGraph using Karp's
 * algorithm.
 * Note that the following assumptions are made about the MCMgraph
 * 1. it is assumed that the edge weights have integer values.
 * 2. it is assumed that all nodes in the graph are 'visible'
 * 3. it is assumed that the node have id's ranging from 0 up to the number of nodes.
 *
 * The critical cycle is only returned if cycle and len are not NULL. Then *cycle points
 * to an array of MCMEdges of the critical cycle and *len indicates the length of the cycle.
 * *cycle is a freshly allocated array and it is the caller's obligation to deallocate it
 * in due time.
 * cycle  - pointers to arcs on maximum ratio cycle,
 *          ordered in array from top to bottom with
 *          respect to subsequent arcs on cycle
 *
 * len    - number of elements of "cycle"
 */
CDouble maximumCycleMeanKarp(MCMgraph *mcmGraph) {
    int k, n;
    long long **d;
    double l, ld;
    MCMnode *u;

    // Allocate memory d[n+1][n]
    n = mcmGraph->nrVisibleNodes();
    d = new long long *[n + 1];
    for (int i = 0; i < n + 1; i++)
        d[i] = new long long[n];

    // Initialize
    // d[k][u], 1<=k<n+1, 0<=u<n with value -inf
    for (k = 1; k < n + 1; k++)
        for (int u = 0; u < n; u++)
            d[k][u] = -INT_MAX;
    // d[0][u], 0<=u<n with value 0
    for (int u = 0; u < n; u++)
        d[0][u] = 0;

    // Compute the distances
    for (k = 1; k < n + 1; k++) {
        for (MCMnodesCIter iter = mcmGraph->getNodes().begin(); iter != mcmGraph->getNodes().end();
             iter++) {
            MCMnode *v = *iter;

            for (MCMedgesIter e = v->in.begin(); e != v->in.end(); e++) {
                MCMnode *u = (*e)->src;

                d[k][v->id] = MAX(d[k][v->id], d[k - 1][u->id] + ((int)(*e)->w));
            }
        }
    }

    // Compute lambda using Karp's theorem
    l = -INT_MAX;
    for (MCMnodesCIter iter = mcmGraph->getNodes().begin(); iter != mcmGraph->getNodes().end();
         iter++) {
        u = *iter;
        ld = INT_MAX;
        for (k = 0; k < n; k++) {
            ld = MIN(ld, (double)(d[n][u->id] - d[k][u->id]) / (double)(n - k));
        }
        l = MAX(l, ld);
    }

    // Cleanup
    for (int i = 0; i < n + 1; i++)
        delete[] d[i];
    delete[] d;

    return l;
}

/**
 * maximumCycleMeanKarpDouble ()
 * The function computes the maximum cycle mean of an MCMgGraph using Karp's
 * algorithm.
 * Note that the following assumptions are made about the MCMgraph
 * 1. it is assumed that all nodes in the graph are 'visible'
 * 2. it is assumed that the node have id's ranging from 0 up to the number of nodes.
 *
 * A critical node is only returned if criticalNode is not NULL.
 */

CDouble maximumCycleMeanKarpDouble(MCMgraph *mcmGraph, MCMnode **criticalNode = NULL) {
    int k, n;
    typedef CDouble *CDoublePtr;
    CDoublePtr *d;
    double l, ld;
    MCMnode *u;

    // Allocate memory d[n+1][n]
    n = mcmGraph->nrVisibleNodes();
    d = new CDoublePtr[n + 1];
    for (int i = 0; i < n + 1; i++)
        d[i] = new CDouble[n];

    // Initialize
    // d[k][u], 1<=k<n+1, 0<=u<n with value -inf
    for (k = 1; k < n + 1; k++)
        for (int u = 0; u < n; u++)
            d[k][u] = -DBL_MAX;
    // d[0][u], 0<=u<n with value 0
    for (int u = 0; u < n; u++)
        d[0][u] = 0.0;

    // Compute the distances
    for (k = 1; k < n + 1; k++) {
        for (MCMnodesCIter iter = mcmGraph->getNodes().begin(); iter != mcmGraph->getNodes().end();
             iter++) {
            MCMnode *v = *iter;

            for (MCMedgesIter e = v->in.begin(); e != v->in.end(); e++) {
                MCMnode *u = (*e)->src;

                d[k][v->id] = MAX(d[k][v->id], d[k - 1][u->id] + ((*e)->w));
            }
        }
    }

    // Compute lambda using Karp's theorem
    l = -DBL_MAX;
    for (MCMnodesCIter iter = mcmGraph->getNodes().begin(); iter != mcmGraph->getNodes().end();
         iter++) {
        u = *iter;
        ld = DBL_MAX;
        for (k = 0; k < n; k++) {
            CDouble nld = (d[n][u->id] - d[k][u->id]) / (double)(n - k);
            if (nld < ld) {
                ld = nld;
            }
        }
        if (ld > l) {
            l = ld;
            if (criticalNode != NULL) {
                *criticalNode = u;
            }
        }
    }

    // Cleanup
    for (int i = 0; i < n + 1; i++)
        delete[] d[i];

    delete[] d;

    return l;
}

} // namespace Graphs

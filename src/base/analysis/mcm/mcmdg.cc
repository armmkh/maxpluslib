/*
 *  Eindhoven University of Technology
 *  Eindhoven, The Netherlands
 *  Dept. of Electrical Engineering
 *  Electronics Systems Group
 *  Model Based Design Lab (https://computationalmodeling.info/)
 *
 *  Name            :   mcmhoward.cc
 *
 *  Author          :   Sander Stuijk (sander@ics.ele.tue.nl)
 *
 *  Date            :   November 8, 2005
 *
 *  Function        :   Compute the MCM for an HSDF graph using Dasdan-Gupta's
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
#include <memory>


namespace Graphs {
/**
 * mcmDG ()
 * The function computes the maximum cycle mean of a HSDF graph using
 * Dasdan-Gupta's algorithm.
 * Note: this algorithm assumes that edge weights are integer valued !
 * @todo
 * check if algorithm can be generalized to float edge weights
 */
CDouble mcmDG(MCMgraph *mcmGraph) {
    int k, n;
    int *level;
    int **pi, **d;
    CDouble l, ld;
    std::shared_ptr<MCMnode> u;
    list<int> Q_k;
    list<std::shared_ptr<MCMnode>> Q_u;

    // Allocate memory
    n = mcmGraph->nrVisibleNodes();
    level = new int[n];
    pi = new int *[n + 1];
    d = new int *[n + 1];
    for (int i = 0; i < n + 1; i++) {
        pi[i] = new int[n];
        d[i] = new int[n];
    }

    // Initialize
    for (int i = 0; i < n; i++)
        level[i] = -1;
    d[0][0] = 0;
    pi[0][0] = -1;
    level[0] = 0;
    Q_k.push_back(0);
    Q_u.push_back(mcmGraph->getNodes().front());

    // Compute the distances
    k = Q_k.front();
    Q_k.pop_front();
    u = Q_u.front();
    Q_u.pop_front();
    do {
        for (auto iter = u->out.begin(); iter != u->out.end(); iter++) {
            std::shared_ptr<MCMedge> e = *iter;
            std::shared_ptr<MCMnode> v = e->dst;

            if (level[v->id] < k + 1) {
                Q_k.push_back(k + 1);
                Q_u.push_back(v);
                pi[k + 1][v->id] = level[v->id];
                level[v->id] = k + 1;
                d[k + 1][v->id] = -INT_MAX;
            }
            d[k + 1][v->id] = MAX(d[k + 1][v->id], d[k][u->id] + ((int)e->w));
        }
        k = Q_k.front();
        Q_k.pop_front();
        u = Q_u.front();
        Q_u.pop_front();
    } while (k < n);

    // Compute lambda using Karp's theorem
    l = -INT_MAX;
    for (auto iter = mcmGraph->getNodes().begin(); iter != mcmGraph->getNodes().end();
         iter++) {
        u = *iter;

        if (level[u->id] == n) {
            ld = INT_MAX;
            k = pi[n][u->id];
            while (k > -1) {
                ld = MIN(ld, (CDouble)(d[n][u->id] - d[k][u->id]) / (CDouble)(n - k));
                k = pi[k][u->id];
            }
            l = MAX(l, ld);
        }
    }

    // Cleanup
    delete[] level;
    for (int i = 0; i < n + 1; i++) {
        delete[] pi[i];
        delete[] d[i];
    }
    delete[] pi;
    delete[] d;

    return l;
}

} // namespace Graphs
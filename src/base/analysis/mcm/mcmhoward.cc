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
 *  Function        :   Compute the MCM for an HSDF graph using Howard's
 *                      algorithm implemented in Max-Plus algebra.
 *
 *  Acknowledgement :   This code is based on 'Howard Policy Iteration Algorithm
 *                      for Max Plus Matrices' written by Stephane Gaubert
 *                      (Stephane.Gaubert@inria.fr). The max-plus version of
 *                      Howard's algorithm is described in the paper:
 *                      'Numerical computation of spectral elements in max-plus
 *                      algebra'. Jean Cochet-Terrasson, Guy Cohen, Stephane
 *                      Gaubert, Michael Mc Gettrick, Jean-Pierre Quadrat
 *                      IFAC Workshop on System Structure and Control,
 *                      Nantes, July 1997.
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
#include "base/exception/exception.h"

#include <math.h>
#include <memory>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>

namespace Graphs {
/*
 * Howard terminates with an error if MAX_NIterations occur.
 */
#define MAX_NIterations 1000

#define EPSILON -HUGE_VAL

static int *ij;
static double *a;
static int nr_nodes;
static int narcs;
static double *chi;
static double *v;
static int *pi;
static int *NIterations;
static int *NComponents;

static int *new_pi = NULL; /*  new policy */

/* The inverse policy is coded by a linearly chained list.
 * pi_inv_idx[i]= pointer to the chain of inverses of node i.
 */
static int *pi_inv_idx = NULL;

/* pi_inv_idx[j]= pointer to the next inverse */
static int *pi_inv_succ = NULL;

/* corresponding node  */
static int *pi_inv_elem = NULL;

/* pi_inv_last[i]= last inverse of i */
static int *pi_inv_last = NULL;

static double *c = NULL;
static double *v_aux = NULL;
static double *new_c = NULL;
static double *new_chi = NULL;
static int *visited = NULL;
static int *component = NULL;
static double lambda = 0;
static double epsilon = 0;
static int color = 1;

/**
 * Epsilon ()
 * The termination tests are performed up to an epsilon constant, which is fixed
 * heuristically by the following routine.
 */
static void Epsilon(double *a, int narcs, double *epsilon) {
    int i;
    double MAX, MIN;

    MAX = a[0];
    MIN = a[0];

    for (i = 1; i < narcs; i++) {
        if (a[i] > MAX)
            MAX = a[i];
        if (a[i] < MIN)
            MIN = a[i];
    }

    *epsilon = (MAX - MIN) * 0.000000001;
}

/**
 * Initial_Policy
 * Build an admissible policy pi and its associated cost vector c from ij and A.
 * Reasonable greedy rule to determine the first policy. pi(node i) = arc with
 * maximal weight starting from i for full random matrices, this choice of
 * initial policy seems to cut the number of iterations by a factor 1.5, by
 * comparison with a random initial policy.
 */
static void Initial_Policy() {
    int i;

    /* we loose a O(nr_nodes) time here ... */
    /* we use the auxiliary variable v_aux to compute the row max of A */
    for (i = 0; i < nr_nodes; i++)
        v_aux[i] = EPSILON;

    for (i = 0; i < narcs; i++) {
        if (v_aux[ij[i * 2]] <= a[i]) {
            pi[ij[i * 2]] = ij[i * 2 + 1];
            c[ij[i * 2]] = a[i];
            v_aux[ij[i * 2]] = a[i];
        }
    }
}

static void New_Build_Inverse() {
    int i, j, locus;
    int ptr = 0;

    for (i = 0; i < nr_nodes; i++) {
        pi_inv_idx[i] = -1;
        pi_inv_last[i] = -1;
    }

    for (i = 0; i < nr_nodes; i++) {
        j = pi[i];
        if (pi_inv_idx[j] == -1) {
            pi_inv_succ[ptr] = -1;
            pi_inv_elem[ptr] = i;
            pi_inv_last[j] = ptr;
            pi_inv_idx[j] = ptr;
            ptr++;
        } else {
            pi_inv_succ[ptr] = -1;
            pi_inv_elem[ptr] = i;
            locus = pi_inv_last[j];
            pi_inv_succ[locus] = ptr;
            pi_inv_last[j] = ptr;
            ptr++;
        };
    }
}

static void Init_Depth_First() {
    int j;

    for (j = 0; j < nr_nodes; j++) {
        visited[j] = 0;
        component[j] = 0;
    }
}

/**
 *
 * Given the value of v at initial point i, we compute v[j] for all predecessor
 * j of i, according to the spectral equation, v[j]+ lambda = A(arc from j to i)
 * v[i] the array visited is changed by side effect.
 */
static void New_Depth_First_Label(int i) {
    int next_i, a;
    a = pi_inv_idx[i];
    while (a != -1 && visited[pi_inv_elem[a]] == 0) {
        next_i = pi_inv_elem[a];
        visited[next_i] = 1;
        v[next_i] = -lambda + c[next_i] + v[i];
        component[next_i] = color;
        chi[next_i] = lambda;
        New_Depth_First_Label(next_i);
        a = pi_inv_succ[a];
    }
}

static void Visit_From(int initial_point, int color) {
    int index, newindex, i;
    double weight;
    int length;

    index = initial_point;
    component[index] = color;
    newindex = pi[index];

    while (component[newindex] == 0) {
        component[newindex] = color;
        index = newindex;
        newindex = pi[index];
    }

    /* a cycle has been detected, since newindex is already visited */
    weight = 0;
    length = 0;
    i = index;
    do {
        weight += c[i];
        length++;
        i = pi[i];
    } while (i != index);

    lambda = weight / length;
    v[i] = v_aux[i]; /* keeping the previous value */
    chi[i] = lambda;
    New_Depth_First_Label(index);
}

/**
 * Value()
 * Computes the value (v,chi) associated with a policy pi.
 */
static void Value() {
    int initial_point;
    color = 1;

    Init_Depth_First();
    initial_point = 0;

    do {
        Visit_From(initial_point, color);
        while ((initial_point < nr_nodes) && (component[initial_point] != 0))
            initial_point++;
        color++;
    } while (initial_point < nr_nodes);

    *NComponents = --color;
}

static void Init_Improve() {
    int i;

    for (i = 0; i < nr_nodes; i++) {
        new_chi[i] = chi[i];
        v_aux[i] = v[i];
        new_pi[i] = pi[i];
        new_c[i] = c[i];
    }
}

static void First_Order_Improvement(int *improved) {
    int i;
    for (i = 0; i < narcs; i++) {
        if (chi[ij[i * 2 + 1]] > new_chi[ij[i * 2]]) {
            *improved = 1;
            new_pi[ij[i * 2]] = ij[i * 2 + 1];
            new_chi[ij[i * 2]] = chi[ij[i * 2 + 1]];
            new_c[ij[i * 2]] = a[i];
        }
    }
}

static void Second_Order_Improvement(int *improved) {
    int i;
    double w;
    if (*NComponents > 1) {
        for (i = 0; i < narcs; i++) {
            /* arc i is critical */
            if (chi[ij[i * 2 + 1]] == new_chi[ij[i * 2]]) {
                w = a[i] + v[ij[i * 2 + 1]] - chi[ij[i * 2 + 1]];
                if (w > v_aux[ij[i * 2]] + epsilon) {
                    *improved = 1;
                    v_aux[ij[i * 2]] = w;
                    new_pi[ij[i * 2]] = ij[i * 2 + 1];
                    new_c[ij[i * 2]] = a[i];
                }
            }
        }
    } else {
        /* we know that all the arcs realize the max in the
        first order improvement */
        for (i = 0; i < narcs; i++) {
            w = a[i] + v[ij[i * 2 + 1]] - chi[ij[i * 2 + 1]];
            if (w > v_aux[ij[i * 2]] + epsilon) {
                *improved = 1;
                v_aux[ij[i * 2]] = w;
                new_pi[ij[i * 2]] = ij[i * 2 + 1];
                new_c[ij[i * 2]] = a[i];
            }
        }
    }
}

static void Improve(int *improved) {
    *improved = 0;
    Init_Improve();

    /* a first order policy improvement may occur */
    if (*NComponents > 1)
        First_Order_Improvement(improved);

    if (*improved == 0)
        Second_Order_Improvement(improved);
}

static void Allocate_Memory() {
    new_pi = (int *)calloc(nr_nodes, sizeof(int));
    pi_inv_idx = (int *)calloc(nr_nodes, sizeof(int));
    pi_inv_succ = (int *)calloc(nr_nodes, sizeof(int));
    pi_inv_elem = (int *)calloc(nr_nodes, sizeof(int));
    pi_inv_last = (int *)calloc(nr_nodes, sizeof(int));
    visited = (int *)calloc(nr_nodes, sizeof(int));
    component = (int *)calloc(nr_nodes, sizeof(int));
    c = (double *)calloc(nr_nodes, sizeof(double));
    new_c = (double *)calloc(nr_nodes, sizeof(double));
    v_aux = (double *)calloc(nr_nodes, sizeof(double));
    new_chi = (double *)calloc(nr_nodes, sizeof(double));

    if ((new_chi == NULL) || (v_aux == NULL) || (new_c == NULL) || (c == NULL)
        || (component == NULL) || (visited == NULL) || (pi_inv_idx == NULL) || (pi_inv_succ == NULL)
        || (pi_inv_elem == NULL) || (pi_inv_last == NULL) || (new_pi == NULL)) {
        throw CException("Failed allocation memory");
    }
}

static void Free_Memory() {
    free(new_pi);
    free(pi_inv_idx);
    free(pi_inv_succ);
    free(pi_inv_elem);
    free(pi_inv_last);
    free(visited);
    free(component);
    free(c);
    free(new_c);
    free(v_aux);
    free(new_chi);
}

static void Check_Rows() {
    int i;
    int *u = NULL;
    u = (int *)calloc(nr_nodes, sizeof(int));

    for (i = 0; i < narcs; i++)
        u[ij[2 * i]] = 1;

    for (i = 0; i < nr_nodes; i++) {
        if (u[i] == 0)
            throw CException("Failed check on rows in Howard's MCM algorithm.");
    }

    free(u);
}

static void Security_Check() {
    if (nr_nodes < 1)
        throw CException("Howard: number of nodes must be a positive integer.");

    if (narcs < 1)
        throw CException("Howard: number of arcs must be a positive integer.");

    Check_Rows();
}

static void Import_Arguments(int *ij,
                             double *A,
                             int nr_nodes,
                             int nr_arc,
                             double *chi,
                             double *v,
                             int *policy,
                             int *nr_iterations,
                             int *nr_components) {
    ij = ij;
    a = A;
    nr_nodes = nr_nodes;
    narcs = nr_arc;
    chi = chi;
    v = v;
    pi = policy;
    NIterations = nr_iterations;
    NComponents = nr_components;
}

static void Update_Policy() {
    int i;
    for (i = 0; i < nr_nodes; i++) {
        pi[i] = new_pi[i];
        c[i] = new_c[i];
        v_aux[i] = v[i]; /* Keep a copy of the current value function */
    }
}

static void End_Message() {
    if (*NIterations == MAX_NIterations)
        throw CException("Howard: exceeded maximum number of iterations.");
}

/**
 * Howard ()
 * Howard Policy Iteration Algorithm for Max Plus Matrices.
 *
 * INPUT of Howard Algorithm:
 *      ij,A,nr_nodes,narcs : sparse description of a matrix.
 *
 * OUTPUT:
 *      chi cycle time vector
 *      v bias
 *      pi optimal policy
 *      NIterations: Number of iterations of the algorithm
 *      NComponents: Number of connected components of the optimal policy
 *
 * REQUIRES: O(nr_nodes) SPACE
 * One iteration requires: O(narcs+nr_nodes) TIME
 *
 * The following variables should be defined in the environment from which the
 * Howard routine is called.
 *
 * INPUT VARIABLES
 * int nr_nodes;  number of nodes of the graph
 * int nr_arcs;   number of arcs of the graph
 * int *ij;       array of integers of size 2*narcs
 *                for (0 <=k <narcs), the arc numbered k  goes from
 *                IJ[k][0] =(IJ[2k]) to IJ[k][1] (=IJ[2k+1])
 * double *A;     array of double of size narcs
 *                A[k]=weight of the arc numbered k
 *
 * OUTPUT VARIABLES
 * double *v;   array of double of size nr_nodes (the bias vector)
 * double *chi; array of double of size nr_nodes (the cycle time vector)
 * int *POLICY; array of integer of size nr_nodes (an optimal policy)
 * int nr_iterations; the number of iterations of the algorithm
 * int nr_components; the number of connected components of the optimal
 *               policy which is returned.
 */
void Howard(int *ij,
            double *A,
            int nr_nodes,
            int nr_arcs,
            double *chi,
            double *v,
            int *policy,
            int *nr_iterations,
            int *nr_components) {
    int improved = 0;
    *nr_iterations = 0;

    Import_Arguments(ij, A, nr_nodes, nr_arcs, chi, v, policy, nr_iterations, nr_components);
    Security_Check();
    Allocate_Memory();
    Epsilon(a, narcs, &epsilon);
    Initial_Policy();
    New_Build_Inverse();

    do {
        Value();
        Improve(&improved);
        Update_Policy();
        New_Build_Inverse();
        (*NIterations)++;
    } while ((improved != 0) && *NIterations < MAX_NIterations);

    End_Message();
    Free_Memory();
}

/**
 * convertMCMgraphToMatrix ()
 * The function converts a weighted directed graph used in the MCM algorithms
 * to a sparse matrix input for Howard's algorithm.
 */
void convertMCMgraphToMatrix(const MCMgraph& g, int *ij, double *A) {
    int k = 0;
    uint i = 0;
    uint j = 0;
    v_uint mapId(g.getNodes().size());

    // Re-map the id of all visible nodes back to the range [0, g->nrNodes())
    for (auto iter = g.getNodes().begin(); iter != g.getNodes().end(); iter++) {
        const std::shared_ptr<MCMnode>& n = *iter;

        if (n->visible) {
            mapId[i] = j;
            j++;
        }
        i++;
    }

    // Create an entry in the matrices for each edge
    for (auto iter = g.getEdges().begin(); iter != g.getEdges().end(); iter++) {
        const std::shared_ptr<MCMedge>& e = *iter;

        // Is the edge a existing edge in the graph?
        if (e->visible) {
            ij[2 * k] = mapId[e->src->id];
            ij[2 * k + 1] = mapId[e->dst->id];
            A[k] = e->w;

            // Next edge
            k++;
        }
    }
}

} // namespace Graphs
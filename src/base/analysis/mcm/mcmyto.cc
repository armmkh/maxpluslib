/*
 *  Eindhoven University of Technology
 *  Eindhoven, The Netherlands
 *  Dept. of Electrical Engineering
 *  Electronics Systems Group
 *  Model Based Design Lab (https://computationalmodeling.info/)
 *
 *  Name            :   mcmyto.cc
 *
 *  Author          :   Sander Stuijk (sander@ics.ele.tue.nl)
 *
 *  Date            :   November 8, 2005
 *
 *  Function        :   Compute the MCM for an HSDF graph using Young-
 *                      Tarjan-Orlin's algorithm.
 *
 *  Acknowledgement :   This code is based the 'mmcycle' program which can be
 *                      found at 'http:// elib.zib.de/pub/Packages/mathprog
 *                      /netopt/mmc-info'. The original code is written by
 *                      Georg Skorobohatyj (bzfskoro@zib.de) and is free
 *                      for redistribution.
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
#include "base/lookup/clookup.h"
#include "base/exception/exception.h"
#include "base/analysis/mcm/mcmyto.h"


#include <float.h>


namespace Graphs
{


#define NILN (node *) NULL
#define  NILA (arc *) NULL

    /*
     * d_heap implementation
     * ---------------------
     * See R.E. Tarjan: Data Structures and Network
     * Algorithms (Society for Industrial and Applied
     * Mathematics)", Philadelphia, PA, 1983) for a
     * decsription.
     *
     * Indices are numbered from 0 to hz-1 , hz =
     * maximum heap size, therefore the parent of
     * node x is int((x-1)/dh), dh = d_heap parameter
     * = number of children per node, and the children
     * of node x are the nodes in the interval
     *
     * [dh*x+1, dh*x+2, ..., min (dh *(x+1), hz-1)].
     */

#define  dh 4L

    typedef  arc  *item;

    typedef struct D_heap
    {
        long max_size;
        long size;
        item *items;
    } d_heap;

#if 0
    static
    void printHeap(d_heap *h)
    {
        cerr << "Heap: ";
        for (long i = 0; i < h->size; i++)
            cerr << h->items[i]->key << ", ";
        cerr << endl;
    }
#endif

    static inline
    long MINCHILD(long x, d_heap *h)
    {
        double min;
        long k, kmin, upb;

        kmin = -1;
        if (x != h->size - 1)
        {
            min = DBL_MAX;
            if (h->size - 1 < dh * (x + 1))
                upb = h->size - 1;
            else
                upb = dh * (x + 1);

            for (k = dh * x + 1; k <= upb; k++)
            {
                if (h->items[k]->key < min)
                {
                    min = h->items[k]->key;
                    kmin = k;
                }
            }
        }

        return (kmin);
    }

	static inline
		long MINCHILD_ROBUST(long x, d_heap* h, double epsilon)
	{
		double min;
		long k, kmin, upb;

		kmin = -1;
		if (x != h->size - 1)
		{
			min = DBL_MAX;
			if (h->size - 1 < dh * (x + 1))
				upb = h->size - 1;
			else
				upb = dh * (x + 1);

			for (k = dh * x + 1; k <= upb; k++)
			{
				if (min - h->items[k]->key > epsilon)
				{
					min = h->items[k]->key;
					kmin = k;
				}
			}
		}

		return (kmin);
	}


    static inline
    void SIFTDOWN(d_heap *h, item i, long x)
    {
        long c, cc;

        c = MINCHILD(x, h);
        while (c >= 0 && h->items[c]->key < i->key)
        {
            cc = MINCHILD(c, h);
            h->items[x] = h->items[c];
            h->items[x]->hpos = x;
            x = c;
            c = cc;
        }
        h->items[x] = i;
        i->hpos = x;
    }

	static inline
		void SIFTDOWN_ROBUST(d_heap* h, item i, long x, double epsilon)
	{
		long c, cc;

		c = MINCHILD_ROBUST(x, h, epsilon);
		while (c >= 0 && i->key - h->items[c]->key > epsilon)
		{
			cc = MINCHILD_ROBUST(c, h, epsilon);
			h->items[x] = h->items[c];
			h->items[x]->hpos = x;
			x = c;
			c = cc;
		}
		h->items[x] = i;
		i->hpos = x;
	}

    static inline
    void SIFTUP(d_heap *h, item i, long x)
    {
        long p;

        p = (x == 0) ? - 1 : (x - 1) / dh;
        while (p >= 0 && h->items[p]->key > i->key)
        {
            h->items[x] = h->items[p];
            (h->items[x])->hpos = x;
            x = p;
            p = (p == 0) ? -1 : (p - 1) / dh;
        }
        h->items[x] = i;
        i->hpos = x;
    }

	static inline
		void SIFTUP_ROBUST(d_heap* h, item i, long x, double epsilon)
	{
		long p;

		p = (x == 0) ? -1 : (x - 1) / dh;
		while (p >= 0 && h->items[p]->key - i->key > epsilon)
		{
			h->items[x] = h->items[p];
			(h->items[x])->hpos = x;
			x = p;
			p = (p == 0) ? -1 : (p - 1) / dh;
		}
		h->items[x] = i;
		i->hpos = x;
	}

    static inline
    void INSERT(d_heap *h, item i)
    {
        SIFTUP(h, i, h->size);
        ++(h->size);
    }

	static inline
		void INSERT_ROBUST(d_heap* h, item i, double epsilon)
	{
		SIFTUP_ROBUST(h, i, h->size, epsilon);
		++(h->size);
	}

#ifdef _MSC_VER
    // Remove define DELETE because of a conflict with a define in winnt.h.
#undef DELETE
#endif

    static inline
    void DELETE(d_heap *h, item i)
    {
        item j;

        j = h->items[h->size - 1];
        --(h->size);
        if (i != j)
        {
            if (j->key <= i->key)
                SIFTUP(h, j, i->hpos);
            else
                SIFTDOWN(h, j, i->hpos);
        }
    }

	static inline
		void DELETE_ROBUST(d_heap* h, item i, double epsilon)
	{
		item j;

		j = h->items[h->size - 1];
		--(h->size);
		if (i != j)
		{
			if (!(j->key - i->key > epsilon))
				SIFTUP_ROBUST(h, j, i->hpos, epsilon);
			else
				SIFTDOWN_ROBUST(h, j, i->hpos, epsilon);
		}
	}

    static inline
    item DELETE_MIN(d_heap *h)
    {
        item i;

        if (h->size == 0)
            return (NILA);
        else
        {
            i = h->items[0];
            DELETE(h, i);
            return (i);
        }
    }

    static inline
    item GET_MIN(d_heap *h)
    {
        if (h->size == 0)
            return (NILA);
        else
            return (h->items[0]);
    }

    static inline
    bool ALLOC_HEAP(d_heap *h, long k)
    {
        h->items = (item *) malloc(k * sizeof(item));
        if (h->items == (item *) NULL)
            return false;
        else
        {
            h->max_size = k;
            h->size = 0;
            return true;
        }
    }

    static inline
    void DEALLOC_HEAP(d_heap *h)
    {
        free(h->items);
    }

    static long update_level;
    static node *upd_nodes;

    /**
     * update_subtree ()
     * recursive subtree traversal function, produces a one-way liked list of nodes
     * contained in subtree updates node levels and costs of paths along sub-tree to
     * nodes contained in it.
     */
    static
    void update_subtree(node *root)
    {
        node *vptr;

        root->level = update_level;
        root->link = upd_nodes;
        upd_nodes = root;
		root->in_list = true;

        if (root->first_child != NILN)
        {
            ++update_level;
            vptr = root->first_child;
            do
            {
                vptr->cost_t = root->cost_t + vptr->parent_in->cost;
                vptr->transit_time_t = root->transit_time_t
                                       + vptr->parent_in->transit_time;
                update_subtree(vptr);
                vptr = vptr->right_sibl;
            }
            while (vptr != root->first_child);
            --update_level;
        }
    }

    /**
     * mmcycle ()
     *
     * Determines maximum ratio cycle in a directed graph
     * G = (V, E, c), c: E->IR a "cost" function on the
     * edges, alternatively called "length" or "weight".
     *
	 * Note that it is assumed that every cycle has transit time  > 0 !
	 * 
	 * Note that it has been observed that the algorithm can behave non-deterministically 
	 * across compilers, when the graph has multiple equally critical cycles
	 * Due to the order in which floating point calculations are scheduled, the
	 * algorithm's control flow may follow different paths and lead to different 
	 * critical cycles as the output.
	 * The function mmcycle_robust implements a version of the algorithm that is
	 * deterministic across compilers. It may be a little slower.
	 *
     * Input parameter:
     * ---------------
     * gr     - graph structure with incidence lists of
     *          incoming and outgoing edges for each node
     *
     * Output parameters:
     * -----------------
     * lambda - maximum cycle ratio
     *
     * cycle  - pointers to arcs on maximum ratio cycle,
     *          ordered in array from top to bottom with
     *          respect to subsequent arcs on cycle
     *
     * len    - number of elements of "cycle"
     *
     * If cycle or len is a NULL-pointer, then these parameters
     * are not assigned a value.
     *
     * Reference
     * ---------
     * N.E. Young, R. E. Tarjan, J. B. Orlin: "Faster Parametric
     * Shortest Path and Minimum-Balance Algorithms", Networks
     * 21 (1991), 205-221
     *
     *
     * Sketch of algorithm:
     * -------------------
     *
     * 1) Introduce an extra node s and an emanating arc with
     * cost zero to each node of the graph, let V' and E'
     * be the extended node and edge sets respectively and
     * G' = (V',E',c).
     *
     * 2) Let lambda = lambda_ini = sum (abs(c(e)) + 1.0 over
     * all edges e of E';
     *
     * 3) Introduce edge costs c_lambda(e) = c(e) - lambda for
     * all edges e of E', let G'_lambda = (V', E', c_lambda),
     * (all edge costs c_lambda(e) are positive);
     *
     * 4) Set up tree T_lambda of shortest paths from s to all
     * other nodes in G'_lambda, i.e. with respect to edge
     * costs c_lambda, this tree consists of all edges ema-
     * nating from node s, for each node v retain cost ct(v)
     * of tree path in G' and the number of edges on the path
     * from s to v, i.e. the tree level of v;
     *
     * 5) For all edges e = (u,v) in G - T_lambda compute edge
     * (pivot) key pk as
     *
     * pk(e) = (ct(u) + c(u,v) - ct(v))/(lev(u) + 1 - lev(v))
     *
     * if lev(u) + 1 > lev(v),
     * otherwise set pk(e) = -infinite
     *
     * For each node v select an incoming edge whose key is
     * maximum among all incoming edges and assign it to v
     * as vertex key.
     *
     * 6) Determine an edge emin = (u',v') such that pk(u',v')
     * is maximum among all edge keys by way of vertex keys -
     * such a key always exists at this point - let lambda =
     * pk(u', v');
     *
     * 7) If v' is an ancestor of u' in the tree, inserting edge
     * (u',v') into it creates a cycle with mean value lambda,
     * this is a maximum mean cycle, return cycle and lambda;
     *
     * 8) Delete edge (u'',v') from the tree and insert edge (u',v')
     * instead, for all nodes v of subtree rooted at v', update
     * lev (v) and ct(v) by adding lev(u') + 1 - lev (v') and
     * ct(u') + c(u',v') - ct(v') respectively, compute edge key
     * for edge (u'',v'), remove edge key of (u',v') and update
     * edge keys for all edges (u,v) with exactly one of its
     * incident nodes u or v in subtree rooted at v' (such edges
     * are the final ones on new shortest paths to destination
     * node v), determine new vertex key for each vertex that
     * has an incoming edge whose key has been changed;
     *
     * 9) goto (6);
     */
    void mmcycle(graph *gr, double *lambda, arc **cycle, long *len)
    {
        double min, infty, akey;
        arc *aptr, *par_aptr, *vmin_aptr, *min_aptr;
        node *sptr, *uptr, *vptr, *wptr;
        bool foundCycle;
        d_heap h;

        // set up initial tree
        sptr = gr->vs;
        sptr->in_list = false;
        sptr->cost_t = 0.0L;
        sptr->transit_time_t = 0.0L;
        sptr->level = 0L;
        sptr->parent_in = NILA;
        sptr->left_sibl = sptr;
        sptr->right_sibl = sptr;
        aptr = sptr->first_arc_out;
        vptr = aptr->head;
        sptr->first_child = vptr;
        while (aptr != NILA)
        {
            wptr = aptr->head;
            wptr->cost_t = 0.0L;
            wptr->transit_time_t = 0.0L;
            wptr->level = 1L;
            wptr->parent_in = aptr;
            wptr->first_child = NILN;
            wptr->left_sibl = vptr;
            vptr->right_sibl = wptr;
            wptr->in_list = false;
            aptr->in_tree = true;
            aptr->hpos = -1;  // arc does not go into heap
            aptr->key = DBL_MAX;
            vptr = wptr;
            aptr = aptr->next_out;
        }
        sptr->first_child->left_sibl = vptr;
        vptr->right_sibl = sptr->first_child;


        // determine upper bound on lambda, can be used as 'infinity'
        // adds up all costs, and divide by smallest non-zero transit time
		// requires that there are no cycles with zero transit time!
		// Also, determine epsilon value on transit times, cost and cost/time ratios
		// as a constant fraction of the smallest observed values
        double total_cost_plus_one = 1.0L;
		double min_transit_time = DBL_MAX;
		for (aptr = &(gr->arcs[gr->n_arcs - 1L]); aptr >= gr->arcs; aptr--)
        {
			// add costs to total cost
			total_cost_plus_one += fabs(aptr->cost);
			// keep min of transit times
			if (aptr->transit_time > 0.0L) {
				if (aptr->transit_time < min_transit_time) {
					min_transit_time = aptr->transit_time;
				}
			}
		}
        infty = total_cost_plus_one / min_transit_time;

        // initial keys of non tree edges are equal to arc costs
        for (aptr = &(gr->arcs[gr->n_arcs - 1L]); aptr >= gr->arcs; aptr--)
        {
			if (aptr->transit_time > 0) {
				aptr->key = aptr->cost / aptr->transit_time;
			}
			else {
				aptr->key = infty;
			}
            aptr->in_tree = false;
        }



        // d-heap used for maintenance of vertex keys
        if (! ALLOC_HEAP(&h, gr->n_nodes))
            throw CException("Failed allocating heap");

        // compute initial vertex keys
        for (vptr = &(gr->nodes[gr->n_nodes - 1L]); vptr >= gr->nodes; vptr--)
        {
            min = DBL_MAX;
            vmin_aptr = NILA;
            aptr = vptr->first_arc_in;
            while (aptr != NILA)
            {
				if (! aptr->in_tree && aptr->key < min)
                {
					min = aptr->key;
                    vmin_aptr = aptr;
                }
                aptr = aptr->next_in;
            }
            vptr->vkey = vmin_aptr;
			if (vmin_aptr != NILA) {
				INSERT(&h, vmin_aptr);
			}
        }
        gr->vs->vkey = NILA;

        while (true)
        {
			min_aptr = GET_MIN(&h);
            ASSERT(min_aptr != NILA, "No element on heap!");

            *lambda = min_aptr->key;
            if (*lambda >= infty)
            {
				min_aptr = NILA;
                break; // input graph is acyclic in this case
            }

            uptr = min_aptr->tail;
            vptr = min_aptr->head;

            /* check if *vptr is an ancestor of *uptr in tree */

            foundCycle = false;
            par_aptr = uptr->parent_in;

            // MG: below is a fix, not in the original algorithm, since the original algorithm
            // does not seem to anticipate the possibility of self-edges in the graph.
            if (uptr == vptr)
            {
                // statement below was added to makje the calculation of the critical cycle work
                // correctly for critical self-edges.
                uptr->parent_in = min_aptr;
                break;
            }

            while (par_aptr != NILA)
            {
				if (par_aptr->tail == vptr)
                {
					foundCycle = true;
                    break;
                }
                else
                    par_aptr = par_aptr->tail->parent_in;
            }
            if (foundCycle) break;

            // it is not, remove edge (parent(v),v) from tree and make edge (u,v) a
            // tree edge instead
            par_aptr = vptr->parent_in;
            par_aptr->in_tree = false;
            min_aptr->in_tree = true;

            vptr->cost_t = uptr->cost_t + min_aptr->cost;
            vptr->transit_time_t = uptr->transit_time_t + min_aptr->transit_time;
            wptr = par_aptr->tail;

            // delete link (wptr,vptr) from tree
			if (vptr->right_sibl == vptr) {
				wptr->first_child = NILN;
			}
            else
            {
				vptr->right_sibl->left_sibl = vptr->left_sibl;
                vptr->left_sibl->right_sibl = vptr->right_sibl;
				if (wptr->first_child == vptr) {
					wptr->first_child = vptr->right_sibl;
				}
            }

            // insert link (uptr,vptr) into tree
            vptr->parent_in = min_aptr;
            if (uptr->first_child == NILN)
            {
				uptr->first_child = vptr;
                vptr->right_sibl = vptr;
                vptr->left_sibl = vptr;
            }
            else
            {
				vptr->right_sibl = uptr->first_child->right_sibl;
                uptr->first_child->right_sibl->left_sibl = vptr;
                vptr->left_sibl = uptr->first_child;
                uptr->first_child->right_sibl = vptr;
            }


            // subtree rooted at v has u as parent node now, update level and cost
            // entries of its nodes accordingly and produce list of nodes contained
            // in subtree
            upd_nodes = NILN;
            update_level = uptr->level + 1L;

            update_subtree(vptr);
            // now compute new keys of arcs into nodes that have acquired a new
            // shortest path, such arcs have head or tail in the subtree rooted at
            // "vptr", update vertex keys at the same time, nodes to be checked are
            // those contained in the subtree and the ones pointed to by arcs
            // emanating from nodes in the subtree
            vptr = upd_nodes;
            while (vptr != NILN)
            {
				if (vptr->vkey != NILA)
                    DELETE(&h, vptr->vkey);
                min = DBL_MAX;
                vmin_aptr = NILA;
                aptr = vptr->first_arc_in;
                while (aptr != NILA)
                {
					if (! aptr->in_tree)
                    {
						uptr = aptr->tail;
                        aptr->key = uptr->transit_time_t + aptr->transit_time
                                    > vptr->transit_time_t ?
                                    (double)(uptr->cost_t + aptr->cost - vptr->cost_t) /
                                    (double)(uptr->transit_time_t + aptr->transit_time
                                             - vptr->transit_time_t) : infty;

						if (aptr->key < min)
                        {
							min = aptr->key;
                            vmin_aptr = aptr;
                        }
                    }
                    aptr = aptr->next_in;
                }
				if (vmin_aptr != NILA) {
					INSERT(&h, vmin_aptr);
				}
                vptr->vkey = vmin_aptr;

                vptr = vptr->link;
            }

            min_aptr->key = DBL_MAX;

            // now update keys of arcs from nodes in subtree to nodes not contained
            // in subtree and update vertex keys for the latter if necessary
            vptr = upd_nodes;
            while (vptr != NILN)
            {
				aptr = vptr->first_arc_out;
                while (aptr != NILA)
                {
					if (! aptr->in_tree && ! aptr->head->in_list)
                    {
						wptr = aptr->head;
                        akey = vptr->transit_time_t + aptr->transit_time
                               > wptr->transit_time_t ?
                               (double)(vptr->cost_t + aptr->cost - wptr->cost_t)
                               / (double)(vptr->transit_time_t + aptr->transit_time
                                          - wptr->transit_time_t) : infty;

						if (akey < wptr->vkey->key)
                        {
							DELETE(&h, wptr->vkey);
                            aptr->key = akey;
                            INSERT(&h, aptr);
                            wptr->vkey = aptr;
                        }
						else {
							aptr->key = akey;
						}
                    }
                    aptr = aptr->next_out;
                }
                vptr = vptr->link;
            }

            vptr = upd_nodes;
            while (vptr != NILN)
            {
				vptr->in_list = false;
                vptr = vptr->link;
            }
        }

        DEALLOC_HEAP(&h);
        if (cycle != NULL && len != NULL)
        {
            *len = 0L;
            if (min_aptr != NILA)
            {
                cycle[(*len)++] = min_aptr;
                aptr = min_aptr->tail->parent_in;
                // MG: adapted the loop to work also for critical self-edges
                // keeping the original in comment for later reference in case of
                // problems or doubts...
                //    do
                //    {
                //        cycle[(*len)++] = aptr;
                //        aptr = aptr->tail->parent_in;
                //    }
                //    while (aptr->head != min_aptr->head);
                while (aptr->head != min_aptr->head)
                {
                    cycle[(*len)++] = aptr;
                    aptr = aptr->tail->parent_in;
                }
            }
        }
    }

	/**
	 * mmcycle_robust ()
	 *
	 * Determines maximum ratio cycle in a directed graph in a robust way
	 * G = (V, E, c), c: E->IR a "cost" function on the
	 * edges, alternatively called "length" or "weight".
	 * It has been observed that the mmcycle algorithm can behave non - deterministically
	 * across compilers, when the graph has multiple equally critical cycles
	 * Due to the order in which floating point calculations are scheduled, the
	 * algorithm's control flow may follow different paths and lead to different
	 * critical cycles as the output.
	 * This function mmcycle_robust implements a version of the algorithm that is
	 * deterministic across compilers. It may be a little slower.
	 * For more info about the algorithm, see the mmcycle function.
	 *
	 * TODO: see if the algoithms can be unified to remove duplicate code
	 **/

	void mmcycle_robust(graph* gr, double* lambda, arc** cycle, long* len)
	{
		const double MCR_EPSILON_RATIO = 1.0e-8L;
		double min, infty, akey;
		arc* aptr, * par_aptr, * vmin_aptr, * min_aptr;
		node* sptr, * uptr, * vptr, * wptr;
		bool foundCycle;
		d_heap h;

		// set up initial tree
		sptr = gr->vs;
		sptr->in_list = false;
		sptr->cost_t = 0.0L;
		sptr->transit_time_t = 0.0L;
		sptr->level = 0L;
		sptr->parent_in = NILA;
		sptr->left_sibl = sptr;
		sptr->right_sibl = sptr;
		aptr = sptr->first_arc_out;
		vptr = aptr->head;
		sptr->first_child = vptr;
		while (aptr != NILA)
		{
			wptr = aptr->head;
			wptr->cost_t = 0.0L;
			wptr->transit_time_t = 0.0L;
			wptr->level = 1L;
			wptr->parent_in = aptr;
			wptr->first_child = NILN;
			wptr->left_sibl = vptr;
			vptr->right_sibl = wptr;
			wptr->in_list = false;
			aptr->in_tree = true;
			aptr->hpos = -1;  // arc does not go into heap
			aptr->key = DBL_MAX;
			vptr = wptr;
			aptr = aptr->next_out;
		}
		sptr->first_child->left_sibl = vptr;
		vptr->right_sibl = sptr->first_child;


		// determine upper bound on lambda, can be used as 'infinity'
		// adds up all costs, and divide by smallest non-zero transit time
		// requires that there are no cycles with zero transit time!
		// Also, determine epsilon value on transit times, cost and cost/time ratios
		// as a constant fraction of the smallest observed values
		double total_cost_plus_one = 1.0L;
		double min_transit_time = DBL_MAX;
		double min_cost = DBL_MAX;
		for (aptr = &(gr->arcs[gr->n_arcs - 1L]); aptr >= gr->arcs; aptr--)
		{
			// add costs to total cost
			total_cost_plus_one += fabs(aptr->cost);
			// keep min of transit times
			if (aptr->transit_time > 0.0L) {
				if (aptr->transit_time < min_transit_time) {
					min_transit_time = aptr->transit_time;
				}
			}
			// keep min of costs
			if (aptr->cost > 0.0L) {
				if (aptr->cost < min_cost) {
					min_cost = aptr->cost;
				}
			}
		}
		infty = total_cost_plus_one / min_transit_time;
		double epsilon_transit_time = MCR_EPSILON_RATIO * min_transit_time;
		double epsilon_cost_time_ratio = MCR_EPSILON_RATIO * (min_cost / min_transit_time);

		// initial keys of non tree edges are equal to arc costs
		for (aptr = &(gr->arcs[gr->n_arcs - 1L]); aptr >= gr->arcs; aptr--)
		{
			if (aptr->transit_time > epsilon_transit_time) {
				aptr->key = aptr->cost / aptr->transit_time;
			}
			else {
				aptr->key = infty;
			}
			aptr->in_tree = false;
		}



		// d-heap used for maintenance of vertex keys
		if (!ALLOC_HEAP(&h, gr->n_nodes))
			throw CException("Failed allocating heap");

		// compute initial vertex keys
		for (vptr = &(gr->nodes[gr->n_nodes - 1L]); vptr >= gr->nodes; vptr--)
		{
			min = DBL_MAX;
			vmin_aptr = NILA;
			aptr = vptr->first_arc_in;
			while (aptr != NILA)
			{
				if (!aptr->in_tree && (min - aptr->key > epsilon_cost_time_ratio))
				{
					min = aptr->key;
					vmin_aptr = aptr;
				}
				aptr = aptr->next_in;
			}
			vptr->vkey = vmin_aptr;
			if (vmin_aptr != NILA) {
				INSERT_ROBUST(&h, vmin_aptr, epsilon_cost_time_ratio);
			}
		}
		gr->vs->vkey = NILA;

		while (true)
		{
			min_aptr = GET_MIN(&h);
			ASSERT(min_aptr != NILA, "No element on heap!");

			*lambda = min_aptr->key;
			if (*lambda >= infty)
			{
				min_aptr = NILA;
				break; // input graph is acyclic in this case
			}

			uptr = min_aptr->tail;
			vptr = min_aptr->head;

			/* check if *vptr is an ancestor of *uptr in tree */

			foundCycle = false;
			par_aptr = uptr->parent_in;

			// MG: below is a fix, not in the original algorithm, since the original algorithm
			// does not seem to anticipate the possibility of self-edges in the graph.
			if (uptr == vptr)
			{
				// statement below was added to makje the calculation of the critical cycle work
				// correctly for critical self-edges.
				uptr->parent_in = min_aptr;
				break;
			}

			while (par_aptr != NILA)
			{
				if (par_aptr->tail == vptr)
				{
					foundCycle = true;
					break;
				}
				else
					par_aptr = par_aptr->tail->parent_in;
			}
			if (foundCycle) break;

			// it is not, remove edge (parent(v),v) from tree and make edge (u,v) a
			// tree edge instead
			par_aptr = vptr->parent_in;
			par_aptr->in_tree = false;
			min_aptr->in_tree = true;

			vptr->cost_t = uptr->cost_t + min_aptr->cost;
			vptr->transit_time_t = uptr->transit_time_t + min_aptr->transit_time;
			wptr = par_aptr->tail;

			// delete link (wptr,vptr) from tree
			if (vptr->right_sibl == vptr) {
				wptr->first_child = NILN;
			}
			else
			{
				vptr->right_sibl->left_sibl = vptr->left_sibl;
				vptr->left_sibl->right_sibl = vptr->right_sibl;
				if (wptr->first_child == vptr) {
					wptr->first_child = vptr->right_sibl;
				}
			}

			// insert link (uptr,vptr) into tree
			vptr->parent_in = min_aptr;
			if (uptr->first_child == NILN)
			{
				uptr->first_child = vptr;
				vptr->right_sibl = vptr;
				vptr->left_sibl = vptr;
			}
			else
			{
				vptr->right_sibl = uptr->first_child->right_sibl;
				uptr->first_child->right_sibl->left_sibl = vptr;
				vptr->left_sibl = uptr->first_child;
				uptr->first_child->right_sibl = vptr;
			}


			// subtree rooted at v has u as parent node now, update level and cost
			// entries of its nodes accordingly and produce list of nodes contained
			// in subtree
			upd_nodes = NILN;
			update_level = uptr->level + 1L;

			update_subtree(vptr);
			// now compute new keys of arcs into nodes that have acquired a new
			// shortest path, such arcs have head or tail in the subtree rooted at
			// "vptr", update vertex keys at the same time, nodes to be checked are
			// those contained in the subtree and the ones pointed to by arcs
			// emanating from nodes in the subtree
			vptr = upd_nodes;
			while (vptr != NILN)
			{
				if (vptr->vkey != NILA)
					DELETE_ROBUST(&h, vptr->vkey, epsilon_cost_time_ratio);
				min = DBL_MAX;
				vmin_aptr = NILA;
				aptr = vptr->first_arc_in;
				while (aptr != NILA)
				{
					if (!aptr->in_tree)
					{
						uptr = aptr->tail;
						if (uptr->transit_time_t + aptr->transit_time - vptr->transit_time_t > epsilon_transit_time) {
							aptr->key = (uptr->cost_t + aptr->cost - vptr->cost_t) / (uptr->transit_time_t + aptr->transit_time - vptr->transit_time_t);
						}
						else {
							aptr->key = infty;
						}

						if (min - aptr->key > epsilon_cost_time_ratio)
						{
							min = aptr->key;
							vmin_aptr = aptr;
						}
					}
					aptr = aptr->next_in;
				}
				if (vmin_aptr != NILA) {
					INSERT_ROBUST(&h, vmin_aptr, epsilon_cost_time_ratio);
				}
				vptr->vkey = vmin_aptr;
				vptr = vptr->link;
			}

			min_aptr->key = DBL_MAX;

			// now update keys of arcs from nodes in subtree to nodes not contained
			// in subtree and update vertex keys for the latter if necessary
			vptr = upd_nodes;
			while (vptr != NILN)
			{
				aptr = vptr->first_arc_out;
				while (aptr != NILA)
				{
					if (!aptr->in_tree && !aptr->head->in_list)
					{
						wptr = aptr->head;
						if (vptr->transit_time_t + aptr->transit_time - wptr->transit_time_t > epsilon_transit_time) {
							akey = (vptr->cost_t + aptr->cost - wptr->cost_t) / (vptr->transit_time_t + aptr->transit_time - wptr->transit_time_t);
						}
						else {
							akey = infty;
						}
						if (wptr->vkey->key - akey > epsilon_cost_time_ratio)
						{
							DELETE_ROBUST(&h, wptr->vkey, epsilon_cost_time_ratio);
							aptr->key = akey;
							INSERT_ROBUST(&h, aptr, epsilon_cost_time_ratio);
							wptr->vkey = aptr;
						}
						else {
							aptr->key = akey;
						}
					}
					aptr = aptr->next_out;
				}
				vptr = vptr->link;
			}

			vptr = upd_nodes;
			while (vptr != NILN)
			{
				vptr->in_list = false;
				vptr = vptr->link;
			}
		}

		DEALLOC_HEAP(&h);
		if (cycle != NULL && len != NULL)
		{
			*len = 0L;
			if (min_aptr != NILA)
			{
				cycle[(*len)++] = min_aptr;
				aptr = min_aptr->tail->parent_in;
				// MG: adapted the loop to work also for critical self-edges
				// keeping the original in comment for later reference in case of
				// problems or doubts...
				//    do
				//    {
				//        cycle[(*len)++] = aptr;
				//        aptr = aptr->tail->parent_in;
				//    }
				//    while (aptr->head != min_aptr->head);
				while (aptr->head != min_aptr->head)
				{
					cycle[(*len)++] = aptr;
					aptr = aptr->tail->parent_in;
				}
			}
		}
	}




    /**
     * convertMCMgraphToYTOgraph ()
     * The function converts a weighted directed graph used in the MCM algorithms
     * to graph input for Young-Tarjan-Orlin's algorithm.
	 * It assumes that the id's of the nodes are 0 <= id < number of nodes
     */
    void convertMCMgraphToYTOgraph(MCMgraph *g,
                                   graph *gr, double(*costFunction)(MCMedge *e),
                                   double(*transit_timeFunction)(MCMedge *e))
    {
        node *x;
        arc *a;

        gr->n_nodes = g->nrVisibleNodes();
        gr->n_arcs = g->nrVisibleEdges();
		// allocate space for the nodes, plus one for the exta source node that will be added
        gr->nodes = (node *)malloc((gr->n_nodes + 1) * sizeof(node));
        gr->arcs = (arc *)malloc((gr->n_arcs + gr->n_nodes) * sizeof(arc));

        // create nodes
		// keep an index of node id's
		CLookupIntInt nodeIndex;
		uint ind = 0;
		for (MCMnodesCIter iter = g->getNodes().begin(); iter != g->getNodes().end(); iter++)
        {
            MCMnode *n = *iter;
			nodeIndex.put(n->id, ind);
			x = &((gr->nodes)[ind]);
            x->id = n->id + 1;
            x->first_arc_out = NULL;
            x->first_arc_in = NULL;

            // Next
            x++;
			ind++;
        }

        // create arcs
        a = gr->arcs;
        for (MCMedgesCIter iter = g->getEdges().begin(); iter != g->getEdges().end(); iter++)
        {
            MCMedge *e = *iter;
            MCMnode *u = e->src;
            MCMnode *v = e->dst;

            a->tail = &(gr->nodes[nodeIndex.get(u->id)]);
            a->head = &(gr->nodes[nodeIndex.get(v->id)]);
            a->cost = (*costFunction)(e);
            a->transit_time = (*transit_timeFunction)(e);
            a->next_out = a->tail->first_arc_out;
            a->tail->first_arc_out = a;
            a->next_in = a->head->first_arc_in;
            a->head->first_arc_in = a;
            a->mcmEdge = e;
            // Next
            a++;
        }

        // Create a source node which has an edge to all nodes
        gr->vs = &(gr->nodes[gr->n_nodes]);
        gr->vs->id = 0;
        gr->vs->first_arc_out = NULL;
        gr->vs->first_arc_in = NULL;
        for (int i = 0; i < gr->n_nodes; i++)
        {
            a->cost = 0;
            a->transit_time = 0.0;
            a->tail = gr->vs;
            a->head = &(gr->nodes[i]);
            a->next_out = gr->vs->first_arc_out;
            gr->vs->first_arc_out = a;
            a->next_in = a->head->first_arc_in;
            a->head->first_arc_in = a;

            // Next
            a++;
        }

#if 0
        // Print the MCM graph
        cerr << "#nodes: " << g->nodes.size() << endl;
        cerr << "#edges: " << g->edges.size() << endl;
        cerr << "edge: (u, v, w, d)" << endl;
        for (MCMedgesIter iter = g->edges.begin();
             iter != g->edges.end(); iter++)
        {
            MCMedge *e = *iter;

            if (!e->visible)
                continue;

            cerr << "(" << e->src->id;
            cerr << ", " << e->dst->id;
            cerr << ", " << e->w;
            cerr << ", " << e->d;
            cerr << ")" << endl;
        }
        cerr << endl;

        // Print the graph used in the YTO algorithm
        x = gr->nodes;
        for (int i = 0; i <= gr->n_nodes; i++)
        {
            cerr << "vertex: " << x->id << endl;

            for (arc *a = x->first_arc_in; a != NULL; a = a->next_in)
            {
                cerr << "   edge from: " << a->tail->id;
                cerr << " (weight: " << a->cost << ")" << endl;
            }
            for (arc *a = x->first_arc_out; a != NULL; a = a->next_out)
            {
                cerr << "   edge to:   " << a->head->id;
                cerr << " (weight: " << a->cost << ")" << endl;
            }

            x++;
        }
#endif
    }

    /**
     * constOne ()
     * The function returns the unit cost associated with an edge.
     */
    static double constOne(MCMedge *e)
    {
        return 1.0;
    }

    /**
     * getWeight ()
     * The function returns the weight associated with an edge.
     */
    double getWeight(MCMedge *e)
    {
        return e->w;
    }

    /**
     * getDelay ()
     * The function returns the delay associated with an edge.
     */
    double getDelay(MCMedge *e)
    {
        return e->d;
    }

    /**
     * maxCycleMeanAndCriticalCycleYoungTarjanOrlin ()
     * The function computes the maximum cycle mean of edge weight of
     * an MCMgraph using Young-Tarjan-Orlin's algorithm.
     * It returns both the MCM and a critical cycle
     * The critical cycle is only returned if cycle and len are not NULL. Then *cycle points
     * to an array of *MCMEdges of the critical cycle and *len indicates the length of the cycle.
     * *cycle is a freshly allocated array and it is the caller's obligation to deallocate it
     * in due time.
	 *
	 * Note that the following assumed are made about the MCMgraph
	 * 1. it is assumed that all nodes in the graph are 'visible'
	 * 2. it is assumed that the node have id's ranging from 0 up to the number of nodes.
	 * 3. it is assumed that cycles have a weight > 0 ! 
	 */
    CDouble maxCycleMeanAndCriticalCycleYoungTarjanOrlin(MCMgraph *mcmGraph, MCMedge *** cycle, uint *len)
    {
        double mincr;
        graph ytoGraph;
        arc **ytoCycle;
        long ytoCycLen;

        // Convert the graph to an input graph for the YTO algorithm
        convertMCMgraphToYTOgraph(mcmGraph, &ytoGraph, constOne, getWeight);

        if (cycle != NULL && len != NULL)
        {
            // allocate space for the critical cycle
            ytoCycle = (arc **) malloc(sizeof(arc *) * mcmGraph->nrVisibleEdges());

            // Find maximum cycle mean
            mmcycle_robust(&ytoGraph, &mincr, ytoCycle, &ytoCycLen);

            *len = ytoCycLen;

            *cycle = (MCMedge **) malloc(sizeof(MCMedge *) * (*len));

            for (uint i = 0; i < *len; i++)
            {
                (*cycle)[i] = ytoCycle[i]->mcmEdge;
            }

            free(ytoCycle);

        }
        else
        {
            // Find maximum cycle mean without cycle
            mmcycle_robust(&ytoGraph, &mincr, NULL, NULL);
        }

        // Cleanup
        free(ytoGraph.nodes);
        free(ytoGraph.arcs);

        return 1.0 / mincr;
    }

    /**
     * mcmYoungTarjanOrlin ()
     * The function computes the maximum cycle mean of edge weight per edge of
     * an MCMgraph using Young-Tarjan-Orlin's algorithm.
     */

    CDouble maxCycleMeanYoungTarjanOrlin(MCMgraph *mcmGraph)
    {
        return maxCycleMeanAndCriticalCycleYoungTarjanOrlin(mcmGraph, NULL, NULL);
    }

    /**
     * maxCycleRatioAndCriticalCycleYoungTarjanOrlin ()
     * The function computes the maximum cycle ratio of edge weight over delay of
     * an MCMgraph using Young-Tarjan-Orlin's algorithm.
     * It returns both the MCR and a critical cycle
     * Since MCM is in C-style code, let's do it the C-way.
     * The critical cycle is only returned if cycle and len are not NULL. Then *cycle points
     * to an array of MCMEdges of the critical cycle and *len indicates the length of the cycle.
     * *cycle is a freshly allocated array and it is the caller's obligation to deallocate it
     * in due time.
     */

    CDouble maxCycleRatioAndCriticalCycleYoungTarjanOrlin(MCMgraph *mcmGraph, MCMedge *** cycle, uint *len)
    {
        double mincr;
        graph ytoGraph;
        arc **ytoCycle;
        long ytoCycLen;

		// catch special case when the graph has no edges
		if (mcmGraph->nrVisibleEdges() == 0) {
            if ( len != NULL ){
                *len = 0;
            }
            if ( cycle != NULL ){
                *cycle = NULL;
            }
			return 0.0;
		}

        // Convert the graph to an input graph for the YTO algorithm
        convertMCMgraphToYTOgraph(mcmGraph, &ytoGraph, getDelay, getWeight);

        if (cycle != NULL && len != NULL)
        {
            // allocate space for the critical cycle
            ytoCycle = (arc **) malloc(sizeof(arc *) * mcmGraph->nrVisibleEdges());

            // Find maximum cycle ratio
            mmcycle_robust(&ytoGraph, &mincr, ytoCycle, &ytoCycLen);

            *len = ytoCycLen;

            *cycle = (MCMedge **) malloc(sizeof(MCMedge *) * (*len));

            // note that mmcycle returns the critical cycle following edges backwards
            // therefore reverse the order of the edges.
            for (uint i = 0; i < *len; i++)
            {
                (*cycle)[i] = ytoCycle[(*len) - 1 - i]->mcmEdge;
            }

            free(ytoCycle);

        }
        else
        {
            // Find maximum cycle ratio without cycle
            mmcycle_robust(&ytoGraph, &mincr, NULL, NULL);
        }

        // Cleanup
        free(ytoGraph.nodes);
        free(ytoGraph.arcs);

        return 1.0 / mincr;
    }

    /**
     * maxCycleRatioYoungTarjanOrlin ()
     * The function computes the maximum cycle ratio of edge weight over delay of
     * an MCMgraph using Young-Tarjan-Orlin's algorithm.
     */

    CDouble maxCycleRatioYoungTarjanOrlin(MCMgraph *mcmGraph)
    {
        return maxCycleRatioAndCriticalCycleYoungTarjanOrlin(mcmGraph, NULL, NULL);
    }

    /**
     * minCycleRatioAndCriticalCycleYoungTarjanOrlin ()
     * The function computes the minimum cycle ratio of edge weight over delay of
     * an MCMgraph using Young-Tarjan-Orlin's algorithm.
     * It returns both the MCR and a critical cycle
     * The critical cycle is only returned if cycle and len are not NULL. Then *cycle points
     * to an array of MCMEdges of the critical cycle and *len indicates the length of the cycle.
     * *cycle is a freshly allocated array and it is the caller's obligation to deallocate it
     * in due time.
     */

    CDouble minCycleRatioAndCriticalCycleYoungTarjanOrlin(MCMgraph *mcmGraph, MCMedge *** cycle, uint *len)
    {
        double mincr;
        graph ytoGraph;
        arc **ytoCycle;
        long ytoCycLen;

        // Convert the graph to an input graph for the YTO algorithm
        convertMCMgraphToYTOgraph(mcmGraph, &ytoGraph, getWeight, getDelay);

        if (cycle != NULL && len != NULL)
        {
            // allocate space for the critical cycle
            ytoCycle = (arc **) malloc(sizeof(arc *) * mcmGraph->nrVisibleEdges());

            // Find minimum cycle ratio
            mmcycle_robust(&ytoGraph, &mincr, ytoCycle, &ytoCycLen);

            *len = ytoCycLen;

            *cycle = (MCMedge **) malloc(sizeof(MCMedge *) * (*len));

            for (uint i = 0; i < *len; i++)
            {
                (*cycle)[i] = ytoCycle[i]->mcmEdge;
            }

            free(ytoCycle);

        }
        else
        {
            // Find minimum cycle ratio without cycle
            mmcycle_robust(&ytoGraph, &mincr, NULL, NULL);
        }

        // Cleanup
        free(ytoGraph.nodes);
        free(ytoGraph.arcs);

        return mincr;
    }

    /**
     * minCycleRatioYoungTarjanOrlin ()
     * The function computes the minimum cycle ratio of edge weight over delay of
     * an MCMgraph using Young-Tarjan-Orlin's algorithm.
     */

    CDouble minCycleRatioYoungTarjanOrlin(MCMgraph *mcmGraph)
    {
        return minCycleRatioAndCriticalCycleYoungTarjanOrlin(mcmGraph, NULL, NULL);
    }


}//namespace SDF

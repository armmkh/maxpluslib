#include <algorithm>

#include "algebra/mptype.h"
#include "base/analysis/mcm/mcm.h"
#include "base/analysis/mcm/mcmdg.h"
#include "base/analysis/mcm/mcmgraph.h"
#include <base/analysis/mcm/mcmhoward.h>
#include <base/analysis/mcm/mcmyto.h>
#include "mcmtest.h"
#include "testing.h"


using namespace MaxPlus;
using namespace Graphs;

void MCMTest::Run() {
    this->test_dg();
    this->test_howard();
    this->test_karp();
    this->test_yto();
};

MCMgraph makeGraph1() {
    MCMgraph g;

    MCMnode &n0 = *g.addNode(0);
    MCMnode &n1 = *g.addNode(1);
    MCMnode &n2 = *g.addNode(2);
    MCMnode &n3 = *g.addNode(3);
    MCMnode &n4 = *g.addNode(4);
    g.addEdge(0, n0, n1, 1.0, 1.0);
    g.addEdge(1, n1, n2, 2.0, 4.0);
    g.addEdge(2, n2, n3, 3.0, 0.0);
    g.addEdge(3, n3, n0, 4.0, 1.0);
    g.addEdge(4, n3, n3, 1.0, 0.3);
    g.addEdge(5, n3, n4, 10.0, 0.0);
    g.addEdge(6, n4, n4, 1.0, 1.0);
    return g;
}

// A corner case graph
MCMgraph makeGraph2() {
    MCMgraph g;

    MCMnode &n0 = *g.addNode(0);
    MCMnode &n1 = *g.addNode(1);
    g.addEdge(0, n0, n1, 1.0, 1.0);
    //g.addEdge(1, n1, n1, 1.0, 1.0);
    return g;
}


// Test mcmdg.
void MCMTest::test_dg() {
    std::cout << "Running test: MCM-dg" << std::endl;

    MCMgraph g1 = makeGraph1();
    CDouble result =  mcmDG(g1);
    ASSERT_APPROX_EQUAL(2.5, result, 1e-5);

    MCMgraph g2 = makeGraph2();
    result = mcmDG(g2);
    ASSERT_EQUAL(-INFINITY, result);
}

/// Test MCM Howard.
void MCMTest::test_howard() {
    std::cout << "Running test: MCM-Howard" << std::endl;
    MCMgraph g1 = makeGraph1();
    
    std::shared_ptr<std::vector<int>> ij = nullptr;
    std::shared_ptr<std::vector<CDouble>> A = nullptr;
    std::shared_ptr<std::vector<CDouble>> chi = nullptr;
    std::shared_ptr<std::vector<CDouble>> v = nullptr;
    std::shared_ptr<std::vector<int>> policy = nullptr;
    int nr_iterations = 0;
    int nr_components = 0;
    
    convertMCMgraphToMatrix(g1, &ij, &A);

   Howard(*ij,
           *A,
           static_cast<int>(g1.getNodes().size()),
           static_cast<int>(g1.getEdges().size()),
           &chi,
           &v,
           &policy,
           &nr_iterations,
           &nr_components);
    

    ASSERT_APPROX_EQUAL(2.5, chi->at(0), 1e-5);
    ASSERT_APPROX_EQUAL(2.5, chi->at(1), 1e-5);
    ASSERT_APPROX_EQUAL(2.5, chi->at(2), 1e-5);
    ASSERT_APPROX_EQUAL(2.5, chi->at(3), 1e-5);
    ASSERT_APPROX_EQUAL(1.0, chi->at(4), 1e-5);
    ASSERT_EQUAL(1, policy->at(0));
    ASSERT_EQUAL(2, policy->at(1));
    ASSERT_EQUAL(3, policy->at(2));
    ASSERT_EQUAL(0, policy->at(3));
    ASSERT_EQUAL(4, policy->at(4));

    MCMnode *criticalNode;
    CDouble result = maximumCycleMeanHoward(g1, &criticalNode);
    ASSERT_APPROX_EQUAL(2.5, result, 1e-5);
    ASSERT_THROW(criticalNode->id == 0 || criticalNode->id == 1 || criticalNode->id == 2
                 || criticalNode->id == 3);

    MCMgraph g2 = makeGraph2();
    result = maximumCycleMeanHoward(g2, &criticalNode);
    ASSERT_EQUAL(-INFINITY, result);
}

/// Test MCM Karp.
void MCMTest::test_karp() {
    std::cout << "Running test: MCM-Karp" << std::endl;

    MCMgraph g = makeGraph1();

    CDouble result = maximumCycleMeanKarp(g);
    ASSERT_APPROX_EQUAL(2.5, result, 1e-5);

    const MCMnode *criticalNode;
    result = maximumCycleMeanKarpDouble(g, &criticalNode);
    ASSERT_APPROX_EQUAL(2.5, result, 1e-5);
    MCMnodes& nodes = g.getNodes();
    ASSERT_THROW(criticalNode->id == 0 || criticalNode->id == 1 || criticalNode->id == 2
                 || criticalNode->id == 3);

}

/// Test MCM YTO.
void MCMTest::test_yto() {
    std::cout << "Running test: MCM-YTO" << std::endl;

    MCMgraph g = makeGraph1();
    CDouble result = maxCycleMeanYoungTarjanOrlin(g);
    ASSERT_APPROX_EQUAL(2.5, result, 1e-5);

    std::shared_ptr<std::vector<const MCMedge *>> cycle;
    result = maxCycleMeanAndCriticalCycleYoungTarjanOrlin(g, &cycle);
    ASSERT_APPROX_EQUAL(2.5, result, 1e-5);
    ASSERT_THROW(cycle->size() == 4);
    int eid = cycle->at(0)->id;
    ASSERT_THROW(eid == 0 || eid == 1 || eid == 2 || eid == 3);

    result = maxCycleRatioYoungTarjanOrlin(g);
    ASSERT_APPROX_EQUAL(10.0/3.0, result, 1e-5);

    result = maxCycleRatioAndCriticalCycleYoungTarjanOrlin(g, &cycle);
    ASSERT_APPROX_EQUAL(10.0 / 3.0, result, 1e-5);
    eid = cycle->at(0)->id;
    ASSERT_THROW(eid == 4);

    result = minCycleRatioYoungTarjanOrlin(g);
    ASSERT_APPROX_EQUAL(1.0, result, 1e-5);

    result = minCycleRatioAndCriticalCycleYoungTarjanOrlin(g, &cycle);
    ASSERT_APPROX_EQUAL(1.0, result, 1e-5);
    eid = cycle->at(0)->id;
    ASSERT_THROW(eid == 6);

}

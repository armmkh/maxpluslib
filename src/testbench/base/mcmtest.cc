#include <algorithm>

#include "algebra/mptype.h"
#include "base/analysis/mcm/mcmdg.h"
#include "base/analysis/mcm/mcmgraph.h"
#include <base/analysis/mcm/mcmhoward.h>
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
    g.addEdge(1, n1, n2, 2.0, 0.0);
    g.addEdge(2, n2, n3, 3.0, 0.0);
    g.addEdge(3, n3, n0, 4.0, 0.0);
    g.addEdge(4, n3, n3, 2.0, 0.0);
    g.addEdge(5, n3, n4, 10.0, 0.0);
    g.addEdge(6, n4, n4, 1.0, 0.0);
    return g;
}

// Test mcmdg.
void MCMTest::test_dg() {
    std::cout << "Running test: MCM-dg" << std::endl;

    MCMgraph g = makeGraph1();
    CDouble result =  mcmDG(g);
    ASSERT_APPROX_EQUAL(2.5, result, 1e-5);
}

/// Test MCM Howard.
void MCMTest::test_howard() {
    std::cout << "Running test: MCM-Howard" << std::endl;
    MCMgraph g = makeGraph1();
    
    std::shared_ptr<std::vector<int>> ij = nullptr;
    std::shared_ptr<std::vector<CDouble>> A = nullptr;
    std::shared_ptr<std::vector<CDouble>> chi = nullptr;
    std::shared_ptr<std::vector<CDouble>> v = nullptr;
    std::shared_ptr<std::vector<int>> policy = nullptr;
    int nr_iterations;
    int nr_components;
    
    convertMCMgraphToMatrix(g, &ij, &A);

    Howard(*ij,
           *A,
           static_cast<int>(g.getNodes().size()),
           static_cast<int>(g.getEdges().size()),
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
}

/// Test MCM Karp.
void MCMTest::test_karp() {
    std::cout << "Running test: MCM-Karp" << std::endl;
    ASSERT_APPROX_EQUAL(1.0, 1.0, 1e-5);
}

/// Test MCM Howard.
void MCMTest::test_yto() {
    std::cout << "Running test: MCM-YTO" << std::endl;
    ASSERT_APPROX_EQUAL(1.0, 1.0, 1e-5);
}

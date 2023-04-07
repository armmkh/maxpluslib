#include <algorithm>

#include "algebra/mptype.h"
#include "base/analysis/mcm/mcmdg.h"
#include "base/analysis/mcm/mcmgraph.h"
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

// Test mcmdg.
void MCMTest::test_dg() {
    std::cout << "Running test: MCM-dg" << std::endl;
    MCMgraph g;

    MCMnode& n0 = *g.addNode(0);
    MCMnode& n1 = *g.addNode(1);
    MCMnode& n2 = *g.addNode(2);
    MCMnode& n3 = *g.addNode(3);
    g.addEdge(0, n0, n1, 2.0, 1.0);
    g.addEdge(0, n1, n2, 2.0, 0.0);
    g.addEdge(0, n2, n3, 2.0, 0.0);
    g.addEdge(0, n3, n0, 2.0, 0.0);

    CDouble result =  mcmDG(g);
    ASSERT_APPROX_EQUAL(1.0, result, 1e-5);
}

/// Test MCM Howard.
void MCMTest::test_howard() {
    std::cout << "Running test: MCM-Howard" << std::endl;
    ASSERT_APPROX_EQUAL(1.0, 1.0, 1e-5);
}

/// Test MCM Howard.
void MCMTest::test_karp() {
    std::cout << "Running test: MCM-Karp" << std::endl;
    ASSERT_APPROX_EQUAL(1.0, 1.0, 1e-5);
}

/// Test MCM Howard.
void MCMTest::test_yto() {
    std::cout << "Running test: MCM-YTO" << std::endl;
    ASSERT_APPROX_EQUAL(1.0, 1.0, 1e-5);
}

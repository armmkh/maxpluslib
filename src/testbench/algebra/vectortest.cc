#include <algorithm>

#include "algebra/mpmatrix.h"
#include "testing.h"
#include "vectortest.h"

using namespace MaxPlus;

void VectorTest::Run() { this->test_Infinity(); };

// Test vector operations.
void VectorTest::test_Infinity(void) {
    std::cout << "Running test: Infinity" << std::endl;

    Vector *vec = new Vector(3);
    vec->put(0, MPTime(4.0));
    vec->put(1, MPTime(5.0));
    vec->put(2, MPTime(6.0));
    ASSERT_EQUAL(MPTime(4.0), vec->get(0));
    ASSERT_EQUAL(MPTime(5.0), vec->get(1));
    ASSERT_EQUAL(MPTime(6.0), vec->get(2));
}
#include <algorithm>

#include "algebra/mpmatrix.h"
#include "testing.h"
#include "vectortest.h"

using namespace MaxPlus;

void VectorTest::Run() { this->test_Infinity(); };

// Test vector operations.
void VectorTest::test_Infinity() {
    std::cout << "Running test: Infinity" << std::endl;

    Vector vec(3);
    vec.put(0, MPTime(4.0));
    vec.put(1, MPTime(5.0));
    vec.put(2, MPTime(6.0));
    ASSERT_EQUAL(4.0, static_cast<CDouble>(vec.get(0)));
    ASSERT_EQUAL(5.0, static_cast<CDouble>(vec.get(1)));
    ASSERT_EQUAL(6.0, static_cast<CDouble>(vec.get(2)));
}
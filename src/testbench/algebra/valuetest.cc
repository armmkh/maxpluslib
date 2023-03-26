#include <algorithm>

#include "algebra/mptype.h"
#include "valuetest.h"
#include "testing.h"

using namespace MaxPlus;


void ValueTest::Run() {
    this->test_mptime();
    this->test_Max();
    this->test_Min();
    this->test_BasicArithmetic();
};

// Test infinity operations.
void ValueTest::test_mptime() {
    std::cout << "Running test: mptime"<< std::endl;
    ASSERT_EQUAL(1.0, MPTime(1.0));
}

/// Test max operator.
void ValueTest::test_Max() {
    std::cout << "Running test: Max"<< std::endl;
    // Test normal double values.
    double a = 5.0;
    double b = 3.14;
    double max_ab = 5;

    ASSERT_EQUAL(MPTime(max_ab), MP_MAX(MPTime(a), MPTime(b)));
    ASSERT_EQUAL(MPTime(max_ab), MP_MAX(MPTime(b), MPTime(a)));
    ASSERT_EQUAL(MPTime(b), MP_MAX(MPTime(-a), MPTime(b)));
    ASSERT_EQUAL(MPTime(a), MP_MAX(MPTime(-b), MPTime(a)));

    // -Infinity case.
    ASSERT_THROW(MP_ISMINUSINFINITY(MP_MAX(MP_MINUSINFINITY, MP_MINUSINFINITY)));
    ASSERT_EQUAL(MPTime(3.14), MP_MAX(MPTime(3.14), MP_MINUSINFINITY));
}

/// Test min operator.
void ValueTest::test_Min() {
    std::cout << "Running test: Min"<< std::endl;
    // Test normal double values.
    ASSERT_EQUAL(MPTime(3.14), MP_MIN(MPTime(5.0), MPTime(3.14)));
    ASSERT_EQUAL(MPTime(3.14), MP_MIN(MPTime(3.14), MPTime(5.0)));
    ASSERT_EQUAL(MPTime(-5.0), MP_MIN(MPTime(-5.0), MPTime(3.14)));
    ASSERT_EQUAL(MPTime(-3.14), MP_MIN(MPTime(-3.14), MPTime(5.0)));

    // -Infinity case.
    ASSERT_THROW(MP_ISMINUSINFINITY(MP_MIN(MP_MINUSINFINITY, MP_MINUSINFINITY)));
    ASSERT_THROW(MP_ISMINUSINFINITY(MP_MIN(MPTime(-3.14), MP_MINUSINFINITY)));
    ASSERT_THROW(MP_ISMINUSINFINITY(MP_MIN(MPTime(3.14), MP_MINUSINFINITY)));
}

/// Test basic arithmetic operations.
void ValueTest::test_BasicArithmetic() {
    std::cout << "Running test: BasicArithmetic"<< std::endl;
    MPTime a(3.14);
    MPTime b(6.0);

    // Addition.
    ASSERT_EQUAL(MPTime(9.14), a + b);
    ASSERT_EQUAL(MPTime(9.14), b + a);
    // Subtraction.
    ASSERT_EQUAL(MPTime(2.86), b - a);
    ASSERT_EQUAL(MPTime(-2.86), a - b);
    // Multiplication.
    ASSERT_EQUAL(MPTime(18.84), a * b);
    // Multiplication.
    ASSERT_EQUAL(MPTime(18.84), b * a);
}


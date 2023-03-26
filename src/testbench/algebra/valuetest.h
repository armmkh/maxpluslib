#pragma once

#include <algorithm>

#include "algebra/mptype.h"
#include "testing.h"

using namespace MaxPlus;

class ValueTest : public ::testing::Test {

public:
    ValueTest() {}
    virtual void Run();
    virtual void SetUp() {};
    virtual void TearDown() {};

    void test_mptime();
    void test_Max();
    void test_Min();
    void test_BasicArithmetic();
};


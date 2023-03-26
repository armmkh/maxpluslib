#pragma once

#include <algorithm>

#include "testing.h"
#include "algebra/mpmatrix.h"



using namespace MaxPlus;

class MatrixTest : public ::testing::Test {

public:
    MatrixTest() {}

    virtual void SetUp() {};
    virtual void TearDown() {};
    int test_SetMPTimeInMatrix(void);
    int test_PasteMatrix(void);
    int test_SubMatrix(void);
    int test_Equality(void);
    int test_Addition(void);
    virtual void Run(void);
};



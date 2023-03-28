#pragma once

#include <algorithm>

#include "algebra/mpsparsematrix.h"
#include "testing.h"

using namespace MaxPlus;

class SparseMatrixTest : public ::testing::Test {

public:
    SparseMatrixTest() {}

    virtual void SetUp(){};
    virtual void TearDown(){};
    virtual void Run(void);

    int test_Vectors();
    int test_StarClosure();
    int test_EigenVectors();
    int test_GetPutMatrix();
    int test_Addition();
    int test_Multiplication();
};

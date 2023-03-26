#pragma once

#include <algorithm>
#include <base/basic_types.h>

#include "testing.h"

#include "graph/mpautomaton.h"

using namespace MaxPlus;

class MPAutomatonTest : public ::testing::Test {

public:
    MPAutomatonTest();

    virtual void SetUp();
    virtual void TearDown();
    virtual void Run();
    void testCreateFSM(void);
};

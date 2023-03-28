#pragma once

#include "game/mpgameautomaton.h"
#include "game/strategyvector.h"
#include <algorithm>

using namespace MaxPlus;

class StrategyVectorTest : public ::testing::Test {

public:
    StrategyVectorTest();

    virtual void Run();
    virtual void SetUp();
    virtual void TearDown();

    void testSimpleTest(void);
};

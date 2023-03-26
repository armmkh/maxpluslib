#pragma once

#include <algorithm>
#include "game/mpgameautomaton.h"
#include "game/strategyvector.h"
#include "game/policyiteration.h"

#include "testing.h"

using namespace MaxPlus;

class PolicyIterationTest : public ::testing::Test {

public:
    PolicyIterationTest();

    virtual void SetUp();
    virtual void Run();
    virtual void TearDown();
    void testPlayer1CycleTest(void);
    void testPlayer1CycleTest2(void);
    void testTwoPlayersTest(void);
    void testSimpleTest(void);
    void testInvalidInputGraphTest(void);
};




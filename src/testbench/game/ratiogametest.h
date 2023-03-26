#pragma once

#include <algorithm>

#include "game/ratiogame.h"

using namespace MaxPlus;

class RatioGameTest : public ::testing::Test {

public:
    virtual void Run();
    virtual void SetUp();
    virtual void TearDown();
    void testInitializeRandomStrategy(void);
};


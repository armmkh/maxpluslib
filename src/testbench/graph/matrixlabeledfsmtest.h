#pragma once

class MatrixLabeledFSMTest : public ::testing::Test {

public:
    MatrixLabeledFSMTest();
    virtual void SetUp();
    virtual void TearDown();
    virtual void Run();

    void testCreateFSM(void);
};

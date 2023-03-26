#include <algorithm>
#include "game/mpgameautomaton.h"

#include "testing.h"
#include "matrixlabeledfsmtest.h"

using namespace FSM;
using namespace MaxPlus;

MatrixLabeledFSMTest::MatrixLabeledFSMTest(){
}

void MatrixLabeledFSMTest::SetUp(){
}

void MatrixLabeledFSMTest::TearDown(){
}

void MatrixLabeledFSMTest::Run(){
    testCreateFSM();
}


void MatrixLabeledFSMTest::testCreateFSM(void) {

    // MG: This test appears to never have worked.
    // class MatrixLabeledFSM does not seem to exist ?

    // MatrixLabeledFSM *fsm = new MatrixLabeledFSM();

    // Matrix *m = new Matrix(1);
    // Value *v = new Value(4.0);
    // m->put(0,0,*v);

    // MLState *s0 = new MLState(0);
    // MLState *s1 = new MLState(1);
    // fsm->addState(s0);
    // fsm->addState(s1);
    // fsm->addEdge(s0,m,s1);
    // fsm->addEdge(s1,m,s0);

    // fsm->setInitialState(s0);

    // ASSERT_EQUAL(s0, fsm->getInitialState());
}
#include <algorithm>
#include <base/basic_types.h>

#include "mpautomatontest.h"

#include "graph/mpautomaton.h"

using namespace MaxPlus;

MPAutomatonTest::MPAutomatonTest(){
}

void MPAutomatonTest::SetUp(){
}

void MPAutomatonTest::TearDown(){
}

void MPAutomatonTest::Run(){
    testCreateFSM();
}

void MPAutomatonTest::testCreateFSM(void) {

    std::cout << "Running test: CreateFSM"<< std::endl;

    MaxPlusAutomatonWithRewards *mpa = new MaxPlusAutomatonWithRewards();
    // One FSM state, three tokens:
    CId fsm_s0 = 0;

    MPARState *s1 = new MPARState(makeMPAStateLabel(fsm_s0, 0));
    MPARState *s2 = new MPARState(makeMPAStateLabel(fsm_s0, 1));
    MPARState *s3 = new MPARState(makeMPAStateLabel(fsm_s0, 2));

    mpa->addState(s1);
    mpa->addState(s2);
    mpa->addState(s3);
    mpa->addEdge(s1, makeRewardEdgeLabel(MPTime(3.0), new CString("A"), 1.0), s2);
    mpa->addEdge(s1, makeRewardEdgeLabel(MPTime(3.0), new CString("A"), 1.0), s3);
    mpa->addEdge(s2, makeRewardEdgeLabel(MPTime(1.0), new CString("A"), 1.0), s1);
    mpa->addEdge(s3, makeRewardEdgeLabel(MPTime(7.0), new CString("A"), 1.0), s1);

    SetOfEdges<MPAStateLabel, MPAREdgeLabel> *es =
            (SetOfEdges<MPAStateLabel, MPAREdgeLabel> *) s1->getOutgoingEdges();
    typename SetOfEdges<MPAStateLabel, MPAREdgeLabel>::CIter i;
    for (i = es->begin(); i != es->end(); i++) {
        Edge<MPAStateLabel, MPAREdgeLabel> *e = (Edge<MPAStateLabel, MPAREdgeLabel> *) *i;

        ASSERT_THROW(e->getDestination() == s2 || e->getDestination() == s3);
    }

    mpa->setInitialState(s1);

    ASSERT_EQUAL_NOPRINT(s1, mpa->getInitialState());
}
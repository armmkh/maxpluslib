#include <algorithm>
#include <base/basic_types.h>
#include <memory>

#include "base/fsm/fsm.h"
#include "mpautomatontest.h"

#include "graph/mpautomaton.h"
#include "testing.h"

using namespace MaxPlus;

MPAutomatonTest::MPAutomatonTest() = default;

void MPAutomatonTest::SetUp() {}

void MPAutomatonTest::TearDown() {}

void MPAutomatonTest::Run() {
    testCreateFSM(); 
    testDeterminizeFSM();
}

void MPAutomatonTest::testCreateFSM() {

    std::cout << "Running test: CreateFSM" << std::endl;

    MaxPlusAutomatonWithRewards mpa;
    // One FSM state, three tokens:
    CId fsm_s0 = 0;

    MPARState *s1 = mpa.addState(makeMPAStateLabel(fsm_s0, 0));
    MPARState *s2 = mpa.addState(makeMPAStateLabel(fsm_s0, 1));
    MPARState *s3 = mpa.addState(makeMPAStateLabel(fsm_s0, 2));

    mpa.addEdge(*s1, makeRewardEdgeLabel(MPTime(3.0), CString("A"), 1.0), *s2);
    mpa.addEdge(*s1, makeRewardEdgeLabel(MPTime(3.0), CString("A"), 1.0), *s3);
    mpa.addEdge(*s2, makeRewardEdgeLabel(MPTime(1.0), CString("A"), 1.0), *s1);
    mpa.addEdge(*s3, makeRewardEdgeLabel(MPTime(7.0), CString("A"), 1.0), *s1);

    auto es = dynamic_cast<const FSM::Abstract::SetOfEdgeRefs &>(
            s1->getOutgoingEdges());
    for (const auto &i : es) {
        auto *e = dynamic_cast<Edge<MPAStateLabel, MPAREdgeLabel> *>(i);

        ASSERT_THROW((&(e->getDestination())) == s2 || (&(e->getDestination())) == s3);
    }

    mpa.setInitialState(*s1);

    ASSERT_EQUAL_NOPRINT(s1, &(mpa.getInitialState()));
}

void MPAutomatonTest::testDeterminizeFSM() {

    std::cout << "Running test: DeterminizeFSM" << std::endl;

    MaxPlusAutomatonWithRewards mpa;
    // One FSM state, three tokens:
    CId fsm_s0 = 0;

    MPARState *s1 = mpa.addState(makeMPAStateLabel(fsm_s0, 0));
    MPARState *s2 = mpa.addState(makeMPAStateLabel(fsm_s0, 1));
    MPARState *s3 = mpa.addState(makeMPAStateLabel(fsm_s0, 2));

    // s1 -- (3,A,1) -> s2 
    mpa.addEdge(*s1, makeRewardEdgeLabel(MPTime(3.0), CString("A"), 1.0), *s2);
    // s1 -- (3,A,1) -> s3 
    mpa.addEdge(*s1, makeRewardEdgeLabel(MPTime(3.0), CString("A"), 1.0), *s3);
    // s2 -- (1,B,1) -> s1
    mpa.addEdge(*s2, makeRewardEdgeLabel(MPTime(1.0), CString("B"), 1.0), *s1);
    // s3 -- (7,A,1) -> s1 
    mpa.addEdge(*s3, makeRewardEdgeLabel(MPTime(7.0), CString("A"), 1.0), *s1);


    auto es = dynamic_cast<const FSM::Abstract::SetOfEdgeRefs &>(
            s1->getOutgoingEdges());
    for (const auto &i : es) {
        auto *e = dynamic_cast<Edge<MPAStateLabel, MPAREdgeLabel> *>(i);

        ASSERT_THROW((&(e->getDestination())) == s2 || (&(e->getDestination())) == s3);
    }

    mpa.setInitialState(*s1);

    std::shared_ptr<MaxPlusAutomatonWithRewards> mpaDeterminized = std::dynamic_pointer_cast<MaxPlusAutomatonWithRewards>(mpa.determinizeEdgeLabels());

    ASSERT_THROW(mpaDeterminized->getStates().size() == 2);
    ASSERT_THROW(mpaDeterminized->getEdges().size()==3);

    // TODO: implement MCR on MaxPlusAutomatonWithRewards and MCM on MaxPlusAutomaton
    // CDouble mcr = mpaDeterminized->calculateMCR();
    // std::cout << "MCR: " << mcr << std::endl;

}
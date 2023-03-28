#include "game/strategyvector.h"
#include "game/mpgameautomaton.h"
#include <algorithm>


#include "strategyvectortest.h"
#include "testing.h"


using namespace MaxPlus;

StrategyVectorTest::StrategyVectorTest(void) {}

void StrategyVectorTest::Run(void) { testSimpleTest(); }

void StrategyVectorTest::SetUp(void) {}
void StrategyVectorTest::TearDown(void) {}

void StrategyVectorTest::testSimpleTest(void) {

    std::cout << "Running test: SimpleTest" << std::endl;

    MaxPlusGameAutomatonWithRewards *mpa = new MaxPlusGameAutomatonWithRewards();

    // One FSM state, three tokens:
    CId fsm_s0 = 0;

    MPARState *s1 = new MPARState(makeMPAStateLabel(fsm_s0, 0));
    MPARState *s2 = new MPARState(makeMPAStateLabel(fsm_s0, 1));
    MPARState *s3 = new MPARState(makeMPAStateLabel(fsm_s0, 2));

    mpa->addState(s1);
    mpa->addState(s2);
    mpa->addState(s3);
    mpa->addEdge(s1, makeRewardEdgeLabel(3.0, new CString("A"), 1.0), s2);
    mpa->addEdge(s1, makeRewardEdgeLabel(3.0, new CString("A"), 1.0), s3);
    mpa->addEdge(s2, makeRewardEdgeLabel(1.0, new CString("A"), 1.0), s1);
    mpa->addEdge(s3, makeRewardEdgeLabel(7.0, new CString("A"), 1.0), s1);
    mpa->addV0(s1);
    mpa->addV1(s2);
    mpa->addV1(s3);

    // Get a strategy vector on the given game.
    StrategyVector<MPAStateLabel, MPAREdgeLabel> *vec =
            new StrategyVector<MPAStateLabel, MPAREdgeLabel>();
    vec->initializeRandomStrategy(mpa);

    MPARState *result = vec->getSuccessor(s1);
    ASSERT_THROW(result == s2 || result == s3);
}

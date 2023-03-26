#include <algorithm>
#include "game/mpgameautomaton.h"
#include "game/strategyvector.h"
#include "game/policyiteration.h"

#include "testing.h"
#include "policyiterationtest.h"

using namespace MaxPlus;

PolicyIterationTest::PolicyIterationTest(void) {
};

void PolicyIterationTest::Run(void) {
    testPlayer1CycleTest();
    testPlayer1CycleTest2();
    testTwoPlayersTest();
    testSimpleTest();
    testInvalidInputGraphTest();
};

void PolicyIterationTest::SetUp(void) {
};

void PolicyIterationTest::TearDown(void) {
};

void PolicyIterationTest::testPlayer1CycleTest(void) {

    std::cout << "Running test: Player1CycleTest"<< std::endl;

    MaxPlusGameAutomatonWithRewards *mpa = new MaxPlusGameAutomatonWithRewards();

    // One FSM state, four tokens:
    CId fsm_s0 = 0;

    MPARState *s1 = new MPARState(make_mpastatelabel(fsm_s0, 0));
    MPARState *s2 = new MPARState(make_mpastatelabel(fsm_s0, 1));
    MPARState *s3 = new MPARState(make_mpastatelabel(fsm_s0, 2));
    MPARState *s4 = new MPARState(make_mpastatelabel(fsm_s0, 3));

    mpa->addV1(s1);
    mpa->addV1(s2);
    mpa->addV1(s3);
    mpa->addV1(s4);
    mpa->addState(s1);
    mpa->addState(s2);
    mpa->addState(s3);
    mpa->addState(s4);

    mpa->addEdge(s1, make_rewedgelabel(1.0, new CString("A"), 1.0), s2);
    mpa->addEdge(s2, make_rewedgelabel(1.0, new CString("A"), 2.0), s3);
    mpa->addEdge(s3, make_rewedgelabel(1.0, new CString("A"), 3.0), s4);
    mpa->addEdge(s4, make_rewedgelabel(1.0, new CString("A"), 4.0), s2);

    PolicyIteration<MPAStateLabel,MPAREdgeLabel> pi;
    PolicyIteration<MPAStateLabel,MPAREdgeLabel>::PolicyIterationResult result = pi.solve(mpa);

    double expected_ratio = 9.0/3.0;
    const std::map<State<MPAStateLabel, MPAREdgeLabel> *, double> & values = result.values;

    ASSERT_DOUBLE_EQUAL(expected_ratio, values.at(s1));
    ASSERT_DOUBLE_EQUAL(expected_ratio, values.at(s2));
    ASSERT_DOUBLE_EQUAL(expected_ratio, values.at(s3));
    ASSERT_DOUBLE_EQUAL(expected_ratio, values.at(s4));
    
}



void PolicyIterationTest::testPlayer1CycleTest2(void) {

    std::cout << "Running test: Player1CycleTest2"<< std::endl;

    MaxPlusGameAutomatonWithRewards *mpa = new MaxPlusGameAutomatonWithRewards();

    // One FSM state, four tokens:
    CId fsm_s0 = 0;

    MPARState *s1 = new MPARState(make_mpastatelabel(fsm_s0, 0));
    MPARState *s2 = new MPARState(make_mpastatelabel(fsm_s0, 1));
    MPARState *s3 = new MPARState(make_mpastatelabel(fsm_s0, 2));
    MPARState *s4 = new MPARState(make_mpastatelabel(fsm_s0, 3));

    mpa->addV1(s1);
    mpa->addV1(s2);
    mpa->addV1(s3);
    mpa->addV1(s4);
    mpa->addState(s1);
    mpa->addState(s2);
    mpa->addState(s3);
    mpa->addState(s4);

    mpa->addEdge(s1, make_rewedgelabel(3.0, new CString("A"), 1.0), s2);
    mpa->addEdge(s2, make_rewedgelabel(3.0, new CString("A"), 2.0), s3);
    mpa->addEdge(s3, make_rewedgelabel(1.0, new CString("A"), 3.0), s4);
    mpa->addEdge(s4, make_rewedgelabel(7.0, new CString("A"), 4.0), s2);

    PolicyIteration<MPAStateLabel,MPAREdgeLabel> *pi = new PolicyIteration<MPAStateLabel,MPAREdgeLabel>();
    PolicyIteration<MPAStateLabel,MPAREdgeLabel>::PolicyIterationResult result = pi->solve(mpa);

    double expected_ratio = 9.0/11.0;
    std::map<State<MPAStateLabel, MPAREdgeLabel> *, double> values = result.values;

    ASSERT_DOUBLE_EQUAL(expected_ratio, values[s1]);
    ASSERT_DOUBLE_EQUAL(expected_ratio, values[s2]);
    ASSERT_DOUBLE_EQUAL(expected_ratio, values[s3]);
    ASSERT_DOUBLE_EQUAL(expected_ratio, values[s4]);
}


void PolicyIterationTest::testTwoPlayersTest(void) {

    std::cout << "Running test: TwoPlayersTest"<< std::endl;

    MaxPlusGameAutomatonWithRewards *mpa = new MaxPlusGameAutomatonWithRewards();

    // One FSM state, four tokens:
    CId fsm_s0 = 0;

    MPARState *s1 = new MPARState(make_mpastatelabel(fsm_s0, 0));
    MPARState *s2 = new MPARState(make_mpastatelabel(fsm_s0, 1));
    MPARState *s3 = new MPARState(make_mpastatelabel(fsm_s0, 2));
    MPARState *s4 = new MPARState(make_mpastatelabel(fsm_s0, 3));
    MPARState *s5 = new MPARState(make_mpastatelabel(fsm_s0, 4));

    mpa->addV0(s1);
    mpa->addV0(s3);
    mpa->addV0(s5);
    mpa->addV1(s2);
    mpa->addV1(s4);
    mpa->addState(s1);
    mpa->addState(s2);
    mpa->addState(s3);
    mpa->addState(s4);
    mpa->addState(s5);

    mpa->addEdge(s1, make_rewedgelabel(1.0, new CString("A"), 0.0), s4);
    mpa->addEdge(s1, make_rewedgelabel(1.0, new CString("A"), 0.0), s2);
    mpa->addEdge(s4, make_rewedgelabel(1.0, new CString("A"), 0.0), s5);
    mpa->addEdge(s5, make_rewedgelabel(2.0, new CString("A"), 1.0), s4);
    mpa->addEdge(s2, make_rewedgelabel(0.0, new CString("A"), 3.0), s3);
    mpa->addEdge(s3, make_rewedgelabel(2.0, new CString("A"), 4.0), s2);

    PolicyIteration<MPAStateLabel,MPAREdgeLabel> *pi = new PolicyIteration<MPAStateLabel,MPAREdgeLabel>();
    PolicyIteration<MPAStateLabel,MPAREdgeLabel>::PolicyIterationResult result = pi->solve(mpa);

    double val1 = 7.0/2.0;
    double val2 = 1.0/3.0;
    std::map<State<MPAStateLabel, MPAREdgeLabel> *, double> values = result.values;

    ASSERT_DOUBLE_EQUAL(val1, values[s1]);
    ASSERT_DOUBLE_EQUAL(val1, values[s2]);
    ASSERT_DOUBLE_EQUAL(val1, values[s3]);
    ASSERT_DOUBLE_EQUAL(val2, values[s4]);
    ASSERT_DOUBLE_EQUAL(val2, values[s5]);
}


void PolicyIterationTest::testSimpleTest(void) {

    std::cout << "Running test: SimpleTest"<< std::endl;

    MaxPlusGameAutomatonWithRewards *mpa = new MaxPlusGameAutomatonWithRewards();

    // One FSM state, three tokens:
    CId fsm_s0 = 0;

    MPARState *s1 = new MPARState(make_mpastatelabel(fsm_s0, 0));
    MPARState *s2 = new MPARState(make_mpastatelabel(fsm_s0, 1));
    MPARState *s3 = new MPARState(make_mpastatelabel(fsm_s0, 2));

    mpa->addState(s1);
    mpa->addState(s2);
    mpa->addState(s3);
    mpa->addEdge(s1, make_rewedgelabel(3.0, new CString("A"), 1.0), s2);
    mpa->addEdge(s1, make_rewedgelabel(3.0, new CString("A"), 1.0), s3);
    mpa->addEdge(s2, make_rewedgelabel(1.0, new CString("A"), 1.0), s1);
    mpa->addEdge(s3, make_rewedgelabel(7.0, new CString("A"), 1.0), s1);
    mpa->addV0(s1);
    mpa->addV1(s2);
    mpa->addV1(s3);

    PolicyIteration<MPAStateLabel,MPAREdgeLabel> *pi = new PolicyIteration<MPAStateLabel,MPAREdgeLabel>();
    PolicyIteration<MPAStateLabel,MPAREdgeLabel>::PolicyIterationResult result = pi->solve(mpa);
}


void PolicyIterationTest::testInvalidInputGraphTest(void) {

    std::cout << "Running test: InvalidInputGraphTest"<< std::endl;

    MaxPlusGameAutomatonWithRewards *mpa = new MaxPlusGameAutomatonWithRewards();

    // One FSM state, three tokens:
    CId fsm_s0 = 0;

    MPARState *s1 = new MPARState(make_mpastatelabel(fsm_s0, 0));
    MPARState *s2 = new MPARState(make_mpastatelabel(fsm_s0, 1));

    mpa->addState(s1);
    mpa->addState(s2);
    mpa->addEdge(s1, make_rewedgelabel(3.0, new CString("A"), 1.0), s2);
    mpa->addV0(s1);
    mpa->addV1(s2);

    try {
        PolicyIteration<MPAStateLabel,MPAREdgeLabel> *pi = new PolicyIteration<MPAStateLabel,MPAREdgeLabel>();
        PolicyIteration<MPAStateLabel,MPAREdgeLabel>::PolicyIterationResult result = pi->solve(mpa);
    } catch (const std::runtime_error& error) {
        // test passes
        //SUCCEED();
    }

}

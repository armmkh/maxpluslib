/*
 *  Eindhoven University of Technology
 *  Eindhoven, The Netherlands
 *  Dept. of Electrical Engineering
 *  Electronics Systems Group
 *  Model Based Design Lab (https://computationalmodeling.info/)
 *
 *  Name            :   fsm.cc
 *
 *  Author          :   Marc Geilen (m.c.w.geilen@tue.nl)
 *
 *  Date            :   March 23, 2009
 *
 *  Function        :   finite state machines
 *
 *  History         :
 *      23-03-09    :   Initial version.
 *
 *
 *  Copyright 2023 Eindhoven University of Technology
 *
 *  Permission is hereby granted, free of charge, to any person obtaining
 *  a copy of this software and associated documentation files (the “Software”),
 *  to deal in the Software without restriction, including without limitation
 *  the rights to use, copy, modify, merge, publish, distribute, sublicense,
 *  and/or sell copies of the Software, and to permit persons to whom the
 *  Software is furnished to do so, subject to the following conditions:
 *
 *  The above copyright notice and this permission notice shall be included
 *  in all copies or substantial portions of the Software.
 *
 *  THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 *  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 *  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 *  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 *  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 *  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 *  SOFTWARE.
 */

#include "base/fsm/fsm.h"

namespace FSM {

CId FSM::Abstract::WithUniqueID::nextID = 0;

/*removes states of the fsm with no outgoing edges.*/
void EdgeLabeledScenarioFSM::removeDanglingStates() {

    // temporary sets of states and edges
    ELSSetOfEdges edgesToBeRemoved;
    ELSSetOfStates statesToBeRemoved;

    ELSSetOfStates *elsStates = this->getStates();
    ELSSetOfEdges *elsEdges = this->getEdges();

    /*go through all edges and find all edges that end in
    dangling states. Also store dangling states.*/
    for (ELSSetOfEdges::iterator it = elsEdges->begin(); it != elsEdges->end(); it++) {

        ELSEdge *e = (ELSEdge *)*it;
        ELSState *s = (ELSState *)e->getDestination();
        ELSSetOfEdges *oEdges = (ELSSetOfEdges *)s->getOutgoingEdges();
        if (oEdges->size() == 0) {
            edgesToBeRemoved.insert(e);
            statesToBeRemoved.insert(s);
        }
    }

    while (edgesToBeRemoved.size() != 0) {

        // remove dangling states
        for (ELSSetOfStates::iterator sr_it = statesToBeRemoved.begin();
             sr_it != statesToBeRemoved.end();
             sr_it++) {
            elsStates->erase(*sr_it);
            delete *sr_it;
        }

        // remove edges ending in dangling states
        // remove edges ending in dangling states from the outgoing edges of their source states
        for (ELSSetOfEdges::iterator er_it = edgesToBeRemoved.begin();
             er_it != edgesToBeRemoved.end();
             er_it++) {
            ELSEdge *e = (ELSEdge *)*er_it;
            elsEdges->erase(e);
            ELSState *s = (ELSState *)e->getSource();
            s->removeOutgoingEdge(e);
            delete e;
        }

        // empty the temporary sets
        edgesToBeRemoved.clear();
        statesToBeRemoved.clear();

        elsStates = this->getStates();
        elsEdges = this->getEdges();

        /*go through all edges and find all edges that end in
        dangling states. Also store dangling states.*/
        for (ELSSetOfEdges::iterator it = elsEdges->begin(); it != elsEdges->end(); it++) {

            ELSEdge *e = (ELSEdge *)*it;
            ELSState *s = (ELSState *)e->getDestination();
            ELSSetOfEdges *oEdges = (ELSSetOfEdges *)s->getOutgoingEdges();
            if (oEdges->size() == 0) {
                edgesToBeRemoved.insert(e);
                statesToBeRemoved.insert(s);
            }
        }
    }
}

namespace StateStringLabeled {

void FiniteStateMachine::addStateLabeled(const CString &sl) {
    StateStringLabeled::State *s = new StateStringLabeled::State(sl);
    this->addState(s);
}

void FiniteStateMachine::addEdgeLabeled(const CString &src, const CString &dst) {
    Labeled::State<CString, char> *s_src = this->getStateLabeled(src);
    Labeled::State<CString, char> *s_dst = this->getStateLabeled(dst);
    Labeled::FiniteStateMachine<CString, char>::addEdge(s_src, 'X', s_dst);
}

SetOfStates *FiniteStateMachine::reachableStates(void) {
    return (SetOfStates *)Labeled::FiniteStateMachine<CString, char>::reachableStates();
}

void FiniteStateMachine::setInitialStateLabeled(const CString &sl) {
    this->setInitialState(this->getStateLabeled(sl));
}

} // namespace StateStringLabeled

namespace Product {

State *FiniteStateMachine::getInitialState() {
    return new State(this->fsm_a->getInitialState(), this->fsm_b->getInitialState(), this);
}

Abstract::SetOfEdges *State::getOutgoingEdges() {
    if (!outgoingEdgesDone) {
        // compute outgoing edges
        const Abstract::SetOfEdges *oea = this->sa->getOutgoingEdges();
        const Abstract::SetOfEdges *oeb = this->sb->getOutgoingEdges();
        Abstract::SetOfEdges::const_iterator i = oea->begin();
        Abstract::SetOfEdges::const_iterator j = oeb->begin();
        while (i != oea->end()) {
            while (j != oeb->end()) {
                if (this->fsm->matchEdges(*i, *j)) {
                    Abstract::Edge *e = this->fsm->ensureEdge(*i, *j);
                    this->outgoingEdges->insert(e);
                }
                j++;
            }
            i++;
        }
        outgoingEdgesDone = true;
    }
    return this->outgoingEdges;
}
} // namespace Product

} // namespace FSM
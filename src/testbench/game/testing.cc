#include "policyiterationtest.h"
#include "ratiogametest.h"
#include "strategyvectortest.h"


int main(void) {
    
    PolicyIterationTest T1;
    T1.Run();

    RatioGameTest T2;
    T2.Run();

    StrategyVectorTest T3;
    T3.Run();

    return 0;
}
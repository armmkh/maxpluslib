#include "mpautomatontest.h"
#include "matrixlabeledfsmtest.h"


int main(void) {
    
    std::cout << "Running Tests Maxplus - Graph"<< std::endl;

    MatrixLabeledFSMTest T1;
    T1.Run();

    MPAutomatonTest T2;
    T2.Run();



    return 0;
}
#include "matrixtest.h"
#include "sparsematrixtest.h"
#include "vectortest.h"
#include "valuetest.h"



int main(void) {
    
    std::cout << "Running Tests Maxplus - Algebra"<< std::endl;

    ValueTest T1;
    T1.Run();

    VectorTest T2;
    T2.Run();

    MatrixTest T3;
    T3.Run();

    SparseMatrixTest T4;
    T4.Run();

    return 0;
}
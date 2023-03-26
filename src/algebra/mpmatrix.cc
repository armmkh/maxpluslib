/*
 *  Eindhoven University of Technology
 *  Eindhoven, The Netherlands
 *  Dept. of Electrical Engineering
 *  Electronics Systems Group
 *  Model Based Design Lab (https://computationalmodeling.info/)
 *
 *  Name            :   mpmatrix.cc
 *
 *  Author          :   Marc Geilen (m.c.w.geilen@tue.nl)
 *
 *  Date            :   March 23, 2009
 *
 *  Function        :   MaxPlus vectors and matrices
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

#include "algebra/mpmatrix.h"
#include "base/exception/exception.h"
#include "base/analysis/mcm/mcmgraph.h"
#include "base/analysis/mcm/mcmyto.h"
#include <stdlib.h>
#include <math.h>
#include "algebra/mpmatrix.h"
#include "base/exception/exception.h"

using namespace Graphs;

namespace MaxPlus
{

    /**
     * Construct a max-plus vector of size
     */
    Vector::Vector(unsigned int size, MPTime value) {
        this->table.resize(size);
        for (unsigned int k = 0; k < size; k++) {
            this->table[k] = value;
        }
    }

    /**
     * Construct a max-plus vector from an std:vector
     */
    Vector::Vector(std::vector <MPTime> *v) {
        this->table.resize(v->size());
        for (unsigned int i = 0; i < v->size(); i++) {
            this->table[i] = v->at(i);
        }
    }

    /**
     * copy constructor
     */
    Vector::Vector(const Vector &other) {
        this->table.resize(other.getSize());
        for (unsigned int row = 0; row < this->getSize(); row++) {
            this->table[row] = other.table[row];
        }
    }

    /**
     * vector assignment
     */
    Vector &Vector::operator=(const Vector &other) {
        if (this->getSize() != other.getSize()) {
            throw CException("Vectors of different size in"
                                     "Vector::operator=");
        }
        for (unsigned int row = 0; row < this->getSize(); row++) {
            this->table[row] = other.table[row];
        }
        return *this;
    }

    /**
     * vector negate
     */
    void Vector::negate() {
        for (unsigned int row = 0; row < this->getSize(); row++) {
            if (this->get(row) == MP_MINUSINFINITY) {
                throw CException("Cannot negate vectors with MP_MINUSINFINITY elements in"
                                         "Vector::negate");
            } else {
                this->put(row, -this->get(row));
            }
        }
    }


    /**
     * calculate vector norm
     */
    MPTime Vector::norm() {
        MPTime maxEl = MP_MINUSINFINITY;
        for (unsigned int row = 0; row < this->getSize(); row++) {
            maxEl = MP_MAX(maxEl, this->get(row));
        }
        return maxEl;
    }

    /**
     * normalize vector
     */
    MPTime Vector::normalize() {
        double maxEl = this->norm();

        if (maxEl == MP_MINUSINFINITY) {
            throw CException("Cannot normalize vector with norm MP_MINUSINFINITY"
                                     "Vector::normalize");
        } else {
            for (unsigned int row = 0; row < this->getSize(); row++) {
                MPTime x_i = this->get(row); // MPTime handles -INF correctly
                x_i = x_i - maxEl;          // overloaded using MP_PLUS
                this->put(row, x_i);
            }
        }
        return maxEl;
    }

    /**
     * add scalar to vector
     */
    Vector *Vector::add(MPTime increase) const {
        unsigned int M = this->getSize();
        Vector *result = new Vector(M);
        this->add(increase, result);
        return result;
    }

    /**
     * add scalar to vector
     */
    void Vector::add(MPTime increase, Vector *result) const {
        unsigned int M = this->getSize();
        assert(result->getSize() == M);

        for (unsigned int pos = 0; pos < M; pos++) {
            result->put(pos, this->get(pos) + increase); // uses MP_PLUS()
        }
    }

    /**
     * max of vectors
     */
    void Vector::maximum(const Vector *vecB, Vector *result) const {
        unsigned int M = this->getSize();
        assert(vecB->getSize() == M);
        assert(result->getSize() == M);

        for (unsigned int pos = 0; pos < M; pos++) {
            result->put(pos, MP_MAX(this->get(pos), vecB->get(pos)));
        }
    }


    /**
     * Destructor of max-plus vector
     */
    Vector::~Vector() {
    }

    /**
     * Put an entry into the vector. Grows vector if necessary
     */
    void Vector::put(unsigned int row, MPTime value) {
        if (this->table.size() <= row)
            this->table.resize(row + 1, 0);
        this->table[row] = value;
    }

    /**
     * String representation of vector
     */
    void Vector::toString(CString &outString, double scale) const {
        outString = "";
        for (unsigned int row = 0; row < this->getSize(); row++) {
            outString += timeToString(scale * this->get(row)) + " ";
        }
    }

    /**
     * Vector::add
     * add vectors
     */
    void Vector::add(const Vector *vecB, Vector *res) const {
        assert(this->getSize() == vecB->getSize());
        assert(this->getSize() == res->getSize());
        for (unsigned int row = 0; row < this->getSize(); row++) {
            res->put(row, this->get(row) + vecB->get(row));
        }
    }

    /**
     * add vectors
     */
    Vector *Vector::add(const Vector *vecB) const {
        Vector *res = new Vector(this->getSize());
        this->add(vecB, res);
        return res;
    }

    /**
     * Get minimal finite element
     * returns the smallest amoung the finite elements in the vector or
     * MP_MINUSINFINITY if no finite elements exist
     * itsPosition returns the index of the (a) smallest finite element is set
     * to a pointer to unsigned int, otherwise set or defaults to NULL
     */
    MPTime Vector::minimalFiniteElement(unsigned int *itsPosition_Ptr) const {
        unsigned int itsPosition_tmp;
        unsigned int *itsPosition;
        if (itsPosition_Ptr != NULL)
            itsPosition = itsPosition_Ptr;
        else
            itsPosition = &itsPosition_tmp;
        *itsPosition = this->getSize() + 1; // arbitrary value, invalid

        MPTime minEl = MP_MINUSINFINITY;
        bool minEl_initialized = false;
        for (unsigned int row = 0; row < this->getSize(); row++) {
            MPTime val = this->get(row);
            if (MaxPlus::MP_ISMINUSINFINITY(val)) continue;
            if (!minEl_initialized || val < minEl) {
                minEl = val;
                minEl_initialized = true;
                *itsPosition = row;
            }
        }
        return minEl;
    }


    /**
     * Compare vectors up to MP_EPSILON
     */
    bool Vector::compare(const Vector &v) {
        if (this->getSize() != v.getSize())
            return false;
        for (unsigned int k = 0; k < this->getSize(); k++) {
            if (fabs(this->get(k) - v.get(k)) > MP_EPSILON)
                return false;
        }
        return true;
    }

    Matrix::Matrix(unsigned int nrows, unsigned int ncols, MatrixFill fill) {
        this->szRows = nrows;
        this->szCols = ncols;
        init(fill);
    }


    /**
     * initialize matrix
     */
    void Matrix::init(MatrixFill fill) {
        // Generate a table given the number of rows and columns.
        unsigned int nels = getRows() * getCols();
        this->table.resize(nels);
        // Fill the matrix according to the given fill pattern.
        MPTime zeroValue = MPTime(0.0);
        switch(fill) {
            case MatrixFill::MinusInfinity:
                for (unsigned int pos = 0; pos < nels; pos++) {
                    this->table[pos] = MP_MINUSINFINITY;
                }
                break;
            case MatrixFill::Zero:
                for (unsigned int pos = 0; pos < nels; pos++) {
                    this->table[pos] = zeroValue;
                }
                break;
            case MatrixFill::Identity:
                for (unsigned int rowIndex = 0; rowIndex < getRows(); rowIndex++) {
                    for (unsigned int colIndex = 0; colIndex < getCols(); colIndex++) {
                        if(rowIndex == colIndex) {
                            put(rowIndex,colIndex,zeroValue);
                        } else {
                            put(rowIndex,colIndex,MP_MINUSINFINITY);
                        }
                    }
                }
            default:
                ;
        }
    }


    /**
     * Initialize matrix
     */
    void Matrix::init() {
        init(MatrixFill::MinusInfinity);
    }


    /**
     * Construct a square max-plus matrix of N by N
     */
    Matrix::Matrix(unsigned int N) {
        this->szRows = N;
        this->szCols = N;
        this->init();
    }
    /**
    * Construct a max-plus matrix of M by N
    
    Matrix::Matrix(unsigned int M, unsigned int N)
    {
        this->szRows = M;
        this->szCols = N;
        this->init(MatrixFill::MinusInfinity);
    }
   */
    /**
    *Creates a matrix with reserved memory
    */
    Matrix::Matrix(unsigned int nrows, unsigned int ncols, unsigned int nrel)
    {
        this->szRows = nrows;
        this->szCols = ncols;
        this->table.reserve(nrel);
        this->init();
    }

    /**
     * Destructor of MaxPlus matrix
     */
    Matrix::~Matrix() {
    }

    /**
    * Increases the number of rows of the matrix by n and fills the new elements with -\infty.  
    */
    void Matrix::addRows(uint n)
    {
        this->szRows = this->szRows + n;
        unsigned int nels = this->getRows() * this->getCols();
        this->table.resize(nels);
        for (unsigned int pos = nels - n * this->getCols(); pos < nels; pos++)
        {
            this->table[pos] = MP_MINUSINFINITY;
        }
    }


    /**
     * Get size of a square matrix
     */
    unsigned int Matrix::getSize(void) const {
        assert(this->getRows() == this->getCols());
        return this->getRows();
    }


    /**
     * Get an entry from the matrix. Row and column index must be between 0 and size-1
     */
    MPTime Matrix::get(unsigned int row, unsigned int column) const {
        if ((row >= this->getRows()) || (column >= this->getCols())) {
            throw CException("Index out of bounds in"
                                     "Matrix::get");
        }
        return this->table[row * this->getCols() + column];
    }

    /**
     * Get a row from the matrix as a vector. Row index must be between 0 and row size-1
     */
    Vector Matrix::getRowVector(unsigned int row) const {
        unsigned int sz = this->getCols();
        Vector result(sz);
        for (unsigned int i = 0; i < sz; i++) {
            result.put(i, this->get(row, i));
        }
        return result;
    }


    /**
     * Put a value in the matrix. Row and column index must be between 0 and size-1
     */
    void Matrix::put(unsigned int row, unsigned int column, MPTime value) {
        if ((row >= this->getRows()) || (column >= this->getCols())) {
            throw CException("Index out of bounds in"
                                     "Matrix::put");
        }
        this->table[row * this->getCols() + column] = value;
    }

    /**
     * Paste submatrix into matrix
     */
    void Matrix::paste(unsigned int top_row, unsigned int left_column, const Matrix *pastedMatrix) {
        const unsigned int p_rsz = pastedMatrix->getRows();
        const unsigned int p_csz = pastedMatrix->getCols();
        const unsigned int bottom_row = top_row + p_rsz - 1;
        const unsigned int right_column = left_column + p_csz - 1;
        assert(bottom_row < this->getRows());
        assert(right_column < this->getCols());
        unsigned int p_row = 0;
        unsigned int p_col = 0;

        for (unsigned int row = top_row; row <= bottom_row; row++, p_row++) {
            p_col = 0;
            for (unsigned int col = left_column; col <= right_column; col++, p_col++) {
                this->put(row, col, pastedMatrix->get(p_row, p_col));
            }
        }
    }

    /**
     * Paste rowvector into matrix
     */
    void Matrix::pasteRowVector(unsigned int top_row, unsigned int left_column, const Vector* pastedVector) {
        const unsigned int p_csz = pastedVector->getSize();
        const unsigned int right_column = left_column + p_csz - 1;
        assert(right_column < this->getCols());
        for (unsigned int col = left_column; col <= right_column; col++) {
            this->put(top_row, col, pastedVector->get(col));
        }
    }


    /**
     * mpmultiply()
     * Matrix-vector multiplication.
     */
    Vector *Matrix::mpmultiply(const Vector &v) const {
        // Check size of the matrix and vector
        if (this->getCols() != v.getSize())
        {
            throw CException("Matrix and vector are of unequal size in "
                             "Matrix::mpmultiply");
        }

        // Allocate space of the resulting vector
        Vector *res = new Vector(this->getRows());

        // Perform point-wise multiplication
        for (unsigned int i = 0; i < this->getRows(); i++) {
            MPTime m = MP_MINUSINFINITY;
            for (unsigned int k = 0; k < this->getCols(); k++) {
                m = MP_MAX(m, MP_PLUS(this->get(i, k), v.get(k)));
            }
            res->put(i, m);
        }
        return res;
    }

	/**
	* mpmultiply()
	* Matrix-matrix multiplication.
	*/
	Matrix *Matrix::mpmultiply(const Matrix &m) const
	{
		// Check sizes of the matrices
		if (this->getCols() != m.getRows())
		{
			throw CException("Matrices are of incompatible size in"
				"Matrix::mpmultiply(Matrix)");
		}

		// Allocate space of the resulting matrix
		Matrix *res = makeMatrix(this->getRows(), m.getCols());

		// Perform point-wise multiplication
		for (unsigned int i = 0; i < this->getRows(); i++)
		{
			for (unsigned int j = 0; j < m.getCols(); j++)
			{
				MPTime mpt = MP_MINUSINFINITY;
				for (unsigned int k = 0; k < this->getCols(); k++)
				{
					mpt = MP_MAX(mpt, MP_PLUS(this->get(i, k), m.get(k, j)));
				}
				res->put(i, j, mpt);
			}
		}
		return res;
	}

	/**
    * mpmultiply()
    * Matrix-matrix subtraction.
    */
    Matrix *Matrix::mpsub(const Matrix &m) const
    {
        // Check sizes of the matrices
        if ((m.getRows() != this->getRows()) ||
            (m.getCols() != this->getCols()))
        {
            throw CException("Matrices are of different size in"
                "Matrix::mpsub(Matrix&");
        }

        // Allocate space of the resulting matrix
        Matrix *res = makeMatrix(this->getRows(), this->getCols());

        // Perform element-wise substraction
        for (unsigned int i = 0; i < this->getRows(); i++)
        {
            for (unsigned int j = 0; j < this->getCols(); j++)
            {
                res->put(i, j, this->get(i, j) - m.get(i, j));
            }
        }
        return res;
    }

    /**
    * mpmultiply()
    * Matrix-matrix maximization.
    */
    Matrix *Matrix::mpmaximum(const Matrix &m) const
    {
        // Check sizes of the matrices
        if ((m.getRows() != this->getRows()) ||
            (m.getCols() != this->getCols()))
        {
            throw CException("Matrices are of different size in"
                "Matrix::maximum(Matrix*, Matrix*");
        }

        // Allocate space of the resulting matrix
        Matrix *res = makeMatrix(this->getRows(), this->getCols());

        // Perform element-wise substraction
        for (unsigned int i = 0; i < this->getRows(); i++)
        {
            for (unsigned int j = 0; j < this->getCols(); j++)
            {
                res->put(i, j, MP_MAX(this->get(i, j) , m.get(i, j)));
            }
        }
        return res;
    }
	/**
	* mppower()
	* Raise matrix to a positive integer power >= 1.
	*/
	Matrix *Matrix::mppower(const unsigned int p) const
	{
		Matrix* res;

		// base case p==1
		if (p == 1)
		{
			return this->createCopy();
		}

		// check if p is odd
		if(p & 1) {
			Matrix* mpow = this->mppower(p - 1);
			res = this->mpmultiply(*mpow);
			delete mpow;
		}
		else { //  p is even
			Matrix* mpow = this->mppower(p>>1);
			res = mpow->mpmultiply(*mpow);
			delete mpow;
		}
		return res;
	}

    /**
     * Matrix copy.
     */
    Matrix *Matrix::createCopy() const
    {
        Matrix *newMatrix = makeMatrix(this->getRows(), this->getCols());
        unsigned int nels = this->getRows() * this->getCols();
        for (unsigned int pos = 0; pos < nels; pos++) {
            newMatrix->table[pos] = this->table[pos];
        }
        return newMatrix;
    }

    /**
     * Matrix transposed copy.
     */
    Matrix *Matrix::getTransposedCopy() const {
        unsigned int MR = this->getCols();
        unsigned int MC = this->getRows();
        Matrix *newMatrix = makeMatrix(MR, MC);
        for (unsigned int col = 0; col < MC; col++)
        {
            for (unsigned int row = 0; row < MR; row++)
            {
                newMatrix->put(row, col, this->get(col, row));
            }
        }
        return newMatrix;
    }

    /**
     * Make submatrix with indices in list.
     */
    Matrix *Matrix::getSubMatrix(const list<unsigned int> &rowIndices, const list<unsigned int> &colIndices) const
    {
        unsigned int NR = (unsigned int) rowIndices.size();
        unsigned int NC = (unsigned int) colIndices.size();
        Matrix *newMatrix = makeMatrix(NR, NC);

        list<unsigned int>::const_iterator rit = rowIndices.begin();
        for (unsigned int r = 0; r < NR; r++, rit++) {
            unsigned int ri = (*rit);
            list<unsigned int>::const_iterator cit = colIndices.begin();
            for (unsigned int c = 0; c < NC; c++, cit++) {
                unsigned int ci = (*cit);
                newMatrix->put(r, c, this->get(ri, ci));
            }
        }
        return newMatrix;
    }

    /**
     * Make submatrix with indices in list from square matrix
     */
    Matrix *Matrix::getSubMatrix(const list<unsigned int> &indices) const {
        assert(this->getRows() == this->getCols());
        return this->getSubMatrix(indices, indices);
    }

    /**
    * Make submatrix with indices in list for non-square matrix. The new matrix only keeps the coloumns of the original matrix with the selected indices.
    */
    Matrix *Matrix::getSubMatrixNonSquare(const list<unsigned int> &colIndices) const
    {
        unsigned int NC = (unsigned int) colIndices.size();
        Matrix *newMatrix = makeMatrix(this->getRows(), NC);

        for (unsigned int r = 0; r < this->getRows(); r++) {
            list<unsigned int>::const_iterator cit = colIndices.begin();
            for (unsigned int c = 0; c < NC; c++, cit++) {
                unsigned int ci = (*cit);
                newMatrix->put(r, c, this->get(r, ci));
            }
        }
        return newMatrix;
    }

    /**
     * Matrix addition of scalar.
     */
    Matrix *Matrix::add(MPTime increase) const {
        unsigned int MR = this->getRows();
        unsigned int MC = this->getCols();
        Matrix *result = makeMatrix(MR, MC);

        this->add(increase, result);
        return result;
    }

    /**
     * Matrix addition of scalar with existing result matrix.
     */
    void Matrix::add(MPTime increase, Matrix *result) const {
        unsigned int MR = this->getRows();
        unsigned int MC = this->getCols();
        if ((MR != result->getRows()) || (MC != result->getCols())) {
            throw CException("Matrices are of different size in"
                                     "Matrix::add(Matrix*, MPTime, Matrix*");
        }
        for (unsigned int r = 0; r < MR; r++) {
            for (unsigned int c = 0; c < MC; c++) {
                result->put(r, c, this->get(r, c) + increase); // uses MP_PLUS()
            }
        }
    }


    /**
     * Matrix maximum with existing result matrix.
     */
    void Matrix::maximum(const Matrix *matB, Matrix *result) {
        unsigned int MR = this->getRows();
        unsigned int MC = this->getCols();
        if ((matB->getRows() != MR) ||
            (matB->getCols() != MC) ||
            (result->getRows() != MR) ||
            (result->getCols() != MC)) {
            throw CException("Matrices are of different size in"
                                     "Matrix::maximum(Matrix*, Matrix*, Matrix*");
        }

        for (unsigned int r = 0; r < MR; r++) {
            for (unsigned int c = 0; c < MC; c++) {
                result->put(r, c, MP_MAX(this->get(r, c), matB->get(r, c)));
            }
        }
    }

    /**
     * Matrix to string.
     */
    void Matrix::toString(CString &outString, double scale) const {
        outString = "";
        unsigned int MR = this->getRows();
        unsigned int MC = this->getCols();
        for (unsigned int i = 0; i < MR; i++) {
            for (unsigned int j = 0; j < MC; j++) {
                outString += timeToString(this->get(i, j) * scale) + " ";
            }
            outString += "\n";
        }
    }

	/**
	 * Matrix to string.
	 */
	void Matrix::toMatlabString(CString& outString, double scale) const {
		outString = "";
		unsigned int MR = this->getRows();
		unsigned int MC = this->getCols();
		outString += "[\n";
		for (unsigned int i = 0; i < MR; i++) {
			for (unsigned int j = 0; j < MC; j++) {
				outString += timeToMatlabString(this->get(i, j) * scale) + " ";
			}
			if (i < MR - 1) {
				outString += ";\n";
			}
			else {
				outString += "\n";
			}
		}
		outString += "]\n";
	}

	/**
	 * Matrix to LaTex string.
	 */
	void Matrix::toLaTeXString(CString& outString, double scale) const {
		outString = "";
		unsigned int MR = this->getRows();
		unsigned int MC = this->getCols();
		outString += "\\begin{bmatrix}\n";
		for (unsigned int i = 0; i < MR; i++) {
			for (unsigned int j = 0; j < MC; j++) {
				outString += timeToLaTeXString(this->get(i, j) * scale) + " ";
				if (j < MC - 1) {
					outString += "&";
				}
				else {
					outString += "\\\\";
				}
			}
			outString += "\n";
		}
		outString += "\\end{bmatrix}\n";
	}



    /**
     * Matrix return largest element.
     */
    MPTime Matrix::largestFiniteElement() const {
        // finite element with the largest absolute value
        // all -INF => 0
        //
        MPTime largestEl = 0, largestMag = 0;
        unsigned int MR = this->getRows();
        unsigned int MC = this->getCols();

        for (unsigned int r = 0; r < MR; r++) {
            for (unsigned int c = 0; c < MC; c++) {
                if (this->get(r, c) == MP_MINUSINFINITY) continue;

                MPTime mag = fabs(this->get(r, c));
                if (mag > largestMag) {
                    largestEl = this->get(r, c);
                    largestMag = mag;
                }
            }
        }

        return largestEl;
    }

    /**
     * Matrix get minimal element
     */
    MPTime Matrix::minimalFiniteElement() const {
        // smallest finite element
        // all -INF => -INF
        unsigned int MR = this->getRows();
        unsigned int MC = this->getCols();

        if (MR == 0 || MC == 0) return MP_MINUSINFINITY;

        MPTime minimalEl = MP_MINUSINFINITY;

        for (unsigned int r = 0; r < MR; r++) {
            for (unsigned int c = 0; c < MC; c++) {
                if (this->get(r, c) == MP_MINUSINFINITY) continue;

                if (minimalEl == MP_MINUSINFINITY) {
                    minimalEl = this->get(r, c);
                } else {
                    minimalEl = MP_MIN(minimalEl, this->get(r, c));
                }
            }
        }

        return minimalEl;
    }

    /**
     * Matrix plus closure.
     */
    Matrix *Matrix::plusClosureMatrix(MPTime posCycleThre) const {
        // notation: A^(+) = max(A, A^2, ...)
        return Matrix::allPairLongestPathMatrix(posCycleThre, false /*implyZeroSelfEdges*/);
    }

    /**
     * Matrix star closure.
     */
    Matrix *Matrix::starClosureMatrix(MPTime posCycleThre) const {
        // notation: A^(*) = max(E,A,A^2,...)
        // E - diagonal matrix with 'e'==0 on the diagonal

        return Matrix::allPairLongestPathMatrix(posCycleThre, true /*implyZeroSelfEdges*/);
    }


    /**
     * Matrix all pair longest path.
     */
    Matrix *Matrix::allPairLongestPathMatrix(MPTime posCycleThre, bool implyZeroSelfEdges) const {
        // Floyd-Warshall algorithm
        if (this->getRows() != this->getCols())
            throw CException("Matrix must be square in Matrix::allPaiLongestPathMatrix.");
        unsigned int N = this->getRows();


        Matrix *distMat = this->createCopy();

        // k - intermediate node
        for (unsigned int k = 0; k < N; k++) {
            for (unsigned int u = 0; u < N; u++) {
                for (unsigned int v = 0; v < N; v++) {
                    MPTime extra = (implyZeroSelfEdges && u == v) ? (MPTime) 0 : MP_MINUSINFINITY;
                    MPTime path_u2v = MP_MAX(distMat->get(v, u), extra);
                    MPTime path_u2k = distMat->get(k, u);
                    MPTime path_k2v = distMat->get(v, k);

                    MPTime path_u2v_candidate = (path_u2k + path_k2v);  // uses MP_PLUS()
                    if (path_u2v_candidate > path_u2v) {
                        path_u2v = path_u2v_candidate;
                    }
                    distMat->put(v, u, path_u2v);
                }
            }
        }

        for (unsigned int k = 0; k < N; k++) {
            if (distMat->get(k, k) > posCycleThre) {
                CString tmp;
                distMat->toString(tmp, 1.0e-06);
                cout << tmp << endl;
                throw CException("Positive cycle!");
            }
        }

        return distMat;
    }


    /**
    * Matrix all pair longest path. Retruns ture if there is a positive cycle. 
	* TODO: unify with the method above to have only one longest path implementation
    */
    bool Matrix::allPairLongestPathMatrix(MPTime posCycleThre, bool implyZeroSelfEdges, Matrix &res) const
    {
        // Floyd-Warshall algorithm
        if (this->getRows() != this->getCols())
            throw CException("Matrix must be square in Matrix::allPaiLongestPathMatrix.");
        unsigned int N = this->getRows();

        if ((N != res.getRows()) || (N != res.getCols()))
        {
            throw CException("The matrix of the longest pathes "
                "should have the same size as the given matirx.");
        }

		// TODO: make / use copy operation
        for (unsigned int u = 0; u < N; u++)
        {
            for (unsigned int v = 0; v < N; v++)
            {
                res.put(u, v, this->get(u, v));
            }
        }

        // k - intermediate node
        for (unsigned int k = 0; k < N; k++)
        {
            for (unsigned int u = 0; u < N; u++)
            {
                for (unsigned int v = 0; v < N; v++)
                {
                    MPTime extra = (implyZeroSelfEdges && u == v) ? (MPTime)0 : MP_MINUSINFINITY;
                    MPTime path_u2v = MP_MAX(res.get(v, u), extra);
                    MPTime path_u2k = res.get(k, u);
                    MPTime path_k2v = res.get(v, k);

                    MPTime path_u2v_candidate = (path_u2k + path_k2v);  // uses MP_PLUS()
                    if (path_u2v_candidate > path_u2v)
                    {
                        path_u2v = path_u2v_candidate;
                    }
                    res.put(v, u, path_u2v);
                }
            }
        }

        for (unsigned int k = 0; k < N; k++)
        {
            if (res.get(k, k) > posCycleThre)
            {
                return true;
            }
        }

        return false;
    }

    // factory methods
    Matrix* Matrix::makeMatrix(unsigned int nrows, unsigned int ncols) const
    {
        return new Matrix(nrows, ncols);
    }
    /**
     * class VectorList
     */

    /**
     * VectorList::toString()
     */
    void VectorList::toString(CString &outString, double scale) const
    {
        outString = "";
        for (unsigned int i = 0; i < this->getSize(); i++)
        {
            CString vec_str;
            const Vector *v = &vectorRefAt(i);
            assert(v->getSize() == this->oneVectorSize);
            v->toString(vec_str, scale);
            outString += vec_str;
            outString += "\n";
        }
    }

  //  /**
  //   * VectorList::findSimilar
  //   * test if list contains a vector which differs to vecX by less than threshold
	 //* implementation incomplete!
  //   */
  //  bool VectorList::findSimilar(const Vector &vecX, double threshold) const
  //  {

  //      // similar - differs by a constant within a threshold
  //      bool found = false;
  //      assert(vecX.getSize() == this->oneVectorSize);
  //      // test all vectors in list
  //      for (unsigned int i = 0; i < this->size(); i++)
  //      {
  //          const Vector &vecY = this->vectorRefAt(i);
  //          assert(vecY.getSize() == this->oneVectorSize);

  //          // determine min and max difference
  //          double minDiff = 0;
  //          double maxDiff = 0;
  //          bool min_def = false;
  //          bool max_def = false;
  //          for (unsigned int j = 0; j < this->oneVectorSize; j++)
  //          {
  //              bool vecX_inf = MP_ISMINUSINFINITY(vecX.get(j));
  //              bool vecY_inf = MP_ISMINUSINFINITY(vecY.get(j));
  //              if (vecX_inf || vecY_inf)
  //              {
  //                  if (vecX_inf && vecY_inf)
  //                  {
  //                      continue;
  //                  }
  //                  else
  //                  {
  //                      // Difference is infinity
  //                      minDiff = MP_MINUSINFINITY;
  //                      maxDiff = 0;
  //                      break;
  //                  }
  //              }
  //          }
  //          return true;
  //      }
		//return false;
  //  }

    CDouble Matrix::mpeigenvalue() const
    {
        // check if matrix is square.
        if (this->getRows() != this->getCols()) {
            throw CException("Matrix is not square in Matrix::mpeigenvalue().");
        }

        uint sz = this->getRows();

        // convert matrix to mcm graph
        
        // Create a new MCM graph
        MCMgraph *mcmGraph = new MCMgraph();

        // store vector of nodes
        vector<MCMnode*> nodes(sz);

        // Generate ids by countings
        CId id = 0;
        for (uint i = 0; i != sz; i++)
        {
            // Create an MCM node for this state
            MCMnode *n = new MCMnode(id, true);
            id++;

            // Add the node to the MCM graph
            mcmGraph->addNode(n);
            nodes[i] = n;
        }

        // Add edges to the MCM graph
        // for all states of the fsm...
        CId edgeId = 0;
        uint row = 0;
        uint col = 0;
        vector<MPTime>::const_iterator i;
        for (i = this->table.begin(); i != this->table.end(); i++)
        {
            MCMedge *e = new MCMedge(edgeId++, true);
            e->src = nodes[col];
            e->dst = nodes[row];
            e->w = (*i);
            e->d = 1.0;

            // Add the edge to the MCM graph and the src and dst node
            mcmGraph->addEdge(e);
            edgeId++;

            col++;
            if (col == sz) {
                col = 0;
                row++;
            }
        }

        // prune the graph
        MCMgraph *pruned = mcmGraph->pruneEdges();
        delete mcmGraph;

        // compute MCM
        CDouble res = pruned->calculateMaximumCycleMeanKarpDouble();
        delete pruned;
        return res;

    }

    MCMgraph Matrix::mpMatrixToPrecedenceGraph() const
    {
        // check if matrix is square.
        if (this->getRows() != this->getCols()) {
            throw CException("Matrix is not square in Matrix::mpMatrixToPrecedenceGraph().");
        }

        uint sz = this->getRows();

        // convert matrix to precedence graph

        // Create a new MCM graph
        MCMgraph precGraph;

        // store vector of nodes
        vector<MCMnode*> nodes(sz);

        // Generate ids by countings
        CId id = 0;
        for (uint i = 0; i != sz; i++)
        {
            // Create an MCM node for state element i
            MCMnode* n = new MCMnode(id, true);
            id++;

            // Add the node to the MCM graph
            precGraph.addNode(n);
            nodes[i] = n;
        }

        // Add edges to the MCM graph
        CId edgeId = 0;
        uint row = 0;
        uint col = 0;
        for (auto i = this->table.begin(); i != this->table.end(); i++)
        {
            // add edge is the value is not minus infinity
            if ((*i) != MP_MINUSINFINITY) {
                MCMedge* e = new MCMedge(edgeId++, true);
                e->src = nodes[col];
                e->dst = nodes[row];
                e->w = (CDouble)(*i); 
                e->d = 1.0;

                // Add the edge to the MCM graph and the src and dst node
                precGraph.addEdge(e);
                edgeId++;
            }
            // update row and col
            col++;
            if (col == sz) {
                col = 0;
                row++;
            }
        }

        return precGraph;
    }


    std::pair<Matrix::EigenvectorList, Matrix::GeneralizedEigenvectorList> Matrix::mpGeneralizedEigenvectors() const
    {
        // check if matrix is square.
        if (this->getRows() != this->getCols()) {
            throw CException("Matrix is not square in Matrix::mpeigenvector().");
        }

        //gr = mpMatrixToPrecendenceGraph(M)
        // compute precedence graph
        MCMgraph precGraph = this->mpMatrixToPrecedenceGraph();
        MCMgraphs sccs;

        // compute SCCs including single nodes without edges.
        stronglyConnectedMCMgraph(&precGraph, sccs, true);

        // map from number of SCC to SCC
        std::map<uint, MCMgraph*> sccMapInv;

        // map from node of precGraph to its SCC
        std::map<MCMnode*, uint> sccMap;

        // map from nodes of precGraph to the cycle mean of its SCC
        std::map<MCMnode*, MPTime> cycleMeansMap;

        // vector such that element k is a node from precGraph in the critical path of SCC k
        std::vector<MCMnode*> criticalNodes;

        // vector such that element k is the maximum cycle mean of SCC k
        std::vector<MPTime> cycleMeans;

        // SCC counter k
        uint k = 0;
        // for each SCC
        for (auto scci = sccs.cbegin(); scci != sccs.cend(); scci++, k++) {
            MCMgraph* scc = *scci;
            sccMapInv[k] = scc;

            // MCM calculation requires node relabelling
            std::map<int, int> sccNodeIdMap;
            scc->relabelNodeIds(sccNodeIdMap);
            
            if (scc->nrVisibleEdges() > 0) {
                // compute MCM mu and critical node n of scc
                MCMnode* n;
                MPTime mu = scc->calculateMaximumCycleMeanKarpDouble(&n);
                criticalNodes.push_back(precGraph.getNode(sccNodeIdMap[n->id]));
                cycleMeans.push_back(mu);
            }
            else {
                // a SCC without edges is a single node and has no cycle mean
                criticalNodes.push_back(NULL);
                cycleMeans.push_back(MP_MINUSINFINITY);
            }
            // for each node in the scc
            for (auto ni = scc->getNodes().cbegin(); ni != scc->getNodes().cend(); ni++) {
                auto nn = *ni;
                auto nnn = precGraph.getNode(sccNodeIdMap[nn->id]);
                // map node to SCC index
                sccMap[nnn] = k;
                // map node to the cycle mean of its SCC
                cycleMeansMap[nnn] = cycleMeans[k];
            }
        }


        // compute the eigenvectors
        Matrix::EigenvectorList eigenVectors;
        Matrix::GeneralizedEigenvectorList genEigenVectors;

        for (size_t k = 0; k < sccs.size(); k++) {
            // there is one for each SCC with a cycle mean larger than -inf
            if (cycleMeans[k] != MP_MINUSINFINITY) {
                // eigenvector is formed by normalized longest paths from critical node to all other nodes 
                // compute transitive cycle means such that all nodes in SCC k and downstream SCCs get a cyclemean
                // that is the maximum if all (reflexive) upstream SCCs
                std::map<CId,CDouble> trCycleMeans;
                // initialize all nodes to undefined (represented by -DBL_MAX), except the root node
                // which is initialized with its own cycle mean
                for (unsigned int n = 0; n != precGraph.getNodes().size(); n++) {
                    trCycleMeans[n] = (n == criticalNodes[k]->id) ? (cycleMeans[k]) : MP_MINUSINFINITY;
                }
                bool change = true;
                while (change) {
                    change = false;
                    // for all edges in precGraph
                    for (auto i = precGraph.getEdges().begin(); i != precGraph.getEdges().end(); i++) {
                        MCMedge& e = **i;
                        if (trCycleMeans[e.src->id] > trCycleMeans[e.dst->id]) {
                            change = true;
                            trCycleMeans[e.dst->id] = MP_MAX(trCycleMeans[e.src->id], cycleMeans[sccMap[e.dst]]);
                        }
                    }
                }

                // compute normalization map replace MP_MINUSINFINITE by -DBL_MAX
                std::map<CId,CDouble> muMap;
                for (unsigned int n = 0; n < trCycleMeans.size(); n++) {
                    muMap[n] = (trCycleMeans[n] == MP_MINUSINFINITY) ? -DBL_MAX : trCycleMeans[n];
                }

                // compute normalized longest paths
                std::map<CId, CDouble> lengths = precGraph.normalizedLongestPaths(criticalNodes[k]->id, muMap);
                // make an eigenvector
                Vector v(this->getCols());
                for (auto i = lengths.cbegin(); i != lengths.cend(); i++) {
                    CId n = (*i).first;
                    MPTime  value = (((MPTime) (*i).second) <= MP_MINUSINFINITY ? MP_MINUSINFINITY : (MPTime) (*i).second);
                    v.put((unsigned int) n, value);
                }

                // check if it is a generalized eigenvalue
                bool isGeneralized = false;
                MPTime lambda = MP_MINUSINFINITY;
                Vector ev(this->getCols());
                int l = 0;
                for (auto i = precGraph.getNodes().cbegin(); i != precGraph.getNodes().cend(); i++, l++) {
                    const MCMnode& n = **i;
                    if (lambda == MP_MINUSINFINITY) {
                        lambda = trCycleMeans[n.id];
                    }
                    else {
                        if ((lambda != MP_MINUSINFINITY) && (trCycleMeans[n.id] != MP_MINUSINFINITY) && (lambda != trCycleMeans[n.id])) isGeneralized = true;
                    }
                    ev.put(l, trCycleMeans[n.id]);
                }

                if (isGeneralized) {
                    genEigenVectors.push_back(std::make_pair(v, ev));
                }
                else {
                    eigenVectors.push_back(std::make_pair(v, lambda));
                }
            }
        }
 
        // the SCC graphs still need to be destroyed
        for (auto i = sccs.begin(); i != sccs.end(); i++) {
            MCMgraph* g = *i;
            delete g;
        }

        return std::make_pair(eigenVectors, genEigenVectors);

    }

    Matrix::EigenvectorList Matrix::mpEigenvectors() const
    {
        auto gev = this->mpGeneralizedEigenvectors();
        return gev.first;
    }



    void ExtendedMatrix::put(unsigned int row, unsigned int column, MPTime value, std::unordered_set<int>& bs)
    {
        this->Matrix::put(row, column, value);
        this->bufferSets[row * this->getCols() + column] = bs;
    }

    Matrix * ExtendedMatrix::getSubMatrix(const list<unsigned int> &rowIndices, const list<unsigned int> &colIndices) const
    {
        unsigned int NR = rowIndices.size();
        unsigned int NC = colIndices.size();
        Matrix *newMatrix = makeMatrix(NR, NC);
        ExtendedMatrix *newExtendedMatrix = dynamic_cast<ExtendedMatrix*> (newMatrix);

        list<unsigned int>::const_iterator rit = rowIndices.begin();
        for (unsigned int r = 0; r < NR; r++, rit++)
        {
            unsigned int ri = (*rit);
            list<unsigned int>::const_iterator cit = colIndices.begin();
            for (unsigned int c = 0; c < NC; c++, cit++)
            {
                unsigned int ci = (*cit);
				std::unordered_set<int> bs = this->getBufferSet(ri, ci);
                newExtendedMatrix->put(r, c, this->get(ri, ci), bs);
            }
        }
        return newExtendedMatrix;
    }

    Matrix* ExtendedMatrix::makeMatrix(unsigned int nrows, unsigned int ncols) const
    {
        return new ExtendedMatrix(nrows, ncols);
    }

    std::unordered_set<int> ExtendedMatrix::getBufferSet(unsigned int row, unsigned int column) const
    {
        if ((row >= this->getRows()) || (column >= this->getCols()))
        {
            throw CException("Index out of bounds in"
                "Matrix::get");
        }
        return this->bufferSets[row * this->getCols() + column];
    }
}


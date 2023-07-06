/*
 *  Eindhoven University of Technology
 *  Eindhoven, The Netherlands
 *  Dept. of Electrical Engineering
 *  Electronics Systems Group
 *  Model Based Design Lab (https://computationalmodeling.info/)
 *
 *  Name            :   mpmatrix.h
 *
 *  Author          :   Marc Geilen (m.c.w.geilen@tue.nl)
 *
 *  Date            :   March 23, 2009
 *
 *  Function        :   MaxPlus matrices
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


#ifndef MAXPLUS_ALGEBRA_MATRIX_H_INCLUDED
#define MAXPLUS_ALGEBRA_MATRIX_H_INCLUDED

#include "mptype.h"
#include "base/analysis/mcm/mcmgraph.h"
#include <vector>
#include <unordered_set>

using namespace Graphs;


class CString;

namespace MaxPlus {

    /**
     * Vector, represents a max-plus column vector.
     */
    class Vector {
    public:
        Vector(unsigned int size = 0, MPTime value = MP_MINUSINFINITY);

        Vector(std::vector <MPTime> *v);

        ~Vector();

        inline unsigned int getSize(void) const {
            return (unsigned int) this->table.size();
        }

        inline MPTime get(unsigned int row) const {
            return this->table[row];
        }
        void toLaTeXString(CString& outString, double scale = 1.0) const;
        void put(unsigned int row, MPTime value);

        void toString(CString &outString, double scale = 1.0) const;

        Vector(const Vector &);
        void randomize(unsigned int max);
        Vector &operator=(const Vector &);

        MPTime norm();

        void negate();
        MPTime normalize();


        Vector *add(MPTime increase) const;

        void add(MPTime increase, Vector *result) const;

        void maximum(const Vector *matB, Vector *result) const;

        Vector *add(const Vector *vecB) const;

        void add(const Vector *vecB, Vector *res) const;

        Vector &operator+=(MPTime increase) {
            this->add(increase, this);
            return *this;
        }

        Vector &operator-=(MPTime decrease) {
            assert(!MP_ISMINUSINFINITY(decrease));
            this->add(-decrease, this);
            return *this;
        }

        bool compare(const Vector &v);

        Vector &incrementalMaximum(const Vector *vec) {
            this->maximum(vec, this);
            return *this;
        }

        /**
         * Get minimal finite element
         * returns the smallest among the finite elements in the vector or
         * MP_MINUSINFINITY if no finite elements exist
         * itsPosition returns the index of the (a) smallest finite element is set
         * to a pointer to unsigned int, otherwise set or defaults to NULL
         */
        MPTime minimalFiniteElement(unsigned int *itsPosition_Ptr = NULL) const;

    private:
        vector <MPTime> table;
    };


    enum class MatrixFill {MinusInfinity, Zero, Identity};

    class Matrix {
    public:

		/**
		 * Construct a max-plus matrix of \p N by \p N using the default fill pattern (minus infinity).
		 * @param N Number of rows and columns
		 * @param ncols Number of columns
		 */
		Matrix(unsigned int N);

        /**
         * Construct a max-plus matrix of M by N using the default fill pattern (minus infinity).
         * @param N Number of rows and columns
         * @param ncols Number of columns
         */


        //Matrix(unsigned int M, unsigned int N);
        /**
         * Construct a max-plus matrix of \p nrows by \p ncols using the given fill pattern.
         * @param nrows Number of rows
         * @param ncols Number of columns
         * @param fill Fill pattern for matrix
         */
        Matrix(unsigned int nrows, unsigned int ncols, MatrixFill fill = MatrixFill::MinusInfinity);

		Matrix(unsigned int nrows, unsigned int ncols, unsigned int nrel);

        virtual ~Matrix();

        inline unsigned int getRows(void) const {
            return this->szRows;
        }

        inline unsigned int getCols(void) const {
            return this->szCols;
        }

        unsigned int getSize(void) const;

        MPTime get(unsigned int row, unsigned int column) const;
        Vector getRowVector(unsigned int row) const;

        void put(unsigned int row, unsigned int column, MPTime value);

        void paste(unsigned int top_row, unsigned int left_column, const Matrix *pastedMatrix);

        void pasteRowVector(unsigned int top_row, unsigned int left_column, const Vector* pastedVector);

        virtual Matrix *createCopy() const;

        Matrix *getTransposedCopy() const;

        Matrix* getSubMatrixNonSquareRows(const list<unsigned int>& rowIndices) const;
        virtual Matrix *getSubMatrix(const list<unsigned int> &rowIndices, const list<unsigned int> &colIndices) const;

        Matrix *getSubMatrix(const list<unsigned int> &indices) const;

        Matrix *getSubMatrixNonSquare(const list<unsigned int> &indices) const;

		/**
		* Increases the number of rows of the matrix by n and fills the new elements with -\infty.
		*/
		void addRows(uint n);

        /**
        * Increases the number of cols of the matrix by n and fills the new elements with -\infty.
        */
        void addCols(uint n);

        void toString(CString &outString, double scale = 1.0) const;
		void toMatlabString(CString& outString, double scale = 1.0) const;
		void toLaTeXString(CString& outString, double scale = 1.0) const;



        // Algebraic operations.
        Matrix *add(MPTime increase) const;

        void add(MPTime increase, Matrix *result) const;

		Matrix* mpsub(const Matrix& m) const;

		Matrix* mpmaximum(const Matrix& m) const;

        void maximum(const Matrix *matB, Matrix *result);

        Vector *mpmultiply(const Vector &v) const;

        Matrix *mpmultiply(const Matrix &m) const;

        Matrix *mppower(const unsigned int p) const;

		CDouble mpeigenvalue() const;
        
        typedef std::list<std::pair<Vector, CDouble>> EigenvectorList;
        typedef std::list<std::pair<Vector, Vector>> GeneralizedEigenvectorList;
        std::pair<EigenvectorList,GeneralizedEigenvectorList> mpGeneralizedEigenvectors() const;
        EigenvectorList mpEigenvectors() const;

        Matrix &operator+=(MPTime increase) {
            this->add(increase, this);
            return *this;
        }

        Matrix &operator-=(MPTime decrease) {
            assert(!MP_ISMINUSINFINITY(decrease));
            this->add(-decrease, this);
            return *this;
        }

        Matrix &incrementalMaximum(const Matrix *matrix) {
            this->maximum(matrix, this);
            return *this;
        }

        bool operator==(const Matrix &other);

        /**
         * Return the element having the largest abs() value.
         * @return element having the largest abs() value
         */
        MPTime largestFiniteElement() const;
        MPTime minimalFiniteElement() const;
        MPTime getMaxOfCol(uint colNumber) const;
        MPTime getMaxOfRow(uint rowNumber) const;
        MPTime  getMaxOfRowUntilCol(uint rowNumber, uint colNumber) const;
        MPTime getMaxOfColUntilRow(uint colNumber,uint rowNumber) const;

        Matrix *plusClosureMatrix(MPTime posCycleThre = MP_EPSILON) const;

        Matrix *starClosureMatrix(MPTime posCycleThre = MP_EPSILON) const;
        // randomizes all elements of an existing matrix up to and including val max
        void randomize(unsigned int max);
        Matrix *allPairLongestPathMatrix(MPTime posCycleThre, bool implyZeroSelfEdges) const;
		bool allPairLongestPathMatrix(MPTime posCycleThre, bool implyZeroSelfEdges, Matrix& res) const;

        MCMgraph mpMatrixToPrecedenceGraph() const;

        // factory methods
        virtual Matrix* makeMatrix(unsigned int nrows, unsigned int ncols) const;

     private:
        // Implicit copying is not allowed
        //  => Intentionally private and not implemented
        Matrix(const Matrix &);

        Matrix &operator=(const Matrix &);

        void init(MatrixFill fill);
        void init();



    private:
        Matrix();

    private:
        vector <MPTime> table;
        unsigned int szRows;
        unsigned int szCols;
    };


    class ExtendedMatrix : public Matrix
    {
    public:
        ExtendedMatrix(unsigned int nrows, unsigned int ncols) : Matrix(nrows, ncols)
        {
            unsigned int nels = this->getRows() * this->getCols();
            this->bufferSets.resize(nels);
        }

        //Creates a matrix with reserved memory for nrel expected entities
        ExtendedMatrix(unsigned int nrows, unsigned int ncols, unsigned int nrel) : Matrix(nrows, ncols, nrel)
        {
            this->bufferSets.reserve(nrel);
        }
        ExtendedMatrix(unsigned int N) : Matrix(N)
        {}

        virtual Matrix *createCopy() const
        {
            Matrix *newMatrix = Matrix::createCopy();
            ExtendedMatrix* newExtendedMatrix = dynamic_cast<ExtendedMatrix*> (newMatrix);

            unsigned int nels = this->getRows() * this->getCols();
            for (unsigned int pos = 0; pos < nels; pos++)
            {
                newExtendedMatrix->bufferSets[pos] = this->bufferSets[pos];
            }

            return newExtendedMatrix;
        }

        void put(unsigned int row, unsigned int column, MPTime value, std::unordered_set<int>&);
        virtual Matrix *getSubMatrix(const list<unsigned int> &rowIndices, const list<unsigned int> &colIndices) const;

        // factory methods
        virtual Matrix* makeMatrix(unsigned int nrows, unsigned int ncols) const;

        std::unordered_set<int> getBufferSet(unsigned int row, unsigned int column) const;

    private:

        vector<std::unordered_set<int>> bufferSets;
    };


	/****************************************************
	* VectorList: usually represents a set of eigenvectors
	* More efficient than vector<MaxPlus::Vector>
	****************************************************/

	class VectorList : private std::vector<Vector* >
	{
	public:
		VectorList(unsigned int oneVectorSizeInit);
		~VectorList();

		const Vector& vectorRefAt(int n) const; // vector at index 'n'
		Vector& vectorRefAt(int n);

		const Vector& lastVectorRef() const; // last vector
		Vector& lastVectorRef();

		unsigned int getSize() const; // vector count
		unsigned int getOneVectorSize() const
		{
			return this->oneVectorSize;
		}

		void grow(); // append one vector place

		void toString(CString& outString, double scale = 1.0) const;

		//bool findSimilar(const Vector& vec, double threshold) const;
		// similar - differs by a constant within a threshold

	private:
		// Implicit copying is not allowed
		//  => Intentionally private and not implemented
		VectorList(const VectorList&);
		VectorList& operator=(const VectorList&);

	private:
		const unsigned int oneVectorSize;
	};

	inline VectorList::VectorList(unsigned int oneVectorSizeInit)
		: oneVectorSize(oneVectorSizeInit)
	{
		assert(oneVectorSize > 0);
	}

	inline VectorList::~VectorList()
	{
		for (unsigned int pos = 0; pos < size(); pos++)
		{
			delete this->at(pos);
		}
	}

	inline const Vector& VectorList::vectorRefAt(int n) const
	{
		return *this->at(n);
	}

	inline Vector& VectorList::vectorRefAt(int n)
	{
		return *this->at(n);
	}

	inline const Vector& VectorList::lastVectorRef() const
	{
		return *this->at(this->size() - 1);
	}

	inline Vector& VectorList::lastVectorRef()
	{
		return *this->at(this->size() - 1);
	}

	inline unsigned int VectorList::getSize() const
	{
		return (unsigned int)vector<Vector* >::size();
	}

	inline void VectorList::grow()
	{
		unsigned int last = (unsigned int)this->size();
		this->resize((size_t)(last + 1));
		this->at(last) = new Vector(oneVectorSize, MP_MINUSINFINITY);
	}

}

#endif
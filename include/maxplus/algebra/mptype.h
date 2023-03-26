/*
 *  Eindhoven University of Technology
 *  Eindhoven, The Netherlands
 *  Dept. of Electrical Engineering
 *  Electronics Systems Group
 *  Model Based Design Lab (https://computationalmodeling.info/)
 *
 *  Name            :   mptype.h
 *
 *  Author          :   Marc Geilen (m.c.w.geilen@tue.nl)
 *
 *  Date            :   March 23, 2009
 *
 *  Function        :   maxplus type definitions and operations
 *
 *  History         :
 *      23-03-09    :   Initial version.
 *      23-12-09    :   Peter Poplavko. Rewritten. Use C++ tools such as
 *                      inline functions, classes and operator overloading instead
 *                      of C tools.
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

#ifndef BASE_MAXPLUS_MPTYPE_H_INCLUDED
#define BASE_MAXPLUS_MPTYPE_H_INCLUDED

#include <assert.h>
#include "base/basic_types.h"
#include "base/string/cstring.h"

namespace MaxPlus
{

    typedef CDouble MPThroughput;

    class MPTime
    {
        public:
            MPTime(CDouble val = 1.0e+30)
            {
                myVal = val;
            }

            operator CDouble() const
            {
                return myVal;
            }
            MPTime &operator-();
            MPTime &operator+=(MPTime val);
            MPTime &operator-=(MPTime val);

        private:
            CDouble myVal;
    };

    typedef MPTime MPDelay;

    //==============================
    // MP_MAX()
    //==============================

    inline MPTime MP_MAX(MPTime a, MPTime b)
    {
        return ((a) > (b)) ? (a) : (b);
    }

    inline MPTime MP_MAX(CDouble a, MPTime b)
    {
        return MP_MAX(MPTime(a), b);
    }

    inline MPTime MP_MAX(MPTime a, CDouble b)
    {
        return MP_MAX(a, MPTime(b));
    }

    inline CDouble MP_MAX(CDouble a, CDouble b)
    {
        return CDouble(MP_MAX(MPTime(a), MPTime(b)));
    }

    //==============================
    // MP_MIN()
    //==============================

    inline MPTime MP_MIN(MPTime a, MPTime b)
    {
        return ((a) < (b)) ? (a) : (b);
    }

    inline MPTime MP_MIN(CDouble a, MPTime b)
    {
        return MP_MIN(MPTime(a), b);
    }

    inline MPTime MP_MIN(MPTime a, CDouble b)
    {
        return MP_MIN(a, MPTime(b));
    }

    inline CDouble MP_MIN(CDouble a, CDouble b)
    {
        return CDouble(MP_MIN(MPTime(a), MPTime(b)));
    }

    //==============================
    // MP_INFINITY
    //==============================

    // the quick and dirty way of representing -infinity
    const MPTime MP_MINUSINFINITY = -1.0e+30;

    inline bool MP_ISMINUSINFINITY(CDouble a)
    {
        return a < -1.0e+20;
    }

    inline const MPTime MP_PLUS(CDouble a, CDouble b)
    {
        return (MP_ISMINUSINFINITY(a) || MP_ISMINUSINFINITY(b)) ? MP_MINUSINFINITY : (MPTime)((CDouble)a + (CDouble)b);
    }


    // MaxPlus epsilon (used to compare floating point numbers for equality)
    const MPTime MP_EPSILON = 1e-10;


    //==============================
    // MPTime operators
    //==============================
    inline const MPTime operator+(MPTime a, MPTime b)
    {
        return MP_PLUS(a, b);
    }

    inline const MPTime operator+(MPTime a, CDouble b)
    {
        return MPTime(a) + MPTime(b);
    }

    inline const MPTime operator+(CDouble a, MPTime b)
    {
        return MPTime(a) + MPTime(b);
    }

    inline const MPTime operator-(MPTime a, MPTime b)
    {
        assert(!MP_ISMINUSINFINITY(b));
        return a + MPTime(-b);
    }

    inline const MPTime operator-(MPTime a, CDouble b)
    {
        return MPTime(a) - MPTime(b);
    }

    inline const MPTime operator-(CDouble a, MPTime b)
    {
        return MPTime(a) - MPTime(b);
    }

    inline MPTime &MPTime::operator-()
    {
        assert(!MP_ISMINUSINFINITY(myVal));
        myVal = -myVal;
        return *this;
    }

    inline const MPTime operator*(MPTime a, MPTime b)
    {
        if (MP_ISMINUSINFINITY(a))
        {
            assert(b > 0);
            return MP_MINUSINFINITY;
        }
        if (MP_ISMINUSINFINITY(b))
        {
            assert(a > 0);
            return MP_MINUSINFINITY;
        }
        return CDouble(a) * CDouble(b);
    }

    inline const MPTime operator*(CDouble a, MPTime b)
    {
        return MPTime(a) * MPTime(b);
    }

    inline const MPTime operator*(MPTime a, CDouble b)
    {
        return MPTime(a) * MPTime(b);
    }

    inline MPTime &MPTime::operator+=(MPTime a)
    {
        *this = *this + a;
        return *this;
    }

    inline MPTime &MPTime::operator-=(MPTime a)
    {
        assert(!MP_ISMINUSINFINITY(a));
        *this = *this + MPTime(CDouble(-a));
        return *this;
    }


    //==============================
    // toString
    //==============================

    inline const CString timeToString(MPTime val)
    {
        // We intentionally dont use MP_ISMINUSINFINITY(val) here,
        // so that we can expose the unwanted "impure" infinities here.
        //
        if (val == MP_MINUSINFINITY)
        {
            return CString("-mp_inf");
        }
        else return CString(val);
    }
	inline const CString timeToMatlabString(MPTime val)
	{
		if (val == MP_MINUSINFINITY)
		{
			return CString("-Inf");
		}
		else return CString(val);
	}
	inline const CString timeToLaTeXString(MPTime val)
	{
		if (val == MP_MINUSINFINITY)
		{
			return CString("-\\infty{}");
		}
		else return CString(val);
	}


}
#endif
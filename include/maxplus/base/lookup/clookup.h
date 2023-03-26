/*
 *  Eindhoven University of Technology
 *  Eindhoven, The Netherlands
 *  Dept. of Electrical Engineering
 *  Electronics Systems Group
 *  Model Based Design Lab (https://computationalmodeling.info/)
 *
 *  Name            :   clookup.h
 *
 *  Author          :   Peter Poplavko (poplavko@ics.ele.tue.nl)
 *
 *  Date            :   September 26, 2003
 *
 *  Function        :   CLookupXXX = handy wrappers to access hash maps
 *                         e.g. from pointer to pointer
 *
 *  History         :
 *      15-08-09    :   Initial version.
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

#ifndef BASE_STRING_CLOOKUP_H
#define BASE_STRING_CLOOKUP_H

/*****************************************************/
/* CLookupXXX classes                                */
/*   Handy wrappers to access hash maps              */
/*****************************************************/


#include  <unordered_map>

struct hash_void
{
    size_t operator()(const void *__x) const
    {
        return reinterpret_cast<size_t>(__x);
    }
};
typedef std::unordered_map<void *, void *, hash_void> CLOOKUP_MAP_VOID;
typedef std::unordered_map<void *, int, hash_void>   CLOOKUP_MAP_INT;
typedef std::unordered_map<int, int>   CLOOKUP_MAP_INT_INT;


/**
 * CLookupPtr<T>
 * Lookup from void pointer to pointer 'T*'
 */

template< class T>
class CLookupPtr
{

    public:
        friend class CIterator;

        /**
        * Iterator over CLookupPtr: visit all <key, value> pairs
        *
        */
        class CIterator
        {
            public:
                CIterator(const CLookupPtr<T>* lookup) : m_lookup(lookup), m_iterator(lookup->m_impl.begin())
                {
                }
                bool next(void*& key, T*& value)
                {
                    if (m_iterator != m_lookup->m_impl.end())
                    {
                        key = (*m_iterator).first;
                        value = static_cast<T *>((*m_iterator).second);
                        m_iterator++;
                        return true;
                    }
                    else
                    {
                        key = NULL;
                        value = NULL;
                        return false;
                    }
                }
                void reset()
                {
                    m_iterator = m_lookup->m_impl.begin();
                }

            private:
                const CLookupPtr<T>* m_lookup;
                CLOOKUP_MAP_VOID::const_iterator m_iterator;
        };

    public:
        //------------------------------------------
        //----- CLookupPtr -------------------------
        //------------------------------------------

        // TODO: void erase(key)

        CLookupPtr(bool freePointers) : m_freePointers(freePointers) {}
        ~CLookupPtr()
        {
            if (m_freePointers)
            {
                CIterator iterator(this);
                void *key;
                T *value;
                while (iterator.next(key, value))
                {
                    delete value;
                }
            }
        }

        T *get(void *key) const
        {
            CLOOKUP_MAP_VOID::const_iterator foundPointer = m_impl.find(key);
            if (foundPointer == m_impl.end())
            {
                // not found
                return NULL;
            }
            return static_cast<T *>((*foundPointer).second);
        }

        T *getOrCreate(void *key)
        {
            CLOOKUP_MAP_VOID::const_iterator foundPointer = m_impl.find(key);
            T *value;
            if (foundPointer == m_impl.end())
            {
                // not found
                value = new T;
                m_impl[key] = value;
            }
            else
            {
                value = static_cast<T *>((*foundPointer).second);
            }
            return value;
        }

        void put(void *key, T *value)
        {
            CLOOKUP_MAP_VOID::iterator foundPointer = m_impl.find(key);
            if (foundPointer == m_impl.end())
            {
                // not found
                m_impl[key] = value;
                return;
            }
            // found
            if (m_freePointers)
            {
                delete static_cast<T *>((*foundPointer).second);
            }
            (*foundPointer).second = value;
        }

        bool exists(void *key) const
        {
            CLOOKUP_MAP_VOID::const_iterator foundPointer = m_impl.find(key);
            if (foundPointer == m_impl.end())
            {
                return false;
            }
            else
            {
                return true;
            }
        }
    private:
        CLOOKUP_MAP_VOID m_impl;
        bool m_freePointers;
};

/**
*  CLookupInt: Lookup from pointer to integer
*
*/

class CLookupInt
{
    public:
        friend class CIterator;

        /**
        * Iterator over CLookupInt: visit all <key_pointer, integer_value> pairs
        *
        */
        class CIterator
        {
            public:
                CIterator(const CLookupInt *lookup) : m_lookup(lookup), m_iterator(lookup->m_impl.begin())
                {
                }
                bool next(void*& key, int &value)
                {
                    if (m_iterator != m_lookup->m_impl.end())
                    {
                        key   = (*m_iterator).first;
                        value = (*m_iterator).second;
                        m_iterator++;
                        return true;
                    }
                    else
                    {
                        key = NULL;
                        value = 0;
                        return false;
                    }
                }
                void reset()
                {
                    m_iterator = m_lookup->m_impl.begin();
                }

            private:
                const CLookupInt *m_lookup;
                CLOOKUP_MAP_INT::const_iterator m_iterator;
        };

    public:
        //------------------------------------------
        //----- CLookupInt -------------------------
        //------------------------------------------

        // TODO: void erase(key)

        int get(void *key) const
        {
            CLOOKUP_MAP_INT::const_iterator foundPointer = m_impl.find(key);
            if (foundPointer == m_impl.end())
            {
                // not found
                return 0;
            }
            return (*foundPointer).second;
        }

        void put(void *key, int value)
        {
            m_impl[key] = value;
        }

        bool exists(void *key) const
        {
            CLOOKUP_MAP_INT::const_iterator foundPointer = m_impl.find(key);
            if (foundPointer == m_impl.end())
            {
                return false;
            }
            else
            {
                return true;
            }
        }

    private:
        CLOOKUP_MAP_INT m_impl;
};

/**
*  CLookupIntInt: Lookup from integer to integer
*
*/

class CLookupIntInt
{
public:
	friend class CIterator;

	/**
	* Iterator over CLookupIntInt: visit all <integer_value, integer_value> pairs
	*
	*/
	class CIterator
	{
	public:
		CIterator(const CLookupIntInt *lookup) : m_lookup(lookup), m_iterator(lookup->m_impl.begin())
		{
		}
		bool next(int& key, int &value)
		{
			if (m_iterator != m_lookup->m_impl.end())
			{
				key = (*m_iterator).first;
				value = (*m_iterator).second;
				m_iterator++;
				return true;
			}
			else
			{
				key = 0;
				value = 0;
				return false;
			}
		}
		void reset()
		{
			m_iterator = m_lookup->m_impl.begin();
		}

	private:
		const CLookupIntInt *m_lookup;
		CLOOKUP_MAP_INT_INT::const_iterator m_iterator;
	};

public:
	//------------------------------------------
	//----- CLookupIntInt -------------------------
	//------------------------------------------

	// TODO: void erase(key)

	int get(int key) const
	{
		CLOOKUP_MAP_INT_INT::const_iterator foundPointer = m_impl.find(key);
		if (foundPointer == m_impl.end())
		{
			// not found
			return -1;
		}
		return (*foundPointer).second;
	}

	void put(int key, int value)
	{
		m_impl[key] = value;
	}

	bool exists(int key) const
	{
		CLOOKUP_MAP_INT_INT::const_iterator foundPointer = m_impl.find(key);
		if (foundPointer == m_impl.end())
		{
			return false;
		}
		else
		{
			return true;
		}
	}

private:
	CLOOKUP_MAP_INT_INT m_impl;
};


/**
 * CLookupAny<T>
 *   Lookup from void pointer to 'T'
 *   T must be a compact data structure that allows copying by the '=' operator
 *
 */

template< class T>
class CLookupAny
{
    public:
        friend class CIterator;

        /**
        * Iterator over CLookupAny: visit all <key_pointer, T_value> pairs
        *
        */
        class CIterator : private CLookupPtr<T >::CIterator
        {
            public:
                CIterator(const CLookupAny *lookup)
                    : CLookupPtr<T >::CIterator(&lookup->m_lookup_ptr), m_default_value(lookup->m_default_value)
                {
                }
                bool next(void*& key, T &value)
                {
                    T *value_ptr;
                    bool have_next = CLookupPtr<T >::CIterator::next(key, value_ptr);

                    if (value_ptr)
                    {
                        value = *value_ptr;
                    }
                    else
                    {
                        value = m_default_value;
                    }

                    return have_next;
                }

                void reset()
                {
                    CLookupPtr<T >::CIterator::reset();
                }

            private:
                T m_default_value;
        };

    public:
        //------------------------------------------
        //------------ CLookupAny ------------------
        //------------------------------------------

        // TODO: void erase(key)

        CLookupAny(T default_value) : m_lookup_ptr(true/*free pointers*/), m_default_value(default_value)
        {
        }

        const T get(void *key) const
        {
            T *ptr = m_lookup_ptr.get(key);
            if (!ptr)
            {
                ptr = new T;
                *ptr = m_default_value;
                CLookupPtr<T>* lookup_impl = const_cast<CLookupPtr<T>* >(&m_lookup_ptr);
                lookup_impl->put(key, ptr);
            }
            return *ptr;
        }

        void put(void *key, T value)
        {
            T *ptr = new T;
            *ptr = value;
            m_lookup_ptr.put(key, ptr);
        }

        bool exists(void *key) const
        {
            return m_lookup_ptr.exists(key);
        }

    private:
        CLookupPtr<T> m_lookup_ptr;
        T m_default_value;
};
#endif
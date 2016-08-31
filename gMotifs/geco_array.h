/**
 * @file geco_array.h
 * @author  Uberto Pozzoli  <uberto.pozzoli@bp.lnf.it>
 * @author  Matteo Cereda   <matteo.cereda@bp.lnf.it>
 * @version 0.1
 *
 * @section LICENSE
 * Copyright (C) 2010 by Uberto Pozzoli and Matteo Cereda
 * uberto.pozzoli@bp.lnf.it
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the
 * Free Software Foundation, Inc.,
 * 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 * @section DESCRIPTION
 * This file contains declarations gArray template class,
 * some classes and template classes that inherits form gArray.
 * It also contains abstract gArraRetriever and gArrayRetrieverImplementation classes:
 * for their use see documentation.
 */


#ifndef __GECO_ARRAY__
#define __GECO_ARRAY__

#include "geco_base.h"
#include <math.h>
#include <string>
#include <vector>
#include <memory.h>

using namespace std;

namespace geco {

/** @brief Array Template Class
 *
 * gArray is a template class used throughout the library to represent arrays of different types. It has been developed for two main reasons: to make subsetting efficient in terms of memory usage and to keep trace of Not Available values (NAs). Efficient subsetting is obtained though a reference counting mechanism so that when you instantiate a subset of an array the data are not duplicated but an internal reference is maintained to the original data with a specification of where the new instance starts and ends in the previously allolcate memory. Only when a value in the new object is changed, a new copy is allocated. NAs tracing is obtained by mean of a specilized gBitsArray class that manage bits arrays.
 * @tparam T the array elements type
 */
template <class T> class gArray {
    template <class U> friend class gArrayRetrieverImplementation;
    template <class U> friend class gArray;
    template <class U> friend class gMatrix;

private:
    gArrayInternal<T> *i_array;
    gPos               i_start;
    gPos               i_end;
    void acquirecontent(gArrayInternal<T> * content);
    void releasecontent();
    gArray(T* data,gSize length,const gArray<gPos> & isNaFlags,bool sorted);

protected:
    T * getData() const;
    gBytes * getNAs() const;
    template <class C> gArray<T> operator () (const gArray<C> & rows,const gArray<C> & cols,gSize ncols) const;
    template <class C> gArray<T> & operator () (gArray<T> & res,const gArray<C> & rows,const gArray<C> & cols,gSize ncols) const;

public:
    /** @name Public constructors & destructors */
    //@{
    gArray();
    gArray(const T* data,gSize length,const gArray<gPos> & NaPositions,bool sorted=false);
    gArray(T value,gSize length=1,bool isna=false,bool init=true);
    gArray(const gArray<T> & array,gPos start=0,gPos end=0);
    ~gArray();
    //@}

    /** @name Array information members */
    //@{
    gSize              NACount(gPos start=0,gPos end=0) const;
    gSize              getSize() const;
    //@}

    /** @name Array value setters */
    //@{
    void setValue(gPos pos, T value,bool isna=false,T fillvalue = (T) 0,bool fillNA=true);
    void setValue(gPos destpos,const gArray<T> array,gPos sourcepos,T fillvalue = (T) 0,bool fillNA=true);
    void setValues(const geco::gArray< geco::gPos >& positions, const gArray< T >& values, T fillvalue = (T) 0, bool fillNA=true);
    void setData(gSize length,const T* data,const gArray<gPos> & napos=gArray<gPos>(),bool sorted=false);
    void ownData(gSize length,T* data,const gArray<gPos> & napos=gArray<gPos>(),bool sorted=false);
    void setNoNa();
    void setAllNa();
    gArray<T> & operator = (const gArray<T> & array);
    gArray<T> & operator = (T value);
    //@}

    /** @name Elements access operators and function members */
    //@{
    const T            operator [] (gPos pos) const;
    gBool               isNA(gPos pos) const;
    gArray<gBool>       isNA() const;    
    template <class C> gArray<T> operator [] (const gArray<C> & positions) const;
    //@}

    /** @name Arithmetic Logic and comparison atomic operators */
    //@{
    gArray<T>      operator + (const gArray<T> & array) const;
    gArray<T>      operator - (const gArray<T> & array) const;
    gArray<T>      operator * (const gArray<T> & array) const;
    gArray<T>      operator / (const gArray<T> & array) const;
    gArray<gBool> operator <  (const gArray<T> & array) const;
    gArray<gBool> operator <= (const gArray<T> & array) const;
    gArray<gBool> operator == (const gArray<T> & array) const;
    gArray<gBool> operator >= (const gArray<T> & array) const;
    gArray<gBool> operator >  (const gArray<T> & array) const;
    gArray<gBool> operator != (const gArray<T> & array) const;
    gArray<gBool> operator || (const gArray<T> & array) const;
    gArray<gBool> operator && (const gArray<T> & array) const;
    gArray<gBool> operator ! () const;
    //@}

    /** @name Arithemtic compound assignment atomic operators */
    //@{
    gArray<T> & operator += (const gArray<T> & array);
    gArray<T> & operator -= (const gArray<T> & array);
    gArray<T> & operator *= (const gArray<T> & array);
    gArray<T> & operator /= (const gArray<T> & array);
    //@}

    /** @name Transformation function members */
    //@{
    gArray<T> & concatenate(const gArray<T> & array,gPos astart=0,gPos aend=0);
    gArray<T> & replace(gPos pos,const gArray<T> & array,gPos astart=0,gPos aend=0);
    gArray<T> & insert(gPos pos,const gArray<T> & array,gPos astart=0,gPos aend=0);
    gArray<T> & remove(gPos pos,gSize npos);
    template <class C> operator gArray<C> () const;
    gArray<T> & revert(gPos start=0, gPos end=0);
    gArray<gPos> sort(bool index_return=false);
    //@}

    /** @name Sorting,searching and counting function members */
    //@{
    gArray<gPos> find(T value,gPos start=0,gPos end=0,bool excludeNA=true) const;
    gArray<gPos> match(const gArray<T> & values,gPos start=0,gPos end=0) const;    
    gArray<T> getUnique(gPos start=0,gPos end=0) const;
    gArray<T> getReverted(gPos start=0,gPos end=0) const;
    gArray<T> & getReverted(gArray<T> & res,gPos start=0,gPos end=0) const;
    gArray<T> getSorted(gPos start=0,gPos end=0) const;
    gArray<T> & getSorted(gArray<T> & res,gPos start=0,gPos end=0) const;
    gArray<gSize> getCounts(const gArray<T> values,bool countall=false,gPos start=0,gPos end=0,bool skipnan=true) const;
    gArray<gSize> & getCounts(gArray<gSize> & res,const gArray<T> values,bool countall=false,gPos start=0,gPos end=0,bool skipnan=true) const;
    gArray<gSize> getCounts(T value,gSize winlength,gSize pos=0,gPos start=0,gPos end=0,bool skipnan=true) const;
    gArray<gSize> & getCounts(gArray<gSize> & res,T value,gSize winlength,gSize pos=0,gPos start=0,gPos end=0,bool skipnan=true) const;
    //@}

    /** @name Whole array comparison members */
    //@{
    bool equals(const gArray<T> & array) const;
    //@}

    /** @name Calculation function members */
    //@{
    gArray<T> getMin(gPos start=0,gPos end=0,bool skipnan=true,gSize winlength=0,gPos winpos=0) const;
    gArray<T> & getMin(gArray<T> & res,gPos start=0,gPos end=0,bool skipnan=true,gSize winlength=0,gPos winpos=0) const;
    gArray<T> getMax(gPos start=0,gPos end=0,bool skipnan=true,gSize winlength=0,gPos winpos=0) const;
    gArray<T> & getMax(gArray<T> & res, gPos start=0,gPos end=0,bool skipnan=true,gSize winlength=0,gPos winpos=0) const;
    gArray<T> getSum(gPos start=0,gPos end=0,bool skipnan=true,gSize winlength=0,gPos winpos=0) const;
    gArray<T> & getSum(gArray<T> & res, gPos start=0,gPos end=0,bool skipnan=true,gSize winlength=0,gPos winpos=0) const;
    gArray<T> getSquareSum(gPos start=0,gPos end=0,bool skipnan=true,gSize winlength=0,gPos winpos=0) const;
    gArray<T> & getSquareSum(gArray<T> & res, gPos start=0,gPos end=0,bool skipnan=true,gSize winlength=0,gPos winpos=0) const;
    gArray<gScore> getMean(gPos start=0,gPos end=0,bool skipnan=true,gSize winlength=0,gPos winpos=0) const;
    gArray<gScore> & getMean(gArray<gScore> & res, gPos start=0,gPos end=0,bool skipnan=true,gSize winlength=0,gPos winpos=0) const;
    gArray<gScore> getStdDev(bool unbiased=true,gPos start=0,gPos end=0,bool skipnan=true,gSize winlength=0,gPos winpos=0) const;
    gArray<gScore> & getStdDev(gArray<gScore> & res, bool unbiased=true,gPos start=0,gPos end=0,bool skipnan=true,gSize winlength=0,gPos winpos=0) const;
    //@}
};

/** @brief Array retriever implementation template class
 *
 * From this class retrievers implementation must be drived to gain access to the
 * internal data and reference counts of arrays
 * It is used when implementing new sequence or feature retrievers implementation.
 */
template <class T> class gArrayRetrieverImplementation:public gRetrieverImplementation {
  friend class gArrayRetriever<T>;
private:
  gSize i_nFeatures;
protected:
    template <class C> const C * getArrayDataPtr(const gArray<C> & array) const;
    template <class C> C * getArrayDataPtr(gArray<C> & array) const;
    unsigned getArrayReferenceCount(gArray<T> & array) const;
public:
    gArrayRetrieverImplementation(gSize nFeatures=1);
    virtual ~gArrayRetrieverImplementation();
};

/** @brief Array retriever template class
 *
 * This class is a container for array retrievers implementation
 * From this you may derive new kind or retrievers
 */
template <class T> class gArrayRetriever:public gRetriever {
public:
    gArrayRetriever();
    gArrayRetriever(const gArrayRetrieverImplementation<T> & implementation);
    gArrayRetriever(const gArrayRetriever<T> & retriever);
    ~gArrayRetriever();
    gSize getFeaturesCount() const;
};

/** @brief Matrix Class
 *
 * Array Class Template specilization to efficiently manage matrices.
 */
template<class T> class gMatrix:public gArray<T> {
private:
    gSize nrows;
    gSize ncols;
public:
    /** @name Contructors */
    //@{
    gMatrix();
    gMatrix( gSize rows,gSize cols,const T *data);
    gMatrix( gSize rows,gSize cols,const gArray<T> & data);
    gMatrix( gSize rows,gSize cols,T fillvalue,bool fillNA=false);
    gMatrix(gSize rows,gSize cols);
    gMatrix( const gMatrix<T> & matrix );
    //@}

    /** @name Matrix values setters */
    //@{
    void setValue( gPos row,gPos col,T value,bool isna=false);
    void setValues( const gArray<gPos> & rows,const gArray<gPos> & cols,const gArray<T> & values);
    void setData(gSize rows, gSize cols, const T* data, const gArray<gPos> & napos);
    void setRow( gPos row,const gArray<T> & values );
    void setCol( gPos col,const gArray<T> & values );
    gMatrix<T> & operator = ( const gMatrix<T> & matrix );
    //@}

    /** @name Matrix information function members */
    //@{
    gSize            getRowsNum() const;
    gSize            getColsNum() const;
    //@}

    /** @name Matrix values access operators and function members */
    //@{
    gBool isNA(gPos row,gPos col) const;
    gArray<T> getRow(gPos rownum) const;
    gArray<T> getCol(gPos colnum) const;
    T operator () ( gPos row,gPos col ) const;
    template <class C> gArray<T> & operator () (gArray<T> & res,const gArray<C> & rows, const gArray<C> & cols) const;
    template <class C> gArray<T> operator () (const gArray<C> & rows, const gArray<C> & cols) const;
    //@}
};

/** @brief String Class
 *
 * Array Class template specialization to efficiently manage strings.
 */
class gString:public gArray<gChar> {
    //friend std::ostream & operator << (std::ostream & sout,const gString & str);
public:
    /** @name Constructors */
    //@{
    gString();
    gString(const gArray<gChar> & str,gPos start=0,gPos end=0);
    gString(const gChar * str,gPos start=0,gPos end=0);
    gString ( const std::string & str,gPos start=0,gPos end=0);
    gString(gSize length,gChar symbol='n');
    //@}

    /** @name Cast operators */
    //@{
    operator const std::string () const;
    //@}

    /** @name String information function members */
    //@{
    gSize getLength() const;
    //@}

    /** @name String manipulation function members (self) */
    //@{
    gString & lower();
    gString & upper();
    gString & operator =  (const std::string & str);
    gString & operator =  (const gChar *str);
    gString & operator += (const gString & str);
    //@}

    /** @name String manipulation function members (instatiate new objects) */
    //@{
    gString   operator +  (const gString & str) const;
    std::vector<gString> &  split(gChar splitChar,std::vector<gString> & parts) const;
    gString getLower() const;
    gString getUpper() const;
    //@}

    /** @name String comparison operators and function members */
    //@{
    bool operator == (const char * str) const;
    bool operator == (const gString & str) const;
    bool operator == (const std::string & str) const;

    bool operator < (const char * str) const;
    bool operator < (const gString & str) const;
    bool operator < (const std::string & str) const;
    //@}

    /** @name String output operators and function members */
    //@{
    std::ostream & out(std::ostream &stream) const;
    //@}
};

/** @brief Interval Class
 *
 * Positions intervals class template to represent array intervals.
 * This is a template but only works (and is intended to) meaningfully with integer types
 */
template <class T> class gPositionInterval:private gArray<T> {
public:
    /** @name Constructors and assignment */
    //@{
    gPositionInterval(T start =  0, T end =  0, bool validStart = true, bool validEnd = true);
    gPositionInterval(const gArray<T> & positions);
    gPositionInterval(const gPositionInterval<T> & positionInterval, T startOffset= 0, T endOffset= 0);
    gPositionInterval<T> & operator = (const gPositionInterval<T> & positionInterval);
    //@}

    /** @name Interval information */
    //@{
    gSize getLength() const;
    T getStart() const;
    T getEnd() const;
    bool validStart() const;
    bool validEnd() const;
    gArray<T> getPositions() const;
    gArray<bool> contains(const gArray<T> & positions) const;
    gBool operator == (const gPositionInterval<T> & interval) const;
    gBool operator > (const gPositionInterval<T> & interval) const;
    gBool operator < (const gPositionInterval<T> & interval) const;
    //@}

    /** @name Interval operations */
    //@{
    gPositionInterval<T> getIntersection(const gPositionInterval<T> & interval) const;
    gPositionInterval<T> getLeftDiff(const gPositionInterval<T>  & interval) const ;
    gPositionInterval<T> getRightDiff(const gPositionInterval<T> & interval) const;
    gPositionInterval<T> getUnion(const gPositionInterval<T> & interval) const;
    //@}
};

/** @typedef gPositionInterval<gPos> gReferenceInterval;
 *  @brief A specialization of gPositionInterval to represent genomic reference intervals. 
 */
typedef gPositionInterval<gPos> gReferenceInterval;

/** @typedef gPositionInterval<gRelativePos> gElementInterval;
 *  @brief A specialization of gPositionInterval to represent genomic element intervals. 
 */
typedef gPositionInterval<gRelativePos> gElementInterval;

/**  @brief Class to manage flag conditions
 *
 * This template class can be 
 */
template<class T>
class gCondition{
 private:
   gArray<unsigned char> i_sbits;
   gArray<unsigned char> i_ubits;
 public:
   gCondition(){
    i_sbits=gArray<unsigned char>((unsigned char) 0);
    i_ubits=gArray<unsigned char>((unsigned char) 0);
   };
   
   gCondition(T v){
    i_sbits=gArray<unsigned char>((unsigned char) v);
    i_ubits=gArray<unsigned char>((unsigned char) 0);
   };

   gCondition operator ~ (){
     gCondition<T> ret=*this;
     gArray<unsigned char> s;
     gArray<unsigned char> u;
     s.setValue(0,i_ubits[0],false);
     u.setValue(0,i_sbits[0],false);
     for(gPos i=1;i<i_sbits.getSize();i++){
        s.setValue(0,s[0] | i_ubits[i],false);
        u.setValue(0,u[0] | i_sbits[i],false);
     }
     ret.i_sbits=s;
     ret.i_ubits=u;
     return ret;
   }

   gCondition operator & (const gCondition & o){
     if (i_sbits.getSize()>=o.i_sbits.getSize()) {
       gCondition<T> ret=*this;
       for (gPos i=0;i<i_sbits.getSize();i++) {
         for (gPos j=0;j<o.i_sbits.getSize();j++) {
           ret.i_sbits.setValue(i,i_sbits[i] | o.i_sbits[j],false);
           ret.i_ubits.setValue(i,i_ubits[i] | o.i_ubits[j],false);
         }
       }
       return ret;
     } else {
       gCondition<T> ret=o;
       for (gPos i=0;i<i_sbits.getSize();i++) {
         for (gPos j=0;j<o.i_sbits.getSize();j++) {
           ret.i_sbits.setValue(i,i_sbits[i] | o.i_sbits[j],false);
           ret.i_ubits.setValue(i,i_ubits[i] | o.i_ubits[j],false);
         }
       }
       return ret;
     }
   }

   gCondition operator | (const gCondition & o){
     gCondition<T> ret=*this;
     ret.i_sbits.concatenate(o.i_sbits);
     ret.i_ubits.concatenate(o.i_ubits);
     return ret;
   }

  gArray<gPos> select(const gArray<unsigned char> & flags) const{
    gArray<gPos> ret;
    for(gPos j=0;j<flags.getSize();j++){
      for(gPos i=0;i<i_sbits.getSize();i++){
        unsigned char bits_required=i_sbits[i]|i_ubits[i];
        if(((flags[j] ^ i_ubits[i]) & bits_required)==bits_required){
          ret.concatenate(gArray<gPos>(j,1,false));
          break;
        }
      }
    }
    return ret;
  }
};



//--------------------------------------------------------------------------
// Template helper functions definitions
//--------------------------------------------------------------------------
/** @brief Monotonic array creation utility.
 *
 * This template function returns an array of ordered equispaced(step) values of type T
 * starting from * firstValue upt to lastValue firstvalue can be greater or lesser than
 * lastvalue in the latter case the array values will be monotonic decreasing
 * @param firstValue first value of the array
 * @param lastValue  last value
 * @param step absolute separation between adjacent values (must be positive)
 * @return a  gArray<T> object
 * @ingroup arrays
 */
template <class T> gArray<T> getArray(T firstValue,T lastValue,T step) {
    gSize nelm;
    T * data;
    if (lastValue>firstValue) {
        nelm = (gSize) floor(((gScore)lastValue-(gScore)firstValue)/(gScore)step)+1;
        data=new T[nelm];
        data[0]=firstValue;
        T *act=data+1;
        T *prev=data;
        for (gPos i=1;i<nelm;i++) {
            *act=*prev + step;
            prev++;
            act++;
        }
    } else {
        nelm= (gSize) floor(((gScore)firstValue-(gScore)lastValue)/(gScore)step)+1;
        data=new T[nelm];
        data[0]=firstValue;
        T *act=data+1;
        T *prev=data;
        for (gPos i=1;i<nelm;i++) {
            *act=*prev - step;
            prev++;
            act++;
        }
    }
    gArray<T> ret;
    ret.ownData(nelm,data,gArray<gPos>(0,0,false),false);
    return ret;
}

/** @brief Atomic array function application
 *
 * This template function apply the provided function to the elements of the passed array
 * @param array an array containing elements to which apply the function
 * @param function the function to be applied
 * @return a gArray<T> containing the results.
 * @ingroup arrays
 */
template <class T> gArray<T> apply(const gArray<T> &array,T (*function)(T val)) {
    gArray<T> ret(array.getSize(),0,true,false);
    T ares;
    for (gSize i=0;i<array.getSize();i++) {
        if (!array.isNA(i)) {
            try {
                ares=function(array[i]);
                if (ares!=ares) ret.setValue(i,function(array[i]),true);
                else ret.setValue(i,function(array[i]),false);
            } catch (...) {
                ret.setValue(i,0,true);
            }
        } else ret.setValue(i,0,true);
    }
    return ret;
}
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
// Non-template helper functions declarations
//--------------------------------------------------------------------------
gArray<gPos> which ( const gArray<gBool> & condition,bool includeNA=false );
std::ostream & operator << (std::ostream & sout,const gString & str);
std::istream & operator >> (std::istream & sin,gString & str);
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
// template class gArray definition
//--------------------------------------------------------------------------
/** @brief Default constructor
 *
 * The default constructor instantiates a zero length array.
 */
template <class T> gArray<T>::gArray() {
    i_array=NULL;
    i_start=0;
    i_end=0;
}

/** @brief Raw data constructor
 *
 * This constructor instantiates an array of the given length initialized using
 * the values prvided in data. Positions provided in NApositions are marked
 * as NA.
 * @param data const pointer used for initialization. It must point to an allocated block of size length.
 * @param length the length of the array.
 * @param NaPositions a gArray containing the positions of NA values (defaults to an empty array)
 * @param sorted to be set if data contains sorted values (defaults to false)
 */
template <class T> gArray<T>::gArray ( const T* data,gSize length,const gArray<gPos> & NaPositions,bool sorted ) {
    i_array=NULL;
    i_start=0;
    i_end=0;
    setData ( length,data,NaPositions,sorted );
}

/** @brief Single value constructor
 * 
 * This constructor instatiates an arry of size length. If init is true all elements are initialized
 * using the value provied and marked as NA if isan is true.
 *
 * @param value initializing value
 * @param length length of the array (defaults to 1)
 * @param isna initializing NA flag (dfaults to false)
 * @param init wheter or not to initialize the array using value and isna (defaults to true)
 */
template <class T> gArray<T>::gArray ( T value,gSize length,bool isna,bool init ) {
    i_array=NULL;
    i_start=0;
    i_end=0;
    if ( init ) acquirecontent ( new gArrayInternal<T> ( length,value,isna ) );
    else acquirecontent ( new gArrayInternal<T> ( length ) );
}

/** @brief Range Copy constructor
 *
 * This constructor instantiate an array as a subarray og the provided one between start and end.
 * @param array array to copy from
 * @param start zero base copy start position
 * @param end one based copy end position
 */
template <class T> gArray<T>::gArray ( const gArray<T> & array,gPos start,gPos end ) {
    gPos aend= ( end==0 ) ? ( array.getSize() ) : ( end );
    i_array=NULL;
    i_start=0;
    i_end=0;
    if ( ( aend<=array.getSize() ) && ( start<aend ) ) {
        acquirecontent ( array.i_array );
        i_start=array.i_start+start;
        i_end=array.i_start+aend;
    } else if ( end>0 ) throw ( gException ( "gArray: invalid range" ) );
}

/** @brief Destructor */
template <class T> gArray<T>::~gArray() {
    releasecontent();
}

/** @brief Number of NA values in the array
 *
 * This function returns the number of elemnts in the specified range foe which NA flag is set to true.
 * @param start the beginning of the range (0 based)
 * @param end the end of the range (1 based)
 * @return teh number of NA
 */
template <class T> gSize gArray<T>::NACount ( gPos start, gPos end ) const {
    if ( !i_array ) throw gException ( "operation on object not yet initialized" );
    gPos aend= ( end==0 ) ? ( getSize() ) : ( end );
    return i_array->NACount ( i_start+start,i_start+aend );
}

/** @brief Size of the array
 *
 * This function returns the size of the array
 * @return the size of the array
 */
template <class T> gSize gArray<T>::getSize() const {
    return i_end-i_start;
}

/** @brief Set a value in the array
 *
 * Sets the value of the array element indexd by pos. if pos is greater than array size
 * expands the array accordingly filling the new elements with the values and na specified
 * @param pos position of the element to set
 * @param value value of the element to set
 * @param isna if the element shoud be NA (defalut to false)
 * @param fillvalue value to use to fill value when array  is resized
 * @param fillNA NA status to use whena rray is resized
 */
template <class T> void gArray<T>::setValue ( gPos pos, T value,bool isna,T fillvalue,bool fillNA ) {
    if ( !i_array ) acquirecontent ( new gArrayInternal<T> ( pos+1,fillvalue,fillNA ) );
    if ( i_array->i_refcount>1 ) acquirecontent ( new gArrayInternal<T> ( *i_array,i_start,i_end,fillvalue,fillNA ) );
    i_array->setValue ( i_start+pos,value,isna,fillvalue,fillNA );
    i_end=max ( i_end,i_start+pos+1 );
}

/** @brief Sets values from another array
 *
 * This method copies the values from an array.
 * @param destpos Destination position to wich copy the values from array.
 * @param array Array from which values are taken.
 * @param sourcepos Source position of the first value to be copied (to the end)
 * @param fillvalue Value used to initialize newly allocated positions when needed (defaults to 0)
 * @param fillNA Value used to initialize newly allocated positions NA status when needed (defaults to true)
 */
template <class T> void gArray<T>::setValue ( gPos destpos,const gArray<T> array,gPos sourcepos,T fillvalue,bool fillNA ) {
    if ( !i_array ) acquirecontent ( new gArrayInternal<T> ( destpos+1,fillvalue,fillNA ) );
    if ( i_array->i_refcount>1 ) acquirecontent ( new gArrayInternal<T> ( *i_array,i_start,i_end,fillvalue,fillNA ) );
    i_array->setValue ( i_start+destpos,array[sourcepos],array.isNA ( sourcepos ),fillvalue,fillNA );
    i_end=max ( i_end,i_start+destpos+1 );
}

/** @brief Sets values at given positions
 *
 * This method copies the values contained in an array to the positions provided. Number of
 * values must match number of positions.
 * @param positions gPos array containig the target positions
 * @param values array containig the target values
 * @param fillvalue Value used to initialize newly allocated positions when needed (defaults to 0)
 * @param fillNA Value used to initialize newly allocated positions NA status when needed (defaults to true)
 */
template <class T> void gArray<T>::setValues ( const gArray<gPos> & positions, const gArray<T> & values, T fillvalue,bool fillNA ) {
    if (positions.getSize()>0) {
        gPos  mpos=positions.getMax ( 0,0,true ) [0];
        if ( !i_array ) acquirecontent ( new gArrayInternal<T> ( mpos+1,fillvalue,fillNA ) );
        if ( i_array->i_refcount>1 ) acquirecontent ( new gArrayInternal<T> ( *i_array,i_start,i_end,fillvalue,fillNA ) );
        i_array->setValues(i_start,positions.i_start,positions.i_end,*positions.i_array,values.i_start,values.i_end,*values.i_array,fillvalue,fillNA );
        i_end=max(i_end,i_start+mpos+1);
    }
}

/** @brief Copy values from raw data
 *
 * This function copy the values in a block of allocated memory pointed by data into this array.
 * @param length  length of the data provided
 * @param data pointer to the data (allocated memory must be at least of length length)
 * @param napos positions to be marked as NA (defaults to an empty set)
 * @param sorted whether raw data are sorted (defaults to false)
 */
template <class T> void gArray<T>::setData ( gSize length,const T* data,const gArray<gPos> & napos,bool sorted ) {
    if ( length>0 ) {
        T *values=new T[length];
        gBitsArray nans ( length,false );
        if ( data!=NULL ) memmove ( values,data,sizeof ( T ) *length );
        else memset ( values,0,sizeof ( T ) *length );
        for ( gPos i=0;i<napos.getSize();i++ ) {
            nans.setValue ( napos[i],true );
        }
        if ( ( !i_array ) || ( i_array->i_refcount>1 ) ) {
            acquirecontent ( new gArrayInternal<T> ( length,values,nans,true,sorted ) );
        } else {
            i_array->setData ( length,values,nans,0,true,true,sorted );
        }
        i_start=0;
        i_end=i_array->i_length;
    }
}

/** @brief Use values from raw data
 *
 * This method is similar to setData but it doesn't copy the values from data
 * instead it "owes" the pointer itself. It will be freed when the object will be destroied
 * any program calling this method suould not use anymore the pointer after this call.
 * @param length  length of the data provided
 * @param data pointer to the data (allocated memory must be at least of length length)
 * @param napos positions to be marked as NA (defaults to an empty set)
 * @param sorted whether raw data are sorted (defaults to false)
 */
template <class T> void gArray<T>::ownData ( gSize length,T* data,const gArray<gPos> & napos,bool sorted ) {
    if ( length>0 ) {
        gBitsArray nans ( length,false );
        for ( gPos i=0;i<napos.getSize();i++ ) nans.setValue ( napos[i],true );
        if ( ( !i_array ) || ( i_array->i_refcount>1 ) ) {
            acquirecontent ( new gArrayInternal<T> ( length,data,nans,true,sorted ) );
        } else {
            i_array->setData ( length,data,nans,0,true,true,sorted );
        }
        i_start=0;
        i_end=i_array->i_length;
    } else {
    }
}

/** @brief Set all values NA flag to false
 *
 * This method set all the NA flags to false
 */
template <class T> void gArray<T>::setNoNa() {
    if (i_array) i_array->i_nans.resetAll();
}

/** @brief Set all values NA flag to true
 *
 * This method set all the NA flags to true
 */
template <class T> void gArray<T>::setAllNa() {
    if (i_array) i_array->i_nans.setAll();
}

/** @brief Array assignment operator
 *
 * Makes this array a copy of the array provided
 * @param array array to be copied
 * @return a reference to this array
 */
template <class T> gArray<T> & gArray<T>::operator = ( const gArray<T> &array ) {
    acquirecontent ( array.i_array );
    i_start=array.i_start;
    i_end=array.i_end;
    return *this;
}

/** @brief Single value assignment operator
 *
 * Makes this array a single valued one, initialized to the value provided marked as !NA
 * @param value
 * @return a reference to this array
 */
template <class T> gArray<T> & gArray<T>::operator = ( T value ) {
    acquirecontent ( new gArrayInternal<T> ( 1,value,false ) );
    i_start=0;
    i_end=1;
    return *this;
}

/** @brief Single value access operator
 *
 * Returns the const value of the element at the specified position.
 * @param pos the required element position (0 based)
 * @return const value of the element
 */
template <class T> const T gArray<T>::operator [] ( gPos pos ) const {
    if ( i_array ) {
        return i_array->i_mem[pos+i_start];
    } else throw gException ( "invalid range" );
}

/** @brief Multiple value access template operator
  *
  * Returns an Array containing the values (and NA status) of the elements in this array at the
  * positions specified by another array. The values in the provided positions array are cast
  * to gPos before indexing. (this is intended to allow for using any type for indexing values)
  * @tparam the type of the array used for indexing
  * @param positions a gArray<C> of values to be interpreted as positions into this
  * @return a newly instatiated gArray containing the values
  */
template <class T> template <class C> gArray<T> gArray<T>::operator [] (const gArray<C> & positions) const {
    gArray<T> ret;
    if (positions.getSize()>0) {
        gArray<gSize> npos=(gArray<gSize>) positions;
        //gArray<gPos> bad=which(npos>=getSize());
        //for (gSize i=0;i<bad.getSize();i++) npos.setValue(bad[i],0,true);
        ret.acquirecontent(i_array->getValues(i_start,i_end,*npos.i_array,npos.i_start,npos.i_end));
    }
    return ret;
}

/** @brief Single value NA status retriever
 *
 * Returns the NA status of the element at a given position
 * @param pos the required element position (0 based)
 * @return true if NA false otherwise
 */
template <class T> gBool gArray<T>::isNA ( gPos pos ) const {
    if ( i_array ) {
        return i_array->isNA ( i_start+pos );
    } else throw gException ( "invalid range" );
}

/** @brief Multiple values NA status retriever
 *
 * Returns the NA status of the elements of the array 
 * @param pos the required element position (0 based)
 * @return a gArray<gbool> true if corresponding element is NA, false otherwise.
 */
template <class T> gArray<gBool> gArray<T>::isNA () const {
    gArray<gBool> ret;  
    if ( i_array ) {
        ret.acquirecontent( i_array->isNA(i_start,i_end ) );
    } 
    return ret;
}



/** @brief Atomic arithmetic operators
 * 
 * These operators apply arithemtic operators between the elements
 * of this array and those of the supplied one.
 * If array is longer than this, only the first this->getSize() are used.
 * Otherwise array is recycled.
 * In particular the division operator returns a new gArray<gScore> object
 * and the values will be cast to gScore before divison.
 * NAs are trated accordingly to the operator.
 * @param array the values.
 * @return a new array object with the results
 */
template <class T> gArray<T> gArray<T>::operator + ( const gArray<T> & array ) const {
    gArray<T> ret;
    if(i_array) ret.acquirecontent ( i_array->add ( *array.i_array,i_start,i_end,array.i_start,array.i_end ) );
    return ret;
}

template <class T> gArray<T> gArray<T>::operator - ( const gArray<T> & array ) const {
    gArray<T> ret;
    if(i_array) ret.acquirecontent ( i_array->subtract ( *array.i_array,i_start,i_end,array.i_start,array.i_end ) );
    return ret;
}

template <class T> gArray<T> gArray<T>::operator * ( const gArray<T> & array ) const {
    gArray<T> ret;
    if(i_array) ret.acquirecontent ( i_array->multiply ( *array.i_array,i_start,i_end,array.i_start,array.i_end ) );
    return ret;
}

template <class T> gArray<T> gArray<T>::operator / ( const gArray<T> & array ) const {
    gArray<T> ret;
    if(i_array) ret.acquirecontent ( i_array->divide ( *array.i_array,i_start,i_end,array.i_start,array.i_end ) );
    return ret;
}

/** @brief Atomic logic operators
 * 
 * These operators apply comparison or logical operators to the elements
 * of this array and those of the supplied one.
 * If array is longer than this, only the first this->getSize() are used.
 * Otherwise array is recycled.
 * They return a gArray<gBool> with values {0,1}
 * In particular the ! operator will return 1 for all those elements of this that are not 0s.
 * NA are trated accordingly to the operator (in particular for the || operator the result
 * will be !NA gicen that ONE of the two operand is !NA)
 * @param array the values.
 * @return a new array<gBool> object with the results
 */
template <class T> gArray<gBool> gArray<T>::operator < ( const gArray<T> & array ) const {
    gArray<gBool> ret;  
    if ( i_array ) ret.acquirecontent ( i_array->lesser ( *array.i_array,i_start,i_end,array.i_start,array.i_end ) );
    return ret;
}

template <class T> gArray<gBool> gArray<T>::operator <= ( const gArray<T> & array ) const {
    gArray<gBool> ret;
    if ( i_array ) ret.acquirecontent ( i_array->lesserEqual ( *array.i_array,i_start,i_end,array.i_start,array.i_end ) );
    return ret;
}

template <class T> gArray<gBool> gArray<T>::operator == ( const gArray<T> & array ) const {
    gArray<gBool> ret;  
    if ( i_array ){
      ret.acquirecontent ( i_array->equal ( *array.i_array,i_start,i_end,array.i_start,array.i_end ) );
    }
    return ret;

}

template <class T> gArray<gBool> gArray<T>::operator >= ( const gArray<T> & array ) const {
    gArray<gBool> ret;
    if ( i_array ) ret.acquirecontent ( i_array->greaterEqual ( *array.i_array,i_start,i_end,array.i_start,array.i_end ) );
    return ret;
}

template <class T> gArray<gBool> gArray<T>::operator > ( const gArray<T> & array ) const {
    gArray<gBool> ret;
    if ( i_array ) ret.acquirecontent ( i_array->greater ( *array.i_array,i_start,i_end,array.i_start,array.i_end ) );
    return ret;
}

template <class T> gArray<gBool> gArray<T>::operator != ( const gArray<T> & array ) const {
    gArray<gBool> ret;
    if ( i_array ) ret.acquirecontent ( i_array->notEqual ( *array.i_array,i_start,i_end,array.i_start,array.i_end ) );
    return ret;
}

template <class T> gArray<gBool> gArray<T>::operator || ( const gArray<T> & array ) const {
    gArray<gBool> ret;
    if ( i_array ) ret.acquirecontent ( i_array->logicalOr ( *array.i_array,i_start,i_end,array.i_start,array.i_end ) );
    return ret;
}

template <class T> gArray<gBool> gArray<T>::operator && ( const gArray<T> & array ) const {
    gArray<gBool> ret;
    if ( i_array ) ret.acquirecontent ( i_array->logicalAnd ( *array.i_array,i_start,i_end,array.i_start,array.i_end ) );
    return ret;
}

/** @brief Negation operator
 * 
 * This operator will return an array of length this->getSize() with values set to 1
 * for all those elements of this that are not 0s.
 * @return a new array<gBool> object with the results
 */
template <class T> gArray<gBool> gArray<T>::operator ! () const {
    gArray<gBool> ret;
    if ( i_array ) ret.acquirecontent ( i_array->negation ( i_start,i_end ) );
    return ret;
}

/** @brief Atomic arithmetic compound operators
 * 
 * These operators apply arithemtic operators on the elements of this array
 * and those of the supplied one. The results are stored in this array
 * If array is longer than this, only the first
 * this->getSize() are used. Otherwise array is recycled.
 * @param array the values to be operated.
 * @return a reference to this array
 */
template <class T> gArray<T> & gArray<T>::operator += ( const gArray<T> & array ) {
    if ( i_array ) {
        if ( array.i_array ) {
            if ( i_array->i_refcount>1 ) acquirecontent ( new gArrayInternal<T> ( *i_array,i_start,i_end,0,false ) );
            i_array->selfAdd ( *array.i_array,i_start,i_end,array.i_start,array.i_end );
        }
    }
    return *this;
}

template <class T> gArray<T> & gArray<T>::operator -= ( const gArray<T> & array ) {
    if ( i_array ) {
        if ( array.i_array ) {
            if ( i_array->i_refcount>1 ) acquirecontent ( new gArrayInternal<T> ( *i_array,i_start,i_end,0,false ) );
            i_array->selfSubtract ( *array.i_array,i_start,i_end,array.i_start,array.i_end );
        }
    }
    return *this;
}

template <class T> gArray<T> & gArray<T>::operator *= ( const gArray<T> & array ) {
    if ( i_array ) {
        if ( array.i_array ) {
            if ( i_array->i_refcount>1 ) acquirecontent ( new gArrayInternal<T> ( *i_array,i_start,i_end,0,false ) );
            i_array->selfMultiply ( *array.i_array,i_start,i_end,array.i_start,array.i_end );
        }
    }
    return *this;
}

template <class T> gArray<T> & gArray<T>::operator /= ( const gArray<T> & array ) {
    if ( i_array ) {
        if ( array.i_array ) {
            if ( i_array->i_refcount>1 ) acquirecontent ( new gArrayInternal<T> ( *i_array,i_start,i_end,0,false ) );
            i_array->selfDivide ( *array.i_array,i_start,i_end,array.i_start,array.i_end );
        }
    }
    return *this;
}

/** @brief Concatenate values
 *
 * Adds to the end of this the values from a range into the provided array
 * @param array gArray containing the values to add
 * @param astart starting position in the provided array (0 based,defaults to 0)
 * @param aend end position in the provided array (1 based, defaults to the last element)
 * @return a reference to this
 */
template <class T> gArray<T> & gArray<T>::concatenate ( const gArray<T> & array,gPos astart,gPos aend ) {
    if (array.getSize()>0) {
        gPos cend= ( aend==0 ) ? ( array.getSize() ) : ( aend );
        if (!i_array) * this=gArray<T>(array,astart,cend);
        else{
         if(getSize()>0) insert ( getSize()-1,array,astart,cend );
         else * this=gArray<T>(array,astart,cend);
        }
    }
    return *this;
}

/** @brief Replace values
 *
 * Replace values in this array starting from a given position with the values
 * from a range in a provided array. This array is resized if needed.
 * @param pos zero based position from wich
 * @param array array containing the values to replace
 * @param astart starting position in the provided array (0 based,defaults to 0)
 * @param aend end position in the provided array (1 based, defaults to the last element)
 * @return a reference to this
 */
template <class T> gArray<T> & gArray<T>::replace ( gPos pos,const gArray<T> & array,gPos astart,gPos aend ) {
    if ( !i_array ) throw gException ( "gArray operator /: undefined first operator" );
    if ( i_array->i_refcount>1 ) acquirecontent ( new gArrayInternal<T> ( *i_array,i_start,i_end,0,true ) );
    gPos cend= ( aend==0 ) ? ( array.getSize() ) : ( aend );
    if ( ( pos+ ( cend-astart ) ) >i_end ) cend=astart+i_end-pos;
    i_array->selfReplace ( i_start+pos,*array.i_array,array.i_start+astart,array.i_start+cend );
    return *this;
}

/** @brief Insert values
 *
 * Insert the values from a range in a provided into a given position in this arraty.
 * This array is resized.
 * @param pos zero based position to which insert the vales
 * @param array array containing the values to replace
 * @param astart starting position in the provided array (0 based,defaults to 0)
 * @param aend end position in the provided array (1 based, defaults to the last element)
 * @return a reference to this
 */
template <class T> gArray<T> & gArray<T>::insert ( gPos pos,const gArray<T> & array,gPos astart,gPos aend ) {
    if ( !i_array ) throw gException ( "gArray operator /: undefined first operator" );
    if ( i_array->i_refcount>1 ) acquirecontent ( new gArrayInternal<T> ( *i_array,i_start,i_end,(T) 0,true ) );
    i_array->selfInsert ( i_start+pos,*array.i_array,array.i_start+astart,array.i_start+ ( ( aend==0 ) ? ( array.getSize() ) : ( aend ) ) );
    i_start=0;
    i_end=i_array->getSize();
    return *this;
}

/** @brief Remove values
 *
 * Removes a given number of values starting from a provided position. The array is resized.
 * @param pos the position starting from which elements have to be removed
 * @param npos the number of elements to remove
 * @return a reference to this
 */
template <class T> gArray<T> & gArray<T>::remove ( gPos pos,gSize npos ) {
    if ( !i_array ) throw gException ( "gArray operator /: undefined first operator" );
    if ( i_array->i_refcount>1 ) acquirecontent ( new gArrayInternal<T> ( *i_array,i_start,i_end,0,true ) );
    i_array->selfRemove ( i_start+pos,npos );
    i_start=0;
    i_end=i_array->getSize();
    return *this;
}

/** @brief Array template type casting operator
 *
 * Allow for array type casting returning an object instatiated as gArray<C>
 * containing original values (C) values.
 * @tparam the type to which array elements must be cast
 * @return a gArray<C> object
 */
template <class T> template <class C> gArray<T>::operator gArray<C> () const {
    if (i_array) {
        gArray<C> ret(0,getSize(),false);
        ret.i_array->castFrom(*i_array,i_start,i_end);
        return ret;
    } else {
        gArray<C> ret;
        return ret;
    }
}

/** @brief Array elements order revert
 *
 * Revert the order of this array elements in the specified range
 * @param start first element to revert (0 based)
 * @param end last element to revert (1 based)
 * @return a reference to this
 */
template <class T> gArray<T> & gArray<T>::revert(gPos start, gPos end) {
    if ( i_array ) {
        if ( i_array->i_refcount>1 ) acquirecontent ( new gArrayInternal<T> ( *i_array,i_start,i_end,0,true ) );
        gPos aend= ( end==0 ) ? ( getSize() ) : ( end );
        i_array->selfRevert(i_start+start,i_start+aend);
    }
    return *this;
}

/** @brief Sort array elements
 *
 * Sorts the elements of this array in ascending order optionally returning their original positions.
 * @param index_return true if original positions of the elemenst should be returned
 * @return if index_return is true an array containing the positions of the elements before sorting, otherwise an empty array
 */
template <class T> gArray<gPos> gArray<T>::sort ( bool index_return ) {
    gArray<gPos> ret;
    if ( i_array->i_refcount>1 ) acquirecontent ( new gArrayInternal<T> ( *i_array,i_start,i_end,0,true ) );
    if (index_return) {
        ret.acquirecontent ( i_array->sort ( index_return ) );
    } else {
        i_array->sort ( index_return );
    }
    return ret;
}

/** @brief Find occurrences of a value
 *
 * Finds occurrences of the provide d value in a subset of this array
 * @param value thh value to search for
 * @param start the subset starting position (0 based, defaults to 0)
 * @param end the subset end position (1 based, defaults to the last element position)
 * @return an array containing occurrence positions
 */
template <class T> gArray<gPos> gArray<T>::find ( T value,gPos start,gPos end, bool excludeNA ) const {
    if ( !i_array ) return gArray<gPos>();
    if( getSize()==0 ) return gArray<gPos>();
    gPos aend= ( end==0 ) ? ( getSize() ) : ( end );
    if ( ( start>=aend ) || ( aend>getSize() ) ) throw gException ( "invalid range" );
    gArray<gPos> ret;
    ret.acquirecontent ( i_array->find ( value,i_start+start,i_start+aend,excludeNA) );
    return ret;
}

/** @brief Find first occurrence of values
 *
 * Finds first occurrences of the values contained in "values"
 * @param value an array of values 
 * @param start the subset starting position (0 based, defaults to 0)
 * @param end the subset end position (1 based, defaults to the last element position)
 * @return an array containing occurrence positions (NA if not present)
 */
template <class T> gArray<gPos> gArray<T>::match(const gArray<T> & values,gPos start,gPos end) const{
    if(values.getSize()==0) return gArray<gPos>();
    if ( (!i_array) || getSize()==0 ) return gArray<gPos>(0,values.getSize(),true);
    gPos aend= ( end==0 ) ? ( getSize() ) : ( end );
    if ( ( start>=aend ) || ( aend>getSize() ) ) throw gException ( "invalid range" );
    gArray<gPos> ret;
    ret.acquirecontent ( i_array->match ( * values.i_array ,i_start+start,i_start+aend ) );
    return ret;
}

/** @brief Array equality comparison
 *
 * Compares this array to the provided one returnintg true if they are identical
 * @param array the array to campare with
 * @return true if this is identical to the procided array, False otherwise
 */
template <class T> bool gArray<T>::equals ( const gArray<T> & array ) const {
    bool ret=false;
    if ( ( getSize() ==0 ) || ( array.getSize() ==0 ) ) {
        ret=getSize() ==array.getSize();
    } else if ( getSize() ==array.getSize() ) {
        gArray<gBool> ee=*this!=array;
        ret= ( ee.getSum ( 0,ee.getSize(),true ) [0]==0 );
    }
    return ret;
}

/** @brief Get unique elements
 *
 * Returns an array containing the a single copy of the
 * different values occurring in a specified subset of this array
 * @param start the subset starting position (0 based, defaults to 0)
 * @param end the subset end position (1 based, defaults to the last element position)
 * @return a newly instatiated array with the results.
 */
template <class T> gArray<T> gArray<T>::getUnique ( gPos start,gPos end ) const {
    if ( !i_array ) throw gException ( "operation on object not yet initialized" );
    gPos aend= ( end==0 ) ? ( getSize() ) : ( end );
    if ( ( start>=aend ) || ( aend>getSize() ) ) throw gException ( "invalid range" );
    gArray<T> ret;
    ret.acquirecontent ( i_array->getUnique ( i_start+start,i_start+aend ) );
    return ret;
}

/**  @brief Revert elements order
 *
 * Returns a reverted copy of the specified subset from the current array
 * @param start the subset starting position (0 based, defaults to 0)
 * @param end the subset end position (1 based, defaults to the last element position)
 * @return a newly instatiated array with the results.
 */
template <class T> gArray<T> gArray<T>::getReverted ( gPos start,gPos end ) const {
    if ( !i_array ) return gArray<T>();
    gArray<T> ret ( *this,start,end );
    ret.revert();
    return ret;
}

/**  @brief Revert elements order
 *
 * Returns a reverted copy of the specified subset from the current array
 * @param[out] res an array to be filled with the results
 * @param start the subset starting position (0 based, defaults to 0)
 * @param end the subset end position (1 based, defaults to the last element position)
 * @return a reference to the passed reults array
 */
template <class T> gArray<T> & gArray<T>::getReverted ( gArray<T> & res,gPos start,gPos end ) const {
    if ( !i_array ) throw gException ( "operation on object not yet initialized" );
    res=gArray<T>( *this,start,end );
    res.revert();
    return res;
}

/**  @brief Sort elements
    *
    * Returns a sorted copy of a specified subset of the current array
    * @param start the subset starting position (0 based, defaults to 0)
    * @param end the subset end position (1 based, defaults to the last element position)
    * @return a newly instatiated array with the results.
    */
template <class T> gArray<T> gArray<T>::getSorted ( gPos start,gPos end ) const {
    if ( !i_array ) throw gException ( "operation on object not yet initialized" );
    gPos aend= ( end==0 ) ? ( getSize() ) : ( end );
    if ( ( start>=aend ) || ( aend>getSize() ) ) throw gException ( "invalid range" );
    gArray<T> ret ( *this,start,aend );
    ret.sort();
    return ret;
}

/**  @brief Sort elements
    *
    * Returns a sorted copy of a specified subset of the current array
    * @param[out] res an array to be filled with the results
    * @param start the subset starting position (0 based, defaults to 0)
    * @param end the subset end position (1 based, defaults to the last element position)
    * @return a reference to the passed results array
    */
template <class T> gArray<T> & gArray<T>::getSorted (gArray<T> & res, gPos start,gPos end ) const {
    if ( !i_array ) throw gException ( "operation on object not yet initialized" );
    gPos aend= ( end==0 ) ? ( getSize() ) : ( end );
    if ( ( start>=aend ) || ( aend>getSize() ) ) throw gException ( "invalid range" );
    res=gArray<T>(*this,start,aend);
    res.sort();
    return res;
}

/** @brief Counts occurrences of the values in an array
    *
    * Counts the occurrences of the values from the provided array in a specified subset of this.
    * If caountall is true a unique count is returned counting all occurences of all elements in teh passed array.
    * Otherwise each elements in values is counted separately.
    * if skipnan is true NA elements are not counted (regardless their value), if it is false the result is NA.
    * @param values The values to count
    * @param countall wheter cont values as a whole or one by one
    * @param start the subset starting position (0 based, defaults to 0)
    * @param end the subset end position (1 based, defaults to the last element position)
    * @param skipnan whether to ignore NA values or not (defaults to true)
    * @return a newly instatiated array with the results.
    */
template <class T> gArray<gSize> gArray<T>::getCounts ( const gArray<T> values,bool countall,gPos start,gPos end,bool skipnan ) const {
    if ( !i_array ) throw gException ( "operation on object not yet initialized" );
    gPos aend= ( end==0 ) ? ( getSize() ) : ( end );
    if ( ( start>=aend ) || ( aend>getSize() ) ) throw gException ( "invalid range" );
    gArray<gSize> ret;
    ret.acquirecontent ( i_array->getCounts ( new gArrayInternal<gSize>((countall)?(1):(values.i_end-values.i_start),0,false),*values.i_array,i_start+start,i_start+aend,values.i_start,values.i_end,skipnan,countall ) );
    return ret;
}

/** @brief Counts occurrences of the values in an array
    *
    * Counts the occurrences of the values from the provided array in a specified subset of this.
    * If caountall is true a unique count is returned counting all occurences of all elements in teh passed array.
    * Otherwise each elements in values is counted separately.
    * if skipnan is true NA elements are not counted (regardless their value), if it is false the result is NA.
    * @param[out] res an array to be filled with the results
    * @param values The values to count
    * @param countall wheter cont values as a whole or one by one
    * @param start the subset starting position (0 based, defaults to 0)
    * @param end the subset end position (1 based, defaults to the last element position)
    * @param skipnan whether to ignore NA values or not (defaults to true)
    * @return a reference to the passed results array
    */
template <class T> gArray<gSize> & gArray<T>::getCounts (gArray<gSize> & res, const gArray<T> values,bool countall,gPos start,gPos end,bool skipnan ) const {
    if ( !i_array ) throw gException ( "operation on object not yet initialized" );
    gPos aend= ( end==0 ) ? ( getSize() ) : ( end );
    if ( ( start>=aend ) || ( aend>getSize() ) ) throw gException ( "invalid range" );
    i_array->getCounts(res.i_array,*values.i_array,i_start+start,i_start+aend,values.i_start,values.i_end,skipnan,countall );
    return res;
}

/** @brief Counts a value occurences in sliding windows
 *
 * Returns occurrence counts of the specified value in windows of
 * specified length sliding through the specified subset of this array.
 * Results are returned in an array of the same length as the subset.
 * The pos parameters specify the offset of the values returned in relation
 * to the window.
 * If skipnan is false a window count is not NA only if there are no NA values in the window.
 * Otherwise NA values in the window are simply ignored.
 * @param value the value to count for
 * @param winlength the length of the sliding windows
 * @param pos the position along the window to which count will be saved
 * @param start the subset starting position (0 based, defaults to 0)
 * @param end the subset end position (1 based, defaults to the last element position)
 * @param skipnan whether to ignore NA values or not (defaults to true)
 * @return an newly instatiated array with the results
 */
template <class T> gArray<gSize> gArray<T>::getCounts ( T value,gSize winlength,gSize pos,gPos start,gPos end,bool skipnan ) const {
    if ( !i_array ) throw gException ( "operation on object not yet initialized" );
    if ( winlength<=0 ) throw gException ( "invalid window size" );
    gPos aend= ( end==0 ) ? ( getSize() ) : ( end );
    if ( ( start>=aend ) || ( aend>getSize() ) ) throw gException ( "invalid range" );
    if ( winlength>getSize() ) throw gException ( "window too large" );
    if ( pos>=getSize() || pos>winlength ) throw gException ( "invalid position" );
    gArray<gSize> ret;
    gSize rlen=(winlength>0)?(aend-start):(1);
    ret.acquirecontent ( i_array->getWinCount (new gArrayInternal<gSize>(rlen,0,true), value,winlength,pos,i_start+start,i_start+aend,skipnan ) );
    return ret;
}

/** @brief Counts a value occurences in sliding windows
 *
 * Returns occurrence counts of the specified value in windows of
 * specified length sliding through the specified subset of this array.
 * Results are returned in an array of the same length as the subset.
 * The pos parameters specify the offset of the values returned in relation
 * to the window.
 * If skipnan is false a window count is not NA only if there are no NA values in the window.
 * Otherwise NA values in the window are simply ignored.
 * @param[out] res an array to be filled with the results
 * @param value the value to count for
 * @param winlength the length of the sliding windows
 * @param pos the position along the window to which count will be saved
 * @param start the subset starting position (0 based, defaults to 0)
 * @param end the subset end position (1 based, defaults to the last element position)
 * @param skipnan whether to ignore NA values or not (defaults to true)
 * @return a reference to the passed results array
 */
template <class T> gArray<gSize> & gArray<T>::getCounts (gArray<gSize> & res, T value,gSize winlength,gSize pos,gPos start,gPos end,bool skipnan ) const {
    if ( !i_array ) throw gException ( "operation on object not yet initialized" );
    if ( winlength<=0 ) throw gException ( "invalid window size" );
    gPos aend= ( end==0 ) ? ( getSize() ) : ( end );
    if ( ( start>=aend ) || ( aend>getSize() ) ) throw gException ( "invalid range" );
    if ( winlength>getSize() ) throw gException ( "window too large" );
    if ( pos>winlength ) throw gException ( "invalid offset" );
    if ((winlength>0)&&(res.getSize()!=(aend-start))) throw gException ( "Result array, invalid size" );
    i_array->getWinCount (res.i_array,value,winlength,pos,i_start+start,i_start+aend,skipnan );
    return res;
}

/** @brief Minimum value(s)
 *
 * Finds the minimum value(s) in a subset of this array or in windows sliding through the susbset.
 * If winlength==0 returns a single value gArray with the result. Otherwise returns
 * results calculated in windows of specified length sliding through the specified subset.
 * In this case results are returned in an array of the same length as the subset.
 * The pos parameters specify the offset of the values returned in relation
 * to the window.
 * If skipnan is false the return values are not NA only if there are no NA values in the subset (or window).
 * Otherwise NA values in the subset are ignored.
 * @param start the subset starting position (0 based, defaults to 0)
 * @param end the subset end position (1 based, defaults to the last element position)
 * @param skipnan whether to ignore NA values or not (defaults to true)
 * @param winlength the length of the sliding windows (defaults to 0, no windowing)
 * @param winpos the position along the window to which count will be saved (ignored if wilength==0,defaults to 0)
 * @return a newly instatiated array with the results.
 */
template <class T> gArray<T> gArray<T>::getMin ( gPos start,gPos end,bool skipnan,gSize winlength,gPos winpos) const {
    if ( !i_array ) throw gException ( "operation on object not yet initialized" );
    gPos aend= ( end==0 ) ? ( getSize() ) : ( end );
    if ( ( start>=aend ) || ( aend>getSize() ) ) throw gException ( "invalid range" );
    gArray<T> ret;
    if (winlength==0) {
        ret.acquirecontent ( i_array->getMin (new gArrayInternal<T>( 1 ), i_start+start,i_start+aend,skipnan ) );
    } else {
        if ( winlength>aend-start) throw gException ( "window too large" );
        if ( winpos>=aend-start) throw gException ( "invalid position" );
        ret.acquirecontent ( i_array->getWinMin (new gArrayInternal<T>(aend-start), winlength,winpos,i_start+start,i_start+aend,skipnan));
    }
    return ret;
}

/** @brief Minimum value(s)
 *
 * Finds the minimum value(s) in a subset of this array or in windows sliding through the susbset.
 * If winlength==0 returns a single value gArray with the result. Otherwise returns
 * results calculated in windows of specified length sliding through the specified subset.
 * In this case results are returned in an array of the same length as the subset.
 * The pos parameters specify the offset of the values returned in relation
 * to the window.
 * If skipnan is false the return values are not NA only if there are no NA values in the subset (or window).
 * Otherwise NA values in the subset are ignored.
 * @param[out] res an array to be filled with the results
 * @param start the subset starting position (0 based, defaults to 0)
 * @param end the subset end position (1 based, defaults to the last element position)
 * @param skipnan whether to ignore NA values or not (defaults to true)
 * @param winlength the length of the sliding windows (defaults to 0, no windowing)
 * @param winpos the position along the window to which count will be saved (ignored if wilength==0,defaults to 0)
 * @return a reference to the passed results array
 */
template <class T> gArray<T> & gArray<T>::getMin (gArray<T> & res, gPos start,gPos end,bool skipnan,gSize winlength,gPos winpos ) const {
    if ( !i_array ) throw gException ( "operation on object not yet initialized" );
    gPos aend= ( end==0 ) ? ( getSize() ) : ( end );
    if ( ( start>=aend ) || ( aend>getSize() ) ) throw gException ( "invalid range" );
    if (winlength==0) {
        i_array->getMin (res.i_array, i_start+start,i_start+aend,skipnan );
    } else {
        if ( winlength>aend-start) throw gException ( "window too large" );
        if ( winpos>=aend-start) throw gException ( "invalid position" );
        i_array->getWinMin (res.i_array, winlength,winpos,i_start+start,i_start+aend,skipnan);
    }
    return res;
}

/** @brief Maximum value(s)
 *
 * Finds the maximum value(s) in a subset of this array or in windows sliding through the susbset.
 * If winlength==0 returns a single value gArray with the result. Otherwise returns
 * results calculated in windows of specified length sliding through the specified subset.
 * In this case results are returned in an array of the same length as the subset.
 * The pos parameters specify the offset of the values returned in relation
 * to the window.
 * If skipnan is false the return values are not NA only if there are no NA values in the subset (or window).
 * Otherwise NA values in the subset are ignored.
 * @param start the subset starting position (0 based, defaults to 0)
 * @param end the subset end position (1 based, defaults to the last element position)
 * @param skipnan whether to ignore NA values or not (defaults to true)
 * @param winlength the length of the sliding windows (defaults to 0, no windowing)
 * @param winpos the position along the window to which count will be saved (ignored if wilength==0,defaults to 0)
 * @return a newly instatiated array with the results.
 */
template <class T> gArray<T> gArray<T>::getMax ( gPos start,gPos end,bool skipnan,gSize winlength,gPos winpos) const {
    if ( i_array ) {
        gPos aend= ( end==0 ) ? ( getSize() ) : ( end );
        if ( ( start>=aend ) || ( aend>getSize() ) ) throw gException ( "invalid range" );
        gArray<T> ret;
        if (winlength==0) {
            ret.acquirecontent ( i_array->getMax (new gArrayInternal<T>( 1 ),  i_start+start,i_start+aend,skipnan ) );
        } else {
            if ( winlength>aend-start) throw gException ( "window too large" );
            if ( winpos>=aend-start) throw gException ( "invalid position" );
            ret.acquirecontent ( i_array->getWinMax (new gArrayInternal<T>(aend-start), winlength,winpos,i_start+start,i_start+aend,skipnan ) );
        }
        return ret;
    } else return gArray<T>();
}

/** @brief Maximum value(s)
 *
 * Finds the maximum value(s) in a subset of this array or in windows sliding through the susbset.
 * If winlength==0 returns a single value gArray with the result. Otherwise returns
 * results calculated in windows of specified length sliding through the specified subset.
 * In this case results are returned in an array of the same length as the subset.
 * The pos parameters specify the offset of the values returned in relation
 * to the window.
 * If skipnan is false the return values are not NA only if there are no NA values in the subset (or window).
 * Otherwise NA values in the subset are ignored.
 * @param[out] res an array to be filled with the results
 * @param start the subset starting position (0 based, defaults to 0)
 * @param end the subset end position (1 based, defaults to the last element position)
 * @param skipnan whether to ignore NA values or not (defaults to true)
 * @param winlength the length of the sliding windows (defaults to 0, no windowing)
 * @param winpos the position along the window to which count will be saved (ignored if wilength==0,defaults to 0)
 * @return a reference to the passed results array
 */
template <class T> gArray<T> & gArray<T>::getMax (gArray<T> & res, gPos start,gPos end,bool skipnan,gSize winlength,gPos winpos ) const {
    if ( !i_array ) throw gException ( "operation on object not yet initialized" );
    gPos aend= ( end==0 ) ? ( getSize() ) : ( end );
    if ( ( start>=aend ) || ( aend>getSize() ) ) throw gException ( "invalid range" );
    if (winlength==0) {
        i_array->getMax (res.i_array,  i_start+start,i_start+aend,skipnan );
    } else {
        if ( winlength>aend-start) throw gException ( "window too large" );
        if ( winpos>=aend-start) throw gException ( "invalid position" );
        i_array->getWinMax (res.i_array, winlength,winpos,i_start+start,i_start+aend,skipnan );
    }
    return res;
}

/** @brief Elements sum
 *
 * Sums the elements  in a subset of this array or in windows sliding through the susbset.
 * If winlength==0 returns a single value gArray with the result. Otherwise returns
 * results calculated in windows of specified length sliding through the specified subset.
 * In this case results are returned in an array of the same length as the subset.
 * The pos parameters specify the offset of the values returned in relation
 * to the window.
 * If skipnan is false the return values are not NA only if there are no NA values in the subset (or window).
 * Otherwise NA values in the subset are ignored.
 * @param start the subset starting position (0 based, defaults to 0)
 * @param end the subset end position (1 based, defaults to the last element position)
 * @param skipnan whether to ignore NA values or not (defaults to true)
 * @param winlength the length of the sliding windows (defaults to 0, no windowing)
 * @param winpos the position along the window to which count will be saved (ignored if wilength==0,defaults to 0)
 * @return a newly instatiated array with the results.
 */
template <class T> gArray<T> gArray<T>::getSum ( gPos start,gPos end,bool skipnan,gSize winlength,gPos winpos) const {
    if ( !i_array ) throw gException ( "operation on object not yet initialized" );
    gPos aend= ( end==0 ) ? ( getSize() ) : ( end );
    if ( ( start>=aend ) || ( aend>getSize() ) ) throw gException ( "invalid range" );
    gArray<T> ret;
    if (winlength==0) {
        ret.acquirecontent ( i_array->getSum (new gArrayInternal<T>( 1 ),i_start+start,i_start+aend,skipnan ) );
    } else {
        if ( winlength>aend-start) throw gException ( "window too large" );
        if ( winpos>=aend-start) throw gException ( "invalid position" );
        ret.acquirecontent ( i_array->getWinSum (new gArrayInternal<T>(aend-start), winlength,winpos,i_start+start,i_start+aend,skipnan ) );
    }
    return ret;
}

/** @brief Elements sum
 *
 * Sums the elements  in a subset of this array or in windows sliding through the susbset.
 * If winlength==0 returns a single value gArray with the result. Otherwise returns
 * results calculated in windows of specified length sliding through the specified subset.
 * In this case results are returned in an array of the same length as the subset.
 * The pos parameters specify the offset of the values returned in relation
 * to the window.
 * If skipnan is false the return values are not NA only if there are no NA values in the subset (or window).
 * Otherwise NA values in the subset are ignored.
 * @param[out] res an array to be filled with the results
 * @param start the subset starting position (0 based, defaults to 0)
 * @param end the subset end position (1 based, defaults to the last element position)
 * @param skipnan whether to ignore NA values or not (defaults to true)
 * @param winlength the length of the sliding windows (defaults to 0, no windowing)
 * @param winpos the position along the window to which count will be saved (ignored if wilength==0,defaults to 0)
 * @return a reference to the passed results array
 */
template <class T> gArray<T> & gArray<T>::getSum ( gArray<T> & res, gPos start,gPos end,bool skipnan,gSize winlength,gPos winpos ) const {
    if ( !i_array ) throw gException ( "operation on object not yet initialized" );
    gPos aend= ( end==0 ) ? ( getSize() ) : ( end );
    if ( ( start>=aend ) || ( aend>getSize() ) ) throw gException ( "invalid range" );
    if (winlength==0) {
        i_array->getSum (res.i_array,i_start+start,i_start+aend,skipnan );
    } else {
        if ( winlength>aend-start) throw gException ( "window too large" );
        if ( winpos>=aend-start) throw gException ( "invalid position" );
        i_array->getWinSum (res.i_array, winlength,winpos,i_start+start,i_start+aend,skipnan );
    }
    return res;
}

/** @brief Squared elements sum
 *
 * Sums the square of elements in a subset of this array or in windows sliding through the susbset.
 * If winlength==0 returns a single value gArray with the result. Otherwise returns
 * results calculated in windows of specified length sliding through the specified subset.
 * In this case results are returned in an array of the same length as the subset.
 * The pos parameters specify the offset of the values returned in relation
 * to the window.
 * If skipnan is false the return values are not NA only if there are no NA values in the subset (or window).
 * Otherwise NA values in the subset are ignored.
 * @param start the subset starting position (0 based, defaults to 0)
 * @param end the subset end position (1 based, defaults to the last element position)
 * @param skipnan whether to ignore NA values or not (defaults to true)
 * @param winlength the length of the sliding windows (defaults to 0, no windowing)
 * @param winpos the position along the window to which count will be saved (ignored if wilength==0,defaults to 0)
 * @return a newly instatiated array with the results.
 */
template <class T> gArray<T> gArray<T>::getSquareSum ( gPos start,gPos end,bool skipnan,gSize winlength,gPos winpos) const {
    if ( !i_array ) throw gException ( "operation on object not yet initialized" );
    gPos aend= ( end==0 ) ? ( getSize() ) : ( end );
    if ( ( start>=aend ) || ( aend>getSize() ) ) throw gException ( "invalid range" );
    gArray<T> ret;
    if (winlength==0) {
        ret.acquirecontent ( i_array->getSquareSum (new gArrayInternal<T>( 1 ),  i_start+start,i_start+aend,skipnan ) );
    } else {
        if ( winlength>aend-start) throw gException ( "window too large" );
        if ( winpos>=aend-start) throw gException ( "invalid position" );
        ret.acquirecontent ( i_array->getWinSquareSum (new gArrayInternal<T>(aend-start), winlength,winpos,i_start+start,i_start+aend,skipnan ) );
    }
    return ret;
}

/** @brief Squared elements sum
 *
 * Sums the square of elements in a subset of this array or in windows sliding through the susbset.
 * If winlength==0 returns a single value gArray with the result. Otherwise returns
 * results calculated in windows of specified length sliding through the specified subset.
 * In this case results are returned in an array of the same length as the subset.
 * The pos parameters specify the offset of the values returned in relation
 * to the window.
 * If skipnan is false the return values are not NA only if there are no NA values in the subset (or window).
 * Otherwise NA values in the subset are ignored.
 * @param[out] res an array to be filled with the results
 * @param start the subset starting position (0 based, defaults to 0)
 * @param end the subset end position (1 based, defaults to the last element position)
 * @param skipnan whether to ignore NA values or not (defaults to true)
 * @param winlength the length of the sliding windows (defaults to 0, no windowing)
 * @param winpos the position along the window to which count will be saved (ignored if wilength==0,defaults to 0)
 * @return a reference to the passed results array
 */
template <class T> gArray<T> & gArray<T>::getSquareSum (gArray<T> & res, gPos start,gPos end,bool skipnan,gSize winlength,gPos winpos ) const {
    if ( !i_array ) throw gException ( "operation on object not yet initialized" );
    gPos aend= ( end==0 ) ? ( getSize() ) : ( end );
    if ( ( start>=aend ) || ( aend>getSize() ) ) throw gException ( "invalid range" );
    if (winlength==0) {
        i_array->getSquareSum (res.i_array,  i_start+start,i_start+aend,skipnan );
    } else {
        if ( winlength>aend-start) throw gException ( "window too large" );
        if ( winpos>=aend-start) throw gException ( "invalid position" );
        i_array->getWinSquareSum (res.i_array, winlength,winpos,i_start+start,i_start+aend,skipnan );
    }
    return res;
}

/** @brief Elements mean
 *
 * Calculates the mean value of the elements in a subset of this array or in windows sliding through the susbset.
 * If winlength==0 returns a single value gArray with the result. Otherwise returns
 * results calculated in windows of specified length sliding through the specified subset.
 * In this case results are returned in an array of the same length as the subset.
 * The pos parameters specify the offset of the values returned in relation
 * to the window.
 * If skipnan is false the return values are not NA only if there are no NA values in the subset (or window).
 * Otherwise NA values in the subset are ignored.
 * @param start the subset starting position (0 based, defaults to 0)
 * @param end the subset end position (1 based, defaults to the last element position)
 * @param skipnan whether to ignore NA values or not (defaults to true)
 * @param winlength the length of the sliding windows (defaults to 0, no windowing)
 * @param winpos the position along the window to which count will be saved (ignored if wilength==0,defaults to 0)
 * @return a newly instatiated array with the results.
 */
template <class T> gArray<gScore> gArray<T>::getMean ( gPos start,gPos end,bool skipnan,gSize winlength,gPos winpos) const {
    if ( !i_array ) throw gException ( "operation on object not yet initialized" );
    gPos aend= ( end==0 ) ? ( getSize() ) : ( end );
    if ( ( start>=aend ) || ( aend>getSize() ) ) throw gException ( "invalid range" );
    gArray<gScore> ret;
    if (winlength==0) {
        ret.acquirecontent ( i_array->getMean (new gArrayInternal<gScore>( 1 ), i_start+start,i_start+aend,skipnan ) );
    } else {
        if ( winlength>aend-start) throw gException ( "window too large" );
        if ( winpos>=aend-start) throw gException ( "invalid position" );
        ret.acquirecontent ( i_array->getWinMean (new gArrayInternal<gScore>(aend-start,0,true), winlength,winpos,i_start+start,i_start+aend,skipnan ) );
    }
    return ret;
}

/** @brief Elements mean
 *
 * Calculates the mean value of the elements in a subset of this array or in windows sliding through the susbset.
 * If winlength==0 returns a single value gArray with the result. Otherwise returns
 * results calculated in windows of specified length sliding through the specified subset.
 * In this case results are returned in an array of the same length as the subset.
 * The pos parameters specify the offset of the values returned in relation
 * to the window.
 * If skipnan is false the return values are not NA only if there are no NA values in the subset (or window).
 * Otherwise NA values in the subset are ignored.
 * @param[out] res an array to be filled with the results
 * @param start the subset starting position (0 based, defaults to 0)
 * @param end the subset end position (1 based, defaults to the last element position)
 * @param skipnan whether to ignore NA values or not (defaults to true)
 * @param winlength the length of the sliding windows (defaults to 0, no windowing)
 * @param winpos the position along the window to which count will be saved (ignored if wilength==0,defaults to 0)
 * @return a reference to the passed results array
 */
template <class T> gArray<gScore> & gArray<T>::getMean (gArray<gScore>  & res, gPos start,gPos end,bool skipnan,gSize winlength,gPos winpos ) const {
    if ( !i_array ) throw gException ( "operation on object not yet initialized" );
    gPos aend= ( end==0 ) ? ( getSize() ) : ( end );
    if ( ( start>=aend ) || ( aend>getSize() ) ) throw gException ( "invalid range" );
    if (winlength==0) {
        i_array->getMean (res.i_array, i_start+start,i_start+aend,skipnan );
    } else {
        if ( winlength>aend-start) throw gException ( "window too large" );
        if ( winpos>=aend-start) throw gException ( "invalid position" );
        i_array->getWinMean ( res.i_array, winlength,winpos,i_start+start,i_start+aend,skipnan );
    }
    return res;
}

/** @brief Elements standard deviation
 *
 * Calculates the standard deviation of the elements from a subset of this array or from windows sliding through the susbset.
 * If winlength==0 returns a single value gArray with the result. Otherwise returns
 * results calculated in windows of specified length sliding through the specified subset.
 * In this case results are returned in an array of the same length as the subset.
 * The pos parameters specify the offset of the values returned in relation
 * to the window.
 * If skipnan is false the return values are not NA only if there are no NA values in the subset (or window).
 * Otherwise NA values in the subset are ignored.
 * @param unbiased whether a corrected (true) or uncorrected sd estimatr will be used (defaults true)
 * @param start the subset starting position (0 based, defaults to 0)
 * @param end the subset end position (1 based, defaults to the last element position)
 * @param skipnan whether to ignore NA values or not (defaults to true)
 * @param winlength the length of the sliding windows (defaults to 0, no windowing)
 * @param winpos the position along the window to which count will be saved (ignored if wilength==0,defaults to 0)
 * @return a newly instatiated array with the results.
 */
template <class T> gArray<gScore> gArray<T>::getStdDev ( bool unbiased, gPos start, gPos end, bool skipnan,gSize winlength,gPos winpos) const {
    if ( !i_array ) throw gException ( "operation on object not yet initialized" );
    gPos aend= ( end==0 ) ? ( getSize() ) : ( end );
    if ( ( start>=aend ) || ( aend>getSize() ) ) throw gException ( "invalid range" );
    gArray<gScore> ret;
    if (winlength==0) {
        ret.acquirecontent ( i_array->getStdDev (new gArrayInternal<gScore>( 1 ), i_start+start,i_start+aend,skipnan,unbiased ) );
    } else {
        if ( winlength>aend-start) throw gException ( "window too large" );
        if ( winpos>=aend-start) throw gException ( "invalid position" );
        ret.acquirecontent ( i_array->getWinStdDev (new gArrayInternal<gScore>(aend-start), winlength,winpos,i_start+start,i_start+aend,skipnan, unbiased ) );
    }
    return ret;
}

/** @brief Elements standard deviation
 *
 * Calculates the standard deviation of the elements from a subset of this array or from windows sliding through the susbset.
 * If winlength==0 returns a single value gArray with the result. Otherwise returns
 * results calculated in windows of specified length sliding through the specified subset.
 * In this case results are returned in an array of the same length as the subset.
 * The pos parameters specify the offset of the values returned in relation
 * to the window.
 * If skipnan is false the return values are not NA only if there are no NA values in the subset (or window).
 * Otherwise NA values in the subset are ignored.
 * @param[out] res an array to be filled with the results
 * @param unbiased whether a corrected (true) or uncorrected sd estimatr will be used (defaults true)
 * @param start the subset starting position (0 based, defaults to 0)
 * @param end the subset end position (1 based, defaults to the last element position)
 * @param skipnan whether to ignore NA values or not (defaults to true)
 * @param winlength the length of the sliding windows (defaults to 0, no windowing)
 * @param winpos the position along the window to which count will be saved (ignored if wilength==0,defaults to 0)
 * @return a reference to the passed results array
 */
template <class T> gArray<gScore> & gArray<T>::getStdDev (gArray<gScore>  & res,  bool unbiased, gPos start,gPos end,bool skipnan,gSize winlength,gPos winpos) const {
    if ( !i_array ) throw gException ( "operation on object not yet initialized" );
    gPos aend= ( end==0 ) ? ( getSize() ) : ( end );
    if ( ( start>=aend ) || ( aend>getSize() ) ) throw gException ( "invalid range" );
    if (winlength==0) {
        i_array->getStdDev (res.i_array, i_start+start,i_start+aend,skipnan,unbiased );
    } else {
        if ( winlength>aend-start) throw gException ( "window too large" );
        if ( winpos>=aend-start) throw gException ( "invalid position" );
        i_array->getWinStdDev (res.i_array, winlength,winpos,i_start+start,i_start+aend,skipnan,unbiased );
    }
    return res;
}

//private & protected
template <class T> gArray<T>::gArray ( T* data,gSize length,const gArray<gPos> & isNaFlags,bool sorted ) {
    i_array=NULL;
    i_start=0;
    i_end=0;
    ownData ( length,data,isNaFlags,sorted );
}

template <class T> void gArray<T>::acquirecontent ( gArrayInternal<T> * content ) {
    releasecontent();
    if ( content ) {
        i_array=content;
        i_array->i_refcount++;
        i_start=0;
        i_end=i_array->i_length;
    }
}

template <class T> void gArray<T>::releasecontent() {
    if ( i_array&& ( i_array->i_refcount>0 ) ) {
        i_array->i_refcount--;
        if ( i_array->i_refcount==0 ) {
            delete i_array;
        }
    }
    i_array=NULL;
    i_start=0;
    i_end=0;
}

template <class T> T * gArray<T>::getData() const {
    T* ret=NULL;
    if ( i_array ) ret=i_array->i_mem+i_start;
    return ret;
}

template <class T> template <class C> gArray<T> gArray<T>::operator () (const gArray<C> & rows,const gArray<C> & cols,gSize ncols) const {
    gArray<T> ret;
    ret.acquirecontent(i_array->getValues(i_start,i_end,*rows.i_array,*cols.i_array,ncols,rows.i_start,rows.i_end,cols.i_start,cols.i_end));
    return ret;
}

template <class T> template <class C> gArray<T> & gArray<T>::operator () (gArray<T> & res,const gArray<C> & rows,const gArray<C> & cols,gSize ncols) const {
    i_array->getValues(*res.i_array,i_start,i_end,*rows.i_array,*cols.i_array,ncols,rows.i_start,rows.i_end,cols.i_start,cols.i_end);
    return res;
}
//--------------------------------------------------------------------------



//--------------------------------------------------------------------------
// template class gArrayRetrieverImplementation implementation
//--------------------------------------------------------------------------
template <class T>  gArrayRetrieverImplementation<T>::gArrayRetrieverImplementation(gSize nFeatures):gRetrieverImplementation() {
  i_nFeatures=nFeatures;
}

template <class T>  gArrayRetrieverImplementation<T>::~gArrayRetrieverImplementation(){
}



template <class T> template<class C> const C * gArrayRetrieverImplementation<T>::getArrayDataPtr ( const gArray<C> & array ) const {
    return array.getData();
}

template <class T> template <class C> C * gArrayRetrieverImplementation<T>::getArrayDataPtr(gArray<C> & array) const{
  return array.getData();
}


template <class T> unsigned gArrayRetrieverImplementation<T>::getArrayReferenceCount(gArray<T> & array) const {
    unsigned ret=0;
    if (array.i_array) ret=array.i_array->i_refCount;
    return ret;
}
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
// template class gArrayRetriever implementation
//--------------------------------------------------------------------------
template <class T> gArrayRetriever<T>::gArrayRetriever():gRetriever() {
}

template <class T> gArrayRetriever<T>::gArrayRetriever(const gArrayRetrieverImplementation<T> & implementation):gRetriever(implementation) {
}

template <class T> gArrayRetriever<T>::gArrayRetriever(const gArrayRetriever<T> & retriever):gRetriever(retriever) {
}

template <class T> gArrayRetriever<T>::~gArrayRetriever() {
}

template <class T> gSize gArrayRetriever<T>::getFeaturesCount() const{
  return ((const gArrayRetrieverImplementation<T> &) getImplementation()).i_nFeatures;
}
//--------------------------------------------------------------------------


//--------------------------------------------------------------------------
// gMatrix implementation
//--------------------------------------------------------------------------
/** @brief Empty constructor
 *
 * Instatiates an empty gMatrix object
 */
template <class T> gMatrix<T>::gMatrix() {
    nrows=0;
    ncols=0;
}

/** @brief Initializing contructor (raw data pointer)
 *
 * Instatiates a gMatrix object with the spcified number of columns and rows.
 * Values are intialized by rows from the pointer passed
 * NA are all set to false
 * @param rows desired number of rows
 * @param cols desired number of columns
 * @param data pointer to data (must be at least nrow*ncolumns)
 */
template <class T> gMatrix<T>::gMatrix ( gSize rows,gSize cols,const T * data ) :gArray<T> ( data,rows*cols ) {
    nrows=rows;
    ncols=cols;
}

/** @brief Initializing contructor (gArray)
 *
 * Instatiates a gMatrix object with the spcified number of columns and rows.
 * Values and NAs are intialized by rows from the gArray passed
 * @param rows desired number of rows
 * @param cols desired number of columns
 * @param data gArray containing initialization data (must be at least rows*columns length)
 */
template <class T> gMatrix<T>::gMatrix ( gSize rows,gSize cols,const gArray<T> & data ) :gArray<T> ( data ) {
    nrows=rows;
    ncols=cols;
    if ( ( rows*cols ) > this->getSize() ) throw gException ( "invalid data array" );
}

/** @brief Initializing contructor (single value)
 *
 * Instatiates a gMatrix object with the specified number of columns and rows.
 * Values and NAs are intialized to the same value specified by parameters
 * @param rows desired number of rows
 * @param cols desired number of columns
 * @param fillvalue initializing value
 * @param fillNA initializing NA status
 */
template <class T> gMatrix<T>::gMatrix ( gSize rows,gSize cols,T fillvalue,bool fillNA ) :gArray<T> ( fillvalue,rows*cols,fillNA ) {
    nrows=rows;
    ncols=cols;
}

/** @brief Non-intializing contructor
 *
 * Instatiates a gMatrix object with the specified number of columns and rows.
 * Values and NAs are not initialized
 * @param rows
 * @param cols
 */
template <class T> gMatrix<T>::gMatrix ( gSize rows,gSize cols ) :gArray<T> ( 0,rows*cols,false,false ) {
    nrows=rows;
    ncols=cols;
}

/** @brief Copy contructor
 *
 * Instatiate a copy of the object passed as const
 * @param matrix the object to copy from
 */
template <class T> gMatrix<T>::gMatrix ( const gMatrix<T> & matrix ) {
    *this=matrix;
}

/** @brief Set a single value
 *
 * Sets a value ata a specified position.
 * @param row the row number (0 based)
 * @param col the column number (0 based)
 * @param value the value to be set
 * @param isna if the values should be marked as NA (default: false)
 */
template <class T> void gMatrix<T>::setValue ( gPos row, gPos col, T value, bool isna) {
  if((nrows==0)||(ncols==0)){
    if(nrows<=row) nrows=row+1;
    if(ncols<=col) ncols=col+1;
	gArray<T>::setValue ( ( nrows*ncols -1  ),0,true,0,true);
    gArray<T>::setValue ( ( col+row*ncols ),value,isna,0,true);
  }else{
	if ( ( col<ncols ) ) {
	  if(nrows<=row){
		nrows=row+1;
		gArray<T>::setValue ( ( nrows*ncols -1  ),0,true,0,true);
	  }
	  gArray<T>::setValue ( ( col+row*ncols ),value,isna,0,true);
	}else{
	  if(nrows<=row) setRow(row,gArray<T>(0,ncols,true));
	  gArray<T> newcol(0,nrows,true);
	  newcol.setValue(row,value,isna);
	  setCol(col,newcol);
    }
  }
}

/** @brief Sets multiple values
 *
 * Sets matrix values at specified positions
 * @param rows an array containing row numbers
 * @param cols an array containing columns numbers
 * @param values the values to be set
 */
template <class T> void gMatrix<T>::setValues ( const gArray<gPos> & rows,const gArray<gPos> & cols, const gArray<T> & values ) {
  gPos rmax=rows.getMax()[0];
  gPos cmax=cols.getMax()[0];
  if ( ( rmax<nrows ) && ( cmax<ncols ) && ( cols.getSize() ==rows.getSize() ) ) {
    gArray<T>::setValues ( cols+ ( rows*ncols ) ,values );
  } else throw ( gException ( "gMatrix: invalid range" ) );
}

/** @brief Sets size and values from row data
 *
 * Re-initializes the matrix specifying dimaensions and data. The raw pointer passed must point to an allocated memory
 * of size rows*cols (or greater). Values are taken by rows.
 * @param rows the number of rows
 * @param cols the number of columns
 * @param data the data to be copied (must be a pointer to an allocated bklock of memory
 * @param napos the positions of the elements to be marked as NA
 */
template <class T> void gMatrix<T>::setData ( gSize rows, gSize cols, const T* data, const gArray<gPos> & napos) {
    nrows=rows;
    ncols=cols;
    gArray<T>::setData ( nrows*ncols,data,napos);
}

/** @brief Sets a row
 *
 * Sets the values of the specified row to those contained in the passed array
 * @param row the row number
 * @param values the values to be set
 */
template <class T> void gMatrix<T>::setRow ( gPos row,const gArray<T> & values ) {
  if(ncols==0) ncols=values.getSize();
  if(nrows<=row) nrows=row+1;
  
    if (  values.getSize() == ncols  ) {
        gArray<T>::setValues ( getArray<gPos> ( 0,ncols-1,1 ) + ( ncols*row ),values );
    } else throw ( gException ( "gMatrix: invalid range" ) );

}

/** @brief Sets a column
 *
 * Sets the values of the specified column to those contained in the passed array
 * @param col the column number
 * @param values the values to be set
 */
template <class T> void gMatrix<T>::setCol ( gPos col,const gArray<T> & values ) {
  if((ncols==0)||(nrows==0)){
   if(ncols==0) ncols=col+1;
   if(nrows==0) nrows=values.getSize();
   gArray<T>::setValues ( getArray<gPos> ( 0,nrows-1,1 ) *ncols+col,values );
  }else{
    if(col<ncols){
      gArray<T>::setValues ( getArray<gPos> ( 0,nrows-1,1 ) *ncols+col,values );
    }else{
      gMatrix<T> newdata(nrows,col+1);
      for(gPos i=0;i<nrows;i++){
        gArray<T> arow=getRow(i);
        arow.setValue(col,values[i],values.isNA(i));
        newdata.setRow(i,arow);
      }
      *this=newdata;
    }
  }
}

/** @brief Assignment operator
 * 
 * Make this object an exacto copy of the passed one
 * @param matrix the object to be copied
 * @return a reference to this
 */
template <class T> gMatrix<T> & gMatrix<T>::operator = ( const gMatrix<T> & matrix ) {
    ( ( gArray<T> & ) *this ) = ( ( gArray<T> & ) matrix );
    nrows=matrix.nrows;
    ncols=matrix.ncols;
    return *this;
}

/** @brief Number of rows
 *
 * Return the number of rows in this matrix
 * @return the number of rows
 */
template <class T> gSize gMatrix<T>::getRowsNum() const {
    return nrows;
}

/** @brief Number of Columns
 *
 * Return the number of rows in this matrix
 * @return the number of columns
 */
template <class T> gSize gMatrix<T>::getColsNum() const {
    return ncols;
}

/** @brief Single value access
 *
 * Returns the value at the specified position
 * @param row the row number
 * @param col the column number
 * @return the value
 */
template <class T> T gMatrix<T>::operator () ( gPos row,gPos col ) const {
    if ( ( row<nrows ) && ( col<ncols ) ) {
        return ( *this ) [col+row*ncols];
    } else throw ( gException ( "gMatrix: invalid range" ) );
}

/** @brief Multiple values access template operator
 *
 * Returns the values at the spcified positions. Positions can be specified using any type of gArray<T>
 * but will be cast to gPos before indexing.
 * @tparam the type of rows and columns numbers
 * @param rows an garray<C> containing the row numbers og the elements to retrieve
 * @param cols an garray<C> containing the column numbers og the elements to retrieve
 * @return a newly instatiated array with the required values.
 */
template <class T> template <class C> gArray<T> gMatrix<T>::operator () (const gArray<C> & rows, const gArray<C> & cols) const {
    return (*(gArray<T>*)this)(rows,cols,ncols);
}

/** @brief Multiple values access template operator
 *
 * Returns the values at the spcified positions. Positions can be specified using any type of gArray<T>
 * but will be cast to gPos before indexing.
 * @tparam the type of rows and columns numbers
 * @param[out] res An array to be filled with teh results
 * @param rows an garray<C> containing the row numbers og the elements to retrieve
 * @param cols an garray<C> containing the column numbers og the elements to retrieve
 * @return a reference to the passed results array
 */
template <class T> template <class C> gArray<T> & gMatrix<T>::operator () (gArray<T> & res,const gArray<C> & rows, const gArray<C> & cols) const {
    return (*(gArray<T>*)this)(res,rows,cols,ncols);
}

/** @brief Single value NA status access
 *
 * Returns the NA status of the value at the specified position
 * @param row the row number
 * @param col the column number
 * @return true if position is NA, false otherwise
 */
template <class T> gBool gMatrix<T>::isNA ( gPos row,gPos col ) const {
    if ( ( row<nrows ) && ( col<ncols ) ) {
        return gArray<T>::isNA ( ( gPos ) ( col+row*ncols ) );
    } else throw ( gException ( "gMatrix: invalid range" ) );
}

/** @brief Row access
 *
 * Returns an arry contining the values of the specified row
 * @param rownum the row number
 * @return an newly allocated array contining the row values
 */
template <class T> gArray<T> gMatrix<T>::getRow ( gPos rownum ) const {
    gArray<T> ret ( 0,ncols,false );
    ret.acquirecontent ( ( *this ).i_array );
    if ( ( rownum>=0 ) && ( rownum<nrows ) ) {
        ret.i_start= ( rownum ) *ncols;
        ret.i_end= ( rownum+1 ) *ncols;
    } else throw ( gException ( "gMatrix: invalid range" ) );
    return ret;
}

/** @brief Column access
 *
 * Returns an arry contining the values of the specified clumn
 * @param colnum the column number
 * @return an newly allocated array contining the column values
 */
template <class T> gArray<T> gMatrix<T>::getCol ( gPos colnum ) const {
    return ( *this ) [getArray<gPos> ( colnum,this->getSize()-1,ncols ) ];
}
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
// template gPositionInterval implementation
//--------------------------------------------------------------------------
/** @brief Constructor from positions
 *
 * Instatiate a new gInterval from absolute boundaries
 * @param start the absolute interval-start position (0 based,default 0)
 * @param end the absolute interval-end position (1 based,default 0)
 * @param validStart wether the starting positionis valid or not (defaults true);
 * @param validEnd wether the ending positionis valid or not (defaults true)
 */
template <class T> gPositionInterval<T>::gPositionInterval(T start,T end,bool validStart,bool validEnd):gArray<T>(0,2,false,false) {
    setValue(0,start,!validStart);
    setValue(1,end,!validEnd);
}

/** @brief Constructor from position array
*
* Instatiate a new gReferenceInterval from absolute boundaries provided as an array of 0 based positions.
* @param positions a gArray<gPos>  containing at least two absolute 0 based positions (only the first two elements will be used)
*/
template <class T> gPositionInterval<T>::gPositionInterval(const gArray<T> & positions):gArray<T>(0,2,false,false)  {
    setValue(0,positions[0],positions.isNA(0));
    setValue(1,positions[1]+1,positions.isNA(1));
}

/** @brief Copy constructor
*
* Instatiate a new gReferenceInterval as a copy of the passed one
* @param positionInterval the interval to copy
* @param startOffset an optional start offset (defaults=0);
* @param endOffset an optional end offset (defaults=0);
*/
template <class T> gPositionInterval<T>::gPositionInterval(const gPositionInterval<T> & positionInterval, T startOffset, T endOffset):gArray<T>(positionInterval) {
    T start=getStart()+startOffset;
    T end=getEnd()+endOffset;
    if (start<=end) {
        setValue(0,start,!validStart());
        setValue(1,end,!validEnd());
    } else throw gException("k_gElementInterval: invalid offsets in copy constructor");
}

/** @brief Assignment operator
*
* Makes this object identical to the passed one
* @param positionInterval the interval to copy
* @return a reference to this object
*/
template <class T> gPositionInterval<T> & gPositionInterval<T>::operator = (const gPositionInterval<T> & positionInterval) {
    setValue(0,positionInterval[0],positionInterval.isNA(0));
    setValue(1,positionInterval[1],positionInterval.isNA(1));
    return *this;
}

/** @brief Length
*
* return this interval's length
* @return a gSize value
*/
template <class T> gSize gPositionInterval<T>::getLength() const {
    return (*this)[1]-(*this)[0];
}

/** @brief Start position
*
* Returns the (absolute, 0 based) start position of this interval
* @return a gPos value
*/
template <class T> T gPositionInterval<T>::getStart() const {
    return (*this)[0];
}

/** @brief End position
*
* Returns the (absolute, 1 based) end position of this interval
* @return a gPos value
*/
template <class T> T gPositionInterval<T>::getEnd() const {
    return (*this)[1];
}

/** @brief Start validity
*
* Returns true if the interval start position is valid
* @return a boolean value
*/
template <class T> bool gPositionInterval<T>::validStart() const {
    return !gArray<T>::isNA(0);
}

/** @brief End validity
*
* Returns true if the interval end position is valid
* @return a boolean value
*/
template <class T> bool gPositionInterval<T>::validEnd() const {
    return !gArray<T>::isNA(1);
}

/** @brief Positions
*
* Returns the boundaries as 0 based absolute positions
* @return a gArray<gPos> object
*/
template <class T> gArray<T> gPositionInterval<T>::getPositions() const {
    gArray<T> ret(*this);
    ret.setValue(1,ret[1]-1,ret.isNA(1));
    return ret;
}

/** @brief Contained positions
*
* Check if "positions" are contained in this interval
* @param positions an array of referencePositions
* @return a gArray<bool>
*/
template <class T> gArray<bool> gPositionInterval<T>::contains(const gArray<T> & positions) const{
    return (positions>=getStart())&&(positions<getEnd());
}

/** @brief Equaliy operator
*
* Check if an interval equal this one
* @param interval an interval
* @return a gBool value
*/
template <class T> gBool gPositionInterval<T>::operator == (const gPositionInterval<T> & interval) const{
  return ((getStart()==interval.getStart())&&(getEnd()==interval.getEnd())&&(validStart()==interval.validEnd())&&(validEnd()==interval.validEnd()));
}

/** @brief Greater than  operator
*
* Check if an interval has the start greater than the other
* @param interval an interval
* @return a gBool value
*/
template <class T> gBool gPositionInterval<T>::operator > (const gPositionInterval<T> & interval) const{
  return validStart() && interval.validStart() && (getStart() > interval.getStart());
}

/** @brief Lesser than  operator
*
* Check if an interval has the start lesser than the other
* @param interval an interval
* @return a gBool value
*/
template <class T> gBool gPositionInterval<T>::operator < (const gPositionInterval<T> & interval) const{
 return validStart() && interval.validStart() && (getStart() < interval.getStart());
}


/** @brief Intersection
*
* Intersection interval
* @param interval an interval to intersecate
* @return a k_gReferenceInterval object
*/
template <class T> gPositionInterval<T> gPositionInterval<T>::getIntersection(const gPositionInterval<T> & interval) const {
    T start,end;
    bool vstart,vend;

    if (getStart()>=interval.getStart()) {
        start=getStart();
        vstart=validStart();
    } else {
        start=interval.getStart();
        vstart=interval.validStart();
    }

    if (getEnd()<=interval.getEnd()) {
        end=getEnd();
        vend=validEnd();
    } else {
        end=interval.getEnd();
        vend=interval.validEnd();
    }
    if (start<end) {
        return gPositionInterval<T>(start,end,vstart,vend);
    } else {
        return gPositionInterval<T>(start,start,vstart,vstart);
    }
}

/** @brief Left difference
*
* Returns the left difference interval
* @param interval an interval to diff
* @return a k_gReferenceInterval object
*/
template <class T> gPositionInterval<T> gPositionInterval<T>::getLeftDiff(const gPositionInterval<T> & interval) const {
    if (interval.getStart()<getStart()) {
        if (interval.getEnd()>getStart()) {
            return gPositionInterval<T>(interval.getStart(),getStart(),interval.validStart(),validStart());
        } else {
            return interval;
        }
    } else {
        return gPositionInterval<T>(getStart(),getStart(),validStart(),validStart());
    }
}

/** @brief Right difference
*
* Returns the right difference interval
* @param interval an interval to diff
* @return a k_gReferenceInterval object
*/
template <class T> gPositionInterval<T> gPositionInterval<T>::getRightDiff(const gPositionInterval<T> & interval) const {
    if (interval.getEnd()>getEnd()) {
        if (getEnd()>interval.getStart()) {
            return gPositionInterval<T>(getEnd(),interval.getEnd(),validEnd(),interval.validStart());
        } else {
            return interval;
        }
    } else {
        return gPositionInterval<T>(getEnd(),getEnd(),validEnd(),validEnd());
    }
}

/** @brief Union
*
* Returns the interval resulting from merging this one to the provided one
* @param interval an interval to merge
* @return a k_gReferenceInterval object
*/
template <class T> gPositionInterval<T> gPositionInterval<T>::getUnion(const gPositionInterval<T> & interval) const {
    T start,end;
    bool vstart,vend;
    if((getIntersection(interval).getLength()>0)||(getStart()==interval.getEnd())||(getEnd()==interval.getStart())){
      if (getStart()<=interval.getStart()) {
          start=getStart();
          vstart=validStart();
      } else {
          start=interval.getStart();
          vstart=interval.validStart();
      }

      if (getEnd()>=interval.getEnd()) {
          end=getEnd();
          vend=validEnd();
      } else {
          end=interval.getEnd();
          vend=interval.validEnd();
      }
      if (start<end) {
          return gPositionInterval<T>(start,end,vstart,vend);
      } else {
          return gPositionInterval<T>(0,0,vstart,vstart);
      }
    }else{
      return gPositionInterval<T>(0,0,vstart,vstart);
    }
}
//--------------------------------------------------------------------------

}//end of geco namespace



#endif

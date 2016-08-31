/**
 * @file geco_base.h
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
 * This file contains the declarations of some base GeCo++ classes.
 *
 */
#ifndef __GECO_BASE__
#define __GECO_BASE__

#include "geco_define.h"
#include <math.h>
#include <exception>
#include <string>
#include <memory.h>

namespace geco {

template <class T> class gArray;
template <class T> class gArrayRetriever;
template <class T> class gMatrix;

/** @brief Exception class.
 *
 * This class is used throughout the library to
 * manage exceptions. It specializes std::exception
 * by adding textual description of what happened.
 * It actually is quite unusefull but has been defined
 * as a placeholder for future improvement.
 */
class gException:public std::exception {
private:
    std::string       i_reason;
public:
    /** @name Constructors and destructor */
    //@{
    gException(const char *reason);
    gException(const gException &e);
    virtual ~gException() throw();
    //@}

    /** @name Copy operators amd function members */
    //@{
    gException & operator = (const gException & e);
    //@}

    /** @name Exception Information function members */
    //@{
    const char* what () const throw();
    //@}
};

/**  @brief Bits Array Class
 *
 * Class to efficiently manage bit arrays.
 */
class gBitsArray {

private:
    gBytes *i_val;
    gSize  i_alength;
    gSize  i_mlength;

protected:

public:
    /** @name Constructors and destructor */
//@{
    gBitsArray();
    gBitsArray(gSize length);
    gBitsArray(gSize length,bool value);
    gBitsArray(const gBitsArray & array,gPos start=0, gPos end=0);
    ~gBitsArray();
//@}

    /** @name Bits array Information operators and access function members */
//@{
    gSize getSize() const;
    gSize getSetBitsCount(gPos start,gPos end) const;///< Set bits counter
    /** @brief Bit value retriever
     *
     * This method returns the number of set bits in the specified range.
     * @param position Position of the bit to retrieve
     * @return the status of the specified bit.
     */
    bool getValue(gPos position) const {
        gBytes b=i_val[position >> mP2];
        if (b==0) return false;
        else return b & MG[position & mP1];
    }
//@}

    /** @name Bits changing function members */
//@{
    void setSize(gSize size,bool initvalues=true,bool fillvalue=false);
    /** @brief Set a bit at a given value
     *
     * This method sets bit at a given position to the given status
     * @param position Position of the bit to set
     * @param value Value to set the bit to.
     */
    void setValue(gPos position,bool value) {
        if (value) i_val[position >> mP2] |=  MG[position & mP1];
        else i_val[position >> mP2] &=  MR[position & mP1];
    }
    /** @brief Copy a bit value
     *
     * This method sets the bit in destpos to the value of the srcpos bit in src array.
     * @param destpos The position of the bit to set in this array (zero based).
     * @param src Source array
     * @param srcpos The position of the bit to copy from (zero based).
     */
    void copyValue(gPos destpos,const gBitsArray & src,gPos srcpos) {
        gBytes dpos=destpos >> mP2;
        gBytes spos=srcpos >> mP2;
        gBytes o1=destpos & mP1;
        gBytes o2=srcpos & mP1;
        gBytes a = i_val[dpos] & MR[o1];
        gBytes b = ((src.i_val[spos] & MG[o2]) << o2) >> o1;
        i_val[dpos]= a & b;
    }
    void copyValues(gPos start,const gBitsArray & array,gPos astart,gPos aend);
    /** @brief Bit setter
     *
     * Set one bit at a given position
     * @param position Zero based bit position
     */
    void setBit(gPos position) {
        i_val[position >> mP2] |= MG[position & mP1];
    }
    /** @brief Bit resetter
     *
     * Reset one bit at a given position
     * @param position Zero based bit position
     */
    void resetBit(gPos position) {
        i_val[position >> mP2] &= MR[position & mP1];
    }
    /** @brief Bit swapper
     *
     * Swap bit value at a given position
     * @param position Zero based bit position
     */
    void swapBit(gPos position) {
        i_val[position >> mP2] ^= MR[position & mP1];
    }
    void setAll();
    void resetAll();
//@}

    /** @name Logical bit atomic operations function members */
//@{
    void Or(gPos start,const gBitsArray & array,gPos astart,gPos aend);
    void And(gPos start,const gBitsArray & array,gPos astart,gPos aend);
//@}
};

/**  @brief Base class for retriever implementations
 *
 * This pure virtual class is a base class for all the retrievers implementation. 
 * It implement a ref counting mechanism (toghether with gRetriever) that allow 
 * for copying retrievers without duplicating memory.It has a virtual method "clone" 
 * and a virtual destructor that must be overriden by derived objects.
 */
class gRetrieverImplementation {
  friend class gRetriever;
private:
    gSize i_refCount;
    /** @brief Cloning method
     *
     * This pure virtual member function should be defined by derived object in order to return a pointer 
     * to an instance copy of the derived object itself. 
     * @return A gRetrieverImplementation pointer to a copy of the object.
     */
    virtual gRetrieverImplementation * clone() const = 0;
public:
  
  gRetrieverImplementation();
  virtual ~gRetrieverImplementation();
};

/**  @brief Base class for retrievers
 *
 * This is a base class for all the retrievers. It implement (toghether with gRetrieverImplementation) 
 * a ref counting mechanism that allow for copying retrievers without memory dupication.
 */
class gRetriever {
private:
  gRetrieverImplementation * i_implementation;
protected:
  const gRetrieverImplementation & getImplementation() const;
public:
  gRetriever();
  gRetriever( const gRetrieverImplementation & implementation);
  gRetriever( const gRetriever & retriever);
  ~gRetriever();
  gRetriever & operator = (const gRetriever & retriever);
};


//This class is used interbnally by gArray. you will never instatiate an object of this type
template <class T> class gArrayInternal {
    template <class U> friend class gArrayRetriever;
    template <class U> friend class gArrayInternal;
    template <class U> friend class gArray;

  private:
    T *      i_mem;
    gBitsArray i_nans;
    unsigned i_refcount;
    bool     i_sorted;
    gSize    i_length;

    void    downHeap(T *a,gPos *b, gPos k, gSize N) const;
    void    internalSort(T *i_sort,gPos *i_map,gSize N) const;
    gSize   internalCount(T value,gPos start,gPos end,gSize nacount) const;
    void    internalSearch(gArrayInternal<gPos> *positions,T value,const gArrayInternal<T> & array,gPos start, gPos end,bool excludeNA) const;
    


    gArrayInternal();
    gArrayInternal(gSize length);
    gArrayInternal(gSize length,T value,bool isna);
    gArrayInternal(gSize length,T* data,const gBitsArray & isNaFlags,bool getdataownership,bool sorted);
    gArrayInternal(const gArrayInternal<T> & array,gPos start,gPos end,T fillvalue,bool fillNA);
    ~gArrayInternal();

    void setLength(gSize length,T fillvalue,bool fillNA);
    void setData(gSize length, T* data,const gBitsArray & isNAFlags,T fillvalue,bool fillNA,bool getdataownership,bool sorted);
    void setValue(gPos pos, T value,bool isna,T fillvalue,bool fillNA);
    void setValues(gPos start,gPos pstart,gPos pend,const gArrayInternal<gPos> & positions,gPos vstart,gPos vend,const gArrayInternal<T> & values,T fillvalue,bool fillNA);
    void datamove(gArrayInternal<T> & dest,gPos dpos,const gArrayInternal<T> & src,gPos spos,gSize npos);

    gArrayInternal<T> & selfRevert(gPos start, gPos end);
    gArrayInternal<T> & selfConcatenate(const gArrayInternal<T> & array,gPos astart,gPos aend);
    gArrayInternal<T> & selfReplace(gPos pos,const gArrayInternal<T> & array,gPos astart,gPos aend);
    gArrayInternal<T> & selfInsert(gPos pos,const gArrayInternal<T> & array,gPos astart,gPos aend);
    gArrayInternal<T> & selfRemove(gPos pos,gSize npos);
    gArrayInternal<T> & selfAdd(const gArrayInternal<T> & array,gPos start,gPos end,gPos astart,gPos aend);
    gArrayInternal<T> & selfSubtract(const gArrayInternal<T> & array,gPos start,gPos end,gPos astart,gPos aend);
    gArrayInternal<T> & selfMultiply(const gArrayInternal<T> & array,gPos start,gPos end,gPos astart,gPos aend);
    gArrayInternal<T> & selfDivide(const gArrayInternal<T> & array,gPos start,gPos end,gPos astart,gPos aend);    

    const T at(gPos pos) const;
    const T operator [] (gPos pos) const;
    template <class C> gArrayInternal<T> * getValues(gPos start,gPos end,const gArrayInternal<C> & pos,gPos pstart,gPos pend) const;
    template <class C> gArrayInternal<T> & getValues(gArrayInternal<T> & res,gPos start,gPos end,const gArrayInternal<C> & rows,const gArrayInternal<C> & cols,gSize ncols,gPos rstart,gPos rend,gPos cstart,gPos cend) const;
    template <class C> gArrayInternal<T> * getValues(gPos start,gPos end,const gArrayInternal<C> & rows,const gArrayInternal<C> & cols,gSize ncols,gPos rstart,gPos rend,gPos cstart,gPos cend) const;
    template <class C> void castFrom(const gArrayInternal<C> & array,gPos start,gPos end);
    gBool  isNA(gPos pos) const;
    gArrayInternal<gBool> * isNA(gPos start, gPos end) const;
    gSize NACount(gPos start=0, gPos end=0) const;
    const T * getData() const;

    // Logical Operations
    //---------------------
    gArrayInternal<gBool> * greater(const gArrayInternal<T> & array,gPos start,gPos end,gPos astart,gPos aend) const;
    gArrayInternal<gBool> * greaterEqual(const gArrayInternal<T> & array,gPos start,gPos end,gPos astart,gPos aend) const;
    gArrayInternal<gBool> * lesser(const gArrayInternal<T> & array,gPos start,gPos end,gPos astart,gPos aend) const;
    gArrayInternal<gBool> * lesserEqual(const gArrayInternal<T> & array,gPos start,gPos end,gPos astart,gPos aend) const;
    gArrayInternal<gBool> * equal(const gArrayInternal<T> & array,gPos start,gPos end,gPos astart,gPos aend) const;
    gArrayInternal<gBool> * notEqual(const gArrayInternal<T> & array,gPos start,gPos end,gPos astart,gPos aend) const;
    gArrayInternal<gBool> * logicalOr(const gArrayInternal<T> & array,gPos start,gPos end,gPos astart,gPos aend) const;
    gArrayInternal<gBool> * logicalAnd(const gArrayInternal<T> & array,gPos start,gPos end,gPos astart,gPos aend) const;
    gArrayInternal<gBool> * negation(gPos start,gPos end) const;


    gArrayInternal<T>      * add(const gArrayInternal<T> & array,gPos start,gPos end,gPos astart,gPos aend) const;
    gArrayInternal<T>      & add(gArrayInternal<T> & res, const gArrayInternal<T> & array,gPos start,gPos end,gPos astart,gPos aend) const;
    gArrayInternal<T>      * subtract(const gArrayInternal<T> & array,gPos start,gPos end,gPos astart,gPos aend) const;
    gArrayInternal<T>      & subtract(gArrayInternal<T> & res, const gArrayInternal<T> & array,gPos start,gPos end,gPos astart,gPos aend) const;
    gArrayInternal<T>      * multiply(const gArrayInternal<T> & array,gPos start,gPos end,gPos astart,gPos aend) const;
    gArrayInternal<T>      & multiply(gArrayInternal<T> & res, const gArrayInternal<T> & array,gPos start,gPos end,gPos astart,gPos aend) const;
    gArrayInternal<T>      * divide(const gArrayInternal<T> & array,gPos start,gPos end,gPos astart,gPos aend) const;
    gArrayInternal<T>      & divide(gArrayInternal<T> & res, const gArrayInternal<T> & array,gPos start,gPos end,gPos astart,gPos aend) const;    
    
    gSize                    getSize() const;
    gSize                    count(T value, gPos start, gPos end,gSize nacount) const;
    gArrayInternal<gPos>   * sort(bool index_return);
    gArrayInternal<gPos>   * find(T value, gPos start, gPos end,bool excludeNA) const;
    gArrayInternal<gPos>   * match(const gArrayInternal<T> & values, gPos start, gPos end) const;
    gArrayInternal<T>      * getUnique(gPos start,gPos end) const;

    gArrayInternal<gSize>  * getCounts(gArrayInternal<gSize> * res,const gArrayInternal<T> & values,gPos start,gPos end,gPos vstart,gPos vend,bool skipnan,bool countall) const;
    gArrayInternal<T>      * getMin(gArrayInternal<T> * res,gPos start,gPos end,bool skipnan) const;
    gArrayInternal<T>      * getMax(gArrayInternal<T> * res,gPos start,gPos end,bool skipnan) const;
    gArrayInternal<T>      * getSum(gArrayInternal<T> * res,gPos start,gPos end,bool skipnan) const;
    gArrayInternal<T>      * getSquareSum(gArrayInternal<T> * res,gPos start,gPos end,bool skipnan) const;
    gArrayInternal<gScore> * getMean(gArrayInternal<gScore> * res,gPos start,gPos end,bool skipnan) const;
    gArrayInternal<gScore> * getStdDev(gArrayInternal<gScore> * res,gPos start,gPos end,bool skipnan,bool unbiased) const;


    gArrayInternal<gSize>  * getWinCount(gArrayInternal<gSize> * res,T value,gSize winlength,gSize pos,gPos start,gPos end,bool skipnan) const;
    gArrayInternal<T>      * getWinMin(gArrayInternal<T> * res,gSize winlength,gSize pos,gPos start,gPos end,bool skipnan) const;
    gArrayInternal<T>      * getWinMax(gArrayInternal<T> * res,gSize winlength,gSize pos,gPos start,gPos end,bool skipnan) const;
    gArrayInternal<T>      * getWinSum(gArrayInternal<T> * res,gSize winlength,gSize pos,gPos start,gPos end,bool skipnan) const;
    gArrayInternal<T>      * getWinSquareSum(gArrayInternal<T> * res,gSize winlength,gSize pos,gPos start,gPos end,bool skipnan) const;
    gArrayInternal<gScore> * getWinMean(gArrayInternal<gScore> * res,gSize winlength,gSize pos,gPos start,gPos end,bool skipnan) const;
    gArrayInternal<gScore> * getWinStdDev(gArrayInternal<gScore> * res,gSize winlength,gSize pos,gPos start,gPos end,bool skipnan,bool unbiased) const;
};

//--------------------------------------------------------------------------
// Template class gArrayInternal definition
//--------------------------------------------------------------------------
using namespace std;

template <class T> gArrayInternal<T>::gArrayInternal() :i_nans() {
    i_mem=NULL;
    i_length=0;
    i_sorted=false;
    i_refcount=0;
}

template <class T> gArrayInternal<T>::gArrayInternal ( gSize length ) :i_nans ( length ) {
    if ( length>0 ) {
        i_mem = new T[length];
        i_length=length;
        i_sorted= ( length==1 );
    } else {
        i_mem=NULL;
        i_length=0;
        i_sorted=false;
    }
    i_refcount=0;
}

template <class T> gArrayInternal<T>::gArrayInternal ( gSize length,T value,bool isna ) :i_nans ( length,isna ) {
    if ( length>0 ) {
        i_mem = new T[length];
        if ( value==0 ) memset ( i_mem,0,sizeof ( T ) *length );
        else for ( gPos pos=0;pos<length;pos++ ) i_mem[pos]=value;
        i_length=length;
        i_sorted=true;
    } else {
        i_mem=NULL;
        i_length=0;
        i_sorted=false;
    }
    i_refcount=0;
}

template <class T> gArrayInternal<T>::gArrayInternal ( gSize length,T* data,const gBitsArray & isNaFlags,bool getdataownership,bool sorted ) :i_nans ( isNaFlags ) {
    if ( ( data!=NULL ) && ( length>0 ) ) {
        if ( getdataownership ) {
            i_mem=data;
        } else {
            i_mem= new T[length];
            memmove ( i_mem,data,sizeof ( T ) *length );
        }
        i_length=length;
        i_sorted=sorted|| ( length==1 );
    } else {
        i_mem=NULL;
        i_length=0;
        i_sorted=false;
    }
    i_refcount=0;
}

template <class T> gArrayInternal<T>::gArrayInternal ( const gArrayInternal<T> & array,gPos start,gPos end,T fillvalue,bool fillNA ) {
    i_mem=NULL;
    i_sorted=false;
    i_length=0;
    i_refcount=0;

    gPos aend= ( end==0 ) ? ( array.getSize() ) : ( end );
    setLength ( aend-start,fillvalue,fillNA );
    datamove ( *this,0,array,start,i_length );
    i_sorted=array.i_sorted;
}

template <class T> gArrayInternal<T>::~gArrayInternal() {
    if ( i_refcount==0 ) {
        if ( i_mem!=NULL ) delete [] i_mem;
        i_mem=NULL;
        i_sorted=false;
        i_length=0;
    }
}

template <class T> void gArrayInternal<T>::downHeap ( T *a,gPos *b, gPos k, gSize N ) const {
    register gPos newidx,child;
    register T newElt;
    if ( b!=NULL ) {
        newElt=a[k];
        newidx=b[k];
        while ( k <= N/2 ) {
            child = 2*k;
            if ( child < N && a[child] < a[child+1] ) child++;
            if ( newElt >= a[child] ) break;
            //else
            a[k] = a[child];
            b[k] = b[child];
            k = child;
        }
        a[k] = newElt;
        b[k] = newidx;
    } else {
        newElt=a[k];
        while ( k <= N/2 ) {
            child = 2*k;
            if ( child < N && a[child] < a[child+1] ) child++;
            if ( newElt >= a[child] ) break;
            //else
            a[k] = a[child];
            k = child;
        }
        a[k] = newElt;
    }
}

template <class T> void gArrayInternal<T>::internalSort ( T *i_sort,gPos *i_map,gSize N ) const {
    register gPos i,tmp2;
    register T tmp1;

    if ( i_map!=NULL ) {
        for ( i=N/2;i>=1;i-- ) {
            downHeap ( i_sort,i_map,i,N );
        }
        for ( i=N; i >  1; i-- ) {
            tmp1 = i_sort[i];
            tmp2 = i_map[i];
            i_sort[i]=i_sort[1];
            i_map[i]=i_map[1];
            i_sort[1]=tmp1;
            i_map[1]=tmp2;
            downHeap ( i_sort,i_map,1,i-1 );
        }
    } else {
        for ( i=N/2;i>=1;i-- ) {
            downHeap ( i_sort,NULL,i,N );
        }
        for ( i=N; i >  1; i-- ) {
            tmp1 = i_sort[i];
            i_sort[i]=i_sort[1];
            i_sort[1]=tmp1;
            downHeap ( i_sort,NULL,1,i-1 );
        }
    }
}

template <class T> gSize gArrayInternal<T>::internalCount ( T value,gPos start,gPos end,gSize nacount ) const {
    gPos aend=min ( (gPos) (i_length-nacount),end );
    register gPos left=start;
    register gPos right=aend;
    register gPos pos;
    bool found=false;
    gPos np=0;
    while ( left < right ) {
        pos = ( left + right ) / 2;
        if ( i_mem[pos]==value ) {
            found=true;
            break;
        } else {
            if ( i_mem[pos] > value ) right = pos;
            else left = pos + 1;
        }
    }
    if ( found ) {
        while ( ( pos<aend ) && ( i_mem[pos]==value ) ) {
            pos++;
            np++;
        }
    }
    return np;
}

template <class T> gSize gArrayInternal<T>::count ( T value, gPos start, gPos end,gSize nacount ) const {
    gSize count=0;
    if ( i_sorted ) {
        count=internalCount ( value,start,end,nacount );
    } else {
        gArrayInternal<T> sub ( *this,start,end,0,false );
        sub.sort ( false );
        count=sub.internalCount ( value,0,sub.i_length,nacount );
    }
    return count;
}

template <class T> void gArrayInternal<T>::internalSearch ( gArrayInternal<gPos> *positions,T value,const gArrayInternal<T> & array,gPos start, gPos end,bool excludeNA) const {
    //gPos aend=min ( array.i_length-nacount,end );
    gPos aend=min ( array.i_length,end );
    register gPos left=start;
    register gPos right=aend;
    register gPos pos;
    register bool found=false;
    while ( left < right ) {
        pos = ( left + right ) / 2;
        if ( array.i_mem[pos]==value ) {
            found=true;
            break;
        } else {
            if ( array.i_mem[pos] > value ) right = pos;
            else left = pos + 1;
        }
    }
    if ( found ) {
        gPos startpos=pos;
        gPos endpos=pos;
        while ( ( startpos>start ) && ( array.i_mem[startpos-1]==value ) ) {
            startpos--;
        }
        while ( ( endpos<aend ) && ( array.i_mem[endpos]==value ) ) {
          endpos++;
        }
        positions->setLength(endpos-startpos,0,false);
        if(excludeNA){
          for(gPos p=startpos;p<endpos;p++){
            if(!array.i_nans.getValue(p)){
              positions->setValue ( p-startpos,p-start,false,0,false );
            }
          }
        }else{
          for(gPos p=startpos;p<endpos;p++){
            if(!array.i_nans.getValue(p)){
              positions->setValue ( p-startpos,p-start,false,0,false );
            }else{
              positions->setValue ( p-startpos,p-start,true,0,false );
            }
          }
        }
    }
}

template <class T> gArrayInternal<gPos> * gArrayInternal<T>::find ( T value, gPos start, gPos end, bool excludeNA ) const {
    gArrayInternal<gPos> *positions=new gArrayInternal<gPos>();
    if ( i_sorted ) {
        internalSearch ( positions,value,*this,start,end,excludeNA );
    } else {
        gArrayInternal<T> sub ( *this,start,end,0,false );
        gArrayInternal<gPos> *map=sub.sort ( true );
        sub.internalSearch ( positions,value,sub,0,sub.i_length,excludeNA );
        for ( gPos i=0;i<positions->i_length;i++ ) positions->i_mem[i]=map->i_mem[positions->i_mem[i]];
        positions->sort(false);
        delete map;
    }
    return positions;
}


template <class T> gArrayInternal<gPos> * gArrayInternal<T>::match(const gArrayInternal<T> & values, gPos start, gPos end) const{
  gArrayInternal<gPos> * res=new gArrayInternal<gPos>(values.getSize(),0,true);  
  gArrayInternal<T> target(*this,start,end,0,false);
  gArrayInternal<T> query(values,0,values.getSize(),0,false);
  gArrayInternal<gPos> * tord=target.sort(true);
  gArrayInternal<gPos> * qord=query.sort(true);
  gPos lp=0;
  gPos low,mid,high;
  bool found;
  T key,val;
  gSize qsize=query.getSize();
  gSize tsize=target.getSize();
  for(gPos i=0;i< qsize;i++){
    found=false;
    low=lp;
    high=target.getSize()-1;
    key=query.i_mem[i];
    while((high >= low)&&(!found)){
      mid = (low + high) / 2;
      val=target.i_mem[mid];
      if(key < val){
        if(mid>0) high=mid-1;
        else break;
      }else if(key>val){
        if(mid<(tsize-1)) low=mid+1;
        else break;
      }else{
        lp=mid;
        res->i_mem[qord->i_mem[i]]=tord->i_mem[lp];
        if((!query.i_nans.getValue(i))&&(!target.i_nans.getValue(mid))){
         res->i_nans.resetBit(qord->i_mem[i]);
        }
        found=true;
        break;
      }
    }
    lp=low;
  }
  delete tord;
  delete qord;
  return res;
}

template <class T> gArrayInternal<gPos> * gArrayInternal<T>::sort ( bool index_return ) {
  gArrayInternal<gPos> *ret=NULL;
    if ( i_mem ) {
        if ( !i_sorted ) {
                T *i_sort=new T[i_length+1];
                gPos * i_map=new gPos[i_length+1];
                i_sort[0]=0;
                i_map[0]=0;
                memcpy ( i_sort+1,i_mem,sizeof ( T ) *i_length );
                for(gPos i=0;i<i_length;i++) i_map[i+1]=i;
                internalSort ( i_sort,i_map,i_length );
                memmove ( i_mem,i_sort+1,i_length*sizeof ( T ) );
                //the following should be made more efficient
                gBitsArray namap(i_length,false);
                for(gPos i=0;i<i_length;i++){
                  if(i_nans.getValue(i_map[i+1])){
                    namap.setBit(i);
                  }
                }
                 i_nans.copyValues (0,namap,0,i_length);
                if(index_return){
                  ret=new gArrayInternal<gPos>(i_length,i_map+1,gBitsArray(i_length,false),false,false);
                }
                delete [] i_sort;
                delete [] i_map;
                i_sorted=true;
        } else if ( index_return ) {
            ret=new gArrayInternal<gPos>();
            for ( gPos i=0;i<i_length;i++ ) ret->setValue ( i,i,false,0,false );
        }
    }
    return ret;
}

template <class T> void gArrayInternal<T>::setLength ( gSize length,T fillvalue,bool fillNA ) {
    if ( length!=i_length ) {
        i_nans.setSize ( length,true,fillNA );
        if ( length>0 ) {
            T       *newmem=new T[length];
            if ( length>i_length ) {
                if ( ( sizeof ( T ) ==1 ) || ( fillvalue==0 ) ) {
                    memset ( newmem+i_length,(_unsigned_char)fillvalue,sizeof ( T ) * ( length-i_length ) );
                } else {
                    for ( gPos i=i_length;i<length;i++ ) newmem[i]=fillvalue;
                }
            }
            if ( i_mem ) {
                memmove ( newmem,i_mem,sizeof ( T ) *min ( i_length,length ) );
                delete [] i_mem;
            }
            i_mem=newmem;
            if ( i_sorted && ( length>i_length ) ) {
                i_sorted= ( i_mem[i_length-1] <= i_mem[i_length] );
            }
        } else {
            if ( i_mem ) {
                delete [] i_mem;
            }
            i_mem=NULL;
        }
        i_length=length;
    }
}

template <class T> void gArrayInternal<T>::setData ( gSize length,T* data,const gBitsArray & isNAFlags,T fillvalue,bool fillNA,bool getdataownership,bool sorted ) {
    if ( length>0 ) {
        if ( data==NULL ) {
            if ( i_length!=length ) {
                if ( i_mem!=NULL ) delete [] i_mem;
                i_mem = new T[length];
            }
            if ( fillvalue==0 ) memset ( i_mem,0,sizeof ( T ) *length );
            else for ( gPos pos=0;pos<length;pos++ ) i_mem[pos]=fillvalue;
        } else {
            if ( getdataownership ) {
                if ( i_mem ) delete [] i_mem;
                i_mem=data;
            } else {
                if ( ( i_length!=length ) && ( i_mem!=NULL ) ) delete i_mem;
                i_mem= new T[length];
                memmove ( i_mem,data,sizeof ( T ) *length );
            }
        }
        if(isNAFlags.getSize()>0){
          i_nans.copyValues ( 0,isNAFlags,0,length );
        }else{
          if(fillNA) i_nans.setAll();
          else i_nans.resetAll();
        }
        i_length=length;
        i_sorted=sorted|| ( i_length==1 );
    } else {
        i_nans.setSize ( 0 );
        if ( i_mem ) {
            delete [] i_mem;
            //delete [] i_val;
            i_length=0;
        }
    }
}

template <class T> void gArrayInternal<T>::setValue ( gPos pos, T value,bool isna,T fillvalue,bool fillNA ) {
    if ( ! ( pos<i_length ) ) setLength ( pos+1,fillvalue,fillNA );
    i_mem[pos]=value;
    i_nans.setValue ( pos,isna );
    if ( i_sorted ) {
        if ( isna ) {
            i_sorted= ( pos==i_length-1 ) || i_nans.getValue ( pos+1 );
        } else {
            i_sorted= ( ( pos>0 ) ? ( i_mem[pos-1]<=value ) : ( true ) ) && ( ( pos< ( i_length-1 ) ) ? ( value<=i_mem[pos+1] ) : ( true ) );
        }
    }
}

template <class T> void gArrayInternal<T>::setValues (gPos start,gPos pstart,gPos pend,const gArrayInternal<gPos> & positions,gPos vstart,gPos vend,const gArrayInternal<T> & values,T fillvalue,bool fillNA ) {
    gPos dp,vp=vstart;
    for ( gPos i=pstart;i<pend;i++ ) {
      if(vp==vend) vp=vstart;
        if ( !positions.i_nans.getValue ( i ) ) {
            dp=start+positions.i_mem[i];
            setValue ( dp,values.i_mem[vp],values.i_nans.getValue ( vp ),fillvalue,fillNA );        
//             if ( dp>=i_length ) {
//                 setValue ( dp,values.i_mem[vp],values.i_nans.getValue ( vp ),fillvalue,fillNA );
//             } else {
//                 i_mem[dp]=values.i_mem[vp];
//                 i_nans.setValue ( dp,values.i_nans.getValue ( vp ) );
//             }
        }
      vp++;
    }
}

template <class T> void gArrayInternal<T>::datamove ( gArrayInternal<T> & dest,gPos dpos,const gArrayInternal<T> & src,gPos spos,gSize npos ) {
    memmove ( dest.i_mem+dpos,src.i_mem+spos,sizeof ( T ) *npos );
    dest.i_nans.copyValues ( dpos,src.i_nans,spos,spos+npos );
}

template <class T> gArrayInternal<T> & gArrayInternal<T>::selfRevert(gPos start, gPos end) {
    T tempVal;
    bool tmpNA;
    gPos dpos;
    gSize len=end-start;
    T *tmem=i_mem+start;
    for ( gPos i=0; i< ( len/2 ); i++ ) {
        dpos=len-i-1;
        tempVal = tmem[i];
        tmpNA=i_nans.getValue ( i+start );
        tmem[i] = tmem[dpos];
        i_nans.setValue ( i+start,i_nans.getValue ( dpos+start ) );
        tmem[dpos] = tempVal;
        i_nans.setValue ( dpos+start,tmpNA );
    }
    i_sorted=false;
    return *this;
}

template <class T> gArrayInternal<T> & gArrayInternal<T>::selfConcatenate ( const gArrayInternal<T> & array,gPos astart,gPos aend ) {
    return selfInsert ( i_length-1,array,astart,aend );
}

template <class T> gArrayInternal<T> & gArrayInternal<T>::selfReplace ( gPos pos,const gArrayInternal<T> & array,gPos astart,gPos aend ) {
    datamove ( *this,pos,array,astart,aend-astart );
    //Possibly this could be changed in order to keep it marked as sorted when possible
    i_sorted=false;
    return *this;
}

template <class T> gArrayInternal<T> & gArrayInternal<T>::selfInsert ( gPos pos,const gArrayInternal<T> & array,gPos astart,gPos aend ) {
    gSize alength=aend-astart;
    setLength ( i_length+alength,(T) 0,true );
    datamove ( *this,pos+1+alength,*this,pos+1,i_length-pos-1-alength );
    datamove ( *this,pos+1,array,astart,alength );
    //Possibly this could be changed in order to keep it marked as sorted when possible
    i_sorted=false;
    return *this;
}

template <class T> gArrayInternal<T> & gArrayInternal<T>::selfRemove ( gPos pos,gSize npos ) {
    datamove ( *this,pos,*this,pos+npos,i_length-pos-npos );
    setLength ( i_length-npos,0,true );
    return *this;
}

template <class T> gArrayInternal<T> & gArrayInternal<T>::selfAdd ( const gArrayInternal<T> & array,gPos start,gPos end,gPos astart,gPos aend ) {
    gSize i,j;
    gSize l1=end-start;
    gSize l2=aend-astart;
    gSize ntimes = l1/l2;
    gSize rem = l1 % l2;
    //gSize nstart=start+l2*ntimes;
    T *tmp1,*tmp2=array.i_mem + astart;

    if ( i_length==0 ) setLength ( l1,0,false );

    for ( i=0;i<ntimes;i++ ) {
        tmp1 = i_mem+start+i*l2;
        for ( j=0;j<l2;j++ ) tmp1[j]+=tmp2[j];
        i_nans.Or ( start+i*l2,array.i_nans,astart,aend );
    }
    tmp1 = i_mem+start+ntimes*l2;
    for ( i=0;i<rem;i++ ) tmp1[i]+=tmp2[i];
    i_nans.Or ( start+ntimes*l2,array.i_nans,astart,astart+rem );
    i_sorted=i_sorted && array.i_sorted;
    return *this;
}

template <class T> gArrayInternal<T> & gArrayInternal<T>::selfSubtract ( const gArrayInternal<T> & array,gPos start,gPos end,gPos astart,gPos aend ) {
    gSize i,j;
    gSize l1=end-start;
    gSize l2=aend-astart;
    gSize ntimes = l1/l2;
    gSize rem = l1 % l2;
    //gSize nstart=start+l2*ntimes;
    T *tmp1,*tmp2=array.i_mem + astart;

    if ( i_length==0 ) setLength ( l1,0,false );

    for ( i=0;i<ntimes;i++ ) {
        tmp1 = i_mem+start+i*l2;
        for ( j=0;j<l2;j++ ) tmp1[j]-=tmp2[j];
        i_nans.Or ( start+i*l2,array.i_nans,astart,aend );
    }
    tmp1 = i_mem+start+ntimes*l2;
    for ( i=0;i<rem;i++ ) tmp1[i]-=tmp2[i];
    i_nans.Or ( start+ntimes*l2,array.i_nans,astart,astart+rem );
    i_sorted=i_sorted && array.i_sorted;
    return *this;
}

template <class T> gArrayInternal<T> & gArrayInternal<T>::selfMultiply ( const gArrayInternal<T> & array,gPos start,gPos end,gPos astart,gPos aend ) {
    gSize i,j;
    gSize l1=end-start;
    gSize l2=aend-astart;
    gSize ntimes = l1/l2;
    gSize rem = l1 % l2;
    //gSize nstart=start+l2*ntimes;

    T *tmp1,*tmp2=array.i_mem + astart;
    if ( i_length==0 ) setLength ( l1,0,false );

    for ( i=0;i<ntimes;i++ ) {
        tmp1 = i_mem+start+i*l2;
        for ( j=0;j<l2;j++ ) tmp1[j]*=tmp2[j];
        i_nans.Or ( start+i*l2,array.i_nans,astart,aend );
    }
    tmp1 = i_mem+start+ntimes*l2;
    for ( i=0;i<rem;i++ ) tmp1[i]*=tmp2[i];
    i_nans.Or ( start+ntimes*l2,array.i_nans,astart,astart+rem );
    //Possibly this could be changesd in order to keep it marked as sorted when possible
    i_sorted=i_sorted && array.i_sorted;
    return *this;
}

template <class T> gArrayInternal<T> & gArrayInternal<T>::selfDivide ( const gArrayInternal<T> & array,gPos start,gPos end,gPos astart,gPos aend ) {
    gSize i,j;
    gSize l1=end-start;
    gSize l2=aend-astart;
    gSize ntimes = l1/l2;
    gSize rem = l1 % l2;
    //gSize nstart=start+l2*ntimes;

    T *tmp1,*tmp2=array.i_mem + astart;
    if ( i_length==0 ) setLength ( l1,0,false );

    for ( i=0;i<ntimes;i++ ) {
        tmp1 = i_mem+start+i*l2;
        for ( j=0;j<l2;j++ ) tmp1[j]/=tmp2[j];
        i_nans.Or ( start+i*l2,array.i_nans,astart,aend );
    }
    tmp1 = i_mem+start+ntimes*l2;
    for ( i=0;i<rem;i++ ) tmp1[i]/=tmp2[i];
    i_nans.Or ( start+ntimes*l2,array.i_nans,astart,astart+rem );
    //Possibly this could be changesd in order to keep it marked as sorted when possible
    i_sorted=i_sorted && array.i_sorted;
    return *this;
}

template <class T> const T gArrayInternal<T>::at ( gPos pos ) const {
    if ( pos<i_length ) {
        return i_mem[pos];
    } else throw ( gException ( "gArrayInternal: invalid range" ) );
}

template <class T> const T gArrayInternal<T>::operator [] ( gPos pos ) const {
    if ( pos<i_length ) {
        return i_mem[pos];
    } else throw ( gException ( "gArrayInternal: invalid range" ) );
}

template <class T> 
template <class C> gArrayInternal<T> * gArrayInternal<T>::getValues(gPos start,gPos end,const gArrayInternal<C> & pos,gPos pstart,gPos pend) const {
        gSize len=pend-pstart;
        gArrayInternal<T> *ret=new gArrayInternal<T>(len);
        gPos p=pstart;
        T *dest=ret->i_mem;
        //This condition will be (hoepfully useful when getSetBitsCount will be optimizrd
        //if (i_nans.getSetBitsCount(start,end)==0) {
          if (false) {
            if (pos.i_nans.getSetBitsCount(pstart,pend)==0) {
                ret->i_nans.resetAll();
                for (gPos i=0;i<len;i++) {
                    *dest=i_mem[start+(gPos)pos.i_mem[p++]];
                    dest++;
                }
            } else {
                ret->i_nans.copyValues(pstart,pos.i_nans,0,len);
                for (gPos i=0;i<len;i++) {
                    *dest=i_mem[start+(gPos)pos.i_mem[p++]];
                    //ret->i_nans.setValue(i,pos.i_nans.getValue(p++));
                    dest++;
                }
            }
        } else {
            //same as before 
            //if (pos.i_nans.getSetBitsCount(pstart,pend)==0) {
              if (false) {
                for (gPos i=0;i<len;i++) {
                    gPos dp=start+(gPos)pos.i_mem[p++];
                    *dest=i_mem[dp];
                    ret->i_nans.setValue(i,i_nans.getValue(dp));
                    dest++;
                }
            } else {
                for (gPos i=0;i<len;i++) {
                    if(pos.i_nans.getValue(p)){
                     *dest=(T)0;
                     ret->i_nans.setValue(i,true);
                    }else{
                     gPos dp=start+(gPos)pos.i_mem[p];
                     *dest=i_mem[dp];
                     ret->i_nans.setValue(i,i_nans.getValue(dp));
                    }
                    p++;
                    dest++;
                }
            }
        }
        ret->i_sorted=(i_sorted && pos.i_sorted);
        return ret;
    }

template <class T> 
template <class C> gArrayInternal<T> & gArrayInternal<T>::getValues(gArrayInternal<T> & res,gPos start,gPos end,const gArrayInternal<C> & rows,const gArrayInternal<C> & cols,gSize ncols,gPos rstart,gPos rend,gPos cstart,gPos cend) const {
        if (res.i_length==rend-rstart) {
            gSize len=rend-rstart;
            gPos p;
            C *row=rows.i_mem+rstart;
            C *col=cols.i_mem+cstart;

            gSize rn=rows.i_nans.getSetBitsCount(0,len);
            gSize cn=cols.i_nans.getSetBitsCount(0,len);
            gSize dn=i_nans.getSetBitsCount(start,end);
            if (dn>0) {
                if ((rn==0)&&(cn==0)) {
                    for (gPos i=0;i<len;i++) {
                        p=(gPos)(row[i]*ncols+col[i]);
                        res.i_mem[i]=i_mem[(gPos)(row[i]*ncols+col[i])];
                        res.i_nans.setValue(i,i_nans.getValue(p));
                    }
                } else if (rn==0) {
                    gSize cpos=cstart;
                    for (gPos i=0;i<len;i++) {
                        p=(gPos)(row[i]*ncols+col[i]);
                        res.i_mem[i]=i_mem[p];
                        res.i_nans.setValue(i,(cols.i_nans.getValue(cpos++))||(i_nans.getValue(p)));
                    }
                } else if (cn==0) {
                    gSize rpos=rstart;
                    for (gPos i=0;i<len;i++) {
                        p=(gPos)(row[i]*ncols+col[i]);
                        res.i_mem[i]=i_mem[p];
                        res.i_nans.setValue(i,(rows.i_nans.getValue(rpos++))||(i_nans.getValue(p)));
                    }
                } else {
                    gSize rpos=rstart,cpos=cstart;
                    for (gPos i=0;i<len;i++) {
                        p=(gPos)(row[i]*ncols+col[i]);
                        res.i_mem[i]=i_mem[p];
                        res.i_nans.setValue(i,(cols.i_nans.getValue(cpos++))||(rows.i_nans.getValue(rpos++))||(i_nans.getValue(p)));
                    }
                }
            } else {
                if ((rn==0)&&(cn==0)) {
                    res.i_nans.resetAll();
                    for (gPos i=0;i<len;i++) {
                        res.i_mem[i]=i_mem[(gPos)(row[i]*ncols+col[i])];
                    }
                } else if (rn==0) {
                    gSize cpos=cstart;
                    for (gPos i=0;i<len;i++) {
                        res.i_mem[i]=i_mem[(gPos)(row[i]*ncols+col[i])];
                        res.i_nans.setValue(i,cols.i_nans.getValue(cpos++));
                    }
                } else if (cn==0) {
                    gSize rpos=rstart;
                    for (gPos i=0;i<len;i++) {
                        res.i_mem[i]=i_mem[(gPos)(row[i]*ncols+col[i])];
                        res.i_nans.setValue(i,rows.i_nans.getValue(rpos++));
                    }
                } else {
                    gSize rpos=rstart,cpos=cstart;
                    for (gPos i=0;i<len;i++) {
                        res.i_mem[i]=i_mem[(gPos)(row[i]*ncols+col[i])];
                        res.i_nans.setValue(i,(cols.i_nans.getValue(cpos++))||(rows.i_nans.getValue(rpos++)));
                    }
                }
            }
        }
        return res;
    }
    
template <class T>
template <class C> gArrayInternal<T> * gArrayInternal<T>::getValues(gPos start,gPos end,const gArrayInternal<C> & rows,const gArrayInternal<C> & cols,gSize ncols,gPos rstart,gPos rend,gPos cstart,gPos cend) const {
        gArrayInternal<T> *ret=new gArrayInternal<T>(rend-rstart);
        getValues(*ret,start,end,rows,cols,ncols,rstart,rend,cstart,cend);
        return ret;
    }
    
template <class T>    
template <class C> void gArrayInternal<T>::castFrom(const gArrayInternal<C> & array,gPos start,gPos end) {
        gPos aend=(end==0)?(array.i_length):(end);
        i_length=aend-start;
        i_sorted=array.i_sorted;
        //i_mem=new T[i_length];
        for (gPos i=0;i<i_length;i++) i_mem[i]=(T)array.i_mem[i+start];
        i_nans.copyValues(0,array.i_nans,start,aend);
    }


template <class T> gBool gArrayInternal<T>::isNA ( gPos pos ) const {
    return i_nans.getValue ( pos );
}

template <class T> gArrayInternal<gBool> * gArrayInternal<T>::isNA(gPos start,gPos end) const {
    gArrayInternal<gBool> *ret=new gArrayInternal<gBool> (end-start,false,false);
    for(gPos i=start;i<end;i++){
      ret->i_mem[i]=i_nans.getValue(i);
    }
    return ret;
}


template <class T> gSize gArrayInternal<T>::NACount ( gPos start, gPos end ) const {
    gPos aend= ( end==0 ) ? ( getSize() ) : ( end );
    return i_nans.getSetBitsCount ( start,aend );
}

template <class T> const T * gArrayInternal<T>::getData() const {
    return i_mem;
}

template <class T> gSize gArrayInternal<T>::getSize() const {
    return i_length;
}

template <class T> gArrayInternal<T> * gArrayInternal<T>::add ( const gArrayInternal<T> & array,gPos start,gPos end,gPos astart,gPos aend ) const {
    gArrayInternal<T> *ret=new gArrayInternal<T> ( *this,start,end,0,false );
    add ( *ret,array,start,end,astart,aend );
    return ret;
}

template <class T> gArrayInternal<T> & gArrayInternal<T>::add ( gArrayInternal & res, const gArrayInternal<T> & array,gPos start,gPos end,gPos astart,gPos aend ) const {
    gSize i,j,k=0;
    gSize l1=end-start;
    gSize l2=aend-astart;
    gSize ntimes = l1/l2;
    gSize rem = l1 % l2;

    T *tmp1=res.i_mem;
    T *tmp2=i_mem+start;
    T *tmp3=array.i_mem + astart;

    for ( i=0;i<ntimes;i++ ) {
        for ( j=0;j<l2;j++ ) tmp1[j]=tmp2[j]+tmp3[j];
        res.i_nans.Or ( k,array.i_nans,astart,aend );
        tmp1+=l2;
        tmp2+=l2;
        k+=l2;
    }
    for ( i=0;i<rem;i++ ) tmp1[i]=tmp2[i]+tmp3[i];
    res.i_nans.Or ( k,array.i_nans,astart,astart+rem );
    return res;
}

template <class T> gArrayInternal<T> * gArrayInternal<T>::subtract ( const gArrayInternal<T> & array,gPos start,gPos end,gPos astart,gPos aend ) const {
    gArrayInternal<T> *ret=new gArrayInternal<T> ( *this,start,end,0,false );
    subtract ( *ret,array,start,end,astart,aend );
    return ret;
}

template <class T> gArrayInternal<T> & gArrayInternal<T>::subtract ( gArrayInternal & res, const gArrayInternal<T> & array,gPos start,gPos end,gPos astart,gPos aend ) const {
    gSize i,j,k=0;
    gSize l1=end-start;
    gSize l2=aend-astart;
    gSize ntimes = l1/l2;
    gSize rem = l1 % l2;

    T *tmp1=res.i_mem;
    T *tmp2=i_mem+start;
    T *tmp3=array.i_mem + astart;

    for ( i=0;i<ntimes;i++ ) {
        for ( j=0;j<l2;j++ ) tmp1[j]=tmp2[j]-tmp3[j];
        res.i_nans.Or ( k,array.i_nans,astart,aend );
        tmp1+=l2;
        tmp2+=l2;
        k+=l2;
    }
    for ( i=0;i<rem;i++ ) tmp1[i]=tmp2[i]-tmp3[i];
    res.i_nans.Or ( k,array.i_nans,astart,astart+rem );
    if(l2==1) res.i_sorted=i_sorted;
    else res.i_sorted=false;
    return res;
}

template <class T> gArrayInternal<T> * gArrayInternal<T>::multiply ( const gArrayInternal<T> & array,gPos start,gPos end,gPos astart,gPos aend ) const {
    gArrayInternal<T> *ret=new gArrayInternal<T> ( *this,start,end,0,false );
    multiply ( *ret,array,start,end,astart,aend );
    return ret;
}

template <class T> gArrayInternal<T> & gArrayInternal<T>::multiply ( gArrayInternal & res, const gArrayInternal<T> & array,gPos start,gPos end,gPos astart,gPos aend ) const {
    gSize i,j,k=0;
    gSize l1=end-start;
    gSize l2=aend-astart;
    gSize ntimes = l1/l2;
    gSize rem = l1 % l2;

    T *tmp1=res.i_mem;
    T *tmp2=i_mem+start;
    T *tmp3=array.i_mem + astart;

    for ( i=0;i<ntimes;i++ ) {
        for ( j=0;j<l2;j++ ) tmp1[j]=tmp2[j]*tmp3[j];
        res.i_nans.Or ( k,array.i_nans,astart,aend );
        tmp1+=l2;
        tmp2+=l2;
        k+=l2;
    }
    for ( i=0;i<rem;i++ ) tmp1[i]=tmp2[i]*tmp3[i];
    res.i_nans.Or ( k,array.i_nans,astart,astart+rem );
    return res;
}

template <class T> gArrayInternal<T> * gArrayInternal<T>::divide(const gArrayInternal<T> & array,gPos start,gPos end,gPos astart,gPos aend) const{
    gArrayInternal<T> *ret=new gArrayInternal<T> ( end-start,0,false );
    //gArrayInternal<T> *ret=new gArrayInternal<T> ( *this );
    divide ( *ret,array,start,end,astart,aend );
    return ret;
}

template <class T> gArrayInternal<T> & gArrayInternal<T>::divide(gArrayInternal<T> & res, const gArrayInternal<T> & array,gPos start,gPos end,gPos astart,gPos aend) const{
    gSize i,j,k=0;
    gSize l1=end-start;
    gSize l2=aend-astart;
    gSize ntimes = l1/l2;
    gSize rem = l1 % l2;
    //res.i_nans=i_nans;

    T *tmp1=res.i_mem;
    T *tmp2=i_mem+start;
    T *tmp3=array.i_mem + astart;

    for ( i=0;i<ntimes;i++ ) {
        for ( j=0;j<l2;j++ ) tmp1[j]=  tmp2[j] / tmp3[j];
        res.i_nans.Or ( k,array.i_nans,astart,aend );
        tmp1+=l2;
        tmp2+=l2;
        k+=l2;
    }
    for ( i=0;i<rem;i++ ) tmp1[i]=  tmp2[i] / tmp3[i];
    res.i_nans.Or ( k,array.i_nans,astart,astart+rem );
    return res;
}

template <class T> gArrayInternal<gBool> * gArrayInternal<T>::greater ( const gArrayInternal<T> & array,gPos start,gPos end,gPos astart,gPos aend ) const {
    gSize i,j;
    gSize l1=end-start;
    gSize l2=aend-astart;
    gArrayInternal<gBool> *ret;
    T * tval=i_mem+start;
    T * aval=array.i_mem+astart;
    if(l2>1){
     ret=new gArrayInternal<gBool> ( l1 ); 
     ret->i_nans.copyValues( 0,i_nans,start,end );     
     if ( l1==l2 ) {
         for ( i=0;i<l1;i++ ) ret->i_mem[i] = tval[i] > aval[i];
         ret->i_nans.Or ( start,array.i_nans,astart,aend );
     } else {
         for ( i=0;i<l1;i+=l2 ) {
             for ( j=0;j<l2;j++ ) ret->i_mem[i+j] = tval[i+j] > aval[j];
             ret->i_nans.Or ( i+j,array.i_nans,astart,aend );
         }
         gSize nstart=l2*(l1/l2);
         for ( i=nstart;i<l1;i++ ) ret->i_mem[i] = i_mem[i] > array.i_mem[i-nstart];
         ret->i_nans.Or ( nstart,array.i_nans,astart,astart+(l1%l2) );
     }
    }else{
      if(i_sorted){
       ret=new gArrayInternal<gBool> ( l1,true,false );
       ret->i_nans.copyValues( 0,i_nans,start,end );
       i=0;
       while((i<l1) && (tval[i] <= aval[0])){
        ret->i_mem[i++] =false;
       }
      }else{
       ret=new gArrayInternal<gBool> ( l1 );
       ret->i_nans.copyValues( 0,i_nans,start,end );
       for ( i=0;i<l1;i++ ){
         ret->i_mem[i] = tval[i] > aval[0];
       }
      }
    }
    return ret;
}

template <class T> gArrayInternal<gBool> * gArrayInternal<T>::greaterEqual ( const gArrayInternal<T> & array,gPos start,gPos end,gPos astart,gPos aend ) const {
    gSize i,j;
    gSize l1=end-start;
    gSize l2=aend-astart;
    T * tval=i_mem+start;
    T * aval=array.i_mem+astart;
    gArrayInternal<gBool> *ret;
    if(l2>1){
      ret=new gArrayInternal<gBool> ( l1 );
      ret->i_nans.copyValues( 0,i_nans,start,end );
      if ( l1==l2 ) {
          for ( i=0;i<l1;i++ ) ret->i_mem[i] = tval[i] >= aval[i];
          ret->i_nans.Or ( start,array.i_nans,astart,aend );
      } else {
          for ( i=0;i<l1;i+=l2 ) {
              for ( j=0;j<l2;j++ ) ret->i_mem[i+j] = tval[i+j] >= aval[j];
              ret->i_nans.Or ( i+j,array.i_nans,astart,aend );
          }
          gSize nstart=l2*(l1/l2);
          for ( i=nstart;i<l1;i++ ) ret->i_mem[i] = i_mem[i] >= array.i_mem[i-nstart];
          ret->i_nans.Or ( nstart,array.i_nans,astart,astart+(l1%l2) );
      }
    }else{
      if(i_sorted){
       ret=new gArrayInternal<gBool> ( l1,true,false );
       ret->i_nans.copyValues( 0,i_nans,start,end );
       i=0;
       while((i<l1) && (tval[i] < aval[0])){
        ret->i_mem[i++] =false;
       }
      }else{
       ret=new gArrayInternal<gBool> ( l1 );
       ret->i_nans.copyValues( 0,i_nans,start,end );
       for ( i=0;i<l1;i++ ){
         ret->i_mem[i] = tval[i] >= aval[0];
       }
      }
    }
    return ret;
}

template <class T> gArrayInternal<gBool> * gArrayInternal<T>::lesser ( const gArrayInternal<T> & array,gPos start,gPos end,gPos astart,gPos aend ) const {
    gSize i,j;
    gSize l1=end-start;
    gSize l2=aend-astart;
    T * tval=i_mem+start;
    T * aval=array.i_mem+astart;
    gArrayInternal<gBool> *ret;    
    if(l2>1){
      ret=new gArrayInternal<gBool> ( l1 );
      ret->i_nans.copyValues( 0,i_nans,start,end );
      if ( l1==l2 ) {
          for ( i=0;i<l1;i++ ) ret->i_mem[i] = tval[i] < aval[i];
          ret->i_nans.Or ( start,array.i_nans,astart,aend );
      } else {
          for ( i=0;i<l1;i+=l2 ) {
              for ( j=0;j<l2;j++ ) ret->i_mem[i+j] = tval[i+j] < aval[j];
              ret->i_nans.Or ( i+j,array.i_nans,astart,aend );
          }
          gSize nstart=l2*(l1/l2);
          for ( i=nstart;i<l1;i++ ) ret->i_mem[i] = i_mem[i] < array.i_mem[i-nstart];
          ret->i_nans.Or ( nstart,array.i_nans,astart,astart+(l1%l2) );
      }
    }else{
      if(i_sorted){
       ret=new gArrayInternal<gBool> ( l1,false,false );
       ret->i_nans.copyValues( 0,i_nans,start,end );
       i=0;
       while((i<l1) && (tval[i] <= aval[0])){
        ret->i_mem[i++] = true;
       }
      }else{
       ret=new gArrayInternal<gBool> ( l1 );
       ret->i_nans.copyValues( 0,i_nans,start,end );
       for ( i=0;i<l1;i++ ){
         ret->i_mem[i] = tval[i] < aval[0];
       }
      }
    }
    return ret;
}

template <class T> gArrayInternal<gBool> * gArrayInternal<T>::lesserEqual ( const gArrayInternal<T> & array,gPos start,gPos end,gPos astart,gPos aend ) const {
    gSize i,j;
    gSize l1=end-start;
    gSize l2=aend-astart;
    T * tval=i_mem+start;
    T * aval=array.i_mem+astart;
    gArrayInternal<gBool> *ret;    
    if(l2>1){
      ret=new gArrayInternal<gBool> ( l1 );
      ret->i_nans.copyValues( 0,i_nans,start,end );      
      if ( l1==l2 ) {
          for ( i=0;i<l1;i++ ) ret->i_mem[i] = tval[i] <= aval[i];
          ret->i_nans.Or ( start,array.i_nans,astart,aend );
      } else {
          for ( i=0;i<l1;i+=l2 ) {
              for ( j=0;j<l2;j++ ) ret->i_mem[i+j] = tval[i+j] <= aval[j];
              ret->i_nans.Or ( i+j,array.i_nans,astart,aend );
          }
          gSize nstart=l2*(l1/l2);
          for ( i=nstart;i<l1;i++ ) ret->i_mem[i] = i_mem[i] <= array.i_mem[i-nstart];
          ret->i_nans.Or ( nstart,array.i_nans,astart,astart+(l1%l2) );
      }
    }else{
      if(i_sorted){
       ret=new gArrayInternal<gBool> ( l1,false,false );
       ret->i_nans.copyValues( 0,i_nans,start,end );
       i=0;
       while((i<l1) && (tval[i] <= aval[0])){
        ret->i_mem[i++] = true;
       }
      }else{
       ret=new gArrayInternal<gBool> ( l1 );
       ret->i_nans.copyValues( 0,i_nans,start,end );
       for ( i=0;i<l1;i++ ){
         ret->i_mem[i] = tval[i] <= aval[0];
       }
      }
    }
    return ret;
}

template <class T> gArrayInternal<gBool> * gArrayInternal<T>::equal ( const gArrayInternal<T> & array,gPos start,gPos end,gPos astart,gPos aend ) const {
    gSize i,j;
    gSize l1=end-start;
    gSize l2=aend-astart;
    gArrayInternal<gBool> *ret=new gArrayInternal<gBool> ( l1 );
    ret->i_nans.copyValues( 0,i_nans,start,end );
    T * tval=i_mem+start;
    T * aval=array.i_mem+astart;
    
    if ( l1==l2 ) {
        for ( i=0;i<l1;i++ ) ret->i_mem[i] = tval[i] == aval[i];
        ret->i_nans.Or ( start,array.i_nans,astart,aend );
    } else {
        for ( i=0;i<l1;i+=l2 ) {
            for ( j=0;j<l2;j++ ) ret->i_mem[i+j] = tval[i+j] == aval[j];
            ret->i_nans.Or ( i+j,array.i_nans,astart,aend );
        }
        gSize nstart=l2*(l1/l2);
        for ( i=nstart;i<l1;i++ ) ret->i_mem[i] = i_mem[i] == array.i_mem[i-nstart];
        ret->i_nans.Or ( nstart,array.i_nans,astart,astart+(l1%l2) );
    }
    return ret;
}

template <class T> gArrayInternal<gBool> * gArrayInternal<T>::notEqual ( const gArrayInternal<T> & array,gPos start,gPos end,gPos astart,gPos aend ) const {
    gSize i,j;
    gSize l1=end-start;
    gSize l2=aend-astart;
    T * tval=i_mem+start;
    T * aval=array.i_mem+astart;
    gArrayInternal<gBool> *ret=new gArrayInternal<gBool> ( l1 );
    ret->i_nans.copyValues ( 0,i_nans,start,end );
    
    if ( l1<=l2 ) {
        for ( i=0;i<l1;i++ ) ret->i_mem[i] = tval[i] != aval[i];
        ret->i_nans.Or(start,array.i_nans,astart,astart+l1);
    } else {
        for ( i=0;i<l1;i+=l2 ) {
            for ( j=0;j<l2;j++ ) ret->i_mem[i+j] = tval[i+j] != aval[j];
            ret->i_nans.Or ( i+j,array.i_nans,astart,aend );
        }
        gSize nstart=l2*(l1/l2);
        for ( i=nstart;i<l1;i++ ) ret->i_mem[i] = tval[i] != aval[i-nstart];
        ret->i_nans.Or (nstart,array.i_nans,astart,astart+(l1 % l2) );
    }
    return ret;
}

template <class T> gArrayInternal<gBool> * gArrayInternal<T>::logicalOr ( const gArrayInternal<T> & array,gPos start,gPos end,gPos astart,gPos aend ) const {
    gSize i,j;
    gSize l1=end-start;
    gSize l2=aend-astart;
    gArrayInternal<gBool> *ret=new gArrayInternal<gBool> ( l1 );
    ret->i_nans.copyValues( 0,i_nans,start,end );
    T * tval=i_mem+start;
    T * aval=array.i_mem+astart;
    
    if ( l1==l2 ) {
        for ( i=0;i<l1;i++ ) ret->i_mem[i] = tval[i] || aval[i];
        ret->i_nans.Or ( start,array.i_nans,astart,aend );
    } else {
        for ( i=0;i<l1;i+=l2 ) {
            for ( j=0;j<l2;j++ ) ret->i_mem[i+j] = tval[i+j] || aval[j];
            ret->i_nans.Or ( i+j,array.i_nans,astart,aend );
        }
        gSize nstart=l2*(l1/l2);
        for ( i=nstart;i<l1;i++ ) ret->i_mem[i] = i_mem[i] || array.i_mem[i-nstart];
        ret->i_nans.Or ( nstart,array.i_nans,astart,astart+(l1%l2) );
    }
    return ret;
}

template <class T> gArrayInternal<gBool> * gArrayInternal<T>::logicalAnd ( const gArrayInternal<T> & array,gPos start,gPos end,gPos astart,gPos aend ) const {
    gSize i,j;
    gSize l1=end-start;
    gSize l2=aend-astart;
    gArrayInternal<gBool> *ret=new gArrayInternal<gBool> ( l1 );
    ret->i_nans.copyValues( 0,i_nans,start,end );
    T * tval=i_mem+start;
    T * aval=array.i_mem+astart;
    
    if ( l1==l2 ) {
        for ( i=0;i<l1;i++ ) ret->i_mem[i] = tval[i] && aval[i];
        ret->i_nans.Or ( start,array.i_nans,astart,aend );
    } else {
        for ( i=0;i<l1;i+=l2 ) {
            for ( j=0;j<l2;j++ ) ret->i_mem[i+j] = tval[i+j] && aval[j];
            ret->i_nans.Or ( i+j,array.i_nans,astart,aend );
        }
        gSize nstart=l2*(l1/l2);
        for ( i=nstart;i<l1;i++ ) ret->i_mem[i] = i_mem[i] && array.i_mem[i-nstart];
        ret->i_nans.Or ( nstart,array.i_nans,astart,astart+(l1%l2) );
    }
    return ret;
}

template <class T> gArrayInternal<gBool> * gArrayInternal<T>::negation ( gPos start,gPos end ) const {
    gArrayInternal<gBool> *ret=new gArrayInternal<gBool> ( end-start );
    ret->i_nans.copyValues ( 0,i_nans,start,end );
    gSize i;
    for ( i=0;i<ret->i_length;i++ ) ret->i_mem[i]= ( gBool ) ( i_mem[i+start]==0 );
    return ret;
}

template <class T> gArrayInternal<T> * gArrayInternal<T>::getUnique (gPos start,gPos end ) const {
    gArrayInternal<T> * ret;
    gPos cpos=start;
    gPos npos=0;
    bool found;
    while ( ( cpos<end ) && ( i_nans.getValue ( cpos ) ) ) {
        cpos++;
    }
    if ( cpos<end ) {
        T *data=new T[end-start];
        data[npos]=i_mem[cpos];
        while ( cpos<end ) {
            found=false;
            if ( !i_nans.getValue ( cpos ) ) {
                for ( gPos ipos=0;ipos<=npos;ipos++ ) {
                    if ( i_mem[cpos]==data[ipos] ) {
                        found=true;
                        break;
                    }
                }
            } else found=true;
            if ( !found ) {
                npos++;
                data[npos]=i_mem[cpos];
            }
            cpos++;
        }
        ret=new gArrayInternal ( npos+1,data,gBitsArray ( npos+1,false ),true,i_sorted );
        //delete [] data;
    } else {
        ret=new gArrayInternal ( 1,0,true );
    }
    return ret;
}

template <class T> gArrayInternal<gSize> * gArrayInternal<T>::getCounts ( gArrayInternal<gSize> * res, const gArrayInternal<T> & values,gPos start,gPos end,gPos vstart,gPos vend,bool skipnan,bool countall ) const {
    //bool isnan;
    gArrayInternal<gPos> * pos;
    gSize nacount=values.NACount ( vstart,vend );
    gSize tnacount=NACount ( start,end );
    if ( skipnan|| ( tnacount==0 ) ) {
        if ( ( !countall ) && ( ( vend-vstart ) >1 ) ) {
            //should be initialized???
            //res=new gArrayInternal<gSize> ( vend-vstart );
            if ( values.i_sorted ) {
                for ( gPos i=start;i<end;i++ ) {
                    if ( !i_nans.getValue ( i ) ) {
                        pos=values.find ( i_mem[i],vstart,vend,true );
                        for ( gPos j=0;j<pos->getSize();j++ ) {
                            res->i_mem[pos->i_mem[j]]++;
                        }
                        delete pos;
                    }
                }
            } else {
                gArrayInternal<T> ovalues ( values,vstart,vend,0,false );
                gArrayInternal<gPos> *map=ovalues.sort ( true );
                for ( gPos i=start;i<end;i++ ) {
                    if ( !i_nans.getValue ( i ) ) {
                        pos=ovalues.find ( i_mem[i],0,ovalues.i_length,true );
                        for ( gPos j=0;j<pos->getSize();j++ ) {
                            res->i_mem[map->i_mem[pos->i_mem[j]]]++;
                        }
                        delete pos;
                    }
                }
                delete map;
            }
        } else {
            //should be initialized???
            //ret=new gArrayInternal<gSize> ( 1 );
            if ( ( vstart-vend ) ==1 ) {
                res->i_mem[0]=count ( values[vstart],start,end,tnacount );
            } else {
                res->i_mem[0]=0;
                if ( values.i_sorted ) {
                    res->i_mem[0]=0;
                    for ( gPos i=start;i<end;i++ ) {
                        if ( !i_nans.getValue ( i ) ) {
                            res->i_mem[0]+=values.count ( i_mem[i],vstart,vend,nacount );
                        }
                    }
                } else {
                    gArrayInternal<T> ovalues ( values,vstart,vend,0,false );
                    ovalues.sort ( false );
                    for ( gPos i=start;i<end;i++ ) {
                        if ( !i_nans.getValue ( i ) ) {
                            res->i_mem[0]+=ovalues.count ( i_mem[i],0,ovalues.i_length,nacount );
                        }
                    }
                }
            }
        }
    } else {
        if ( countall ) res->setValue(0,0,true,0,true);
        else res->setValue(vend-vstart-1,0,true,0,true);
    }
    return res;
}

template <class T> gArrayInternal<T> * gArrayInternal<T>::getMin ( gArrayInternal<T> * res, gPos start,gPos end,bool skipnan ) const {
    gPos pos=start;
    bool isnan = i_nans.getValue ( pos );

    while ( skipnan && isnan && ( pos<end ) ) {
        pos++;
        if ( pos<end ) isnan=i_nans.getValue ( pos );
    }
    if ( pos<end ) {
        if ( !isnan ) {
            if ( skipnan ) {
                for ( gPos i=pos+1; i<end; i++ ) {
                    if ( !i_nans.getValue ( i ) && i_mem[i]<i_mem[pos] ) {
                        pos=i;
                        isnan=i_nans.getValue ( i );
                    }
                }
            } else {
                for ( gPos i=pos+1; i<end; i++ ) {
                    isnan=i_nans.getValue ( i );
                    if ( isnan ) break;
                    else {
                        if ( i_mem[i]<i_mem[pos] ) pos=i;
                    }
                }
            }
        } else pos--;
    }
    if ( res->i_length!=1 ) res->setLength ( 1,i_mem[pos],isnan );
    else {
        res->i_mem[0]=i_mem[pos];
        res->i_nans.setValue ( 0,isnan );
    }
    return res;

}

template <class T> gArrayInternal<T> * gArrayInternal<T>::getMax ( gArrayInternal<T> * res, gPos start,gPos end,bool skipnan ) const {
    gPos pos=start;
    bool isnan=i_nans.getValue ( pos );
    while ( skipnan && isnan && ( pos<end ) ) {
        pos++;
        if ( pos<end ) isnan=i_nans.getValue ( pos );
    }
    if ( pos<end ) {
        if ( !isnan ) {
            if ( skipnan ) {
                for ( gPos i=pos+1; i<end; i++ ) {
                    if ( !i_nans.getValue ( i ) && i_mem[i]>i_mem[pos] ) {
                        pos=i;
                        isnan=i_nans.getValue ( i );
                    }
                }
            } else {
                for ( gPos i=pos+1; i<end; i++ ) {
                    isnan=i_nans.getValue ( i );
                    if ( isnan ) break;
                    else {
                        if ( i_mem[i]>i_mem[pos] ) pos=i;
                    }
                }
            }
        } else pos--;
    }
    if ( res->i_length!=1 ) res->setLength ( 1,i_mem[pos],isnan );
    else {
        res->i_mem[0]=i_mem[pos];
        res->i_nans.setValue ( 0,isnan );
    }
    return res;
}

template <class T> gArrayInternal<T> * gArrayInternal<T>::getSum ( gArrayInternal<T> * res,gPos start,gPos end,bool skipnan ) const {
    T sum=0;
    bool isnan=false;
    gSize nacounts=i_nans.getSetBitsCount ( start,end );
    if ( !skipnan ) {
        isnan= ( nacounts>0 );
        if ( !isnan ) for ( gPos i=start;i<end;i++ ) sum+=i_mem[i];
    } else {
        if ( nacounts==0 ) {
            for ( gPos i=start;i<end;i++ ) sum+=i_mem[i];
        } else if ( nacounts== ( end-start ) ) {
            isnan=true;
        } else {
            //this could be emproved (try finding na positions)
            for ( gPos i=start;i<end;i++ ) if ( !i_nans.getValue ( i ) ) sum+=i_mem[i];
        }
    }
    if ( res->i_length!=1 ){
      res->setLength(1,sum,isnan);
    } else {
        res->i_mem[0]=sum;
        res->i_nans.setValue ( 0,isnan );
    }
    return res;
}

template <class T> gArrayInternal<T> * gArrayInternal<T>::getSquareSum ( gArrayInternal<T> * res, gPos start,gPos end,bool skipnan ) const {
    T ssum=0;
    bool isnan=false;
    gSize nacounts=i_nans.getSetBitsCount ( start,end );
    if ( !skipnan ) {
        isnan= ( nacounts==0 );
        if ( !isnan ) for ( gPos i=start;i<end;i++ ) ssum+=i_mem[i]*i_mem[i];
    } else {
        if ( nacounts==0 ) {
            for ( gPos i=start;i<end;i++ ) ssum+=i_mem[i]*i_mem[i];
        } else if ( nacounts== ( end-start ) ) {
            isnan=true;
        } else {
            //this could be improved (try finding na positions??)
            for ( gPos i=start;i<end;i++ ) if ( !i_nans.getValue ( i ) ) ssum+=i_mem[i]*i_mem[i];
        }
    }
    if ( res->i_length!=1 ) res->setLength ( 1,ssum,isnan );
    else {
        res->i_mem[0]=ssum;
        res->i_nans.setValue ( 0,isnan );
    }
    return res;
}

template <class T> gArrayInternal<gScore> * gArrayInternal<T>::getMean ( gArrayInternal<gScore> * res, gPos start,gPos end,bool skipnan ) const {
    bool isnan=true;
    gScore mean=0;
    gArrayInternal<T> *sum=getSum (new gArrayInternal<T>( 1 ), start,end,skipnan );
    gSize nans=NACount( start,end );
    gScore n= ( gScore ) ( end-start- ( skipnan?nans:0 ) );
    if ( ( !sum->isNA ( 0 ) ) && ( n>0 ) ) {
        mean=(gScore) sum->at ( 0 ) / (gScore) n;
        isnan=false;
    }
    delete sum;
    if ( res->i_length!=1 ) res->setLength ( 1,mean,isnan );
    else {
        res->i_mem[0]=mean;
        res->i_nans.setValue ( 0,isnan );
    }
    return res;
}

template <class T> gArrayInternal<gScore> * gArrayInternal<T>::getStdDev ( gArrayInternal<gScore> * res, gPos start,gPos end,bool skipnan,bool unbiased ) const {
    bool isnan=true;
    gScore stddev=0;
    gArrayInternal<T> * sum=getSum ( new gArrayInternal<T> ( 1 ),start,end,skipnan );
    gArrayInternal<T> * ssum=getSquareSum ( new gArrayInternal<T> ( 1 ),start,end,skipnan );
    gScore n = ( gScore ) end-start- ( skipnan?NACount ( start,end ):0 );
    if ( !sum->isNA ( 0 ) ) {
        stddev=((gScore)1/(gScore)n) * sqrt ( (gScore) n * (gScore) ssum->at ( 0 ) - (gScore) sum->at(0) * (gScore) sum->at(0) );
        if (unbiased) {
            if (n>1) {
                stddev*=(gScore) (n-1) / (gScore) n;
                isnan=false;
            }
        } else isnan=false;
    }
    delete ssum;
    delete sum;
    if ( res->i_length!=1 ) res->setLength ( 1,stddev,isnan );
    else {
        res->i_mem[0]=stddev;
        res->i_nans.setValue ( 0,isnan );
    }
    return res;
}

template <class T> gArrayInternal<gSize> * gArrayInternal<T>::getWinCount ( gArrayInternal<gSize> * res,  T value,gSize winlength,gSize pos,gPos start,gPos end,bool skipnan ) const {
    bool ok;
    register gPos first=start,last=start;
    register gSize count=0;
    //gArrayInternal<gSize> *ret=new gArrayInternal<gSize> ( end-start );
    if ( !skipnan ) {
        while ( last<end ) {
            ok=!i_nans.getValue ( last );
            while ( ( last-first< ( winlength-1 ) ) ) {
                if ( !ok ) break;
                count+= ( i_mem[last]==value );
                last++;
                ok=!i_nans.getValue ( last );
            }
            if ( !ok ) {
                res->i_nans.setValue ( first+pos,true );
                first=last+1;
                last=first;
                if ( winlength> ( end-first ) ) break;
                count=0;
            } else {
                count+= ( i_mem[last]==value );
                res->i_mem[first+pos]=count;
                res->i_nans.setValue ( first+pos,false );
                count-= ( i_mem[first]==value );
                last++;
                first++;
            }
        }
    } else {
        while ( last<end ) {
            ok=!i_nans.getValue ( last );
            while ( last-first< ( winlength-1 ) ) {
                if ( ok ) count+= ( i_mem[last]==value );
                last++;
                ok=!i_nans.getValue ( last );
            }
            if ( ok ) count+= ( i_mem[last]==value );
            res->i_mem[first+pos]=count;
            res->i_nans.setValue ( first+pos,false );
            if ( !i_nans.getValue ( first ) )
                count-= ( i_mem[first]==value );
            last++;
            first++;
        }
    }
    return res;
}

template <class T> gArrayInternal<T> * gArrayInternal<T>::getWinMin ( gArrayInternal<T> * res, gSize winlength,gSize pos,gPos start,gPos end,bool skipnan ) const {
    gArrayInternal<T> * mymin = new  gArrayInternal<T> ( 1 );
    for ( gPos i= start; i<start+pos;i++ ) res->i_nans.setValue (i,true );
    for ( gPos i=start; i< end-winlength+1; i++ ) {
        getMin ( mymin, i,i+winlength,skipnan );
        res->setValue(i+pos, mymin->at ( 0 ), mymin->isNA ( 0 ),0,false );
    }
    delete mymin;
    for ( gPos i= end-winlength+pos+1; i<end;i++ ) res->i_nans.setValue (i,true );
    return res;
}

template <class T> gArrayInternal<T> * gArrayInternal<T>::getWinMax ( gArrayInternal<T> * res,  gSize winlength,gSize pos,gPos start,gPos end,bool skipnan ) const {
    gArrayInternal<T> * mymax = new gArrayInternal<T>( 1 );
    for ( gPos i= start; i<start+pos;i++ ) res->i_nans.setValue (i,true );
    for ( gPos i=start; i< ( end-winlength+1 );i++ ) {
        getMax ( mymax,i,i+winlength,skipnan );
        res->setValue ( i+pos, mymax->at ( 0 ), mymax->isNA ( 0 ),0,false );
    }
    delete mymax;
    for ( gPos i= end-winlength+pos+1; i<end;i++ ) res->i_nans.setValue ( i,true );
    return res;
}

template <class T> gArrayInternal<T> * gArrayInternal<T>::getWinSum ( gArrayInternal<T> * res,  gSize winlength,gSize pos,gPos start,gPos end,bool skipnan ) const {
    gArrayInternal<T> * mysum= new gArrayInternal<T> ( 1 );
    for ( gPos i= start; i<start+pos;i++ ) res->i_nans.setValue (i,true );
    for ( gPos i=start; i< ( end-winlength+1 );i++ ) {
        getSum (mysum, i,i+winlength,skipnan );
        res->setValue ( i+pos, mysum->at ( 0 ), mysum->isNA ( 0 ),0,false );
    }
    for ( gPos i= end-winlength+pos+1; i<end;i++ ) res->i_nans.setValue ( i,true );
    delete mysum;
    return res;
}

template <class T> gArrayInternal<T> * gArrayInternal<T>::getWinSquareSum ( gArrayInternal<T> * res,  gSize winlength,gSize pos,gPos start,gPos end,bool skipnan ) const {
    gArrayInternal<T> * mySsum = new gArrayInternal<T> ( 1 );
    for ( gPos i= start; i<start+pos;i++ ) res->i_nans.setValue (i,true );
    for ( gPos i=start; i< ( end-winlength+1 );i++ ) {
        getSquareSum (mySsum, i,i+winlength,skipnan );
        res->setValue ( i+pos, mySsum->at ( 0 ), mySsum->isNA ( 0 ),0,false );
    }
    delete mySsum;
    for ( gPos i=end-winlength+pos+1; i<end;i++ ) res->i_nans.setValue ( i,true );
    return res;
}

template <class T> gArrayInternal<gScore> * gArrayInternal<T>::getWinMean ( gArrayInternal<gScore> * res,  gSize winlength,gSize pos,gPos start,gPos end,bool skipnan ) const {
    bool ok;
    register gPos first=start,last=start;
    register gScore sum=0,wsize=0;
    for ( gPos i= start; i<start+pos;i++ ) res->i_nans.setValue (i,true );
    if ( !skipnan ) {
        while ( last<end ) {
            ok=!i_nans.getValue ( last );
            while ( ( last-first< ( winlength-1 ) ) ) {
                if ( !ok ) break;
                sum+=i_mem[last];
                wsize++;
                last++;
                ok=!i_nans.getValue ( last );
            }
            if ( !ok ) {
                res->i_nans.setValue ( first+pos,true );
                first=last+1;
                last=first;
                if ( winlength> ( end-first ) ) break;
                sum=0;
                wsize=0;
            } else {
                sum+=i_mem[last];
                wsize++;
                res->i_mem[first+pos]= ( gScore ) ( sum/wsize );
                res->i_nans.setValue ( first+pos,false );
                //resetBit(ret->i_val,first+pos);
                sum-=i_mem[first];
                wsize--;
                last++;
                first++;
            }
        }
    } else {
        while ( last<end ) {
            ok=!i_nans.getValue ( last );
            while ( last-first< ( winlength-1 ) ) {
                if ( ok ) {
                    sum+=i_mem[last];
                    wsize++;
                }
                last++;
                ok=!i_nans.getValue ( last );
            }
            if ( ok ) {
                sum+=i_mem[last];
                wsize++;
            }
            res->i_mem[first+pos]= ( gScore ) ( sum/wsize );
            res->i_nans.setValue ( first+pos,false );
            //resetBit(ret->i_val,first+pos);
            if ( !i_nans.getValue ( first ) ) {
                sum-=i_mem[first];
                wsize--;
            }
            last++;
            first++;
        }
    }
    for ( gPos i=end-winlength+pos+1; i<end;i++ ) res->i_nans.setValue ( i,true );
    return res;
}

template <class T> gArrayInternal<gScore> * gArrayInternal<T>::getWinStdDev ( gArrayInternal<gScore> * res, gSize winlength,gSize pos,gPos start,gPos end,bool skipnan,bool unbiased ) const {
    bool ok;
    register gPos first=start,last=start;
    register gScore sum=0,ssum=0,n= ( unbiased ) ? ( -1 ) : ( 0 );
    //gArrayInternal<gScore> *ret=new gArrayInternal<gScore> ( end-start );
    for ( gPos i= start; i<start+pos;i++ ) res->i_nans.setValue (i,true );
    if ( !skipnan ) {
        while ( last<end ) {
            ok=!i_nans.getValue ( last );
            while ( ( last-first< ( winlength-1 ) ) ) {
                if ( !ok ) break;
                sum+=i_mem[last];
                ssum+= ( i_mem[last]*i_mem[last] );
                n++;
                last++;
                ok=!i_nans.getValue ( last );
            }
            if ( !ok ) {
                first=last+1;
                last=first;
                if ( winlength> ( end-first ) ) break;
                sum=0;
                ssum=0;
                n= ( unbiased ) ? ( -1 ) : ( 0 );
            } else {
                sum+=i_mem[last];
                ssum+= ( i_mem[last]*i_mem[last] );
                n++;
                res->i_mem[first+pos]= ( gScore ) sqrt ( ssum/n- ( sum*sum ) / ( n* ( n+1 ) ) );
                res->i_nans.setValue ( first+pos,false );
                //resetBit(ret->i_val,first+pos);
                sum-=i_mem[first];
                ssum-= ( i_mem[first]*i_mem[first] );
                n--;
                last++;
                first++;
            }
        }
    } else {
        while ( last<end ) {
            ok=!i_nans.getValue ( last );
            while ( last-first< ( winlength-1 ) ) {
                if ( ok ) {
                    sum+=i_mem[last];
                    ssum+= ( i_mem[last]*i_mem[last] );
                    n++;
                }
                last++;
                ok=!i_nans.getValue ( last );
            }
            if ( ok ) {
                sum+=i_mem[last];
                ssum+= ( i_mem[last]*i_mem[last] );
                n++;
                res->i_mem[first+pos]= ( gScore ) sqrt ( ssum/n- ( sum*sum ) / ( n* ( n+1 ) ) );
                res->i_nans.setValue ( first+pos,false );
            }
            if ( !i_nans.getValue ( first ) ) {
                sum-=i_mem[first];
                ssum-= ( i_mem[first]*i_mem[first] );
                n--;
            }
            last++;
            first++;
        }
    }
    for ( gPos i=end-winlength+pos+1; i<end;i++ ) res->i_nans.setValue ( i,true );
    return res;
}
//--------------------------------------------------------------------------

}

extern "C" {
    /** @brief Function for library testing */
    void geco_running();
}

#endif

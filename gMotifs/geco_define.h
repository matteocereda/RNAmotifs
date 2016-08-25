/**
 * @file geco_define.h
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
 * This file contains library spcific typedefs and enum definitions. It also contains the
 * definition of the type used to allocate bit arrays.
 */
#ifndef __GECO_DEFINE__
#define __GECO_DEFINE__

/*!@mainpage GeCo++ library Documentation
 * @image html geco.png
 * @image latex geco.eps
 * @version 0.1
 * @author Matteo Cereda
 * @author Uberto Pozzoli
 * @section intro_sec Introduction
 *
 * GeCo++ (<span style="font-weight:bold;">Genomic Computation C++ Library</span>) is a library developed to 
 * manage genomic elements annotation, sequences and positional genomic features and to help users to keep 
 * The library has been developed starting from the idea to represent and manage the numeric results of computational 
 * algorithms keeping them tied to annotations of genomic elements (transcripts, binding sites, conserved regions,
 * transposable elements etc.), to their sequences and to genomic variations. Given this level of abstraction 
 * and the inherent complexity, we considered an object oriented software model to be the most appropriate. 
 * The choice of ISO C++ guaranteed  speed, portability and most importantly for users, the access to a great 
 * number of other efficient and specialized computational biology libraries.
 * GeCo++ defines two fundamental classes: gArray and gElement. The first one is a general purpose template 
 * array class. The development of such a class was justified by the need to provide both tracking of 
 * undefined/invalid elements (NA: Not Available) and array subsetting in a memory efficient way.  
 * NA tracking is obtained through a speed optimized bits array class (gBitsArray) while memory efficient 
 * subsetting by mean of an internal reference counting mechanism that allows instantiation of array subsets 
 * without data duplication. A number of optimized methods have been added to perform typical array operations 
 * needed in bioinformatics  such as sorting, values finding, counting, reversion, descriptive statistics and 
 * others. Windowed versions of the same methods are also provided. Furthermore a gArray object can be indexed 
 * through both a scalar value or by another gArray for multiple indexing. Type casting capabilities (between 
 * different template specializations) are also provided. Three additional classes have been derived from 
 * gArray: gMatrix still a template, with additional methods to treat matrices; gString that specializes gArray 
 * to manage character strings, and gSequence that inherits from gString adding methods for basic sequence management.<br/>
 * A genomic element is defined as an interval of a given reference sequence along a given strand. Reference 
 * positions are defined as absolute (unsigned) positions along a reference while element ones are relative 
 * (signed) to the beginning of an element along its strand. Sites are defined as particular positions along 
 * the element (for example transcription start sites, splice sites or protein binding sites) while a connection 
 * represent a directed relation between two sites (introns, exons). Positional feature (PFs) are defined as 
 * a property that varies along the element. While no assumption is made on the biological meaning of sites 
 * connection and features, this model is general enough to represent the majority of real-world genomic 
 * elements and their features. 
 * The GeCo++ library defines the class gElement as an implementation of this model: it allows users to 
 * instantiate objects representing genomic element which can contain sequence as well as sites, connections 
 * and features information.  Element positions can be converted to reference ones and vice-versa as well as 
 * one element positions can be mapped to another one. Sequence and features are maintained as gArray objects 
 * avoiding unnecessary data-duplication. Their retrieval or calculation has been kept independent from the 
 * gElement object itself by defining a hierarchy of gArrayRetriever objects from which users can easily derive 
 * new classes implementing sequence retrieval as well as feature calculation algorithms. This mechanism make 
 * easier to develop application that are independent from the specific computational algorithms.
 * The most important characteristic of gElement objects is that they can be instantiated as a sub-intervals 
 * of another one considering a strand and the presence of genomic variations. Sequence, sites, connections 
 * and features are inherited by the new object consistently with the interval, the strand and the variations. 
 * Features recalculation and sequence retrieval are kept to a minimum avoiding unnecessary recalculation at 
 * unaffected positions. This make very straightforward to evaluate the effect of genomic variations on 
 * features, positions and sequence.
 *
 * @section requirement_sec Requirements
 * The library is written in standard ISO C++ and should compile on most platforms given 
 * an ISO C++ compliant compiler is provided. It doesn't depend on other libraries except 
 * than on stl (vector,string and ostream). It comes with a cmake CMakeList file that 
 * requires cmake version 2.8 or higher.
 *
 * @section install_sec Library installation
 * To install the library dowload the package from: http://bioinformatics.emedea.it/geco/geco-0.1.tar.gz. 
 * ungzip it in a directory of your choice and follow the instruction contained in the INSTALL file
 * 
 * @section examples_sec Examples
 * At the address http://bioinformatics.emedea.it/geco you can find a complete tutorial 
 * with code examplestha tcovers many of the features of this library.
 *
 * @section licese_sec License
 * Copyright (C) 2010 by Uberto Pozzoli and Matteo Cereda
 * (uberto.pozzoli@bp.lnf.it)
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
 */

/**
 * @namespace geco
 * @brief The GeCo++ library namespace
 *
 * Namespace to identify all the elements of the library
 */

namespace geco {

//--------------------------------------------------------------------------
// Base type definitions
//--------------------------------------------------------------------------

/** @typedef unsigned char   _unsigned_char;
 * @brief typedef for unsigned char type 
 */
typedef unsigned char   _unsigned_char;

/** @typedef  char            _char;
 *  @brief typedef for char type */
typedef char            _char;

/** @typedef unsigned short  _unsigned_short;
 *  @brief typedef for unsigned short type */
typedef unsigned short  _unsigned_short;

/** @typedef short           _short;
 *  @brief typedef for short type */
typedef short           _short;

/** @typedef  unsigned int    _unsigned_int;
 *  @brief typedef for unsigned int type */
typedef unsigned int    _unsigned_int;

/** @typedef int             _int;
 *  @brief typedef for char type */
typedef int             _int;

/** @typedef unsigned long   _unsigned_long;
 *  @brief typedef for unsigned long type */
typedef unsigned long   _unsigned_long;

/** @typedef long            _long;
 *  @brief typedef for long type */
typedef long            _long;

/** @typedef double          _double;
 *  @brief typedef for double type */
typedef double          _double;

/** @typedef long double     _long_double;
 *  @brief typedef for long double type */
typedef long double     _long_double;
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
// Specialized type definitions
//--------------------------------------------------------------------------
/** @typedef _char gChar;
 *  @brief Typedef used to represent characters */
typedef _char gChar;

/** @typedef _unsigned_long gSize;
 *  @brief Typedef used to represent array size */
typedef _unsigned_long gSize;

/** @typedef _unsigned_int gShortUnsigned;
 *  @brief Typedef used to represent short unsigend integers */
typedef _unsigned_int gShortUnsigned;

/** @typedef _unsigned_short gArrayIndex;
 *  @brief Typedef used for (short) array indexing */
typedef _unsigned_short gArrayIndex;

/** @typedef _unsigned_long gPos;
 *  @brief Typedef used to represent absolute positions into (long) arryas */
typedef _unsigned_long gPos;

/** @typedef _long gRelativePos;
 *  @brief Typedef used to represent relative positions into arrays */
typedef _long gRelativePos;

/** @typedef _unsigned_short gBool;
 *  @brief Typedef used to represent boolean values */
typedef _unsigned_short gBool;

/** @typedef _double gScore;
 *  @brief Typedef used to represent floating point double precision score values */
typedef _double gScore;

/** @typedef _long_double gDoubleScore;
 *  @brief Typedef used to represent floating point quadruple precision scores */
typedef _long_double gDoubleScore;
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
// Specialized enums definition
//--------------------------------------------------------------------------
/** @brief Genomic variation types
 *
 * Values from this enum indicate at which level (element, region or reference) a gElement object holds the sequence
 * @li gSubstitution: a certain number of BPs that are substituted with an equal number of new ones. These latter do nothhave to be all different from the originals matching BPs
 * @li gInsertion: a situatiion wehreby a certain number of BPs is inserted at a given position
 * @li gDeletion: a situation whereby a number of adjacent BPs are removed
 * @ingroup elements
 */
typedef enum {gSubstitution,gInsertion,gDeletion} gVariationType;

/** @brief Elements sequence allocation modes
 *
 * Sequence allocation mode indicates at which level (element, region or reference) a gElement object holds the sequence
 * @li gRef: sequence is maintained at the reference level.
 * This means that sub elements possibly with variations realtive to the original can be instantiated sharing the same sequence data.
 * In this case actual sequence (strand, variations) must be calculated every time it's needed
 * @li gReg: sequence is maintained at the region level.
 * This means that sub elements without variations relative to the original can be instantiated sharing the same sequence.
 * In this case actual sequence strand must be calculated every time it's needed.
 * @li gElm: sequence is maintained at the element level.
 * This means that sub elements with the same strand and without variations relative to the original can be instantiated sharing the same sequence.
 * In this case actual sequence needs no calculation.
 * @ingroup elements
 */
typedef enum {gRet,gRef,gReg,gElm} gElementSequenceMode;

/** @brief Transcript site types
*
* Values from this enum indicate at which level (element, region or reference) a gElement object holds the sequence
* @li tss: Transcription start site
* @li donor: Splice-donor site (5' ss)
* @li t5: Truncated 5' site
* @li acceptor: Splice-acceptor site (3')
* @li t3 Truncated 3'
* @li pas: Polyadenilation site
* @ingroup elements
*/
typedef enum {tss,donor,tdonor,acceptor,tacceptor,pas} gTranscriptSiteType;

//--------------------------------------------------------------------------


//--------------------------------------------------------------------------
// bitsarray type definition & constants
//--------------------------------------------------------------------------
//BitsAarray types
#define __GECO_GBYTES_UCHAR__
#undef __GECO_GBYTES_USHORT__
#undef __GECO_GBYTES_UINT__
#undef __GECO_GBYTES_ULONG__

#ifdef __GECO_GBYTES_UCHAR__
typedef _unsigned_char gBytes;
#endif

#ifdef __GECO_GBYTES_USHORT__
typedef _unsigned_short gBytes;
#endif

#ifdef __GECO_GBYTES_UINT__
typedef _unsigned_int gBytes;
#endif

#ifdef __GECO_GBYTES_ULONG__
typedef _unsigned_long gBytes;
#endif

extern const _unsigned_char mLogNBits[64];

extern const _unsigned_char mNBytes;

extern const _unsigned_char mNBits;

extern const _unsigned_char mP1;

extern const _unsigned_char mP2;

extern const gBytes mP3;

#ifdef __GECO_GBYTES_UCHAR__
extern const gBytes MG[8];
extern const gBytes MR[8];
#endif

#ifdef __GECO_GBYTES_USHORT__
extern const gBytes MG[16];
extern const gBytes MR[16];
#endif

#ifdef __GECO_GBYTES_UINT__
extern const gBytes MG[32];
extern const gBytes MR[32];
#endif

#ifdef __GECO_GBYTES_ULONG__
extern const gBytes MG[64];
extern const gBytes MR[64];
#endif

extern  const _unsigned_char mBitsSetTable[256];
//--------------------------------------------------------------------------

}
//end of namespace

extern "C"{
  void geco_running();
}

#endif

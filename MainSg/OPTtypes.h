/*--------------------------------------------------------------------------*/
/*--------------------------- File OPTtypes.h ------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Standard types and constants definitions: by changing the definitions
 * here, the dimension/precision of all the numbers used within the
 * programs can be customized.
 *
 * It also includes OPTUtils.h for various stuff related to reading the
 * time of a code, generating random numbers and safely reading parameters
 * out of a stream; see the comments in OPTUtils.h.
 *
 * \version 5.00 beta
 *
 * \date 08 - 05 - 2012
 *
 * \author Antonio Frangioni \n
 *         Operations Research Group \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * Copyright &copy 1994 - 2012 by Antonio Frangioni
 */

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __OPTtypes
 #define __OPTtypes  /* self-identification: #endif at the end of the file */

/*@} -----------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "OPTUtils.h"

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

#include <limits.h>

/*--------------------------------------------------------------------------*/
/*--------------------------- NAMESPACE ------------------------------------*/
/*--------------------------------------------------------------------------*/

#if( OPT_USE_NAMESPACES )
namespace OPTtypes_di_unipi_it
{
 /** @namespace OPTtypes_di_unipi_it
     The namespace OPTtypes_di_unipi_it is defined to hold all the data
     types, constants, classes and functions defined here. It also
     comprises the namespace std. */
#endif

/*@} -----------------------------------------------------------------------*/
/*---------------------------- TYPES ---------------------------------------*/
/*--------------------------------------------------------------------------*/
/** @defgroup OPTTYPES_GENPUR General-purpose type definitions
    @{ */

typedef bool           *Bool_Vec;        ///< vector of booleans
typedef Bool_Vec       *Bool_Mat;        ///< matrix of booleans

typedef const bool     cBOOL;            ///< a read-only boolean
typedef cBOOL         *cBool_Vec;        ///< read-only array

/*--------------------------------------------------------------------------*/

typedef unsigned int    Index;           ///< Index in a vector ( >= 0 )
typedef Index          *Index_Set;       ///< set (array) of indices
typedef Index_Set      *Index_Mat;       ///< set of set of indices

typedef const Index    cIndex;           ///< a read-only Index
typedef cIndex        *cIndex_Set;       ///< read-only array

/*--------------------------------------------------------------------------*/

typedef int             SIndex;         ///< A Signed Index, for when an index
                                        ///< may have a "direction"
typedef SIndex         *SIndex_Set;     ///< set (array) of s. indices
typedef SIndex_Set     *SIndex_Mat;     ///< set of set (matrix) of s. indices

typedef const SIndex   cSIndex;         ///< a read-only Signed Index
typedef cSIndex       *cSIndex_Set;     ///< read-only array

/*--------------------------------------------------------------------------*/

typedef long int        INum;           ///< integer numbers
typedef INum           *IRow;           ///< vectors of integers
typedef IRow           *IMat;           ///< matrices (vectors of vectors)
                                        ///< of integers

typedef const INum     cINum;           ///< a read-only integer
typedef cINum         *cIRow;           ///< read-only array

/*--------------------------------------------------------------------------*/

typedef double          Number;         ///< "normal" floating point numbers
typedef Number         *Row;            ///< "normal" (fp) array
typedef Row            *Mat;            ///< "normal" (fp) matrix

typedef const Number   cNumber;         ///< a read-only Number
typedef cNumber       *cRow;            ///< read-only array

/*--------------------------------------------------------------------------*/

typedef double          HpNum;          ///< "finer" floating point numbers
typedef HpNum          *HpRow;          ///< "finer" (fp) array
typedef HpRow          *HpMat;          ///< "finer" (fp) matrix

typedef const HpNum    cHpNum;          ///< a read-only HpNum
typedef cHpNum        *cHpRow;          ///< read-only array

/*@} -----------------------------------------------------------------------*/
/*----------- Type definitions for subgradient-based algorithms ------------*/
/*--------------------------------------------------------------------------*/
/** @defgroup OPTTYPES_SUBGT Type definitions for subgradient-based algorithms
    @{ */

typedef double          SgNum;          ///< subgradient entries
typedef SgNum          *SgRow;          ///< a subgradient
typedef SgRow          *SgMat;          ///< a bundle (set of subgradients)

typedef const SgNum    cSgNum;          ///< a read-only subgradient entry
typedef cSgNum        *cSgRow;          ///< a read-only subgradient

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

typedef double          QuNum;          ///< numbers in Q ( G{i}{T} * G{j} )
typedef QuNum          *QuRow;          ///< row of Q
typedef QuRow          *QuMat;          ///< Q (itself)

typedef const QuNum    cQuNum;          ///< a read-only number in Q
typedef cQuNum        *cQuRow;          ///< a read-only row of Q

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

typedef double          LMNum;          ///< a Lagrangean Multiplier
typedef LMNum          *LMRow;          ///< a vector of Lagrangean Multipliers
typedef LMRow          *LMMat;          ///< a matrix of Lagrangean Multipliers

typedef const LMNum    cLMNum;          ///< a read-only Lagrangean Multiplier
typedef cLMNum        *cLMRow;          ///< a read-only vector of LMs

/*@} -----------------------------------------------------------------------*/

/* @} end( group( OPTTYPES_TYPES ) ) */

/*--------------------------------------------------------------------------*/

#if( OPT_USE_NAMESPACES )
 };  // end( namespace OPTtypes_di_unipi_it )
#endif

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif  /* OPTtypes.h included */

/*--------------------------------------------------------------------------*/
/*---------------------- End File OPTtypes.h -------------------------------*/
/*--------------------------------------------------------------------------*/

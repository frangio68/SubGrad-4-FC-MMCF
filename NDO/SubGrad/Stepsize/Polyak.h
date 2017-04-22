/*--------------------------------------------------------------------------*/
/*--------------------------- File Polyak.h ------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * The class Polyak implements the easiest target value stepsize rule, with a
 * constant level and the factor beta. The class conforms to the interface
 * defined by the class Deflection [see Deflection.h]. For more details we
 * refer to:
 *
 *  B.T. Polyak. <em>Minimization of unsmooth functionals. Zh. Vychisl.</em>
 *  Mat. Fiz. 9(3):509--521, (1969)
 *
 * \date 27 - 10 - 2015
 *
 * \author Antonio Frangioni \n
 *         Operations Research Group \n
 *         Department of Informatics  \n
 *         University of Pisa \n \n
 *
 * \author Enrico Gorgone \n
 *         Graphs and Mathematical Optimization (GOM) Group \n
 *         Department of Informatics \n
 *         Free University of Brussels \n \n
 *
 * Copyright &copy 2001 - 2015 by Antonio Frangioni, Enrico Gorgone.
 */

/*--------------------------------------------------------------------------*/
/*-------------------- DEFINITIONS & IMPLEMENTATION ------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __Polyak
 #define __Polyak  /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "Stepsize.h"

/*--------------------------------------------------------------------------*/
/*------------------------------- MACROS  ----------------------------------*/
/*--------------------------------------------------------------------------*/

#define LOG_STP 0

/* If LOG_STP > 0, the Polyak class produces a log of its activities on the
   ostream object and at the "level of verbosity" set with the method
   SetSTPLog() [see below]. */

#if( LOG_STP )
 #define STPLOG( l , x ) if( STPLLvl > l ) *STPLog << x
 #define STPLOG2( l , c , x ) if( ( STPLLvl > l ) && c ) *STPLog << x
#else
 #define STPLOG( l , x )
 #define STPLOG2( l , c , x )
#endif

/*--------------------------------------------------------------------------*/
/*----------------------------- NAMESPACE ----------------------------------*/
/*--------------------------------------------------------------------------*/

#if( OPT_USE_NAMESPACES )
namespace NDO_di_unipi_it
{
#endif

/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/** Definition of the class Polyak. This class implements  a target value
    stepsize rule, whereby \f$\beta_i\f$ and \f$f^{lev}_i\f$ are constant
    (they do not depend on $i$).                                          */

class Polyak : public Stepsize
{
/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*--------------------------------------------------------------------------*/
/*--------------------- PUBLIC METHODS OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/
/*---------------------------- CONSTRUCTOR ---------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructor
   @{ */

   inline Polyak( SubGrad *slvr, istream *iStrm )
          :
          Stepsize( slvr , iStrm )

 /**<  Constructor of the class. Since the constructor of Polyak is executed
    after the one of Stepsize, the following parameters specific for the
    Polyak have to be found in the stream after those of the base class [see
    the comments in the constructor of Stepsize].

       -# HpNum Beta [ 1e-2 ] beta value. */
   {
    DfltdSfInpt( iStrm , Beta , HpNum( 1e-2 ) );
    }

/*@} -----------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
   @{ */

   inline void SetSTPLog( ostream *outs = 0 , const char lvl = 0 );

/*--------------------------------------------------------------------------*/

   void Format( void ) {
    MaxBeta = Inf<HpNum>();
    FiLev = - Inf<HpNum>();
    }

 /*@}------------------------------------------------------------------------*/
 /*-------------------- METHODS FOR STEPSIZE COMPUTATION --------------------*/
 /*--------------------------------------------------------------------------*/
 /** @name Computing the stepsize
    @{ */

   inline void NewStep( void );

/*@} -----------------------------------------------------------------------*/
/*---------------------- PROTECTED PART OF THE CLASS -----------------------*/
/*--------------------------------------------------------------------------*/

 protected:

    inline bool UpdateTargetLevel( void );

/*@}------------------------------------------------------------------------*/
/*--------------------- PRIVATE PART OF THE CLASS --------------------------*/
/*--------------------------------------------------------------------------*/

 private:

/*--------------------------------------------------------------------------*/
/*----------------------- PRIVATE DATA STRUCTURES  -------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

 };  // end( class Polyak )

/*--------------------------------------------------------------------------*/
/*------------------- inline methods implementation ------------------------*/
/*--------------------------------------------------------------------------*/

inline void Polyak::SetSTPLog( ostream *outs , const char lvl ) {

 Stepsize::SetSTPLog( outs , lvl );

 #if( LOG_STP )
  if( STPLLvl > 1 ) {
   *STPLog << endl << "Polyak: Beta = "<< Beta << " ~ Chi = " <<  LpsFct
	   << endl;
   }
 #endif
 }

/*--------------------------------------------------------------------------*/

inline void Polyak::NewStep( void )
{
 // update the target value - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 UpdateTargetLevel();

 // some exception  - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( FiLev == - Inf<HpNum>() )
  throw( NDOException( "Polyak::GetStepsize: no lower bound is given" ) );

 } // end ( Polyak::NewStep )

/*--------------------------------------------------------------------------*/

inline bool Polyak::UpdateTargetLevel( void )
{
 cHpNum LwrBnd = GetOracle()->GetLowerBound();
 if( LwrBnd > FiLev ) {
  FiLev = LwrBnd;
  return( true );
  }

 return( false );
 }

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#if( OPT_USE_NAMESPACES )
 };  // end( namespace NDO_di_unipi_it )
#endif

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif  /* Polyak.h included */

/*--------------------------------------------------------------------------*/
/*------------------------- End File Polyak.h ------------------------------*/
/*--------------------------------------------------------------------------*/

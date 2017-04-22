/*--------------------------------------------------------------------------*/
/*--------------------------- File FumeroTV.h ------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Definition and implementation of the class FumeroTV, which is a stepsizes
 *  rule (SR) to be used within the subgradient method (SM) [see SubGrad.h].
 * The class conforms to the interface defined by the class Stepsize
 * [see Stepsize.h]. For more details we refer to:
 *
 *  F. Fumero. <em>A modified subgradient algorithm for Lagrangean
 *  relaxation</em> Comput. Oper. Res. 28(1):33--52, (2001)
 *
 * \version 1.00
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

#ifndef __FumeroTV
 #define __FumeroTV  /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "Stepsize.h"

/*--------------------------------------------------------------------------*/
/*------------------------------- MACROS  ----------------------------------*/
/*--------------------------------------------------------------------------*/

#define LOG_STP 0

/* If LOG_STP > 0, the FumeroTV class produces a log of its activities on the
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
/** Definition of the class FumeroTV. This class implements a target value
    stepsize. At the beginning \f$ f^{lev} \f$ is set to the given lower bound
    and exponentially goes towards \f$ f^{rec} \f$. The method works in two
    phases. In the first one \f$ \beta_i \f$ is decremented, while in the
    second one \f$ \beta_i \f$ is still decremented but also incremented
    after a pre-set number of consecutive improving iterations. */

class FumeroTV : public Stepsize
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

   inline FumeroTV( SubGrad *slvr , istream *iStrm = NULL );

/**< Constructor of the class. Since the constructor of FumeroTV is executed
   after the one of Stepsize, the following parameters specific for the
   FumeroTV have to be found in the stream <em>after</em> those of the base
   class [see the comments to the constructor of Stepsize]:

     -# HpNum sigmaMin  [ 1e-4 ] tolerance of the function \f$ sigma_r \f$

     -# Index r1        [ 2 ]    sigma_{r1} = 1/2

     -# HpNum beta0     [ 1 ]    initial value of beta

     -# Index eta1      [ 10 ]   threshold on number of failures (phase I)

     -# index eta2      [ 10 ]   threshold on number of failures (phase II)
     */

/*@} -----------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
   @{ */

   inline void SetSTPLog( ostream *outs = 0 , const char lvl = 0 );

/*--------------------------------------------------------------------------*/

   inline void Format( void );

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
/*-------------------------- PRIVATE METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/

   inline HpNum ExpFun( Index ri );

/*--------------------------------------------------------------------------*/
/*----------------------- PRIVATE DATA STRUCTURES  -------------------------*/
/*--------------------------------------------------------------------------*/

   HpNum FiLambda;   ///< FiLambda

   HpNum LwrBnd;     ///< lower bound
   HpNum FiRef;      ///< the best solution found so far

   HpNum sigma;      ///< current value of \f$ \sigma(r)\ f$
   HpNum sigmaMin;   ///< minimum threshold of \f$ \sigma(r) \f$

   HpNum BetaZero;

   Index rc;         ///< current value of r
   Index r1;         ///<
   Index r2;

   Index etac;       ///< last improvement of FiLev
   Index eta1;       ///< thresholds for max number of failures
   Index eta2;

 }; // end( class FumeroTV )

/*--------------------------------------------------------------------------*/
/*------------------- inline methods implementation ------------------------*/
/*--------------------------------------------------------------------------*/

inline FumeroTV::FumeroTV( SubGrad *slvr , istream *iStrm )
                 :
                 Stepsize( slvr , iStrm )
{
 DfltdSfInpt( iStrm , sigmaMin , HpNum( 1e-4 ) );
 DfltdSfInpt( iStrm , r1 , Index( 5 ) );
 DfltdSfInpt( iStrm , BetaZero , HpNum( 1 ) );
 DfltdSfInpt( iStrm , eta1 , Index( 5 ) );
 DfltdSfInpt( iStrm , eta2 , Index( 5 ) );

 }  // end( FumeroTV::FumeroTV )

/*--------------------------------------------------------------------------*/

inline void FumeroTV::Format( void )
{
 Stepsize::Format( );
 Beta = BetaZero;

 FiLambda = FiRef = Inf<HpNum>();
 // FiLambda and the reference value are not known
 LwrBnd = -Inf<HpNum>();  // no lower bound or target level is known

 etac = rc = 0;
 r2 = ceil( pow( ( - log( sigmaMin ) / 0.6933 ) , ( 1.0 / 3.26 ) ) * r1 );

 sigma = 1;

 }  // end ( FumeroTV::Format() )

/*--------------------------------------------------------------------------*/

inline void FumeroTV::NewStep( void )
{
 // get the function value   - - - - - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 FiLambda = ReadFkVal();

 if( FiLambda == Inf<HpNum>() )
  throw( NDOException( "FumeroTV::GetStepsize: this should not happen" ) );

 // try to improve the target level - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 bool TLchgd = UpdateTargetLevel();

 // update target level- - - - - - - - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( TLchgd )
  FiLev = sigma * LwrBnd + ( 1 - sigma ) * FiRef;

 // fix beta - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( Beta < 1e-6 )
  Beta = 1e-6;

 // print the level - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 STPLOG( 1 , endl << "           " << " level  = " << -FiLev << " ~  beta = "
	     << Beta << endl << "           " );

 }  // end( FumeroTV::NewStep() )

/*--------------------------------------------------------------------------*/

inline void FumeroTV::SetSTPLog( ostream *outs , const char lvl )
{
 Stepsize::SetSTPLog( outs , lvl );

 #if( LOG_STP )
  if( STPLLvl > 1 ) {
   *STPLog << endl << "FumeroTV: ~ BetaZero = " << BetaZero
	   << " ~ (Eps0 , r1) = (" << sigma << " , " << r1 << ")"
	   << " ~  Eta = {" << eta1 << " , " << eta2 << "}" << endl;
   }
 #endif
 }  // end( FumeroTV::SetSTPLog() )

/*--------------------------------------------------------------------------*/

inline bool FumeroTV::UpdateTargetLevel( void )
{
 bool TLchgd = false;

 // the target level is updated if any changes of the lower bound occurs - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 cHpNum LwrBnd_ = GetOracle()->GetLowerBound();
 if( LwrBnd_ > LwrBnd ) {
  LwrBnd = LwrBnd_;
  TLchgd = true;
  }

 // some exception - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( LwrBnd == - Inf<HpNum>() )
  throw( NDOException(
	           "FumeroTV::UpdateTargetLevel: no lower bound is given" ) );
 HpNum FiRef_;
 if( FiRef == Inf<HpNum>() )  // not initialized yet
  FiRef_ = Inf<HpNum>();
 else
  FiRef_ = FiRef - ( 1e-6 * max( 1.0 , ABS( FiRef ) ) );

 // beta updating- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( FiLambda <= FiRef_ ) {
  FiRef = FiLambda;
  if( rc >= r2 )
   Beta *= 2.0;
  etac = 0;
  TLchgd = true;
  }
 else  // FiLambda > FiRef_
  if( rc >= r2 ) {
   if( ++etac >= eta2 ) {
    Beta /= 2;
    etac = 0;
    }
   }
  else
   if( ++etac >= eta1 ) {
    Beta = Beta / (2.0 * Beta + 1.0);
    sigma = ExpFun( ++rc );
    TLchgd = true;
    etac = 0;
    }

 // failure of the target level updating - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( FiLev < LwrBnd ) {  // change the level
  FiLev = LwrBnd;
  TLchgd = true;
  }

 return( TLchgd );

 }  // end( FumeroTV::UpdateTargetLevel() )

/*--------------------------------------------------------------------------*/

inline HpNum FumeroTV::ExpFun( Index ri )
{
 return( exp( -0.6933 * pow( HpNum(ri) / HpNum( r1 ), 3.26 ) ) );

 }  // end( FumeroTV::ExpFun() )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#if( OPT_USE_NAMESPACES )
 };  // end( namespace NDO_di_unipi_it )
#endif

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif  /* FumeroTV.h included */

/*--------------------------------------------------------------------------*/
/*------------------------- End File FumeroTV.h ----------------------------*/
/*--------------------------------------------------------------------------*/

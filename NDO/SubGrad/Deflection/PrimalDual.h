/*--------------------------------------------------------------------------*/
/*----------------------- File PrimalDual.h --------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * The class PrimalDual implements two variants of the Primal-Dual subgradient
 * method (PDSM): the simple and weighted averages. These boils down to the
 * <em>simultaneous</em> choice of a deflection coefficient and a stepsize in
 * a SM. For more details we refer to:
 *
 * Y. Nesterov <em>Primal-dual subgradient methods for convex
 * optimization</em>. Math. Prog. 120, 221-259, 2009
 *
 * Because this class choses <em>both</rm> the deflection coefficient and the
 * stepsize, it conforms to <em>both</em> the interfaces defined by the class
 * Deflection [see Deflection.h] and Stepsize [see Stepsize.h].
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

#ifndef __PrimalDual
 #define __PrimalDual /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------- MACROS  ----------------------------------*/
/*--------------------------------------------------------------------------*/

#define LOG_VOL 0

/* If LOG_VOL > 0, the PrimalDual class produces a log of its activities on the
   ostream object and at the "level of verbosity" set with the method
   SetVOLLog() [see below]. */

#if( LOG_VOL )
 #define VOLLOG( l , x ) if( VOLLLvl > l ) *VOLLog << x
 #define VOLLOG2( l , c , x ) if( ( VOLLLvl > l ) && c ) *VOLLog << x
#else
 #define VOLLOG( l , x )
 #define VOLLOG2( l , c , x )
#endif

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "Deflection.h"
#include "Stepsize.h"

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
/** Definition of the class PrimalDual. This class implements two variants of
    the Primal-Dual subgradient method (PDSM): simple and weighted averages,
    under a unified SM scheme. The method has a fixed stability center, and
    the direction \f$ d_i \f$ can be regarded as the average of all the
    previous weighted subgradients, i.e.:
    \f[
     d_i =  ( \sum_{k=1}^i \upsilon_k g_k) / \Delta_i \quad,\quad
     \Delta_i = \sum_{k=1}^i \upsilon_k
    \f]
    It is worth mentioning that the subgradient weights \f$\upsilon_k \f$
    can be non-vanishing as \f$k\rightarrow \infty\f$. The stepsize rule (SR)
    and the deflection rule (DR), respectively, are:
    \f[
    \alpha_i = \upsilon_i / \Delta_i \quad,\quad
    \nu_i = \Delta_i / \omega_i
    \f]
    where \f$\omega_i\f$ is chosen in a suitable way. */

class PrimalDual : public Deflection , public Stepsize
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

  inline PrimalDual( SubGrad *slvr , istream *iStrm = 0 );

/**< Constructor of the class, which derives from both Deflection and
   Stepsize. Since the constructor of PrimalDual is executed <em>after</em>
   the one of Stepsize, the following parameters specific for PrimalDual have
   to be found in the stream \a after those of the base class Stepsize [see the
   comments to the constructor of Stepsize]:

     -# bool average [true] if true simple averages are used, otherwise
                            weighted averages are used
 
   For this class, the general parameter LpsFct -defined in Stepsize.h- represents
   the gamma factor \a F used to set \f$\omega_i\f$. */

 /*@} -----------------------------------------------------------------------*/
 /*-------------------------- OTHER INITIALIZATIONS -------------------------*/
 /*--------------------------------------------------------------------------*/
 /** @name Other initializations
   @{ */

   inline void SetVOLLog( ostream *outs = 0 , const char lvl = 0 );

   inline void Format( void );

/*@} -----------------------------------------------------------------------*/
/*-----------------------  METHODS FOR DR and SR COMPUTATION  --------------*/
/*--------------------------------------------------------------------------*/
/** @name Reading the solution
     @{ */

/** The variant requires thats DR must be computed before SR. Hence, the call
    to NewDEF comes before the one to NewStep. */

   inline void NewDEF( void );

   inline void NewStep( void );

/*--------------------------------------------------------------------------*/

   inline cHpNum GetDFLCoeff( void );

   inline HpNum GetStepsize( bool StepIsIncr = false );

/*--------------------------------------------------------------------------*/

  HpNum GetLev( void )
  {
   throw( NDOException( "PrimalDual::GetLev: this call is not allowed" ) );
   }

/*--------------------------------------------------------------------------*/

  HpNum GetBeta( void )
  {
   throw( NDOException( "PrimalDual::GetBeta: this call is not allowed" ) );
   }

/*--------------------------------------------------------------------------*/

  void SetMaxBeta( const HpNum alpha )
  {
   throw( NDOException( "PrimalDual::SetMaxBeta: this call is not allowed" ) );
   }

/*@} -----------------------------------------------------------------------*/
/*--------------------- PRIVATE PART OF THE CLASS --------------------------*/
/*--------------------------------------------------------------------------*/

 private:

/*--------------------------------------------------------------------------*/
/*--------------------------- PRIVATE METHODS ------------------------------*/
/*--------------------------------------------------------------------------*/

   inline void UpdateGamma( void );

/*--------------------------------------------------------------------------*/
/*----------------------- PRIVATE DATA STRUCTURES  -------------------------*/
/*--------------------------------------------------------------------------*/

   HpNum delta;
   HpNum Tmpdelta;
   HpNum vk;

   HpNum StepSize;
   HpNum gamma;

   HpNum beta1;
   bool average;

   HpNum LipCnst;      // Lipschitz constant

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

 };  // end( class PrimalDual )

/*--------------------------------------------------------------------------*/
/*------------------- inline methods implementation ------------------------*/
/*--------------------------------------------------------------------------*/

inline PrimalDual::PrimalDual( SubGrad *slvr , istream *iStrm )
                   :
                   Deflection( slvr ) , Stepsize( slvr , iStrm )
 {
  DfltdSfInpt( iStrm , average , bool( true ) );

  }  // end( PrimalDual::PrimalDual() )

/*--------------------------------------------------------------------------*/

inline void PrimalDual::SetVOLLog( ostream *outs , const char lvl )
{
 Deflection::SetVOLLog( outs , lvl );

 #if( LOG_VOL )
  if( VOLLLvl > 1 )
   if( average )
    *VOLLog << endl << "PrimalDual: Simple averages ";
   else
    *VOLLog << endl << "PrimalDual: Weighted averages ";

  *VOLLog << "~ GammaFact = " <<  LpsFct << endl;
 #endif

} // end( PrimalDual::SetVOLLog() )

/*--------------------------------------------------------------------------*/

inline void PrimalDual::Format( void )
{
 Deflection::Format();

 HpNum LipTmp = Deflection::GetOracle()->GetGlobalLipschitz();
 if( LipTmp != LipCnst )
  LipCnst = LipTmp;

 beta1 = 1;
 delta = 0;

 }  // end( PrimalDual::Format() )

/*--------------------------------------------------------------------------*/

inline void PrimalDual::NewStep( void )
{
 // update the Lipschitz constant - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 HpNum LipTmp = Deflection::GetOracle()->GetGlobalLipschitz();
 if( LipTmp != LipCnst )
  LipCnst = LipTmp;

 UpdateGamma( );

 // SR computation - - -  - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 StepSize = Tmpdelta / ( gamma * beta1 );

 VOLLOG( 1, endl << "            Delta = " << Tmpdelta << " ~ Beta = "
	    << beta1 << " ~ gamma = " << gamma );

 delta = Tmpdelta;
 beta1 = ( beta1 + 1/beta1 );

 }  // end( PrimalDual::NewStep() )

/*--------------------------------------------------------------------------*/

inline HpNum PrimalDual::GetStepsize( bool StepIsIncr )
{
 if( StepIsIncr )
  throw( NDOException( "PrimalDual::GetStepsize: this is not allowed" ) );

 return( StepSize );

 }  // end( PrimalDual::GetStepsize() )

/*--------------------------------------------------------------------------*/

inline void PrimalDual::NewDEF( void )
{
 if( average )
  vk = 1;
 else
  vk = 1 / sqrt( Deflection::GetGiNorm() );

 Tmpdelta = delta + vk;

 }  // end( PrimalDual::NewDEF() )

/*--------------------------------------------------------------------------*/

inline cHpNum PrimalDual::GetDFLCoeff( void )
{
 return( vk / Tmpdelta );
 }

/*--------------------------------------------------------------------------*/

inline void PrimalDual::UpdateGamma( void )
{
 HpNum tStar;
 Deflection::Solver->GetPar( NDOSolver::ktStar , tStar );

 if( average )
  gamma = LpsFct * ( sqrt( LipCnst / 2.0 )  / tStar );
 else
  gamma = LpsFct / ( sqrt( 2.0 * LipCnst ) * tStar );

 }  // end( PrimalDual::UpdateGamma() )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#if( OPT_USE_NAMESPACES )
 };  // end( namespace NDO_di_unipi_it )
#endif

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif  /* PrimalDual.h included */

/*--------------------------------------------------------------------------*/
/*-------------------- End File PrimalDual.h -------------------------------*/
/*--------------------------------------------------------------------------*/

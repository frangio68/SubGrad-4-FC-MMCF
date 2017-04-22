/*--------------------------------------------------------------------------*/
/*--------------------------- File Volume.h --------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * The class Volume implements the deflection rule of the modified Volume
 * algorithm, i.e., that of the "poorman's" Bundle method. For more details
 * we refer to:
 *
 * L. Bahiense, N. Maculan, C. Sagastizabal. <em>The volume algorithm
 * revisited: relation with bundle methods. </em> Math. Prog. 94, 41-69, 2002
 *
 * The class conforms to the interface defined by the class Deflection [see
 * Deflection.h].
 *
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
 * Copyright &copy 2001 - 2015 by Antonio Frangioni, Enrico Gorgone
 */

/*--------------------------------------------------------------------------*/
/*-------------------- DEFINITIONS & IMPLEMENTATION ------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __Volume
 #define __Volume  /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------- MACROS  ----------------------------------*/
/*--------------------------------------------------------------------------*/

#define LOG_VOL 0

/* If LOG_VOL > 0, the Volume class produces a log of its activities on the
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
#include "OPTvect.h"

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
/** Definition of the class Volume. This class implements the Volume algorithm.
    The method is revisited and incorporated in the SubGrad solver. The
    deflection coefficient is found solving a quadratic problem that
    involves the subgradient \f$ g_i \f$, the direction \f$ d_{i-1} \f$, and
    their respective linearization errors \f$ \sigma_i \f$  and
    \f$ \epsilon_{i-1} \f$ at the stability center \f$ \bar{\lambda}_i\f$,
    namely
    \f[
     \tau_i = 
     \arg\min \left\{ \; \nu_{i-1} \left\| \tau g_i 
                         + (1 - \tau) d_{i-1} \right\|^2/2
                         +  \tau \sigma_i(\bar{\lambda}_i)
                         + (1 - \tau)\epsilon_{i-1}(\bar{\lambda}_i)
        \;:\; \tau \in [0, 1]\; \right\}
    \f]  */

class Volume : public Deflection
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

  Volume( SubGrad *slvr , istream *iStrm = 0 ) : Deflection( slvr )

/**< Constructor of the class. The parameter `iStrm', if provided, is taken as
   a pointer to a istream from which the algorithmic parameters for the
   Volume are sequentially read in the following order. Each parameter must be
   placed at the beginning of a separate line, max 255 characters long, with
   all the rest of the line up to the first newline character '\n' (apart from
   a separating whitespace) being available for comments. Any line whose first
   character is '#' and any blank line is ignored. If 0 is passed, the file
   ends before reaching a given parameter, or some parameter is in the wrong
   format, each non-specified parameter is given a default value, shown in []
   below.

    -#   HpNum tauInit   [0.05]  \f$ \tau_0 \f$

    -#   HpNum tauMin    [1e-4]  safety threshold of \f$ \tau:~ \tau_{\min}\f$

    -#   HpNum tauFactor [0.5]   factor \f$ \tau_{f} \f$

    -#   HpNum tauIter   [50]    $\tau_p$

    -#   HpNum m         [1e-3]  descent parameter */
   {
    DfltdSfInpt( iStrm , tauInit , HpNum( 1 ) );
    DfltdSfInpt( iStrm , tauMin , HpNum( 1e-4 ) );
    DfltdSfInpt( iStrm , tauFactor , HpNum( 0.5 ) );
    DfltdSfInpt( iStrm , tauIter , Index( 100 ) );
    DfltdSfInpt( iStrm , m , HpNum( 0.1 ) );
    }

 /*@} -----------------------------------------------------------------------*/
 /*-------------------------- OTHER INITIALIZATIONS -------------------------*/
 /*--------------------------------------------------------------------------*/
 /** @name Other initializations
   @{ */

   inline void SetVOLLog( ostream *outs = 0 , const char lvl = 0 );

/*--------------------------------------------------------------------------*/

   void Format( void )
   {
    Deflection::Format( );
    alpha = lastvalue = Inf<HpNum>();  // deflection coefficient is unknown
    tau = tauInit;
    }

/*@} -----------------------------------------------------------------------*/
/*-------------------- METHODS FOR DEFLECTION COMPUTATION ------------------*/
/*--------------------------------------------------------------------------*/
/** @name Reading the solution
     @{ */

   inline void NewDEF( void );

/*@} -----------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/
/** @name Reading the solution
    @{ */

   inline cHpNum GetDFLCoeff( void );

   inline const bool DoSS( void );

   inline cHpNum Delta( void );

/*@} -----------------------------------------------------------------------*/
/*--------------------- PRIVATE PART OF THE CLASS --------------------------*/
/*--------------------------------------------------------------------------*/

 private:

/*--------------------------------------------------------------------------*/
/*--------------------------- PRIVATE METHODS ------------------------------*/
/*--------------------------------------------------------------------------*/

   inline void Solve( void );

/*--------------------------------------------------------------------------*/
/*----------------------- PRIVATE DATA STRUCTURES  -------------------------*/
/*--------------------------------------------------------------------------*/

   HpNum tauInit;     // initial tau
   HpNum tauMin;      // the minimum value of tau
   HpNum tauFactor;   // the decreasing factor
   Index tauIter;     // the number of iterations after a decrease
   HpNum m;           // descent parameter

   HpNum tau;         // the safeguarded parameter tau
   HpNum lastvalue;   // last function value Fi()
   HpNum alpha;       // deflection coefficient

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

 };  // end( class Volume )

/*--------------------------------------------------------------------------*/
/*------------------- inline methods implementation ------------------------*/
/*--------------------------------------------------------------------------*/

inline void Volume::SetVOLLog( ostream *outs , const char lvl )
{
 Deflection::SetVOLLog( outs , lvl );

 #if( LOG_VOL )
  if( VOLLLvl > 1 )
   *VOLLog << endl <<  "Volume: tau = " << tauInit << " in ["<< tauMin
	   << " , -] with ( " << tauFactor << " , " << tauIter << " ) ~ m = "
	   << m << endl;
 #endif
 }

/*--------------------------------------------------------------------------*/

inline void Volume::NewDEF(  void ) {
 // change tau if needed - - - - - - - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( Solver->NrIter() <= 1 ) {
  alpha = 1;
  return;
  }

 HpNum FiRec = Solver->ReadBestFiVal();

 // tau should be decreased if after tauIter iterations Fi() is decreased
 // less than 1% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( ! ( Solver->NrIter() % tauIter ) )
  if( tau > tauMin ) {
   cHpNum x = ( lastvalue - FiRec ) / ABS( FiRec );
   lastvalue = FiRec;
   if( x <= 0.01 ) {
    tau *= tauFactor;
    if( tau < tauMin )
     tau = tauMin;
    }
   }

 // solve a quadratic problem and get the deflection coefficient - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 Solve();

 }  // end( Volume::NewDEF() )

/*--------------------------------------------------------------------------*/

inline cHpNum Volume::GetDFLCoeff( void )
{
 return( alpha );
 }

/*--------------------------------------------------------------------------*/

inline const bool Volume::DoSS( void )
{
 bool serious;
 HpNum Delta;
 Solver->GetPar( NDOSolver::kEpsLin , Delta );
 Delta *= m * max( 1.0 , ABS( Solver->ReadBestFiVal() ) );

 // if there is an improvement of Fi(), a SS occurs; otherwise, a NS occurs
 if( ( Solver->ReadFiVal() - ReadFVal() ) >= Delta )
  serious = true;
 else
  serious = false;

 return( serious );

 }  // end( Volume::DoSS() )

/*--------------------------------------------------------------------------*/

inline const HpNum Volume::Delta( void )
{
 HpNum EpsLin;
 Solver->GetPar( NDOSolver::kEpsLin , EpsLin );
 return( m * EpsLin * max( 1.0 , ABS( Solver->ReadBestFiVal() ) ) );
 }

/*--------------------------------------------------------------------------*/

inline void Volume::Solve( void )
{
 HpNum GiNorm = GetGiNorm();
 HpNum DNorm = GetDNorm();
 HpNum GixD = GetdGk();
 HpNum Sigma = GetSigma();
 HpNum Epsilon = GetEpsilon();

 // compute the proximity parameter  - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 HpNum proximity = GetStepsize() == Inf<HpNum>() ? 1 : GetStepsize();

 // solving the cutting plane model  - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 HpNum den = proximity * ( GiNorm + DNorm - 2.0 * GixD );

 if( den <= 1e-8 )
  alpha = 1;
 else
  alpha = ( Epsilon - Sigma - proximity * ( GixD - DNorm ) ) / den;

 // tau should be projected over [0,1] - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( alpha > 1 )
  alpha = min( tau , 1.0 );
 else
  if( alpha < 1e-8 )
   alpha = tau / 10;

 }  // end ( void Solve )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#if( OPT_USE_NAMESPACES )
 };  // end( namespace NDO_di_unipi_it )
#endif

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif  /* Volume.h included */

/*--------------------------------------------------------------------------*/
/*------------------------- End File Volume.h ------------------------------*/
/*--------------------------------------------------------------------------*/

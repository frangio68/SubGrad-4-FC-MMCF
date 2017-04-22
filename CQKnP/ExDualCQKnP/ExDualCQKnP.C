/*--------------------------------------------------------------------------*/
/*-------------------------- File ExDualCQKnP.C ----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * 
 * Continuous Quadratic Knapsack Problems (CQKnP) solver based on the
 * standard dual-ascent approach, which extends the DualCQKnP class to
 * non-negative quadratic costs and extended real bounds.
 *
 * \version 1.07
 *
 * \date 30 - 04 - 2012
 *
 * \author Antonio Frangioni \n
 *         Operations Research Group \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Enrico Gorgone \n
 *         Dipartimento di Elettronica Informatica e Sistemistica \n
 *         Universita' della Calabria \n
 *
 * Copyright &copy 2011 - 2012 by Antonio Frangioni, Enrico Gorgone.
 */

/*--------------------------------------------------------------------------*/
/*--------------------------- IMPLEMENTATION -------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "ExDualCQKnP.h"

#if( CQKnPClass_LOG )
 #include <iomanip>
#endif

#include <algorithm>

/*--------------------------------------------------------------------------*/
/*-------------------------------- USING -----------------------------------*/
/*--------------------------------------------------------------------------*/

using namespace CQKnPClass_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*------------------------------ LOCAL MACROS ------------------------------*/
/*--------------------------------------------------------------------------*/

#if( CQKnPClass_LOG )
 #define KLOG( l , x ) if( KNPLLvl > l ) *KNPLog << x
#else
 #define KLOG( l , x )
#endif

/*--------------------------------------------------------------------------*/

#if DualCQKnP_SANITY_CHECKS == 0
 #define SanityCheckB()

 #define SanityCheckC()
#endif

/*--------------------------------------------------------------------------*/
/*------------------------------- CONSTANTS --------------------------------*/
/*--------------------------------------------------------------------------*/
/* Only the first three bits of the status flag are used by the base class
   (the values from 0 to 5), so the other bits can be used to signal various
   conditions. */

static const int StatMsk =   7;  // mask for the "official" status bits
static const int Hv2Sort =   8;  // if we need to sort
static const int Hv2ChkP =  16;  // if we need to check primal feasibility
static const int Hv2ChkD =  32;  // if we need to check dual feasibility
static const int Hv2CstI =  64;  // if we need to construct I
static const int HvWrtX  = 128;  // if we know the primal solution

/*--------------------------------------------------------------------------*/
/*-------------------- IMPLEMENTATION OF ExDualCQKnP -----------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- PUBLIC METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/

const double *ExDualCQKnP::KNPGetX( void )
{
 if( ( ( status & StatMsk ) == kOK ) && ( ! ( status & HvWrtX ) ) ) {
  status |= HvWrtX;

  double beta = McB; // derivative of Lagrangian function
  int* up = new int[ n ];
  int* tup = up;
  int lb = Inf<int>();

  for( int i = 0 ; i < n ; i++ ) {
   if( D[ i ] > 0 ) {
    const double muh1 = A[ i ] > - Inf<double>() ? OV[ i ]: - Inf<double>();
    if( muStar < muh1 )
     XSol[ i ] = A[ i ];
    else {
     const double muh2 = B[ i ] < Inf<double>() ? OV[ i + n ]: Inf<double>();
     if( muStar < muh2 )
      XSol[ i ] = 0.5 * ( muStar - C[ i ] ) / D[ i ];
     else
      XSol[ i ] = B[ i ];
     }
    }
   else {  // D[ i ] == 0
    const double muh1 = C[ i ];
    if( muStar > muh1 )
     XSol[ i ] = B[ i ];
    else
     if( muStar < muh1 )
      XSol[ i ] = A[ i ];
     else {  // muStar == muh1
      if( B[ i ] < Inf<double>() ) {
       *(tup++) = i;
       XSol[ i ] = B[ i ];
       }
      else {
       lb = i;
       if( A[ i ] > -Inf<double>() )
	XSol[ i ] = A[ i ];
       else
	XSol[ i ] = 0;
       }
      }
    }

   beta -= XSol[ i ];

   }  // end( for )

  *tup = Inf<double>();

  // lead beta to zero - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  tup = up;
  if( beta > DefEps ) {  // we need to decrease beta
   if( sense | muStar != 0 )
    XSol[ lb ] += beta;
   }
  else
   if( beta < -DefEps ) {  // we need to increase beta
    if( lb < Inf<int>() && A[ lb ] == -Inf<double>() )
     XSol[ lb ] += beta;
    else
     for( int i ; ( i = *(tup++) ) < Inf<int>(); )
      if( A[ i ] <= -Inf<double>() ) {
       XSol[ i ] += beta;
       break;
       }
      else
       if( (-beta ) >  B[ i ] - A[ i ] ) {
	XSol[ i ] -= B[ i ] - A[ i ];
	beta += B[ i ] - A[ i ];
        }
       else {
	XSol[ i ] += beta;
	break;
        }
    }

  delete[] up;
  }

 return( XSol );

 } // end( ExDualCQKnP::KNPGetX )

/*--------------------------------------------------------------------------*/
/*-------------------------- PROTECTED METHODS -----------------------------*/
/*--------------------------------------------------------------------------*/

bool ExDualCQKnP::CheckPFsb( void )
{
 double sumA = 0;
 double sumB = 0;

 for( int k = 0 ; k < n ; k++ ) {
  if( A[ k ] > B[ k ] )
   return( false );

  if( sumA > - Inf<double>() )
   if( A[ k ] > - Inf<double>() )
    sumA += A[ k ];
   else
    sumA = -Inf<double>();

  if( sumB < Inf<double>() )
   if( B[ k ] < Inf<double>() )
    sumB += B[ k ];
   else
    sumB = Inf<double>();
  }

 return( ( sumA <= McB ) && ( !sense || ( sumB >= McB ) ) );

 }  // end( ExDualCQKnP::CheckPFsb )

/*--------------------------------------------------------------------------*/

bool ExDualCQKnP::CheckDFsb( void )
{
 UB = Inf<double>();
 LB = - Inf<double>();

 for( int k = 0 ; k < n ; k++ )
  if( D[ k ] == 0 ) {
   if( ( A[ k ] == - Inf<double>() ) && ( C[ k ] > LB ) )
    LB = C[ k ];

   if( ( B[ k ] == Inf<double>() ) && ( C[ k ] < UB ) )
    UB = C[ k ];
   }

 if( ! sense )
  UB = std::min( UB , double( 0 ) );

 return( LB <= UB );

 }  // end( ExDualCQKnP::CheckDFsb )

/*--------------------------------------------------------------------------*/

#if DualCQKnP_SANITY_CHECKS

void ExDualCQKnP::SanityCheckB( void )
{
 for( int k = 0 ; k < n ; k++ ) {
  if( A[ k ] == Inf<double>() )
   throw( CQKException(
	            "ExDualCQKnP::LoadSet(): lower bounds must be < INF" ) );

  if( B[ k ] == -Inf<double>() )
   throw( CQKException(
	          "ExDualCQKnP::LoadSet(): upper bounds must be < - INF" ) );
  }
 }

/*--------------------------------------------------------------------------*/

void ExDualCQKnP::SanityCheckC( void )
{
 for( int k = 0 ; k < n ; k++ ) {
  if( D[ k ] < 0 )
   throw( CQKException(
                   "ExDualCQKnP::LoadSet(): the instance is not convex" ) );

  if( D[ k ] == Inf<double>() )
   throw( CQKException(
	         "ExDualCQKnP::LoadSet(): quadratic costs must be finite" ) );

  if( ( C[ k ] == Inf<double>() ) || ( C[ k ] == - Inf<double>() ) )
   throw( CQKException(
	            "ExDualCQKnP::LoadSet(): linear costs must be finite" ) );
  }
 }

#endif

/*--------------------------------------------------------------------------*/

void ExDualCQKnP::SetName( void )
{
 // assign the name of items we have to order  - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // if j < 2 * n then D[ i ] > 0 and
 //
 // 0 <= j <  n     , j --->  2 A[ j ] * D[ j ] + C[ j ]
 // n <= j <  2 * n , j --->  2 B[ j - n ] * D[ j - n ] + C[ j - n ]
 //
 // however, these only exist when well-defined (A[ j ] > - INF for the first,
 // B[ j - n ] < INF for the second)
 //
 // if 2 * n <= j < 3 * n then D[ i ] == 0 instead, and
 // j ---> C[ j - 2 * n ]
 //
 // but only when - INF < A[ j - 2 * n ] <= B[ j - 2 * n ] < INF

 if( status & Hv2CstI ) {
  status &= ~Hv2CstI;
  int *tI = I;
  const int n2 = n + n;
  for( int k = 0 ; k < n ; k++ )
   if( D[ k ] > 0 ) {
    if( A[ k ] > - Inf<double>() )
     *(tI++) = k;

    if( B[ k ] < Inf<double>() )
     *(tI++) = k + n;
    }
   else
    if( ( A[ k ] > - Inf<double>() ) && ( B[ k ] < Inf<double>() ) )
     *(tI++) = k + n2;

  *tI = Inf<int>();
  nSort = tI - I;       // how many elements we have to sort
  status |= Hv2Sort;
  }
 }  // end( ExDualCQKnP::SetName )

/*--------------------------------------------------------------------------*/

void ExDualCQKnP::PreSort( void )
{
 // compute values once and for all  - - - - - - - - - - - - - - - - - - - - -

 const int *tI = I;
 for( int k ; ( k = *(tI++) ) < Inf<int>() ; ) {
  if( k < n )
   OV[ k ] = 2 * A[ k ] * D[ k ] + C[ k ];
  else {
   k -= n;
   if( k < n )
    OV[ k + n ] =  2 * B[ k ] * D[ k ] + C[ k ];
   else
    OV[ k - n ] = C[ k - n ];
   }
  }
 }  // end( ExDualCQKnP::PreSort )

/*--------------------------------------------------------------------------*/

void ExDualCQKnP::FindDualSol( void )
{
 status = kOK;
 OptVal = Inf<double>();  // not computed yet

 // initialize the starting point: \mu   - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 const int n2 = n + n;

 double mu;   // current point
 double muk;  // threshold value of k-th item
 muk = nSort ? OV[ I[ 0 ] % n2 ] : Inf<double>();

 if( LB > - Inf<double>() )
  mu = std::min( LB , muk );
 else
  if( nSort )
   mu = std::min( UB , muk );
  else
   if( UB < Inf<double>() )
    mu = UB;
   else
    mu = - Inf<double>();

 KLOG( 1 , std::endl << "muInit = " << mu << std::endl );

 // compute the initial values of variables  - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 double beta = McB;        // derivative of Lagrangian function: phi'( mu )
 double gamma = 0;         // gamma, rate of derivative phi'( mu )

 for( int k = 0 ; k < n ; k++ ) {
  if( A[ k ] > - Inf<double>() )
   beta -= A[ k ];
  else
   if( D[ k ] > 0 ) {
    gamma += 0.5 / D[ k ];
    beta += 0.5 * C[ k ] / D[ k ];
    }
   else
    if( B[ k ] < Inf<double>() )
     beta -= B[ k ];
  }

 // backtracking phase - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( mu == - Inf<double>() ) {
  muStar = beta / gamma;
  return;
  }
 else
  if( LB == -Inf<double>() && ( ( beta - gamma * mu ) < 0 ) ) {
   muStar = beta / gamma; // the case gamma = 0 has been previously considered
   return;
   }

 beta -= gamma * mu; // update derivative of Lagrangian function
 KLOG( 1 , std::endl << " phi'(mu) = " << beta << " - " << gamma << " * mu"
                     << std::endl );

 // move the current point up to lower bound LB  - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 int kIND = 0;
 while( ( mu != LB ) && ( LB != -Inf<double>() ) ) {
  if( I[ kIND ] < n )
   gamma += 0.5 / D[ I[ kIND ] ];
  else
   if( I[ kIND ] >= n && I[ kIND ] < n2 )
    gamma -= 0.5 / D[ I[ kIND ] - n ];
   else
    beta -= ( B[ I[ kIND ] - n2 ] - A[ I[ kIND ] - n2 ] );

  if( kIND == ( nSort - 1 ) ) {
   beta -= gamma * ( LB - muk );
   mu = LB;
   break;
   }
  else {
   kIND++;
   double muk1 = OV[ I[ kIND ] % n2 ];
   beta -= gamma * ( std::min( muk1 , LB ) - muk );
   mu = std::min( muk = muk1 , LB );
   }
  }

 if( beta <= 0 ) {
  muStar = mu;
  return;
  }

 // updating mu up to hat{mu}  - - - - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 for( ; ; ) {
  if( mu == UB ) { // the stationarity point of Lagrangian function is at
   muStar = mu;    // right hand side of the upper bound UB
   break;
   }

  if( ! nSort || ( nSort && mu > muk ) )  {     // move up to UB
   if( ( beta - gamma * ( UB - mu ) ) <= 0 ) {  // the case UB == Inf has
    muStar = mu + ( beta / gamma );             // been previously considered
    break;
    }
   else {
    mu = UB;
    continue;
    }
   }

  if( mu < muk ) {
   double beta1 = beta - gamma * ( std::min( UB , muk ) - mu );
   if( beta1 <= 0 ) {
    muStar = mu + ( beta / gamma );
    break;
    }
   else {
    beta = beta1;
    mu = std::min( muk , UB );
    continue;
    }
   }
  else {
   double beta1 = beta;
   if( I[ kIND ] < n )
    gamma += 0.5 / D[ I[ kIND ] ];
   else
    if( I[ kIND ] >= n && I[ kIND ] < n2 )
     gamma -= 0.5 / D[ I[ kIND ] - n ];
    else
     beta1 -= ( B[ I[ kIND ] - n2 ] - A[ I[ kIND ] - n2 ] );

   KLOG( 1 , std::endl << " phi'(mu) = " << beta << " - " << gamma
	               << " * ( mu - " << mu << " ) " << std::endl );

   if( beta1 <= 0 ) {
    muStar = mu;
    break;
    }

   if( ( nSort - 1 ) == kIND ) {
    if( UB == Inf<double>()  ) {       // the case gamma = 0
     muStar = mu + ( beta1 / gamma );  // has been previously considered
     break;
     }
    mu = UB;
    continue;
    }

   kIND++;
   double muk1 = OV[ I[ kIND ] % n2 ];

   beta = beta1 - gamma * ( std::min( muk1 , UB ) - mu );
   if( beta  <= 0 ) {
    muStar = mu + ( beta1 / gamma );
    break;
    }
   else {
    mu = std::min( muk = muk1 , UB );
    continue;
    }
   }  // end( if( mu == muk ) )
  }  // end( for( ever ) )
 }  // end( ExDualCQKnP::FindDualSol )

/*--------------------------------------------------------------------------*/
/*---------------------- End File ExDualCQKnP.C ----------------------------*/
/*--------------------------------------------------------------------------*/

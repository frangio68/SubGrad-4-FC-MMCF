/*--------------------------------------------------------------------------*/
/*-------------------------- File DualCQKnP.C ------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * 
 * Basis Continuous Quadratic Knapsack Problems (CQKnP) solver based on the
 * standard dual-ascent approach. This class is restricted to instances where
 * *all* items have *strictly positive* quadratic costs and *finite* bounds
 * (both lower and upper).
 *
 * \version 1.08
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

#include "DualCQKnP.h"

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

 #define Log1()

 #define Log2()

 #define Log3()

 #define Log4()
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
static const int Hv2CstI =  64;  // if we need to construct I[]
static const int HvWrtX  = 128;  // if we know the primal solution

/*--------------------------------------------------------------------------*/
/*----------------------------- FUNCTIONS ----------------------------------*/
/*--------------------------------------------------------------------------*/
/*--                                                                      --*/
/*--  These functions are not implemented as methods of the class, since  --*/
/*--  they don't use directly its data structures.                        --*/
/*--                                                                      --*/
/*--------------------------------------------------------------------------*/

template<class T>
inline void VectAssign( T *const g , const T x , const int n )
{
 // g[ i ] = x for each i = 0 .. n - 1

 for( T *tg = g + n ; tg > g ; )
  *(--tg) = x;
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

template<class T1, class T2>
inline void VectAssign( T1 *g1 , const T2 *g2 , int n )
{
 // g1 := g2

 for( ; n-- ; )
  *(g1++) = *(g2++);
 }

/*--------------------------------------------------------------------------*/
/*------------------------- AUXILIARY CLASSES ------------------------------*/
/*--------------------------------------------------------------------------*/

#if ! DualCQKnP_WHCH_QSORT

struct myLess {
 // comparison operator for ordering items: the value is stored in a vector

 myLess( const int tn , const double *ov ) { n = tn; OV = ov; }

 bool operator()( const int x , const int y ) const
 {
  return( OV[ x % (2 * n) ] < OV[ y % (2 * n ) ] );
  }

 int n;
 const double *OV;
 };

#endif

/*--------------------------------------------------------------------------*/
/*-------------------- IMPLEMENTATION OF DualCQKnP -------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- PUBLIC METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/
/*---------------------------- CONSTRUCTOR ---------------------------------*/
/*--------------------------------------------------------------------------*/

DualCQKnP::DualCQKnP( const bool sort , const double Eps )
           :
           CQKnPClass()
{
 WSort = sort;  // use QuickSort by default
 DefEps = Eps;
 #if DualCQKnP_WHCH_QSORT
  InstCntr++;
 #endif

 }  // end( DualCQKnP )

/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

void DualCQKnP::LoadSet( const int pn ,
			 const double *pC , const double *pD ,
			 const double *pA , const double *pB ,
			 const double pV , const bool sns )
{
 // allocating and deallocating memory- - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( n != pn ) {
  if( n )
   MemDeAlloc( );

  n = pn;
  if( n )
   MemAlloc();
  else
   return;      // just sit down in the corner and wait
  }

 // load instance  - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 McB = pV;     // setting volume
 sense = sns;  // setting sense

 // assign the costs and bounds  - - - - - - - - - - - - - - - - - - - - - -

 if( pC )
  VectAssign( C , pC , n );
 else
  VectAssign( C , double( 0 ) , n );

 if( pD )
  VectAssign( D , pD , n );
 else
  VectAssign( D , double( 0 ) , n );

 if( pA )
  VectAssign( A , pA , n );
 else
  VectAssign( A , - Inf<double>() , n );

 if( pB )
  VectAssign( B , pB , n );
 else
  VectAssign( B , + Inf<double>() , n );

 // initialize variables - - - - - - - - - - - - - - - - - - - - - - - - - -

 status = kUnSolved | Hv2Sort | Hv2ChkP | Hv2ChkD | Hv2CstI;

 }  // end( DualCQKnP::LoadSet )

/*--------------------------------------------------------------------------*/

#if CQKnPClass_LOG

void DualCQKnP::SetKNPLog( std::ostream *outs , const char lvl )
{
 CQKnPClass::SetKNPLog( outs , lvl );

 if( KNPLLvl > 1 ) {
  *KNPLog << std::endl << "Sort algorithm: ";
  if( WSort )
   *KNPLog << "Quick";
  else
   *KNPLog << "Bubble";
  *KNPLog << " Sort " << std::endl;
  }
 }

#endif

/*--------------------------------------------------------------------------*/
/*-------------------- METHODS FOR SOLVING THE PROBLEM ---------------------*/
/*--------------------------------------------------------------------------*/

DualCQKnP::CQKStatus DualCQKnP::SolveKNP( void )
{
 if( ! n )
  throw( CQKException( "DualCQKnP::SolveKNP: no instance loaded yet" ) );

 if( ( status & StatMsk ) < kUnSolved )
  return( CQKStatus( status & StatMsk ) );

 SetName();

 // pre-process the solution: meanwhile, find the overall smallest item- - -
 // and put it in the first position - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 Log3();

 // check primal feasibility - - - - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( status & Hv2ChkP ) {
  status &= ~Hv2ChkP;
  SanityCheckB();
  if( ! CheckPFsb() ) {
   KLOG( 1, std::endl << "Primal infeasible / Dual unbounded" << std::endl );
   status &= ~StatMsk;
   status |= kUnfeasible;
   return( kUnfeasible );
   }
  }

 // check dual feasibility - - - - - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( status & Hv2ChkD ) {
  status &= ~Hv2ChkD;
  SanityCheckC();
  if( ! CheckDFsb() ) {
   KLOG( 1, std::endl << "Dual infeasible / Primal unbounded" << std::endl );
   status &= ~StatMsk;
   status |= kUnbounded;
   return( kUnbounded );
   }
  }

 Log4();

 // find the optimal solution, by solving dual problem  - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( nSort && ( status & Hv2Sort ) ) {
  // sorting phase  - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  status &= ~Hv2Sort;

  KLOG( 2 , std::endl << "Order the following vector of " << nSort <<
		  " items: " << std::endl );
  PreSort();
  Log2();  // show the vector we have to order

  // find the smallest element- - - - - - - - - - - - - - - - - - - - - - - -

  KLOG( 2 , std::endl << "Finding the smallest element ..." << std::endl );
  const int n2 = n + n;

  int *InMin = I;
  double Min = Inf<double>();
  for( int *tI = I ; *tI < Inf<int>() ; tI++ ) {
   const double tCh = OV[ *tI % n2 ];
   if( tCh < Min ) {
    Min = tCh;
    InMin = tI;
    }
   }

  std::swap( *InMin , *I );
  Log2();

  // sort the rest (if any)- - - - - - - - - - - - - - - - - - - - - - - - -

  if( nSort > 2 )         // two-elements vectors are already sorted
   if( nSort > 3 )
    if( WSort ) {
     KLOG( 2 , std::endl << "Sort algoritm: Quick Sort" << std::endl );
     qsort();
     }
    else {
     KLOG( 2 , std::endl << "Sort algoritm: Bubble Sort" << std::endl );
     bsort();
     }
   else {                 // special treatment for the case l == 3
    const double p1 = OV[ I[ 1 ] % n2 ];
    const double p2 = OV[ I[ 2 ] % n2 ];
    if( p1 > p2 )
     std::swap( I[ 1 ] , I[ 2 ] );
    }

  }  // end ( sorting )

 FindDualSol();

 KLOG( 1 , std::endl << "Opt. dual sol.: " << muStar << std::endl );

 return( CQKStatus( status ) );

 }  // end( DualCQKnP::SolveKNP )

/*--------------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/

const double *DualCQKnP::KNPGetX( void )
{
 if( ( ( status & StatMsk ) == kOK ) && ( ! ( status & HvWrtX ) ) ) {
  status |= HvWrtX;

  for( int i = 0 ; i < n ; i++ ) {
   const double muh1 =  OV[ i ];
   if( muStar < muh1 )
    XSol[ i ] = A[ i ];
   else {
    const double muh2 = OV[ i + n ];
    if( muStar < muh2 )
     XSol[ i ] = 0.5 * ( muStar - C[ i ] ) / D[ i ];
    else
     XSol[ i ] = B[ i ];
    }
   } // end( for )
  }

 return( XSol );

 } // end( DualCQKnP::KNPGetX )

/*--------------------------------------------------------------------------*/

double DualCQKnP::KNPGetPi( void )
{
 return( muStar );

 }  // end( DualCQKnP::KNPGetPi )

/*--------------------------------------------------------------------------*/

double DualCQKnP::KNPGetFO( void )
{ 
 if( ( status & StatMsk ) == kUnfeasible )
  return( Inf<double>() );

 if( ( status & StatMsk ) == kUnbounded )
  return( - Inf<double>() );

 if( ( status & StatMsk ) != kOK )
  throw( CQKException( "DualCQKnP::KNPGetFO: no solution available" ) );

 if( ! ( status & HvWrtX ) )
  KNPGetX();

 if( OptVal == Inf<double>() ) {  // not computed yet
  OptVal = 0;

  for( int i = 0 ; i < n ; i++ )
   OptVal +=  C[ i ] * XSol[ i ] + D[ i ] * XSol[ i ] * XSol[ i ];

  Log1();
  KLOG( 1 , std::endl << "Opt. value: " << OptVal << std::endl );
  }

 return( OptVal );

 } // end( DualCQKnP::KNPGetFO )

/*--------------------------------------------------------------------------*/
/*-------------- METHODS FOR READING THE DATA OF THE PROBLEM ---------------*/
/*--------------------------------------------------------------------------*/

void DualCQKnP::KNPLCosts( double *csts , const int *nms , int strt , int stp )
{ 
 if( stp > n )
  stp = n;

 if( strt < 0 )
  strt = 0;

 if( nms ) {  
  while( *nms < strt )
   nms++;

  for( int i ; ( i = *( nms++ ) ) < stp ; )
   *csts++ = C[ i ];
  }
 else
  for( int i = strt ; i < stp ; i++ )
   *csts++ = C[ i ];

 } // end( DualCQKnP::KNPLCosts )

/*--------------------------------------------------------------------------*/

void DualCQKnP::KNPQCosts( double *csts , const int *nms , int strt , int stp )
{
 if( stp > n )
  stp = n;

 if( strt < 0 )
  strt = 0;

 if( nms ) {
  while( *nms < strt )
   nms++;

  for( int i ; ( i = *( nms++ ) ) < stp ; )
   *csts++ = D[ i ];
  }
 else
  for( int i = strt ;  i < stp ; i++ )
   *csts++ = D[ i ];

 } // end( DualCQKnP::KNPQCosts )

/*--------------------------------------------------------------------------*/

void DualCQKnP::KNPLBnds( double *bnds , const int *nms , int strt , int stp )
{
 if( stp > n )
  stp = n;

 if( strt < 0 )
  strt = 0;

 if( nms ) {
  while( *nms < strt )
   nms++;

  for( int i ; ( i = *( nms++ ) ) < stp ; )
   *bnds++ = A[ i ];
  }
 else
  for( int i = strt ; i < stp ; i++ )
   *bnds++ = A[ i ];

 }  // end( DualCQKnP::KNPLBnds )

/*--------------------------------------------------------------------------*/

void DualCQKnP::KNPUBnds( double *bnds , const int *nms , int strt , int stp )
{
 if( stp > n )
  stp = n;

 if( strt < 0 )
  strt = 0;

 if( nms ) {
  while( *nms < strt )
   nms++;

  for( int i ; ( i = *( nms++ ) ) < stp ; )
   *bnds++ = B[ i ];
  }
 else
  for( int i = strt ; i < stp ; i++ )
   *bnds++ = B[ i ];

 }  // end( DualCQKnP::KNPUBnds )

/*--------------------------------------------------------------------------*/
/*------------- METHODS FOR ADDING / REMOVING / CHANGING DATA --------------*/
/*--------------------------------------------------------------------------*/

void DualCQKnP::ChgLCosts( const double *csts , const int *nms ,
		      int strt , int stp )
{
 if( stp > n )
  stp = n;

 if( strt < 0 )
  strt = 0;

 if( nms ) {
  while( *nms < strt )
   nms++;

  for( int i ; ( i = *( nms++ ) ) < stp ; )
   C[ i ] = *csts++;
  }
 else
  for( int i = strt ; i < stp ; i++ )
   C[ i ] = *csts++;

 status |= ( Hv2ChkD | Hv2Sort );
 if( ! ( ( status & StatMsk ) == kUnfeasible ) || ( status & Hv2ChkP ) ) {
  status &= ~StatMsk;
  status |= kUnSolved;
  }
 }  // end( DualCQKnP::ChgLCosts )

/*--------------------------------------------------------------------------*/

void DualCQKnP::ChgQCosts( const double *csts , const int *nms ,
			  int strt , int stp )
{
 if( stp > n )
  stp = n;

 if( strt < 0 )
  strt = 0;

 if( nms ) {
  while( *nms < strt )
   nms++;

  for( int i ; ( i = *( nms++ ) ) < stp ; )
   D[ i ] = *csts++;
  }
 else
  for( int i = strt ; i < stp ; i++ )
   D[ i ] = *csts++;

 status |= ( Hv2ChkD | Hv2CstI );
 if( ! ( ( status & StatMsk ) == kUnfeasible ) || ( status & Hv2ChkP ) ) {
  status &= ~StatMsk;
  status |= kUnSolved;
  }
 }  // end( DualCQKnP::ChgQCosts )

/*--------------------------------------------------------------------------*/

void DualCQKnP::ChgLCost( int item , const double cst )
{
 if( ( item >= 0 ) && ( item < n ) ) {
  C[ item ] = cst;

  status |= ( Hv2ChkD | Hv2Sort );
  if( ! ( ( status & StatMsk ) == kUnfeasible ) || ( status & Hv2ChkP ) ) {
   status &= ~StatMsk;
   status |= kUnSolved;
   }
  }
 } // end( DualCQKnP::ChgLCost )

/*--------------------------------------------------------------------------*/

void DualCQKnP::ChgQCost( int item , const double cst )
{
 if( ( item >= 0 ) && ( item < n ) ) {
  D[ item ] = cst;

  status |= ( Hv2ChkD | Hv2CstI );
  if( ! ( ( status & StatMsk ) == kUnfeasible ) || ( status & Hv2ChkP ) ) {
   status &= ~StatMsk;
   status |= kUnSolved;
   }
  }
 } // end( DualCQKnP::ChgQCost )

/*--------------------------------------------------------------------------*/

void DualCQKnP::ChgLBnds( const double *bnds , const int *nms ,
			  int strt , int stp  )
{
 if( stp > n )
  stp = n;

 if( strt < 0 )
  strt = 0;

 if( nms ) {
  while( *nms < strt )
   nms++;

  for( int i ; ( i = *( nms++ ) ) < stp ; )
   A[ i ] = *bnds++;
  }
 else
  for( int i = strt ; i < stp ; i++ )
   A[ i ] = *bnds++;

 status |= ( Hv2ChkP | Hv2ChkD | Hv2CstI );
 status &= ~StatMsk;
 status |= kUnSolved;

 }  // end( DualCQKnP::ChgLBnds )

/*--------------------------------------------------------------------------*/

void DualCQKnP::ChgUBnds( const double *bnds , const int *nms ,
			  int strt , int stp  )
{
 if( stp > n )
  stp = n;

 if( strt < 0 )
  strt = 0;

 if( nms ) {
  while( *nms < strt )
   nms++;

  for( int i ; ( i = *( nms++ ) ) < stp ; )
   B[ i ] = *bnds++;
  }
 else
  for( int i = strt ; i < stp ; i++ )
   B[ i ] = *bnds++;

 status |= ( Hv2ChkP | Hv2ChkD | Hv2CstI );
 status &= ~StatMsk;
 status |= kUnSolved;

 }  // end( DualCQKnP::ChgUBnds )

/*--------------------------------------------------------------------------*/

void DualCQKnP::ChgLBnd( int item , const double bnd )
{
 if( ( item >= 0 ) && ( item < n ) ) {
  A[ item ] = bnd;

  status |= ( Hv2ChkP | Hv2ChkD | Hv2CstI );
  status &= ~StatMsk;
  status |= kUnSolved;
  }
 } // end( DualCQKnP::ChgLBnd )

/*--------------------------------------------------------------------------*/

void DualCQKnP::ChgUBnd( int item , const double bnd )
{
 if( ( item >= 0 ) && ( item < n ) ) {
  B[ item ] = bnd;

  status |= ( Hv2ChkP | Hv2ChkD | Hv2CstI );
  status &= ~StatMsk;
  status |= kUnSolved;
  }
 } // end( DualCQKnP::ChgUBnd )

/*--------------------------------------------------------------------------*/

void DualCQKnP::ChgVlm( const double NVlm )
{
 McB = NVlm;

 status |= Hv2ChkP;
 if( ! ( ( status & StatMsk ) == kUnbounded ) || ( status & Hv2ChkD ) ) {
  status &= ~StatMsk;
  status |= kUnSolved;
  }
 } // end( DualCQKnP::ChgVlm() )

/*--------------------------------------------------------------------------*/
/*------------------------------ DESTRUCTOR --------------------------------*/
/*--------------------------------------------------------------------------*/

DualCQKnP::~DualCQKnP()
{
 #if DualCQKnP_WHCH_QSORT
  InstCntr--;
 #endif
 if( n )
  MemDeAlloc();

 } // end ( ~DualCQKnP )

/*--------------------------------------------------------------------------*/
/*-------------------------- PROTECTED METHODS -----------------------------*/
/*--------------------------------------------------------------------------*/

bool DualCQKnP::CheckPFsb( void )
{
 double sumA = 0;
 double sumB = 0;

 for( int k = 0 ; k < n ; k++ ) {
  if( A[ k ] > B[ k ] )
   return( false );

  sumA += A[ k ];
  sumB += B[ k ];

  }

 return( ( sumA <= McB ) && ( !sense || ( sumB >= McB ) ) );

 }  // end( DualCQKnP::CheckPFsb )

/*--------------------------------------------------------------------------*/

bool DualCQKnP::CheckDFsb( void )
{
 UB = Inf<double>();
 LB = - Inf<double>();

 if( ! sense )
  UB = std::min( UB , double( 0 ) );

 return( true );

 }  // end( DualCQKnP::CheckDFsb )

/*--------------------------------------------------------------------------*/

#if DualCQKnP_SANITY_CHECKS

void DualCQKnP::SanityCheckB( void )
{
 for( int k = 0 ; k < n ; k++ ) {
  if( A[ k ] == -Inf<double>() ||  A[ k ] == Inf<double>() )
   throw( CQKException( "DualCQKnP::LoadSet(): lower bounds must be finite" )
	  );

  if( B[ k ] == -Inf<double>() ||  B[ k ] == Inf<double>() )
   throw( CQKException( "DualCQKnP::LoadSet(): upper bounds must be finite" )
	  );
  }
 }

/*--------------------------------------------------------------------------*/

void DualCQKnP::SanityCheckC( void )
{
 for( int k = 0 ; k < n ; k++ ) {
  if( D[ k ] <= 0 )
   throw( CQKException(
             "DualCQKnP::LoadSet(): the instance is not strictly convex" ) );

  if( D[ k ] == Inf<double>() )
    throw( CQKException(
 	          "DualCQKnP::LoadSet(): quadratic costs must be finite" ) );

   if( ( C[ k ] == Inf<double>() ) || ( C[ k ] == - Inf<double>() ) )
    throw( CQKException( "DualCQKnP::LoadSet(): linear costs must be finite" )
	   );
  }
 }

#endif

/*--------------------------------------------------------------------------*/

void DualCQKnP::SetName( void )
{
 // assign the name of items we have to order  - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // 0 <= j < n     , j  --->  2 A[ j ] * D[ j ] + C[ j ]
 // n <= j < 2 * n , j  --->    2 B[ j - n ] * D[ j - n ] + C[ j - n ]

 if( status & Hv2CstI ) {
  status &= ~Hv2CstI;
  int *tI = I;
  for( int k = 0 ; k < n ; k++ ) {
   *(tI++) = k;
   *(tI++) = k + n;
   }

  *tI = Inf<int>();
  nSort = tI - I;       // how many elements we have to sort
  status |= Hv2Sort;
  }
 } // end ( DualCQKnP::SetName )

/*--------------------------------------------------------------------------*/

void DualCQKnP::PreSort( void )
{
 // compute values once and for all  - - - - - - - - - - - - - - - - - - - - -

 const int *tI = I;
 for( int k ; ( k = *(tI++) ) < Inf<int>() ; ) {
  if( k < n )
   OV[ k ] = 2 * A[ k ] * D[ k ] + C[ k ];
  else
   OV[ k ] =  2 * B[ k - n ] * D[ k - n ] + C[ k - n ];
  }
 } // end ( DualCQKnP::PrSort )

/*--------------------------------------------------------------------------*/

void DualCQKnP::FindDualSol( void )
{
 status = kOK;
 OptVal = Inf<double>();  // not computed yet

 // initialize the starting point: \mu - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 double mu;   // current point
 double muk;  // threshold value of k-th item

 muk = OV[ I[ 0 ] ];
 mu = std::min( UB , muk );

 KLOG( 1 , std::endl << "muInit = " << mu << std::endl );

 // compute the initial values of variables  - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 double beta = McB;        // derivative of Lagrangian function: phi'( mu )
 double gamma = 0;         // gamma, rate of derivative phi'( mu )

 for( int k = 0 ; k < n ; k++ )
  beta -= A[ k ];

 // backtracking phase - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( ( beta - gamma * mu ) < 0  ) {
  muStar = beta / gamma; // the case gamma = 0 has been previously considered
  return;
  }

 beta -= gamma * mu; // update derivative of Lagrangian function
 KLOG( 1 , std::endl << " phi'(mu) = " << beta << " - " << gamma << " * mu"
                     << std::endl );

 // updating mu up to hat{mu}  - - - - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 int kIND = 0;
 for( ; ; ) {
  if( mu == UB ) { // the stationarity point of Lagrangian function is at
   muStar = mu;    // right hand side of the upper bound UB
   break;
   }

  if( mu == muk ) {
   if( I[ kIND ] < n )
    gamma += 0.5 / D[ I[ kIND ] ];
   else
    gamma -= 0.5 / D[ I[ kIND ] - n ];

   KLOG( 1 , std::endl << " phi'(mu) = " << beta << " - " << gamma
	               << " * ( mu - " << mu << " ) " << std::endl );

   if( ( nSort - 1 ) == kIND ) {
    if( UB == Inf<double>()  ) {      // the case gamma = 0
     muStar = mu + ( beta / gamma );  // has been previously considered
     break;
     }
    mu = UB;
    continue;
    }

   kIND++;
   double muk1 = OV[ I[ kIND ] ];
   double beta1 = beta;
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
 } // end ( DualCQKnP::FindDualSol )

/*--------------------------------------------------------------------------*/

#if DualCQKnP_WHCH_QSORT

void DualCQKnP::qsort( void )
{
 // a non-recursive QuickSort implementation- - - - - - - - - - - - - - - - -
 // note that the first element, already known to be the smallest, is skipped

 int *Top = QSStck;
 int l = 0;
 int u = nSort - (++l);

 double pvt;  // pivot element: at the beginning, it is the smallest one

 int lwr;
 int upr;

 for(;;) {
  lwr = l + 1;
  upr = u;

  // set the pivot element - - - - - - - - - - - - - - - - - - - - - - - - - -
  pvt = OV[ I[ l ] % ( 2 * n ) ];

  // put the elements smaller thanpivot at left of pivot, while the bigger
  // ones at right of pivot

  while( lwr <= upr ) {
   while( lwr <= upr )
    if( OV[ I[ lwr ] % ( 2 * n ) ] <= pvt )
     lwr++;
    else
     break;

   while( upr >= lwr )
    if( OV[ I[ upr ] % ( 2 * n ) ] > pvt )
     upr--;
    else
     break;

   if( lwr < upr )
    std::swap( I[ lwr++ ] , I[ upr-- ] );
   }

  std::swap( I[ l ] , I[ upr ] );    // swap the pivot at middle point

  // the following is the "recursive" part, by using the stack- - - - - - - - -

  if( --upr > l )    // left recursion
   if( upr - 1 > l ) {
    *(Top++) = l;
    *(Top++) = upr;
    }
   else
    if( OV[ I[ upr ] % ( 2 * n ) ] < OV[ I[ l ] % ( 2 * n ) ] )
     std::swap( I[ l ] , I[ upr ] );

  upr++;

  if( u > ++upr )    // right recursion
   if( u > upr + 1 ) {
    *(Top++) = upr;
    *(Top++) = u;
    }
   else
    if( OV[ I[ u ] % ( 2 * n ) ] < OV[ I[ upr ] % ( 2 * n ) ] )
     std::swap( I[ u ] , I[ upr ] );

  Log2();

  if( Top == QSStck )
   break;

  u = *(--Top);
  l = *(--Top);

  }  // end for( ever )
 }  // end( DualCQKnP::qsort )

#else

void DualCQKnP::qsort( void )
{
 // quick sort (or whatever sort it is) using the STL function- - - - - - - -
 // note that the first element, already known to be the smallest, is skipped

 std::sort( I + 1 , I + nSort , myLess( n , OV ) );
 Log2();
 }

#endif

/*--------------------------------------------------------------------------*/

void DualCQKnP::bsort( void )
{
 // a Bubble Sort implementation- - - - - - - - - - - - - - - - - - - - - - -
 // note that the first element, already known to be the smallest, is skipped

 int *tI = I;
 int *u = tI + nSort - 1;

 for( tI++ ; u > tI ; ) {
  int *j = tI;
  double Cj = OV[ *j % ( 2 * n ) ];

  for( int *i = j ; i < u ; ) {
   int *h = i++;
   const double Ci = OV[ *i % ( 2 * n ) ];

   if( Ci < Cj )
    std::swap( *i , *( j = h ) );
   else
    Cj = Ci;
   }

  u = j;

  Log2();
  }
 }  // end( DualCQKnP::bsort )

/*--------------------------------------------------------------------------*/

inline void DualCQKnP::MemAlloc( void )
{
 C = new double[ n ];
 D = new double[ n ];
 A = new double[ n ];
 B = new double[ n ];
 I = new int[ 2 * n + 1 ];  // I is INF-terminated

 XSol = new double[ n ];
 OV = new double[ 2 * n ];

 #if DualCQKnP_WHCH_QSORT
  if( n > maxvl ) {
   delete[] QSStck;
   maxvl = n;
   QSStck = new int[ 2 * n ];
   }
 #endif
 }

/*--------------------------------------------------------------------------*/

inline void DualCQKnP::MemDeAlloc( void )
{
 #if DualCQKnP_WHCH_QSORT
  if( ! InstCntr ) {
   delete[] QSStck;
   QSStck = 0;
   maxvl = 0;
   }
 #endif

 delete[] OV;
 delete[] XSol;

 delete[] I;
 delete[] B;
 delete[] A;
 delete[] D;
 delete[] C;
 }

/*--------------------------------------------------------------------------*/

#if CQKnPClass_LOG

void DualCQKnP::Log1( void )
{
 if( KNPLLvl > 1 ) {
  KNPLog->precision( 4 );
  *KNPLog << std::endl << "Opt. primal sol.: ";
  for( int i = 0 ; i < n ; i++ ) {
   if( !( i  % 4 ) )
    *KNPLog << std::endl;

   *KNPLog << "x( " <<  i  << " ) [ " << XSol[ i ] << " ] ";
   }

  *KNPLog << std::endl;
  }
 }  // end( DualCQKnP::Log1 )

/*--------------------------------------------------------------------------*/

void DualCQKnP::Log2( void )
{
 if( KNPLLvl > 2 ) {
  KNPLog->precision( 4 );
  int count = 0;
  int *tI = I;
  for( int h ;  ( h = *(tI++) ) < Inf<int>() ; count++ ) {
   if( ! ( count % 4 ) )
    *KNPLog << std::endl;
   if( h < n || h >= ( 2 * n ) )
    *KNPLog << " u1 ( " <<  h % n << " ) ";
   else
    *KNPLog << " u2 ( " <<  h % n << " ) ";
   *KNPLog << " [ " << OV[ h % ( 2 * n ) ] << " ] ~ ";
   }
  *KNPLog << std::endl;
  }
 } // end( DualCQKnP::Log2 )

/*--------------------------------------------------------------------------*/

void DualCQKnP::Log3( void )
{
 if( KNPLLvl > 2 ) {
   KNPLog->precision( 4 );
  *KNPLog << std::endl << "Var: " << n << std::endl;

  *KNPLog << std::endl << "C: ";
  for( int h = 0 ;  h < n ; h++ ) {
   if( ! ( h % 10 ) )
     *KNPLog << std::endl;
   *KNPLog << "C [ " << h << " ] ( " << C[ h ] << " ) ~ ";
   }

  *KNPLog << std::endl << "D: ";
  for( int h = 0 ;  h < n ; h++ ) {
   if( ! ( h % 10 ) )
    *KNPLog << std::endl;
   *KNPLog << "D [ " << h << " ] ( " << D[ h ] << " ) ~ ";
   }

  *KNPLog << std::endl << "A: ";
  for( int h = 0 ;  h < n ; h++ ) {
   if( ! ( h % 10 ) )
    *KNPLog << std::endl;

   if( A[ h ] == -Inf<double>() )
    *KNPLog << "A[ " << h << " ] ( " << "-INF" << " ) ~ ";
   else
    *KNPLog << "A[ " << h << " ] ( " << A[ h ] << " ) ~ ";
   }

  *KNPLog << std::endl << "B: ";
  for( int h = 0 ;  h < n ; h++ ) {
   if( ! ( h % 10 ) )
    *KNPLog << std::endl;

   if( B[ h ] == Inf<double>() )
    *KNPLog << "B[ " << h << " ] ( " << "+INF" << " ) ~ ";
   else
    *KNPLog << "B[ " << h << " ] ( " << B[ h ] << " ) ~ ";
   }

  *KNPLog << std::endl << "Volume: " << McB
	  << " ~ Knapsack constraint sense: ";
  if( sense )
   *KNPLog <<" E "<< std::endl;
  else
   *KNPLog <<" L "<< std::endl;
  }
 }  // end( DualCQKnP::Log3 )

/*--------------------------------------------------------------------------*/

void DualCQKnP::Log4( void )
{
 if( KNPLLvl > 2 ) {
  KNPLog->precision( 4 );
  *KNPLog << std::endl;

  if( LB > -Inf<double>() )
   *KNPLog << "LBmu = " << LB << " ~ UBmu = ";
  else
   *KNPLog << "LBmu = -INF ~ UBmu = ";

  if( UB < Inf<double>() )
   *KNPLog << UB << std::endl;
  else
   *KNPLog << "INF" << std::endl;
  }
 }  // end( DualCQKnP::Log4 )

#endif

/*-------------------------------------------------------------------------*/
/*------------------ initialize static members ----------------------------*/
/*-------------------------------------------------------------------------*/

#if DualCQKnP_WHCH_QSORT
 int DualCQKnP::InstCntr = 0;

 int *DualCQKnP::QSStck = 0;

 int DualCQKnP::maxvl = 0;
#endif

/*--------------------------------------------------------------------------*/
/*---------------------- End File DualCQKnP.C ------------------------------*/
/*--------------------------------------------------------------------------*/

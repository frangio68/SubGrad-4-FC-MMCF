/*--------------------------------------------------------------------------*/
/*--------------------------- File DualCQKnP.h -----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * 
 * Basis Continuous Quadratic Knapsack Problems (CQKnP) solver based on the
 * standard dual-ascent approach. This class is restricted to instances where
 * \e all items have <em>strictly positive</em> quadratic costs and
 * <em>finite</em> bounds (both lower and upper).
 *
 * Conforms to the standard interface for CQKnP solver defined by the abstract
 * base class CQKnpClass.
 *
 * \version 1.08
 *
 * \date 21 - 12 - 2012
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
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __DualCQKnP
 #define __DualCQKnP  /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*--------------------------------- MACROS ---------------------------------*/
/*--------------------------------------------------------------------------*/
/*--                                                                      --*/
/*--  The following macros control many important details of the actual   --*/
/*--  implementation of the algorithm.                                    --*/
/*--                                                                      --*/
/*--------------------------------------------------------------------------*/
/** @defgroup DualCQKnP_MACROS Compile-time switches in DualCQKnP.h
    These macros control some important details of the class interface.
    Although using macros for activating features of the interface is not
    very C++, switching off some unused features may allow some
    implementation to be more efficient in running time or memory.
    @{ */

#define DualCQKnP_WHCH_QSORT 1

/**< If DualCQKnP_WHCH_QSORT == 0, the sort() function of the STL is used,
   otherwise a hand-made non-recursive quick-sort implementation is used.

   \note Other than the performance impact, this choice has a consequence on
         the thread-safety of the code. If DualCQKnP_WHCH_QSORT == 1, some
	 temporary data structures are created that are shared among all
	 "active" instances of DualCQKnP to save on space in the (frequent)
	 case where many instances are simultaneously in memory. This is
	 *not* thread-safe, while DualCQKnP_WHCH_QSORT == 0 is. */

#define DualCQKnP_SANITY_CHECKS 0

/**< If DualCQKnP_SANITY_CHECKS == 1, sanity checks are done each time the
   data of the instance changes to pick up clearly bad values. Otherwise,
   the user will have to be extra careful to avoid them. */

/*@}  end( group( DualCQKnP_MACROS ) ) */ 
/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "CQKnPClass.h"

/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE and USING ----------------------------*/
/*--------------------------------------------------------------------------*/

namespace CQKnPClass_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*---------------------------- CLASSES -------------------------------------*/
/*--------------------------------------------------------------------------*/
/** Solver of Continuous Quadratic Knapsack Problems (CQKnP) based on the
    standard formulation of the dual problem as a piecewise-convex problem
    in the unique multiplier of the knapsack constraint and the corresponding
    obvious dual-ascent approach. This class is restricted to instances where
    \e all items have <em>strictly positive</em> quadratic costs and
    <em>finite</em> bounds (both lower and upper).
 
    Derives from CQKnpClass and therefore it \e mostly conforms to its
    interface, except for refusing to solve instances without the required
    characteristics. */

class DualCQKnP : public CQKnPClass {

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*--------------------------------------------------------------------------*/
/*--------------------------- PUBLIC METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/
/*---------------------------- CONSTRUCTOR ---------------------------------*/
/*--------------------------------------------------------------------------*/

   DualCQKnP( const bool sort = true , const double eps = 1e-6 );

/**< The most important operation for solving the CQKnP with a dual method is
   the sorting of the items for nondecreasing elements

     2 * A[ i ] * D[ i ] + C[ i ] and 2 * B[ i ] * D[ i ] + C[ i ].

   If the knapsack is a large one, this can be (relatively) time-consuming.
   Different sort procedures can be better in different situations, and the
   parameter 'sort' allows to decide which among the available sorting
   procedures has to be used. Possible values of this parameter are:

   false  Bubble Sort: this is O( n^2 ) on average, but it can be very fast -
          O( n ) - if reoptimizing from a previous problem where the order
	  was not very different (e.g., only few costs have changed).

   true   [default] Quick Sort: this is O( n lg n ) on average and pretty
          efficient in practice, but can be very slow - O( n^2 ) - if the
          vector is already (almost) ordered, e.g. when reoptimizing from a
          previous problem where only few costs have changed.

   The parameter Eps defines the precision required to construct the
   solution [default value is 1e-6].

   The choices can be changed at any time with SetSort() and SetEps(),
   respectively [see below]. */

/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

   void LoadSet( const int pn = 0 ,
		 const double *pC = 0 , const double *pD = 0 ,
		 const double *pA = 0 , const double *pB = 0 ,
		 const double pV = 0 , const bool sns = true );

/*--------------------------------------------------------------------------*/

#if CQKnPClass_LOG

   void SetKNPLog( std::ostream *outs = 0 , const char lvl = 0 );

/**< lvl controls the "level of verbosity" of the code. The first four bits
   of lvl have the following meaning:

    0  =>  no log at all (also assumed if log = 0);

    1  =>  "basic" log: only the errors are reported;

    2  =>  the solution is displayed;

    3  =>  a detailed step-by-step log of the algorithm is displayed. */

#endif

/*--------------------------------------------------------------------------*/

   inline void SetSort( const bool WhchSrt = false );

/**< Allows to change the sorting procedures to be used in the next calls to
   SolveKNP(); see the comments to the constructor for details. */

/*--------------------------------------------------------------------------*/

   inline void SetEps(const double Eps = 1e-6 );

/*--------------------------------------------------------------------------*/
/*---------------------- METHODS FOR SOLVING THE PROBLEM -------------------*/
/*--------------------------------------------------------------------------*/

   CQKnPClass::CQKStatus SolveKNP( void );

/*--------------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/

   const double *KNPGetX( void );

   double KNPGetPi( void );

   double KNPGetFO( void );

/*--------------------------------------------------------------------------*/
/*-------------- METHODS FOR READING THE DATA OF THE PROBLEM ---------------*/
/*--------------------------------------------------------------------------*/

   void KNPLCosts( double *csts , const int *nms = 0 ,
		   int strt = 0 , int stp = Inf<int>() );

   void KNPQCosts( double *csts , const int *nms = 0 ,
		   int strt = 0 , int stp = Inf<int>() );

   inline double KNPLCost( const int i );

   inline double KNPQCost( const int i );

   void KNPLBnds( double *bnds , const int *nms = 0 ,
		   int strt = 0 , int stp = Inf<int>() );

   void KNPUBnds( double *bnds , const int *nms = 0 ,
		  int strt = 0 , int stp = Inf<int>() );

   inline double KNPLBnd( const int i );

   inline double KNPUBnd( const int i );

   inline double KNPVlm( void );

/*--------------------------------------------------------------------------*/
/*------------- METHODS FOR ADDING / REMOVING / CHANGING DATA --------------*/
/*--------------------------------------------------------------------------*/

   void ChgLCosts( const double *csts , const int *nms = 0 ,
		   int strt = 0 , int stp = Inf<int>() );

   void ChgQCosts( const double *csts , const int *nms = 0 ,
		   int strt = 0 , int stp = Inf<int>() );

   void ChgLBnds( const double *bnds , const int *nms = 0 ,
		  int strt = 0 , int stp = Inf<int>() );

   void ChgUBnds( const double *bnds , const int *nms = 0 ,
		  int strt = 0 , int stp = Inf<int>() );

/*--------------------------------------------------------------------------*/

   void ChgLCost( int item , const double cst );

   void ChgQCost( int item , const double cst );

   void ChgLBnd( int item , const double bnd );

   void ChgUBnd( int item , const double bnd );

   void ChgVlm( const double NVlm );

/*--------------------------------------------------------------------------*/
/*------------------------------ DESTRUCTOR --------------------------------*/
/*--------------------------------------------------------------------------*/

   ~DualCQKnP();

/*--------------------------------------------------------------------------*/
/*--------------------- PROTECTED PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/
/*--                                                                      --*/
/*-- Nobody should ever look at this part: everything that is under this  --*/
/*-- advice may be changed without notice in any new release of the code. --*/
/*--                                                                      --*/
/*--------------------------------------------------------------------------*/

 protected:

/*--------------------------------------------------------------------------*/
/*-------------------------- PROTECTED METHODS -----------------------------*/
/*--------------------------------------------------------------------------*/

   virtual bool CheckPFsb( void );

   virtual bool CheckDFsb( void );

/*--------------------------------------------------------------------------*/

  #if DualCQKnP_SANITY_CHECKS
   virtual void SanityCheckB( void );  // bounds

   virtual void SanityCheckC( void );  // objective function
  #endif

/*--------------------------------------------------------------------------*/

   virtual void SetName( void );

   virtual void PreSort( void );

   virtual void FindDualSol ( void );

/* These methods are virtual to allow derived classes to extend the base one
   to handling of more "complicated" instances, such as with zero quadratic
   costs and/or infinite bounds. */

/*--------------------------------------------------------------------------*/
/*----------------------- PROTECTED DATA STRUCTURES  -----------------------*/
/*--------------------------------------------------------------------------*/

  double *A;         ///< vector of lower bounds
  double *B;         ///< vector of upper bounds
  double *C;         ///< vector of linear costs
  double *D;         ///< vector of quadratic costs
  double  McB;       ///< volume value
  bool sense;        ///< sense of knapsack constraint

  double LB;         ///< lower bound on dual variable
  double UB;         ///< upper bound on dual variable

  int *I;            ///< optimal ordering
  int nSort;         ///< how many elements we have to sort

  double *OV;        ///< values upon which to order

  double *XSol;      ///< primal solution
  double muStar;     ///< optimal dual solution

  bool WSort;        ///< which sorting procedure is used
  double OptVal;     ///< The Optimal Value

  double DefEps;   ///< precision required to construct the solution

  #if DualCQKnP_WHCH_QSORT
   static int *QSStck;   ///< the stack to simulate recursive calls in QS
   static int InstCntr;  ///< number of active instances

   static int maxvl;     ///< max value of items
  #endif

/*--------------------------------------------------------------------------*/
/*--------------------- PRIVATE PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

 private:

/*--------------------------------------------------------------------------*/
/*-------------------------- PRIVATE METHODS -----------------------------*/
/*--------------------------------------------------------------------------*/

  void qsort( void );  // the Quick Sort

  void bsort( void );  // the Bubble Sort

/*--------------------------------------------------------------------------*/

  inline void MemAlloc( void );

  inline void MemDeAlloc( void );

/*--------------------------------------------------------------------------*/

 #if CQKnPClass_LOG
  void Log1( void );

  void Log2( void );

  void Log3( void );

  void Log4( void );
 #endif

 };  // end( class DualCQKnP )

/*--------------------------------------------------------------------------*/
/*------------------- inline methods implementation ------------------------*/
/*--------------------------------------------------------------------------*/

inline void DualCQKnP::SetSort( const bool WhchSrt )
{
 WSort = WhchSrt;
 }

/*--------------------------------------------------------------------------*/

inline void DualCQKnP::SetEps( const double eps ) {

 DefEps = eps ;
 }

/*--------------------------------------------------------------------------*/

inline double DualCQKnP::KNPLCost( const int i )
{
 return( C[ i ] );
 }

/*--------------------------------------------------------------------------*/

inline double DualCQKnP::KNPQCost( const int i )
{
 return( D[ i ] );
 }

/*--------------------------------------------------------------------------*/

inline double DualCQKnP::KNPLBnd( const int i )
{
 return( A[ i ] );
 }

/*--------------------------------------------------------------------------*/

inline double DualCQKnP::KNPUBnd( const int i )
{
 return( B[ i ] );
 }

/*--------------------------------------------------------------------------*/

inline double DualCQKnP::KNPVlm( void )
{
 return( McB );
 }

/*--------------------------------------------------------------------------*/

 };  // end( namespace KNPClass_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif  /* DualCQKnP.h included */

/*--------------------------------------------------------------------------*/
/*------------------------- End File DualCQKnP.h ---------------------------*/
/*--------------------------------------------------------------------------*/

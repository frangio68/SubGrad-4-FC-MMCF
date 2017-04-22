/*--------------------------------------------------------------------------*/
/*-------------------------- File CQKnPClass.h -----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * 
 * Definition of the abstract base class CQKnPClass, which implements a
 * standard interface for Continuous Quadratic Knapsack Problem (CQKnP)
 * solvers to be implemented as derived classes. 
 *
 * \version 1.05
 *
 * \date 22 - 12 - 2012
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

#ifndef __CQKnPClass
 #define __CQKnPClass /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*--------------------------------- MACROS ---------------------------------*/
/*--------------------------------------------------------------------------*/
/** @defgroup CQKnPClass_MACROS Compile-time switches in CQKnPClass.h
    These macros control some important details of the class interface.
    Although using macros for activating features of the interface is not
    very C++, switching off some unused features may allow some
    implementation to be more efficient in running time or memory.
    @{ */

#define CQKnPClass_LOG 0

/**< If CQKnPClass_LOG > 0, data structures and methods to log the activities
   of the (actual) solver are added to the (abstract) interface. */

/*@}  end( group( CQKnPClass_MACROS ) ) */ 
/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <limits>
#include <exception>

#include <iostream>

/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE and USING ----------------------------*/
/*--------------------------------------------------------------------------*/

namespace CQKnPClass_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*-------------------------- CLASS CQKnPClass ------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/** The class ContQKsk defines an interface for Continuous Quadratic Knapsack
    problems (CQKnP) solvers. The data of the problem consists of a set of
    items \E = ( 0 , ... , n - 1 ), where each item i has:

    -  a linear cost C[ i ];

    - a quadratic cost D[ i ] \f$\in R_+\f$;

    - a lower bound  A[ i ] \f$\in R \cup\f$ -INF, and

    - an upper bound B[ i ] \f$\in R \cup\f$ +INF.

    The problem requires to find the most valuable set of items which fit in
    a knapsack of given volume V, but <em>partly accepting items is
    allowed</em>. The formulation of the problem is therefore:
    \f[
      \min \sum_{i \in E} C[ i ] * X[ i ] + D[ i ] * X[ i ]^2
    \f]
    \f[
      \sum_{i \in E} X[ i ] \leq V  ( = V )
    \f]
    \f[
       A[ i ] \leq  X[ i ] \leq B[ i ]  \hspace{1cm}  i \in E
    \f]
    This is a convex quadratic problem, hence a "easy" one. However it is
    repeatedly solved as a subproblem in many applications, so a fast solver
    may ultimately be useful. */

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

class CQKnPClass {

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/
/*--                                                                      --*/
/*--  The following methods and data are the actual interface of the      --*/
/*--  class: the standard user should use these methods and data only.    --*/
/*--                                                                      --*/
/*--------------------------------------------------------------------------*/

 public:

/*--------------------------------------------------------------------------*/
/*---------------------------- PUBLIC TYPES --------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Public types
    @{ */

/** Small class for exceptions. Derives from std::exception implementing the
    virtual method what() -- and since what is virtual, remember to always
    catch it by reference (catch exception &e) if you want the thing to work.
    CQKException class are thought to be of the "fatal" type, i.e., problems
    for which no solutions exists apart from aborting the program. Other kinds
    of exceptions can be properly handled by defining derived classes with
    more information. */

  class CQKException : public std::exception {
   public:
    CQKException( const char *const msg = 0 ) { errmsg = msg; }

    const char* what( void ) const throw () { return( errmsg ); }

   private:
    const char *errmsg;
   };

/*--------------------------------------------------------------------------*/
/** Public enum describing the possible status of the MCF solver. */

  enum CQKStatus { kOK         = 0 ,  ///< optimal solution found
		   kStopped    = 1 ,  ///< optimization stopped
		   kUnfeasible = 2 ,  ///< problem is unfeasible
		   kUnbounded  = 3 ,  ///< problem is unbounded
		   kError      = 4 ,  ///< error in the solver
		   kUnSolved   = 5    ///< no solution available yet
                   };

/*--------------------------------------------------------------------------*/
/** Small class using std::numeric_limits to extract the "infinity" value of
    a basic type (just use Inf<type>()). */

   template<typename T>
   class Inf {
    public:
     Inf() {}
     operator T() { return( std::numeric_limits<T>::max() ); }
    };

/*@} -----------------------------------------------------------------------*/
/*--------------------------- PUBLIC METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/
/*---------------------------- CONSTRUCTOR ---------------------------------*/
/*--------------------------------------------------------------------------*/

  CQKnPClass( void )

/**< Constructor of the class: builds an instance of a Continuous Quadratic
   Knapsack Problem as an object of the class CQKnPClass.

   After that an object has been constructed, no problem is loaded; this has
   to be done with LoadSet() [see below]. Thus, it is an error to invoke any
   method which requires the presence of a problem (all except those in the
   initializations part). The base class provide a protected field n for the
   current number of items, that is set to 0 in the constructor precisely to
   indicate that no instance is currently loaded. */
   {
    n = 0;
    status = kUnSolved;

    }  // end( CQKnPClass )

/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

  virtual void LoadSet( const int pn = 0 ,
			const double *pC = 0 , const double *pD = 0 ,
			const double *pA = 0 , const double *pB = 0 ,
			const double pV = 0 , const bool sns = true ) = 0;

/**< Inputs a new CQKnP instance.

   The meaning of the parameters is the following:

   - pn   is the current number of items of the set: passing pn == 0 is
          intended as a signal to the solver to deallocate everything and
          wait for new orders; in this case, all the other parameters are
          ignored.

   - pC   is the n-vector of the linear part of item costs; the linear costs
          must be finite (- INF < pC[ i ] < INF); pC == 0 means that all
	  linear costs are 0.

   - pD   is the n-vector of the quadratic part of item costs; the quadratic
          costs must be non-negative and finite (0 <= pD[ i ] < INF);
	  pD == 0 means that all quadratic costs are 0.

   - pA   is the n-vector of the item lower bounds; lower bound must be
          < + INF; pA == 0 means that all lower bound are - INF.

   - pB   is the n-vector of the item upper bounds; upper bounds must be
          > -INF; pB == 0 means that all upper bounds are + INF.

   - pV   is the knapsack volume, which must be finite (- INF < pV < + INF).

   - sns  is the sense of knapsack constraint: true is for an equality
          constraint (default), false for an inequality (<=) constraint.

   This method *must* be called prior to invoking any other method of the
   class, with the exception of SetKNPLog() and SetEps() (not to mention the
   constructor, of course). */

/*--------------------------------------------------------------------------*/

   virtual inline void ReadInstance( std::istream &inFile , bool RBV = false );

/**< Read the instance from file. While virtual, the method is implemented in
   the base class: the instance data is read in temporary data structures and
   then LoadSet() [see above] is called using these. Thus, the method will work
   for any derived class with no modification, since it uses the class-provided
   implmentation of LoadSet(); however, being virtual it can be re-implemented
   for maximum efficiency if desired.

   The method supports two formats of the file: a "simplified" one and a
   "complete" one. The simplified format, assumed when RBV is false, is:

      <number of items (n)>

      for i = 0 to n - 1
          <linear cost of item i>

      for i = 0 to n - 1
          <quadratic cost of item i>

   In this case all lower bounds are zero, all upper bounds are INF, and the
   volume is 1. Otherwise (RBV == true) the file must continue with

      for i = 0 to n - 1
          <lower bound of item i>

      for i = 0 to n - 1
          <upper bound of item i>

      <volume of knapsack>

   Note that this method is "alternative" to LoadSet() (in the sense that the
   latter is invoked inside), so the object is "ready to use" once this method
   returns. */

/*--------------------------------------------------------------------------*/

#if CQKnPClass_LOG

   virtual void SetKNPLog( std::ostream *outs = 0 , const char lvl = 0 )

/**< The class ouputs "log" information onto the ostream pointed by outs.
   lvl controls the "level of verbosity" of the code; lvl == 0 means that
   nothing at all is printed, and values larger than 0 mean increasing
   amounts of information, the specific effect of each value being derived-
   class-dependent. outs == 0 implies lvl == 0. */
   {
    if( ( KNPLog = outs ) )
     KNPLLvl = lvl;
    else
     KNPLLvl = 0;
    }

#endif

/*--------------------------------------------------------------------------*/

   virtual void SetEps( const double eps = 1e-6 ) = 0;

/**< Defines the precision required to construct the solution of the
    knapsack problem. */

/*--------------------------------------------------------------------------*/
/*---------------------- METHODS FOR SOLVING THE PROBLEM -------------------*/
/*--------------------------------------------------------------------------*/

   virtual CQKStatus SolveKNP( void ) = 0;

/**< SolveKNP() attempts to solve the current CQKnP instance, returning the
   status indicating the outcome of the optimization. Possible values are:

   - kUnSolved    SolveKNP() has not been called yet, or the data of the
                  problem has been changed since the last call;

   - kOK          optimization has been carried out succesfully;

   - kStopped     optimization have been stopped before that the stopping
                  conditions of the solver applied, e.g. because of the
		  maximum allowed number of "iterations" have been reached;
		  this is not necessarily an error, as it might just be
		  required to re-call SolveKNP() giving it more "resources"
		  in order to solve the problem;

   - kUnfeasible  the current CQKnP instance is (primal) unfeasible;

   - kUnbounded   if the current CQKnP instance is (primal) unbounded: this
                  can only happen if some of the items have +/-INF bounds
		  and zero quadratic cost.

   - kError       there was an error during the optimization; this typically
                  indicates that computation cannot be resumed, although
		  solver-dependent ways of dealing with solver-dependent
		  errors may exhist. */

/*--------------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/

   virtual const double *KNPGetX( void ) = 0;

/**< Return the best feasible solution of the CQKnP found so far (assuming
   CheckFsb() == true, see above), which is optimal if SolveKNP() has
   returned kOK (see above). */

/*--------------------------------------------------------------------------*/

   virtual double KNPGetPi( void ) = 0;

/**< Returns the optimal dual multiplier of the knapsack constraint; the
   returned value is dependable only if SolveKNP() has  returned kOK (see
   above).*/

/*--------------------------------------------------------------------------*/

   virtual double KNPGetFO( void ) = 0;

/**< Returns the optimal value of the objective function of the problem.

   The method typically returns INF if SolveKNP() == kUnfeasible. It returns 
   a finite feasible value if SolveKNP() == kOK or SolveKNP() == kStopped,
   but in the latter case it depends on the solver whether this is a lower
   or an upper bound on the optimal value. The return value is undefined in
   all other cases. */

/*--------------------------------------------------------------------------*/
/*-------------- METHODS FOR READING THE DATA OF THE PROBLEM ---------------*/
/*--------------------------------------------------------------------------*/
/*--                                                                      --*/
/*--  KNPn()      =  total number of items in the current set             --*/
/*--                                                                      --*/
/*--  KNPLCosts() =  vectors of the linear costs on the items             --*/
/*--  KNPQCosts() =  vectors of the quadratic costs on the items          --*/
/*--  KNPLBnds()  =  vectors of the lower bounds on the items             --*/
/*--  KNPUBnds()  =  vectors of the upper bounds on the items             --*/
/*--  KNPVlm()    =  value of the volume of the knapsack                  --*/
/*--                                                                      --*/
/*--------------------------------------------------------------------------*/
/*------------------------ "scalar" information ----------------------------*/
/*--------------------------------------------------------------------------*/

   virtual inline int KNPn( void )
   {
    return( n );
    }

/**< Returns the current number of items. The class has its protected field
   for storing this information, so this method is *not* pure virtual, but it
   is indeed virtual to allow re-definition if needed. */

/*--------------------------------------------------------------------------*/

   virtual void KNPLCosts( double *csts , const int *nms = 0 ,
			   int strt = 0 , int stp = Inf<int>() ) = 0;

/**< The linear costs of the items are written into csts[]. If nms == 0 then
   all the costs are written into csts[], otherwise cst[ i ] contains the
   informations relative to item nms[ i ] (nms must be Inf<int>()-terminated).

   The parameters `strt' and `stp' allow to restrict the output of the method
   to all and only the items `i' with strt <= i < min( KNPn() , stp ). `strt'
   and `stp' work in "&&" with nms; that is, if nms != 0 then only the values
   corresponding to items which are *both* in nms and whose index is in the
   correct range are returned. */

/*--------------------------------------------------------------------------*/

   virtual void KNPQCosts( double *csts , const int *nms = 0 ,
			   int strt = 0 , int stp = Inf<int>() ) = 0;

/**< The quadratic costs of the items are written into csts[]. If nms == 0
   then all the costs are written into csts[], otherwise csts[ i ] contains
   the informations relative to item nms[ i ] (nms must be
   Inf<int>()-terminated).

   The parameters `strt' and `stp' allow to restrict the output of the method
   to all and only the items `i' with strt <= i < min( KNPn() , stp ). `strt'
   and `stp' work in "&&" with nms; that is, if nms != 0 then only the values
   corresponding to items which are *both* in nms and whose index is in the
   correct range are returned (see above). */

/*--------------------------------------------------------------------------*/

   virtual double KNPLCost( const int i ) = 0;

/**< Return the linear cost of the i-th item (i = 0 .. n - 1). */

/*--------------------------------------------------------------------------*/

   virtual double KNPQCost( const int i ) = 0;

/**< Return the quadratic cost of the i-th item (i = 0 .. n - 1). */

/*--------------------------------------------------------------------------*/

   virtual void KNPLBnds( double *bnds , const int *nms = 0 ,
			  int strt = 0 , int stp = Inf<int>() ) = 0;

/**< The lower bounds of the items are written into bnds[]. If nms == 0 then
   all the bounds are written, otherwise bnds[ i ] contains the information
   relative to item nms[ i ] (nms must be Inf<int>()-terminated).

   The parameters `strt' and `stp' allow to restrict the output of the method
   to all and only the items `i' with strt <= i < min( KNPn() , stp ). `strt'
   and `stp' work in "&&" with nms; that is, if nms != 0 then only the values
   corresponding to items which are *both* in nms and whose index is in the
   correct range are returned. */

/*--------------------------------------------------------------------------*/

   virtual void KNPUBnds( double *bnds , const int *nms = 0 ,
			  int strt = 0 , int stp = Inf<int>() ) = 0;

/**< The upper bounds of the items are written into bnds[]. If nms == 0 then
   all the bounds are written, otherwise bnds[ i ] contains the information
   relative to item nms[ i ] (nms must be Inf<int>()-terminated).

   The parameters `strt' and `stp' allow to restrict the output of the method
   to all and only the items `i' with strt <= i < min( KNPn() , stp ). `strt'
   and `stp' work in "&&" with nms; that is, if nms != 0 then only the values
   corresponding to items which are *both* in nms and whose index is in the
   correct range are returned (see above). */

/*--------------------------------------------------------------------------*/

   virtual double KNPLBnd( const int i ) = 0;

/**< Return the lower bound of the i-th item (i = 0 .. n - 1). */

/*--------------------------------------------------------------------------*/

   virtual double KNPUBnd( const int i ) = 0;

/**< Return the upper bound of the i-th item (i = 0 .. n - 1). */

/*--------------------------------------------------------------------------*/

   virtual double KNPVlm( void ) = 0;

/**< Returns the volume `V' of the knapsack. */

/*--------------------------------------------------------------------------*/
/*-------------- METHODS FOR WRITING THE DATA OF THE PROBLEM ---------------*/
/*--------------------------------------------------------------------------*/

   virtual inline void WriteInstance( std::ostream &oFile ,
				      const int precc = 16 ,
				      const int precv = 16 );

/**< Write the instance to the provided ostream in the "complete" format read
   by ReadInstance() [see above]. The parameter "precc" and "precv" allow to
   set the precision (number of decimal digits) of the costs (linear and
   quadratic) and the bounds (upper and lower) and volume, respectively, when
   printed to the ostream. Rhe default value is "all digits of a double",
   which results in pretty large files, but smaller precisions are possible if
   one e.g. knows that the instance data is integer. Since this can typically
   be different for costs-related and volume-related information, two
   different parameters are provided.

   While virtual, the method is implemented in the base class: the instance
   data is read from the object using the class-provided implmentation of
   KNPLCost(), KNPQCost() and so on. Thus, the method will work for any
   derived class with no mdification; however, being virtual it can be
   re-implemented for maximum efficiency if desired. */

/*--------------------------------------------------------------------------*/
/*------------- METHODS FOR CHANGING THE DATA OF THE PROBLEM ---------------*/
/*--------------------------------------------------------------------------*/
/*--                                                                      --*/
/*-- Changing the costs, weights and volume of the Knapsack.              --*/
/*--                                                                      --*/
/*--  ChgLCost[s]() change the linear cost(s)                             --*/
/*--  ChgQCost[s]() change the quadratic cost(s)                          --*/
/*--  ChgLBnd[s]()  change the lower bound(s)                             --*/
/*--  ChgUBnd[s]()  change the upper bound(s)                             --*/
/*--  ChgVlm()      change the volume                                     --*/
/*--                                                                      --*/
/*-- The two forms of the methods allow to change either all/a given      --*/
/*-- subset of the entries of the corresponding vector, or to change one  --*/
/*-- given entry.                                                         --*/
/*--                                                                      --*/
/*--------------------------------------------------------------------------*/

   virtual void ChgLCosts( const double *csts , const int *nms = 0 ,
			   int strt = 0 , int stp = Inf<int>() ) = 0;

/**< Change the linear cost coefficients that are:

   - listed in into the vector of indices nms (ordered in increasing sense
     and Inf<int>()-terminated),

   - *and* whose name belongs to the interval [ strt , min( stp , KNPn() ) ).

   That is, if strt <= nms[ i ] < stp, then the coefficient of the nms[ i ]-th
   item will be changed reading from csts[ i ]. If nms == 0 (as the default),
   all the entries in the given range will be changed. */

/*--------------------------------------------------------------------------*/

   virtual void ChgQCosts( const double *csts , const int *nms = 0 ,
			   int strt = 0 , int stp = Inf<int>() ) = 0;

/**< Change the quadratic cost coefficients that are:

   - listed in into the vector of indices nms (ordered in increasing sense
     and Inf<int>()-terminated),

   - *and* whose name belongs to the interval [ strt , min( stp , KNPn() ) ).

   That is, if strt <= nms[ i ] < stp, then the coefficient of the nms[ i ]-th
   item will be changed reading from csts[ i ]. If nms == 0 (as the default),
   all the entries in the given range will be changed (see above). */

/*--------------------------------------------------------------------------*/

   virtual void ChgLCost( int i , const double cst ) = 0;

/**< Change the linear cost coefficient of item i to cst. */

/*--------------------------------------------------------------------------*/

   virtual void ChgQCost( int i , const double cst ) = 0;

/**< Change the quadratic cost coefficient of item i to cst. */

/*--------------------------------------------------------------------------*/

   virtual void ChgLBnds( const double *bnds , const int *nms = 0 ,
			  int strt = 0 , int stp = Inf<int>() ) = 0;

/**< Change the item lower bounds that are:

   - listed in into the vector of indices nms (ordered in increasing sense
     and Inf<int>()-terminated),

   - *and* whose name belongs to the interval [ strt , min( stp , KNPn() ) ).

   That is, if strt <= nms[ i ] < stp, then the bounds of the nms[ i ]-th
   item will be changed reading from bnds[ i ]. If nms == 0 (as the default),
   all the entries in the given range will be changed. */

/*--------------------------------------------------------------------------*/

   virtual void ChgUBnds( const double *bnds , const int *nms = 0 ,
			  int strt = 0 , int stp = Inf<int>() ) = 0;

/**< Change the item upper bounds that are:

   - listed in into the vector of indices nms (ordered in increasing sense
     and Inf<int>()-terminated),

   - *and* whose name belongs to the interval [ strt , min( stp , KNPn() ) ).

   That is, if strt <= nms[ i ] < stp, then the bounds of the nms[ i ]-th
   item will be changed reading from bnds[ i ]. If nms == 0 (as the default),
   all the entries in the given range will be changed. */

/*--------------------------------------------------------------------------*/

   virtual void ChgLBnd( int i , const double bnd ) = 0;

/**< Change the lower bound of item i to bnd. */

/*--------------------------------------------------------------------------*/

   virtual void ChgUBnd( int i , const double bnd ) = 0;

/**< Change the upper bound of item i to bnd. */

/*--------------------------------------------------------------------------*/

   virtual void ChgVlm( const double NVlm ) = 0;

/**< Change the volume. */

/*--------------------------------------------------------------------------*/
/*------------------------------ DESTRUCTOR --------------------------------*/
/*--------------------------------------------------------------------------*/

   virtual ~CQKnPClass() { /* nothing to do for the base class */ }

/*--------------------------------------------------------------------------*/
/*-------------------- PROTECTED PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/
/*--                                                                      --*/
/*--  The following fields are believed to be general enough to make it   --*/
/*--  worth adding them to the abstract base class.                       --*/
/*--                                                                      --*/
/*--------------------------------------------------------------------------*/

 protected:

/*--------------------------------------------------------------------------*/
/*---------------------- PROTECTED DATA STRUCTURES  ------------------------*/
/*--------------------------------------------------------------------------*/

 int n;           ///< total number of items
 int status;      /**< status of the algorithm: it is an int so that derived
		     classes can use it to encode other information apart
		     from the return value of SolveKNP(). */

 #if CQKnPClass_LOG
  std::ostream *KNPLog;   ///< the output stream object for log purposes
  char KNPLLvl;           ///< the "level of verbosity" of the log
 #endif

 };   // end( class KNPClass )

/*--------------------------------------------------------------------------*/
/*------------------- inline methods implementation ------------------------*/
/*--------------------------------------------------------------------------*/

inline void CQKnPClass::ReadInstance( std::istream &inFile , bool RBV )
{
 // reading all data- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 unsigned int length;
 inFile >> length;

 if( length <= 1 )
  throw( CQKnPClass::CQKException( "CQKnPClass::ReadInstance: wrong length" )
	 );

 // allocating memory- - - - - - - - -  - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 double *cC = new double[ length ];  // allocate linear costs

 for( unsigned int i = 0 ; i < length ; i++ ) {
  inFile >> cC[ i ];
  if( ( cC[ i ] <= - Inf<double>() ) || ( cC[ i ] >= Inf<double>() ) )
   throw( CQKnPClass::CQKException(
		       "CQKnPClass::ReadInstance: invalid linear cost" ) );
  }

 double *cD = new double[ length ];  // allocate quadratic costs

 for( unsigned int i = 0 ; i < length ; i++ ) {
  inFile >> cD[ i ];
  if( ( cD[ i ] < 0 ) || ( cC[ i ] >= Inf<double>() ) )
   throw( CQKnPClass::CQKException(
                      "CQKnPClass::ReadInstance: Invalid quadratic cost" ) );
  }

 // upper and lower bounds- - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 double Vol = 1;
 double *bL = new double[ length ];  // allocate lower bounds
 double *bU = new double[ length ];  // allocate upper bounds

 if( RBV ) {
  for( unsigned int i = 0 ; i < length ; i++ )
   inFile >> bL[ i ];

  for( unsigned int i = 0 ; i < length ; i++ ) {
   inFile >> bU[ i ];
   if( bL[ i ] > bU[ i ] )
    throw( CQKnPClass::CQKException(
                          "CQKnPClass::ReadInstance: Invalid bounds" ) );
   }

  inFile >> Vol;
  }
 else
  for( unsigned int i = 0 ; i < length ; i++ ) {
   bL[ i ] = 0;
   bU[ i ] = + Inf<double>();
   }

 LoadSet( length , cC , cD , bL , bU , Vol , true );

 delete[] bU;
 delete[] bL;
 delete[] cD;
 delete[] cC;

 }  // end( ReadInstance )

/*--------------------------------------------------------------------------*/

inline void CQKnPClass::WriteInstance( std::ostream &oFile ,
				       const int precc , const int precv )
{
 oFile << n << std::endl;

 oFile.precision( precc );

 for( unsigned int i = 0 ; i < n ; i++ )
  oFile << KNPLCost( i ) << "\t";
 oFile << std::endl;

 for( unsigned int i = 0 ; i < n ; i++ )
  oFile << KNPQCost( i ) << "\t";
 oFile << std::endl;

 oFile.precision( precv );

 for( unsigned int i = 0 ; i < n ; i++ )
  oFile << KNPLBnd( i ) << "\t";
 oFile << std::endl;

 for( unsigned int i = 0 ; i < n ; i++ )
  oFile << KNPUBnd( i ) << "\t";
 oFile << std::endl;

 oFile << KNPVlm() << std::endl;

 }  // end( WriteInstance )

/*--------------------------------------------------------------------------*/

 };  // end( namespace CQKnPClass_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif  /* CQKnPClass.h included */

/*--------------------------------------------------------------------------*/
/*------------------------- End File CQKnPClass.h --------------------------*/
/*--------------------------------------------------------------------------*/

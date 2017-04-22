/*--------------------------------------------------------------------------*/
/*--------------------------- File MMCFCple.h ------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Multicommodity Min Cost Flow (MMCF) Problems solver, based on calls calls
 * to the Cplex Callable Libraries for solution of generic LPs.
 *
 * \version 3.03
 *
 * \date 20 - 05 - 2012
 *
 * \author Antonio Frangioni \n
 *         Operations Research Group \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Paola Cappanera \n
 *         Operations Research Group \n
 *         Dipartimento di Sistemi e Informatica \n
 *         Universita' di Firenze \n
 *
 * Copyright(C) 1996 - 2012 Antonio Frangioni
 */

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __MMCFCplex
 #define __MMCFCplex

/*--------------------------------------------------------------------------*/
/*----------------------------- INCLUDES -----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "Graph.h"

#include "MMCFClas.h"

#include "cplex.h"

/*--------------------------------------------------------------------------*/
/*------------------------------- MACROS -----------------------------------*/
/*--------------------------------------------------------------------------*/

#define LOG_MMCF 0
/* If LOG_MMCF > 0, the MMCFCple class produces a log of its activities on the
   ostream object and at the "level of verbosity" set with the method
   SetMMCFLog() [see below]. */

/*--------------------------------------------------------------------------*/
/*------------------------ NAMESPACE and USINGS ----------------------------*/
/*--------------------------------------------------------------------------*/

#if( OPT_USE_NAMESPACES )
namespace MMCFClass_di_unipi_it
{
 using namespace MMCFGraph_di_unipi_it;
#endif

/*--------------------------------------------------------------------------*/
/*--------------------------- CLASS MMCFCplex ------------------------------*/
/*--------------------------------------------------------------------------*/
/*---------------------------- GENERAL NOTES -------------------------------*/
/*--------------------------------------------------------------------------*/
/** The class MMCFCplex implements the generic interface for Multicommodity
 *  Min-Cost Flow solvers defined by the class MMCFClass [see MMCFClas.h].
 *  It uses the service class Graph [see Graph.h] for reading the data of the
 *  instance to be solved.
 *
 *  MMCFCplex uses the Cplex Callable Library (version 5.0 or higher,
 *  currently version 12.3) for solving the problem. It therefore accepts
 *  "extra" variobles and "extra" (linear) constraints, if any, in the
 *  instance.
 *
 *  WARNING: this class requires FNumber = MFNumber = CNumber = double and
 *  perhaps Index = int for working properly. */

/*--------------------------------------------------------------------------*/

class MMCFCplex : public virtual MMCFClass
{

/*--------------------------------------------------------------------------*/
/*------------------- PUBLIC PART OF THE CLASS  ----------------------------*/
/*--------------------------------------------------------------------------*/

public:

/*--------------------------------------------------------------------------*/
/*---------------------------- PUBLIC TYPES --------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Public types
   @{ */

   enum MMCFCplexA{ kPPrim , kPDual , kBarrier , kBrCrss , kMIP ,
			    kNPrim , kNDual };

/*--------------------------------------------------------------------------*/
/*-------------------------  END Class CplexState  -------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*---------------------------- PUBLIC METHODS ------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*----------------------------- CONSTRUCTOR --------------------------------*/
/*--------------------------------------------------------------------------*/

   MMCFCplex( Graph* Gh , istream *iStrm = NULL , CPXENVptr extenv = 0 );

/* Constructor of the class.

   Takes the data of the instance from an object of class Graph.

   If extenv != 0, it is taken as a pointer to a valid Cplex environment
   [see GetCplexEnv() below] that will be used for all the lifetime of the
   object. Otherwise, a Cplex environment is possibly initialized in the
   constructor. The environment is shared among all the active instances of
   MmCFCplex objects; this is done in order to save on the number of (costly)
   Cplex licenses required to have multiple MMCF solvers active at the same
   time, since any environment consumes a licence. Thus, the environment is
   actually initialized only when the first instance is constructed, and it
   is released when the last destructor (of an instance using it) is invoked.

   The possibility of passing an external environment is let because having
   all the instances to share the same environment has the drawback that
   all changes made to optimization parameters of the environment [see
   SetCplexParam() below] are "seen" be all active instances, even though the
   changes are invoked for one specific object. If for some reason this is
   unacceptable, the user should provide each "sensitive" instance with its
   own private environment, letting all the others to share the "static" one.
   Of course, it is then user's responsibility to initialize and free the
   environment.

   The parameter `iStrm', if provided, is taken as a pointer to a istream
   from which the algorithmic parameters for the MMCF
   algorithm are sequentially read in the following order. Each parameter
   must be placed at the beginning of a separate line, max 255 carachters
   long, with all the rest of the line up to the first newline carachter
   '\n' (apart from a separating whitespace) being available for comments.
   Any line whose first carachter is '#' and any blank line is ignored.
   If 0 is passed, the file ends before reaching a given parameter, or
   some parameter is in the wrong format, each non-specified parameter is
   given a default value, shown in [] below.

     char LP algorithm ['b']      the LP solver: a = primal , b = dual , c = barrier
                                  d = barrier + crossover , e = primal network,
                                  f = dual network

     bool Populate [false]        true if it needs to populate the solution pool

     bool Lazy [false]            true if a group of constraints can be kept smaller
                                  when these constraints are not included

     double time limit [1000]     the time limit for solving the problem

     int CPX_PARAM_PREIND [1]          The pre-processing parameters by Cplex
     int CPX_PARAM_REPEATPRESOLVE [-1]
     int CPX_PARAM_PRESLVND [0]
     int CPX_PARAM_REDUCE [3]

     The pointer to Graph object is used to update the formulation of the
     problem when the Lazy constraints have to be excluded.     */

/*--------------------------------------------------------------------------*/
/*---------------------- OTHER INITIALIZATIONS -----------------------------*/
/*--------------------------------------------------------------------------*/

   void SetMMCFLog( ostream *outs = 0 , const char lvl = 0 );

/**< The class ouputs "log" information onto the ostream pointed by outs.
   lvl controls the "level of verbosity" of the code: lvl == 0 means that
   nothing at all is printed, and values larger than 0 mean increasing
   amounts of information, the specific effect of each value being derived-
   class-dependent. outs == 0 implies lvl == 0. */

/*--------------------------------------------------------------------------*/

   inline void SetOptEps( const double OE = 0 );

/*--------------------------------------------------------------------------*/

   inline void SetCplexA( MMCFCplexA nCA );

/* Allow the user to choose which kind of algorithm has to be used for solving
   the instance. Available algorithms are:

    kPPrim    Primal Simplex;
    kPDual    Dual Simplex;
    kBarrier  "pure" Barrier (Primal-Dual Interior Point) algorithm;
    kBrCrss   Barrier algorithm followed by crossover to a vertex solution;
    kMIP      Branch & Bound (for Integer and Mixed-Integer problems);
    kNPrim    Network crash-start followed by Primal Simplex;
    kNDual    Network crash-start followed by Dual Simplex.

   The default is kNDual. */

/*--------------------------------------------------------------------------*/

   inline void SetCplexParam( int whichparam , const char * value );

   inline void SetCplexParam( int whichparam , int value );

   inline void SetCplexParam( int whichparam , double value );

   inline CPXENVptr GetCplexEnv( void );

/* The first two methods allow to set the many algorithmic parameters of
   Cplex; see the documentation of CPXsetintparam() and CPXsetdblparam() in
   the Cplex manual for details.

   An alternative is to get the pointer to the internal Cplex environment with
   GetCplexEnv() and then call CPXset***param() directly. This also allows to
   perform any other operation with the environment, such as reading the value
   of the parameters with CPXgetintparam() and CPXgetdblparam(), so care must
   be taken.

   The returned pointer is the same passed to the constructor [see above], if
   any; otherwise it is the "static" environment shared by all the active
   MCFCplex instances. In the latter case, any change in the environment
   simultaneously affect *all* the existing (and future) MCFCplex instances
   which have not been given a "private" environment. */

/*--------------------------------------------------------------------------*/
/*-------------------- METHODS FOR SOLVING THE PROBLEM ---------------------*/
/*--------------------------------------------------------------------------*/

   MMCFStatus SolveMMCF( void );

/*--------------------------------------------------------------------------*/
/*------------------- METHODS FOR READING RESULTS  -------------------------*/
/*--------------------------------------------------------------------------*/

   FONumber GetPVal( void );

   FONumber GetDVal( void );

/*---------------------------------------------------------------------------*/

   FONumber GetUpprBnd( bool &HvSol );

/**< Return true if the primal solution is integer. */

/*---------------------------------------------------------------------------*/

   bool GetPSol( void );

   bool GetDSol( void );

   FONumber CostOf( void );

/*---------------------------------------------------------------------------*/
/*---------- METHODS FOR READING THE DATA OF THE PROBLEM  -------------------*/
/*---------------------------------------------------------------------------*/

   void GetCosts( CRow Csts, cIndex_Set nms = 0 , cIndex strt = 0 ,
		      Index stp = Inf<Index>() );

   void GetXtrCsts( CRow XtrCs, cIndex_Set nms = 0 , cIndex strt = 0 ,
              Index stp = Inf<Index>() );

   void GetICaps( FRow ICps, cIndex_Set nms = 0 , cIndex strt = 0 ,
              Index stp = Inf<Index>() );

   void GetMCaps( FRow MCps, cIndex_Set nms = 0 ,
		  cIndex strt = 0 , Index stp = Inf<Index>() );

/*--------------------------------------------------------------------------*/

   inline Index NrArcsK( cIndex k );

/* Returns the number of defined arcs of commodity k. */

/*--------------------------------------------------------------------------*/

   inline Index NrNodesK( cIndex k );

/* Returns the number of defined nodes of commodity k. */

/*--------------------------------------------------------------------------*/

   inline Index TotArcs( void );

/* Returns the total number of defined arcs relative to all commodities; for
   undirected graphs, TotArcs() is the total number of *undirected* arcs. */

/*--------------------------------------------------------------------------*/

   inline Index TotNodes( void );

/* Returns the total number of defined nodes relative to all commodities. */

/*--------------------------------------------------------------------------*/

   inline Index ArcPosKJ( cIndex k , cIndex j );

/* Returns the ordinal number of arc j for commodity k in the sequence of
   existent arcs (which is also the column index of the arc in the Cplex
   matrix); returns Inf<Index>() if arc j is not defined for commodity k. */

/*--------------------------------------------------------------------------*/

   inline Index NodePosKJ( cIndex k , cIndex j );

/* Returns the ordinal number of node j commodity k in the sequence of
   existent nodes (which is also the row index of the node in the Cplex
   matrix); returns Inf<Index>() if node j is not defined for commodity k. */

/*--------------------------------------------------------------------------*/

   inline cIndex_Set GetActvArcs( void );

   inline Index NumActvArcs( void );

/* ActvArcs() returns a pointer to a vector (ordered in increasing sense and
   Inf<Index>()-terminated) containing the indices of the arcs which have an
   associated mutual capacity constraint; it returns 0 if every arc has
   an associated mutual capacity constraint.

   NumActvArcs() returns the number of the arcs which have an associated
   mutual capacity constraint. */

/*--------------------------------------------------------------------------*/

   inline bool Directed( void );

/* Returns true if the graph is directed */

/*--------------------------------------------------------------------------*/
/*------------- METHODS FOR ADDING / REMOVING / CHANGING DATA --------------*/
/*--------------------------------------------------------------------------*/

   void ChgCosts( cCRow NwCsts, cIndex_Set nms = 0 ,
		  cIndex strt = 0 , Index stp = Inf<Index>() );

   void ChgXtrCsts( cCRow nwXtrCs, cIndex_Set nms = 0 ,
		    cIndex strt = 0 , Index stp = Inf<Index>() );

   void ChgXtrBnds( cRow XLr = 0 , cRow XUr = 0 , cIndex_Set nms = 0 ,
        cIndex strt = 0 , Index stp = Inf<Index>() );

   void ChgICaps( cFRow NwICps, cIndex_Set nms = 0 ,
		  cIndex strt = 0 , Index stp = Inf<Index>() );

   void ChgMCaps( cFRow NwMCps, cIndex_Set nms = 0 ,
		  cIndex strt = 0 , Index stp = Inf<Index>() );

/*--------------------------------------------------------------------------*/

   void ChgIntVar( cIndex k = Inf<Index>() , bool IntVld = true ,
		   cIndex_Set nms = 0 ,
		   Index strt = 0 , Index stp = Inf<Index>() );

/*--------------------------------------------------------------------------*/

#if( CHGARCS_MMCF )

   void CloseArcs( cIndex_Set whch );

   void OpenArcs( cIndex_Set whch );

 #if( CHGARCS_MMCF > 1 )

   void CloseArcs( cIndex k , cIndex_Set whch );

   void OpenArcs( cIndex k , cIndex_Set whch );

 #endif
#endif

/*--------------------------------------------------------------------------*/

   void AddExtraVars( cIndex NXV , cCRow XCst = 0 , cRow XUb = 0 , cRow XLb = 0 ,
		   bool IntVar = false , cIndex_Set WIntVar = 0 );

   void AddExtraConstr( cIndex NXC , int *IBeg , int *Indx , double *Vals ,
			const double *XLr = 0 , const double *XUr = 0 );

/*--------------------------------------------------------------------------*/

   Index GetNode( void );
/**< Returns the number of nodes processed so far in the active branch-and-cut search.*/


/*--------------------------------------------------------------------------*/
/*---------------------------- DESTRUCTOR ----------------------------------*/
/*--------------------------------------------------------------------------*/

   ~MMCFCplex();


 protected:

/*--------------------------------------------------------------------------*/
/*-------------------------- PROTECTED METHODS -----------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*----------------------- PROTECTED DATA STRUCTURES  -----------------------*/
/*--------------------------------------------------------------------------*/
   CPXENVptr    env;        // Cplex environment
   CPXLPptr     lp;         // Cplex LP pointer

   CPXFILEptr   logfile;

   Index_Set    StrtCols;   // 1st column for commodity k
   bool         Drctd;      // true if the problem is directed

   bool         Populate;   // true if cplex generates multiple solutions
                            // whenever deal with MIP problem
   bool         Lazy;

   int preind;              // pre-processing parameters
   int represolve;
   int prenode;
   int reduce;

/*--------------------------------------------------------------------------*/
/*----------------- PRIVATE PART OF THE CLASS ------------------------------*/
/*--------------------------------------------------------------------------*/

 private:

/*--------------------------------------------------------------------------*/
/*-------------------------- PRIVATE TYPES ---------------------------------*/
/*--------------------------------------------------------------------------*/

 typedef Index_Set      *Index_Mat;       ///< set of set of indices

/*--------------------------------------------------------------------------*/
/*------------------------- PRIVATE METHODS --------------------------------*/
/*--------------------------------------------------------------------------*/

   template<class T>
   inline void TranslateK( T *const vect , const cIndex_Set Dict ,
			   Index n , cIndex strt , const T Dflt = 0 );

/* Takes the less-than-n-vector vect[] in the internal representation 
   specified by the (row of the) dictionary Dict[] and translates it into the
   external n-vector representation; vect[] is relative to a commodity, of
   which strt is the starting index in the dictionary.
   Dflt gives the default value for elements not in the dictionary (absent in
   the internal representation). */

/*---------------------------------------------------------------------------*/

   template<class T>
   inline void TranslateF( T *const vect , Index_Mat Dict , cIndex_Set Strt ,
			   Index n , const T Dflt = 0 );

/* Takes the full solution vect[] in the internal representation specified by
   the dictionary Dict[] and Strt[] and translates it into the external
   n-vector representation.
   Dflt gives the default value for elements not in the dictionary (absent in
   the internal representation). */

/*---------------------------------------------------------------------------*/

   inline void GetExtraSol( void );

/* Implements a part of GetPSol() (guess which). */

/*--------------------------------------------------------------------------*/
/*----------------------- PRIVATE DATA STRUCTURE ---------------------------*/
/*--------------------------------------------------------------------------*/

  MMCFCplexA   CA;         // algorithm to be used

  bool         PFeas;      // true if a primal feasible solution is known
  bool         DFeas;      // true if a dual feasible solution is known

  Index_Mat    ADict;      // dictionary of arc names -> columns
  Index_Mat    NDict;      // dictionary of node names -> rows


  Index_Set    StrtRows;   // 1st row of node balance constraints for
                           // commodity k

  Index        NActvA;     // number of arcs with a mutual capacity
  Index_Set    ActvArcs;   // constraint, and their indices

  Index        NumSol;     // number of primal solutions in the pool
  Index        IDSol;      // current primal solution of the pool


  // static members - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  static CPXENVptr genv;   // "global" Cplex environment pointer
  static Index EnvICntr;   // counter of active instances using genv


/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

  };  // end( class MMCFCplex )

/*--------------------------------------------------------------------------*/
/*----------------- INLINE METHODS IMPLEMENTATION --------------------------*/
/*--------------------------------------------------------------------------*/


inline void MMCFCplex::SetOptEps( const double OE )
{
 MMCFClass::SetOptEps( OE );
 SetCplexParam( CPX_PARAM_EPGAP , OE );
 }

/*--------------------------------------------------------------------------*/

inline void MMCFCplex::SetCplexA( MMCFCplexA nCA )
{
 CA = nCA;
 }

inline void MMCFCplex::SetCplexParam( int whichparam , const char * value ) {
 CPXsetstrparam( env , whichparam , value );
 }

/*--------------------------------------------------------------------------*/

inline void MMCFCplex::SetCplexParam( int whichparam , int value )
{
 CPXsetintparam( env , whichparam , value );
 }

/*---------------------------------------------------------------------------*/

inline void MMCFCplex::SetCplexParam( int whichparam , double value )
{
 CPXsetdblparam( env , whichparam , value );
 }

/*--------------------------------------------------------------------------*/

inline CPXENVptr MMCFCplex::GetCplexEnv( void )
{
 return( env );
 }

/*--------------------------------------------------------------------------*/

inline MMCFClass::Index MMCFCplex::NrArcsK( cIndex k )
{
 return( StrtCols[ k + 1 ] - StrtCols[ k ] );
 }

/*--------------------------------------------------------------------------*/

inline MMCFClass::Index MMCFCplex::NrNodesK( cIndex k )
{
 return( StrtRows ? StrtRows[ k + 1 ] - StrtRows[ k ] : NNodes );
 }

/*--------------------------------------------------------------------------*/

inline MMCFClass::Index MMCFCplex::TotArcs( void )
{
 return( StrtCols[ NComm ] );
 }

/*--------------------------------------------------------------------------*/

inline MMCFClass::Index MMCFCplex::TotNodes( void )
{
 return( StrtRows ? StrtRows[ NComm ] : NNodes * NComm );
 }

/*--------------------------------------------------------------------------*/

inline MMCFClass::Index MMCFCplex::ArcPosKJ( cIndex k , cIndex j )
{
 if( ADict[ k ] )
  return( ADict[ k ][ j ] );
 else
  return( StrtCols[ k ] + j );
 }

/*--------------------------------------------------------------------------*/

inline MMCFClass::Index MMCFCplex::NodePosKJ( cIndex k , cIndex j )
{
 if( NDict[ k ] )
  return( NDict[ k ][ j ] );
 else
  if( StrtRows )
   return( StrtRows[ k ] + j );
  else
   return( NNodes * k + j );
 }

/*--------------------------------------------------------------------------*/

inline MMCFClass::cIndex_Set MMCFCplex::GetActvArcs( void )
{
 return( ActvArcs );
 }

/*--------------------------------------------------------------------------*/

inline MMCFClass::Index MMCFCplex::NumActvArcs( void )
{
 return( NActvA );
 } 

/*--------------------------------------------------------------------------*/

inline bool MMCFCplex::Directed( void )
{
 return( Drctd );
 }

/*--------------------------------------------------------------------------*/

#if( OPT_USE_NAMESPACES )
 };  // end( namespace MMCFClass_di_unipi_it )
#endif

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif  /* MMCFCple.h included */

/*--------------------------------------------------------------------------*/
/*-------------------- End File MMCFCple.h ---------------------------------*/
/*--------------------------------------------------------------------------*/

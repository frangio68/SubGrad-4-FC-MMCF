/*--------------------------------------------------------------------------*/
/*--------------------------- File MMCFClas.h ------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the abstract (pure virtual) base class MMCFClass, which
 * defines a standard interface for solvers of Multicommodity Min Cost Flow
 * problems, or, more in general, problems with the MMCF structure "plus or
 * minus" some variables and/or constraints. The class defines the methods
 * that allow the "user" of a "solver object" to solve an instance of MMCF
 * (...), read and change the data of the problem and read the results of the
 * optimization. Note that the way in which the data is initially read
 * depends on the specific solver used, and it is not fixed in the base
 * class, although an independent Graph class [see Graph.h] has also been
 * developed for the purpose.
 *
 * \version 4.02
 *
 * \date 21 - 11 - 2013
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
 * Copyright &copy 1996 - 2013 by Antonio Frangioni
 */
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __MMCFClass
 #define __MMCFClass  /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- MACROS ----------------------------------*/
/*--------------------------------------------------------------------------*/
/** @defgroup MMCFCLASS_MACROS Compile-time switches in MMCFClass.h
    These macros control some important details of the class interface.
    Although using macros for activating features of the interface is not
    very C++, switching off some unused features may allow some
    implementation to be more efficient in running time or memory.
    @{ */

/*---------------------------- CHGARCS_MMCF --------------------------------*/

#define CHGARCS_MMCF 1

/**< If CHGARCS_MMCF > 0, methods for "opening" and "closing" arcs in the
   MMCF are added to the interface of the class and must therefore be
   implemented by the derived classes [see [Open/Close]Arc() below.
   If CHGARCS_MMCF == 1, the arcs can be opened or closed only simultaneously
   for all the commodities, while if CHGARCS_MMCF > 1 then arcs can be opened
   and closed for each commodity individually. */

/* @}-----------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "OPTUtils.h"

/* OPTUtils.h defines standard interfaces for timing and random routines, as
   well as the namespace OPTtypes_di_unipi_it and the macro
   OPT_USE_NAMESPACES, useful for switching off all namespaces in one blow
   for those strange cases where they create problems. */

#include <iostream>

/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE and USING ----------------------------*/
/*--------------------------------------------------------------------------*/

#if( OPT_USE_NAMESPACES )
namespace MMCFClass_di_unipi_it
{
 /** @namespace MMCFClass_di_unipi_it
     The namespace MMCFClass_di_unipi_it is defined to hold the MMCFClass
     class and all the relative stuff. It comprises the namespace
     OPTtypes_di_unipi_it. */

 using namespace OPTtypes_di_unipi_it;
#endif

/*--------------------------------------------------------------------------*/
/*--------------------------- CLASS MMCFClass ------------------------------*/
/*--------------------------------------------------------------------------*/
/*---------------------------- GENERAL NOTES -------------------------------*/
/*--------------------------------------------------------------------------*/

/** This class defines a standard abstract interface for solvers of
  (generalized) Multicommodity Min-Cost Flow problems, i.e., problems of the
  form

   min Sum{k = 0 .. K - 1} C[ k ] * X[ k ] + C[ K ] * Y   s.t.        

   (1.k)   E * X[ k ] = b[ k ]        k = 0 .. K - 1

   (2.k)   0 <= X[ k ] <= U[ k ]      k = 0 .. K - 1

   (3)     Sum{k = 0 .. K - 1} X[ i ] <= U

   (4)     b[ K ] <= Sum{k = 0 .. K - 1} A[ k ] * X[ k ]
                     + A[ K ] * Y <= b[ K + 1 ]

   (5)     U[ K ] <= Y <= U[ K + 1 ]

  X[ k ] are the flow variables, one for each commodity; Y[] are the "extra"
  variables. (1.k) and (2.k) are the Flow Conservation and Upper Bound
  constraints for commodity k, respectively. E is the node-arc incidence
  matrix of an underlying graph G(N, A), which is common to all commodities;
  however, each commodity can in principle flow only on a subset A[ k ] of
  the arcs A, and therefore it is in fact defined only on a subgraph
  G[ k ](N[ k ], A[ k ]), where N[ k ] is the subset of N touched by the
  arcs in A[ k ]. (3) are the Mutual Capacity constraints, linking the
  (otherwise disjoint) different commodities; (4), that may be empty, are
  other "extra" constraints, which may be:
  - commodity-separable, i.e., for some constraint i there exists a
    commodity k such that A[ h ][ i ] == 0 for all h != k, but not
    network-type;
  - other linking constraints between differnet commodities;
  - constraints linking some (or all) of the commodities to the extra
    (non-flow) variables Y;
  - constraints involving only the extra (non-flow) variables Y.
  (5) are bound constraints on the extra variables.                    

 Actually, the cost function and/or the extra constraints may also be
 nonlinear and/or nonseparable; however, the interface does not set any
 standard for this, assuming (temporarly) that the nonlinearities can
 remain "hidden".                                                 

 Both (subsets of) the flow variables and (subsets of) the extra variables
 can be constrained to be integer-valued.

 Examples of problems in this class are:
 - the Linear Continuous Multicommodity Min Cost Flow Problem;        
 - relaxations of the above, such as the Flow Relaxation (obtained
   relaxing the Mutual Capacity Constraints) or the Knapsack Relaxation
   (obtained relaxing the Flow Conservation constraints);
 - extensions of the above, like the Integer MMCF or the Fixed-Charge MMCF.

 The class provides support for the notion that the MMCF-related problem
 at hand may be separable in a certain number of independent subproblems,
 which typically happens in the "relaxations" cases above.

 In all the comments below, m has to be understood as the number of arcs in
 the graph G (#A), n as the number of nodes in the graph G (#N) and K as the
 number of commodities. */

class MMCFClass
{

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*--------------------------------------------------------------------------*/
/*---------------------------- PUBLIC TYPES --------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Public types
    MMCFClass defines five main public types:

    - Index, the type of arc and node indices;

    - FNumber, the type of flow variables, arc capacities, and node deficits;

    - CNumber, the type of flow costs, node potentials, and arc reduced costs;

    - FONumber, the type of objective function value;

    - MFNumber, the type a Multicommodity flow variable;

    - Number, the type an "extra" variable.

    By re-defining the types in this section, one MMCFSolver could be adapted
    to work with "smaller" data types than the obvioous ones (double for all
    but the indices, int for the latter). This may be relevant e.g. for cases
    where flow types are integer but multicommodity flows are not. However,
    *it is the user's responsibility to ensure that these types are set to
    reasonable values*, since *not all solution algorithms will work with
    "restricted" data*. Hence, only the the experienced user may want to
    experiment with changing this, and only if memory footprint and/or speed
    is really a primary concern and it is likely that changing these will
    improve.
    @{ */

/*--------------------------------------------------------------------------*/

 typedef unsigned int    Index;          ///< index of a node or arc ( >= 0 )
 typedef Index          *Index_Set;      ///< set (array) of indices
 typedef const Index    cIndex;          ///< a read-only index
 typedef cIndex        *cIndex_Set;      ///< read-only index array

/*--------------------------------------------------------------------------*/

 typedef double          FNumber;        ///< type of arc flow
 typedef FNumber        *FRow;           ///< vector of flows
 typedef const FNumber  cFNumber;        ///< a read-only flow
 typedef cFNumber      *cFRow;           ///< read-only flow array

/*--------------------------------------------------------------------------*/

 typedef double          CNumber;        ///< type of arc flow cost
 typedef CNumber        *CRow;           ///< vector of costs
 typedef const CNumber  cCNumber;        ///< a read-only cost
 typedef cCNumber      *cCRow;           ///< read-only cost array

/*--------------------------------------------------------------------------*/

 typedef double          MFNumber;       ///< a Multicommodity flow variable
 typedef MFNumber       *MFRow;          ///< vector of Multicommodity flows
 typedef const MFNumber cMFNumber;       ///< a read-only Multicommodity flow
 typedef cMFNumber     *cMFRow;          ///< read-only Multicommodity array

/*--------------------------------------------------------------------------*/

 typedef double          FONumber; 
 /**< type of the objective function: has to hold sums of products of
    MFNumber(s) by CNumber(s) */

 typedef const FONumber cFONumber;       ///< a read-only o.f. value

/*--------------------------------------------------------------------------*/

 typedef double          Number;         ///< an "extra" variable
 typedef Number         *Row;            ///< an "extra" variable array
 typedef const Number   cNumber;         ///< a read-only Number
 typedef cNumber       *cRow;            ///< read-only Number array

/*--------------------------------------------------------------------------*/
/** Very small class to simplify extracting the "+ infinity" value for a
    basic type (FNumber, CNumber, Index); just use Inf<type>().

   template <typename T>
   class Inf {
    public:
     Inf() {}
     operator T() { return( std::numeric_limits<T>::max() ); }
    };
 */
/*--------------------------------------------------------------------------*/
/** Very small class to simplify extracting the "machine epsilon" for a
    basic type (FNumber, CNumber); just use Eps<type>().

   template <typename T>
   class Eps {
    public:
     Eps() {}
     operator T() { return( std::numeric_limits<T>::epsilon() ); }
    };
 */
/*--------------------------------------------------------------------------*/
/** Small class for exceptions. Derives from std::exception implementing the
   virtual method what() -- and since what is virtual, remember to always
   catch it by reference (catch exception &e) if you want the thing to work.
   MMCFException class are thought to be of the "fatal" type, i.e., problems
   for which no solutions exists apart from aborting the program. Other kinds
   of exceptions can be properly handled by defining derived classes with
   more information. */

 class MMCFException : public exception {
 public:
  MMCFException( const char *const msg = 0 ) { errmsg = msg; }

  const char* what( void ) const throw () { return( errmsg ); }

 private:
  const char *errmsg;
  };

/*--------------------------------------------------------------------------*/
/** Public enum describing the possible status of the MMCF solver. */

  enum MMCFStatus { kOK = 0 ,        ///< optimal solution found
		    kStopped ,       ///< optimization stopped
		    kUnfeasible ,    ///< problem is unfeasible
		    kUnbounded ,     ///< problem is unbounded
		    kError           ///< error in the solver
                    };

/*@} -----------------------------------------------------------------------*/
/*--------------------------- PUBLIC METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/
/*---------------------------- CONSTRUCTOR ---------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructors
    @{ */

   MMCFClass( void )

/**< Constructor of the base class: gives some default values to the few data
   structure of the base class, that can be changed later in the constructors
   of the derived classes.

   Note that this is a default constructor (with no arguments), which
   therefore overrides the - surely wrong - default copy constructor. Also,
   by having no arguments one does not need to care to initialize the base
   class in case of multiple inheritance with MMCFClass as a virtual base
   class. */
   {
    NComm = NNodes = NArcs = XtrVrs = XtrCnst = 0;
    NSubPr = 1;

    OptEps = FsbEps = 0;

    FSol = 0;
    MFSol = 0;
    XSol = 0;
    FBse = XBse = 0;
    WhchFS = WhchSP = Inf<Index>();

    NPot = RCst = MCst = 0;
    WhchNP = WhchRC = Inf<Index>();
    XtRC = XtDV = 0;

    MMCFLog = 0;
    MMCFLLvl = 0;

    MMCFt = 0;

    }  // end( MMCFClass )

/* @}-----------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
    @{ */

   virtual void SetMMCFLog( ostream *outs = 0 , const char lvl = 0 )

/**< The class ouputs "log" information onto the ostream pointed by outs.
   lvl controls the "level of verbosity" of the code: lvl == 0 means that
   nothing at all is printed, and values larger than 0 mean increasing
   amounts of information, the specific effect of each value being derived-
   class-dependent. outs == 0 implies lvl == 0. */
   {
    if( ( MMCFLog = outs ) )
     MMCFLLvl = lvl;
    else
     MMCFLLvl = 0;
    }

/*--------------------------------------------------------------------------*/

   void SetMMCFTime( bool TimeIt = true )

/**< SetMMCFTime() allocates an OPTtimers object [see OPTUtils.h] that should
   be used for timing the calls to relevant methods of the class. The time can
   be read with MMCFTime() [see below]. By default, or if SetMMCFTime( false )
   is called, no timing is done. Note that, since all the relevant methods ot
   the class are pure virtual, MMCFClass can only manage the OPTtimers object,
   but it is due to derived classes to actually implement the timing.

   Note that time accumulates over the calls: calling SetMMCFTime(), however,
   resets the counters, allowing to time specific groups of calls. */
   {
    if( TimeIt )
     if( MMCFt )
      MMCFt->ReSet();
     else
      MMCFt = new OPTtimers();
    else
     delete MMCFt;
    }

/*--------------------------------------------------------------------------*/

   virtual void SetOptEps( const double OE = 0 )
   {
    OptEps = OE;
    }

   virtual void SetFsbEps( const double FE = 0 )
   {
    FsbEps = FE;
    }

/**< In many cases, only an "approximate" solution of the problem is possible;
   alternatively, only an "approximate" solution may be required for the
   purposes of the caller (in order to save time).
   The exact meaning of "approximate" is solver-dependent, but the more
   common ways in which this happens are

   - either the value of the solution is not exactly optimal;

   - or the constraints are not exactly satisfied.

   SetOptEps() tells that any solution that is OE-optimal w.r.t. the value of
   the objective function can be considered optimal. SetFsbEps() tells that
   any solution where the violation of the constraints is not larger than FE
   can be considered feasible.

   The exact meaning of OE and FE is not strictly imposed in the interface,
   but the idea is that they are relative accuracies (that is, OE == 1e-6
   should roughly mean that the objectve function value must agree with the
   optimal value for at least about 6 decimal digits ...). */

/* @}-----------------------------------------------------------------------*/
/*-------------------- METHODS FOR SOLVING THE PROBLEM ---------------------*/
/*--------------------------------------------------------------------------*/
/** @name Solving the problem
    @{ */

   virtual void SetSubP( cIndex ws = Inf<Index>() )
   {
    WhchSP = ws;
    }

   virtual MMCFStatus SolveMMCF( void ) = 0;

/**< SolveMMCF() try to solve the current instance of the MMCF-related problem.
   Note that the meaning of "solution" may vary depending on the solver used.
   Some issues, such as the approximations in the value of the objective
   function and the satisfacion of constraints, are (partly) dealt with in
   the Set***Eps() methods [see above].

   The MMCF-related problem may be decomposable into NrSubP() [see below]
   independent SubProblems; furthermore, it may have some "constant term"
   (such as those given by RHS-related terms in a Lagrangean relaxation),
   that is, a subproblem with constant solution. Each of these subproblems
   can be solved, and its results be queried, independently from the others.
   SetSubP( ws ) with tells that all the calls to methods that do something
   subproblem-specific refere to:

   - the 0-th "constant" subproblem if ws == 0;

   - the ws-th subproblem if 1 <= ws <= NrSubP();

   - all the subproblems simultaneously *but* the 0-th "constant" one if
     ws == NrSubP() + 1;

   - all the subproblems simultaneously, *comprised* the 0-th "constant" one,
     if ws > NrSubP() + 1 (e.g., ws == Inf<Index>()).

   If SetSubP() is *never* called, ws == Inf<Index>() is assumed.

   For SolveMMCF(), for instance, 1 <= ws <= NrSubP() only solves the ws-th
   subproblem, ws == 0 only calculates the "constant", while any ws > NrSubP()
   solves all the subproblems, calculating also the constant unless ws is
   exactly NrSubP() + 1.

   SolveMMCF() returns the following codes:

    kOK           if optimization has been carried out succesfully, whatever
                  this means for the solver;

    kStopped      if optimization have been stopped before that the stopping
                  conditions of the solver applied, e.g. because of the
		  maximum allowed number of "iterations" have been reached;
		  this is not necessarily an error, as it might just be
		  required to re-call SolveMMCF() giving it more "resources"
		  in order to solve the problem;

    kUnfeasible   if the current MMCF instance is (primal) unfeasible;

    kUnbounded    if the current MMCF instance is (primal) unbounded;

    kError        if there was an error during the optimization; this
                  typically indicates that computation cannot be resumed,
		  although solver-dependent ways of dealing with
		  solver-dependent errors may exhist. */

/* @}-----------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/
/** @name Reading the solution solution

  In principle, whatever problem is solved by the actual class which
  implements the interface can have two "kinds" of solutions, i.e., primal
  and dual ones. If the problem is an LP (or, more in general, a problem with
  no duality gap), then the corresponding primal and dual objective function
  values will be the same, but in general this may not happen.

  The solvers may be able to provide multiple (almost-)optimal primal and/or
  dual solutions.

  If the problem is decomposable, then the set of (primal and dual) variables
  can be partitioned into disjoint subsets. Hence, each (primal and dual)
  solution can be seen as the composition of solutions for the subproblems,
  and the user can ask for each sub-solution separatedly. However, no
  assumption here is done on how the flow and extra variables are partitioned,
  so each solution for one subproblem is treated as a "full" primal/dual
  solution (although it will be a "very sparse" one). Obviously, querying
  solution information for the ws-th subproblem (after a call to
  SetSubP( ws )) can be done *only* after that SolveMMCF() has been called
  with that setting in effect, or with SetSubP( Inf<Index>() ) in effect.

  Often, as a by-product of the solution of a MMCF-related problem, an
  Upper/Lower Bound on the value of the solution can be obtained. Typically,
  any primal/dual feasible solution provides such a bound: however, note that
  the idea here is that "feasible" means "for the whole problem", as opposed
  to the fact that the prima/dual solutions returned by the solvers can be
  feasible only for a relaxation of the "whole problem". For instance, when
  solving a flow decomposition of MMCF, the unfeasible flow solution might be
  heuristically turned into a feasible multicommodity flow, while when solving
  a Fixed Charge MMCF a feasible solution is easily obtained from a feasible
  flow by just rounding up the fractional design veriables).
  @{ */

   virtual FONumber GetPVal( void ) = 0;

   virtual FONumber GetDVal( void ) = 0;

/**< These methods must return the objective function value, respectively, of
   the primal and dual "optimal" (whatever this means for the actual solver)
   solutions returned by Get[P/D]Sol() [see below].

   The ws set by SetSubP() has the same meaning as in SolveMMCF() above, i.e.,
   for 1 <= ws <= NrSubP() the value of the ws-th subproblem is returned, for
   ws == 0 the value of the "constant" is returned and for any ws > NrSubP()
   the total value is returned.

   If GetPVal() returns +INF, then no feasible primal solution has been
   found (the primal problem is a minimization problem). Analogously, if
   GetDVal() returns -INF, then no feasible dual solution has been found.

   The solver may provide multiple primal/dual solutions, not all of which
   optimal and therefore with the same objective function value. Thus,
   Get[P/D]Val() has to be called for each different primal/dual solution
   extracted with Get[P/D]Sol() [see below]. In particular, since each call
   to Get[P/D]Sol() tells whether or not a new solution exists, and thus it
   may "change the solution that is currently kept in memory", the safe
   order is to Get[P/D]Val() *before* Get[P/D]Sol(), since when Get[P/D]Sol()
   is called again the value returned by Get[P/D]Val() might change, being
   the one corresponding to the solution that Get[P/D]Sol() will return next
   time rather than the value of the one it has just returned. Of course,
   this only applies if Get[P/D]Sol() returns true (but in principle the
   caller has no way to know beforehand whether or not this will happen). */

/*--------------------------------------------------------------------------*/

   virtual FONumber GetUpprBnd( bool &HvSol )
   {
    HvSol = false;              // by default, no UB solution is available
    return( Inf<FONumber>() );  // by default, no UB is known
    }

   virtual FONumber GetLwrBnd( bool &HvSol )
   {
    HvSol = false;                // by default, no LB solution is available
    return( - Inf<FONumber>() );  // by default, no LB is known
    }

/**< These methods have to return any known Upper/Lower Bound on the optimal
   value of the "whole" problem (+/-INF are clearly always a possibility).

   They also have to set HvSol to true if a primal/dual solution of the
   "whole" problem attaining that value is available, and false otherwise:
   in the former case, these solutions can be obtained with a proper
   combinations of Set*****() and Get[U/L]BSol() [see below]. If multiple
   primal/dual solutions are available, each one possibly gives a different
   Upper/Lower Bound: thus, Get[Uppr/Lwr]Bnd() has to be called for each of
   them. */

/*--------------------------------------------------------------------------*/

   virtual void SetFlwSol( FRow Flw = 0 , Index_Set Bse = 0 ,
			   cIndex wf = Inf<Index>() )
   {
    FBse = ( FSol = Flw ) ? Bse : 0;
    WhchFS = wf;
    }

   virtual void SetMFlwSol( MFRow Flw = 0 , Index_Set Bse = 0 ,
			    cIndex wf = Inf<Index>() )
   {
    FBse = ( MFSol = Flw ) ? Bse : 0;
    WhchFS = wf;
    }

   virtual void SetXtrSol( Row Xtr = 0 , Index_Set Bse = 0 )
   {
    XBse = ( XSol = Xtr ) ? Bse : 0;
    }

/**< Set[M]FlwSol() and SetXtrSol() are meant to pass to the object pointers
   to the memory where (respectively) a flow solution and an extra variables
   solution (the "primal information") have to be, along with "instructions"
   about their "format". They are used in pair with GetPSol() and GetUBSol()
   [see below] to access solutions, but they may also be used with CostOf()
   [see below] as a mean for inputting a solution. This double use may
   generate confusion: read carefully the comments to the above methods.

   Set[[M]Flw/Xtr]Sol() must be called before the methods that use the
   information provided therein: for the "output" methods Get[P/UB]Sol(), it
   is also a good practice to call them before MMCFSolve(). Thus, the solver
   will know if the solution are required while it is solving the problem,
   and may for instance use that memory to write the solution rather than
   allocating its own. Whether or not this is required is solver-dependent.

   Passing 0 as the first argument means that that set of variables is
   not required/provided: this is useful for accessing/passing only one of
   the two sets at a time.

   Passing 0 as the second argument means that the corresponding solution
   is in "dense" format, i.e. the i-th variable must be found in the i-th
   position of the vector. For flow variables, the format is commodity-wise

    first arc of the first commodity .. m-th arc of the first commodity ..
      :                                  : 
    first arc of the K-th commodity  .. m-th arc of the K-th commodity.

   while for extra variables it clearly depends on the application. If a
   non-0 is passed instead, it is taken as the pointer to a vector where
   the indices of the *nonzero* variables have to be: hence, the solution is
   in a "sparse" format, with e.g. Flw[ i ] containing the (nonzero) flow of
   the variable Bse[ i ]. Bse must be ordered in increasing sense, without
   duplications and Inf<Index>()-terminated, i.e. an Inf<Index>() must be
   placed right after the last significative entry.

   For Set[M]FlwSol(), wf means that what is required is

   - only the flow of the wf-th commodity if 0 <= wf < K;

   - the *aggregate* flow (i.e. the *sum* of the flows for *all* the
     commodities) if wf == K;

   - the *disaggregate* flow of *all* the commodities if wf > K.

   Important note about "names": for all wf <= K, if a `Bse' has been set,
   then the "names" written in Bse go from 0 to m - 1 because they refere to
   an arc of the graph; the commodity is clear if 0 <= ws < K, and there is
   no commodity at all if ws == K. Conversely, for ws > NComm the names go
   from 0 to m * K - 1, as an item with name `h' referes to arc h % m of
   commodity h / m.

   Note that wf here and ws in SetSubP() refere to two distinct concepts, 
   although decomposition by commodity is clearly one possibility.

   The two similar methods concerning flow solutions, SetFlwSol() and
   SetMFlwSol(), are provided as alternatives in order to take into account
   one (unpleasant) carachteristic of Multicommodity flows, i.e. that the
   optimal Multicommodity flow solution of a problem with all-integer
   capacities and deficits may *not* be integral. Thus, the type `MFNumber'
   is explicitly defined as to be distinct from `FNumber', as the former
   typically has to be a float type, while the latter is often an integer
   type. However, there are cases where the solution found by SolveMMCF()
   is of the `FNumber' type: think e.g. to the Flow Relaxation or to the
   Integer MMCF problem. Thus, the interface provides support for both
   cases - see also **SolIsFNumber() and Get**Sol() below. */

/*--------------------------------------------------------------------------*/

   virtual bool GetPSol( void ) = 0;

   virtual bool GetUBSol( void )
   {
    if( FSol ) throw(
     MMCFException( "MMCFClass::GetUBSol() called with non-0 FSol" ) );

    if( MFSol ) throw(
     MMCFException( "MMCFClass::GetUBSol() called with non-0 MFSol" ) );

    if( XSol ) throw(
     MMCFException( "MMCFClass::GetUBSol() called with non-0 XSol" ) );

    return( false );    // no other solutions are available
    }

/**< Called *after* SolveMMCF() and Set[[M]Flw/Xtr]Sol(), these methods write
   in the memory provided for the purpose respectively the "optimal" primal
   solution of the problem and the primal solution of the "whole" problem
   corresponding to the Upper Bound, according to the required "style"
   [see Set[[M]Flw/Xtr]Sol()].

   The ws set in SetSubP() tells if the solution is relative to one specific
   subproblem or to the whole problem (ws == 0 is clearly not allowed here).
   Note that, although each subproblem only has one subset of the variables,
   each solution is in principle a "full" one, i.e., no assumption is done on
   how the partition intersects the commodities, although decomposition by
   commodity is clearly one possibility.

   The solver may be capable of providing alternative "good" solutions, other
   than the one that it is currently returning; if it is capable (and willing)
   to do that, it should return true, and false otherwise. However, a return
   value of true does not imply that the method *must* be called again, the
   choice being left to the caller. Yet what it does mean is that the value
   returned by Get[PVal/UpprBnd]() [see above] after the call is different
   from that returned before, since this is intended to be the value of
   "the next solution that Get[P/UB]Sol() will return if called", and since
   not all returned solutions are optimal (not even the first, necessarily)
   then these values could be different.

   If no feasible [primal] solution has been found then GetPVal() returns
   +INF [see above] GetPSol() should not be called, as there is no such
   thing as a solution to be got.

   Once these methods have been called, they (like all those "consuming"
   solutions) should "reset" some internal state to indicate that the FSol
   and XSol vectors provided by the previous call[s] to Set[[M]Flw/Xtr]Sol()
   are no more valid: this is to avoid that, for instance, a subsequent call
   to CostOf() think that FSol and XSol be its "input" parameters. A nice way
   for doing that is just to set FSol / MFSol / XSol to 0. Therefore, if
   (these or others) FSol / MFSol and XSol are needed again by some method,
   they must be passed again with Set[[M]Flw/Xtr]Sol() after a call to
   Get[P/UB]Sol(). */

/*--------------------------------------------------------------------------*/

   virtual bool PSolIsFNumber( void )
   {
    return( false );
    }

   virtual bool UBSolIsFNumber( void )
   {
    return( false );
    }

/**< These methods tell if the primal solutions (to be retrieved with
   GetPSol()) and the upper bound solutions (to be retrieved with GetUBSol())
   are of type `FNumber' rather than `MFNumber'.
   If PSolIsFNumber() returns true (which is clearly solver-dependent), then
   SetFlwSol() has to be used when dealing with primal solutions, otherwise
   SetMFlwSol() has to be used; the same holds for upper bound solutions. */

/*--------------------------------------------------------------------------*/

   virtual void SetNPot( CRow NPt = 0 , cIndex wd = Inf<Index>() )
   {
    NPot = NPt;
    WhchNP = wd;
    }

   virtual void SetRCst( CRow RCs = 0 , cIndex wd = Inf<Index>() )
   {
    RCst = RCs;
    WhchRC = wd;
    }

   virtual void SetMCCst( CRow MCs = 0 )
   {
    MCst = MCs;
    }

   virtual void SetXtrRC( Row XRC = 0 )
   {
    XtRC = XRC;
    }

   virtual void SetXtrDV( Row XDV = 0 )
   {
    XtDV = XDV;
    }

/**< SetNPot(), SetRCst(), SetMCCst(), SetXtrRC() and SetXtrDV() are meant to
   pass to the object pointers to the memory where, respectively,

   - the Node Potentials, i.e. the dual costs of the Flow Conservation
     constraints (1.k),

   - the Reduced Costs of the flow variables, i.e. essentially the dual
     costs of the Upper Bound constraints (2.k),

   - the dual costs of Mutual Capacity constraints (3),

   - the Reduced Costs of the "extra" (non-flow) variables, i.e. essentially
     the dual costs of the range constraints (5),

   - and the dual costs of the "extra" constraints (4)

   (the "dual information") have to be, along with "instructions" about their
   "format". They are used in pair with GetDSol() and GetLBSol() [see below]
   to access solutions, but they may also be used with CostOf() and GetSubG()
   [see below] as a mean for inputting a solution.

   Set[NPot/RCst/MCCst/XtrRC/XtrDV]() must be called before the methods that
   use the information provided therein: for the "output" methods
   Get[D/LB]Sol(), it is also a good practice to call them before MMCFSolve(),
   so that the solver will know if the solution are required while it is
   solving the problem, and may for instance use that memory to write the
   solution rather than allocating its own. Whether or not this is required
   is solver-dependent.

   Passing 0 as the first argument means that that set of variables is
   not required/provided: this is useful for accessing/passing only one of
   the three sets at a time.

   For SetNPot() and SetRCst(), wd means that what is required is

   - only the components of the wd-th commodity if 0 <= wd < K;

   - the "full" solution, i.e., the components of *all* the commodities if
     wd >= K.

   If wd >= K, the solution is in "dense" format: for Node Potentials and
   flow Reduced Costs the format is commodity-wise, i.e. respectively

    first node of the first commodity .. n-th node of the first commodity ..
      :                                    : 
    first node of the K-th commodity  .. n-th node of the K-th commodity

   and

    first arc of the first commodity .. m-th arc of the first commodity ..
      :                                  : 
    first arc of the K-th commodity  .. m-th arc of the K-th commodity

   while for the dual Costs of Mutual Capacities it is arc-wise, i.e.

    first arc .. m-th arc.

   Note that wf here and ws in SetSubP() refere to two distinct concepts, 
   although decomposition by commodity is clearly one possibility.

   Note that, in general, the reduced costs of a bounded variable is the
   opposite of the dual cost of the corresponding bound constraint; in fact,
   the dual cost can be > 0 only if the primal variable attains its bound,
   hence its reduced cost is <= 0. */

/*--------------------------------------------------------------------------*/

   virtual bool GetDSol( void ) = 0;

   virtual bool GetLBSol( void )
   {
    if( NPot ) throw(
     MMCFException( "MMCFClass::GetLBSol() called with non-0 NPot" ) );

    if( RCst ) throw(
     MMCFException( "MMCFClass::GetLBSol() called with non-0 RCst" ) );

    if( MCst ) throw(
     MMCFException( "MMCFClass::GetLBSol() called with non-0 MCst" ) );

    if( XtRC ) throw(
     MMCFException( "MMCFClass::GetLBSol() called with non-0 XtRC" ) );

    if( XtDV ) throw(
     MMCFException( "MMCFClass::GetLBSol() called with non-0 XtDV" ) );

    return( false );   // no other solutions are available
    }

/**< Called *after* SolveMMCF() and Set[NPot/RCst/MCCst/XtrRC/XtrDV](), these
   methods write in the memory provided for the purpose respectively the
   "optimal" dual solution of the problem and the dual solution of the
   "whole" problem corresponding to the Lower Bound, according to the
   required "style".

   The ws set in SetSubP() tells if the solution is relative to one specific
   subproblem or to the whole problem (ws == 0 is clearly not allowed here).
   Note that if the problem is separable into independent subproblems then
   each subproblem may have only a subset of the constraints; therefore, for
   each subproblem only a subset of the dual variables is significative, and
   for each dual variable there is (at most) one subproblem which it belongs
   to. Nonethless, each dual solution is in principle a "full" one, i.e., no
   assumption is done on how the dual variables are partitioned among the
   subproblems.

   The solver may be capable of providing alternative "good" solutions, other
   than the one that it is currently returning; if it is capable (and willing)
   to do that, it should return true, and false otherwise. However, a return
   value of true does not imply that the method *must* be called again, the
   choice being left to the caller. Yet what it does mean is that the value
   returned by Get[FVal/LwrBnd]() [see above] after the call is different
   from that returned before, since this is intended to be the value of
   "the next solution that Get[D/LB]Sol() will return if called", and since
   not all returned solutions are optimal (not even the first, necessarily)
   then these values could be different.

   If no feasible [dual] solution has been found then GetDVal() returns
   -INF [see above] GetDSol() should not be called, as there is no such
   thing as a solution to be got.

   These methods (like all those "consuming" solutions) should do something
   equivalent to setting NPot / RCst / MCst / XtRC / XtDV to 0 for avoiding
   them to be wrongly used (e.g. as "input" parameters) after the call: see
   the comments to Get[P/UB]Sol(). */

/*--------------------------------------------------------------------------*/

   virtual FONumber CostOf( void ) = 0;

/**< Evaluates and return the cost of the primal or dual solutions provided by
   the Set****() methods above. The ws set in SetSubP() [see above] has the
   same meaning as in Get[P/D]Val() [see above], i.e., tells the cost of
   which subproblem is required to be calculated.

   Clearly, the "full" primal/dual solution is required in order to compute
   the "full cost" (for ws >  NrSubP()); however, only the "partial" solution
   corresponding to the ws-th subproblem is required if 1 <= ws <= NrSubP().
   Of course, the "full" solution is also fine in this case.

   Note that only one method is provided for both primal and dual costs, so
   only a primal *or* a dual solution should have been given with Set****()
   before the call to CostOf(): otherwise, the return value is undefined.

   This method (like all those "consuming" solutions) should do something
   equivalent to setting FSol / XSol to 0 for avaiding them to be wrongly
   used after the call: see the comments to Get***Sol(). */

/*--------------------------------------------------------------------------*/

   void TimeMMCF( double &t_us , double &t_ss )
   {
    t_us = t_ss = 0;
    if( MMCFt ) MMCFt->Read( t_us , t_ss ); 
    }

   double TimeMMCF( void )
   {
    return( MMCFt ? MMCFt->Read() : 0 );
    }

/**< If these methods are called within any of the methods of the class that
   are "actively timed" (this depends on the subclasses), they return
   respectively the user and sistem time and the total time (in seconds) since
   the start of that method. If methods that are actively timed call other
   methods that are actively timed, these methods return the (...) time since
   the beginning of the *outer* actively timed method.
   If these methods are called outside of any actively timed method, they
   return the (...) time spent in all the previous executions of all the
   actively timed methods of the class.

   Implementing the proper calls to MMCFt->Start() and MMCFt->Stop() is due
   to derived classes; these should at least be placed at the beginning and
   at the end, respectively, of SolveMMCF() and presumably the Chg***()
   methods, that is, at least these methods should be "actively timed". */

/* @}-----------------------------------------------------------------------*/
/*-------------- METHODS FOR READING THE DATA OF THE PROBLEM ---------------*/
/*--------------------------------------------------------------------------*/
/** @name Reading the data of the problem
    @{ */

   Index NrComm( void )
   {
    return( NComm );
    }

/*--------------------------------------------------------------------------*/

   Index NrNodes( void )
   {
    return( NNodes );
    }

/*--------------------------------------------------------------------------*/

   Index NrArcs( void )
   {
    return( NArcs );
    }

/*--------------------------------------------------------------------------*/

   virtual Index NrXtrVrs( void )

/**< This method returns how many "extra" (non-flow) variables are there other
   than the m * K Multicommodity Flow variables; this is the lenght of the
   vectors to be passed to SetXtrSol() and SetXtrRC(). */
   {
    return( XtrVrs );
    }

/*--------------------------------------------------------------------------*/

   virtual Index NrXtrCnst( void )

/**< This method returns how many "extra" constrainst (5) are there; this is
   the lenght of the vectors to be passed to SetXtrDV(). */
   {
    return( XtrCnst );
    }

/*--------------------------------------------------------------------------*/

   virtual Index NrSubP( void )

/**< The class supports the notion that the MMCF-related problem at hand might
   be decomposable into a certain number of independent subproblems; this
   can be exploited by some solution approaches. NrSubP() returns the number
   of independent subproblems into which the MMCF-related problem can be
   decomposed, which must clearly be > 0. */
   {
    return( NSubPr );
    }

/*--------------------------------------------------------------------------*/

   virtual void GetCosts( CRow Csts , cIndex_Set nms = 0 ,
			  cIndex strt = 0 , Index stp = Inf<Index>() ) = 0;

/**< Reads the costs of the flow variables in the MMCF instance.

   Flow variables are numbered from 0 to m * K - 1 according to the format of
   FF in GetSol() [see above]. The flow variables whose costs are written in
   Csts[] are all and only those whose name is
   - comprised between strt (included) and min( stp , m * K ) (excluded);
   - contained in nms[] if nms != 0 (in this case, it has to be a vector of
     indices in [ 0 .. m * K ), with no duplicated elements, ordered in
     increasing sense and Inf<Index>()-terminated). */

/*--------------------------------------------------------------------------*/

   virtual void GetXtrCsts( CRow XtrCs , cIndex_Set nms = 0 ,
			    cIndex strt = 0 , Index stp = Inf<Index>() )

/**< Reads the costs of the "extra" variables in the MMCF instance; these
   variables are numbered from 0 to NrXtrVrs() - 1, according to the format
   of XV in GetSol() [see above].

   The extra variables whose costs are written in XtrCs[] are all and only
   those whose name is
   - comprised between strt (included) and min( stp , NrXtrVrs() ) (excluded);
   - contained in nms[] if nms != 0 (in this case, it has to be a vector of
     indices in [ 0 .. NrXtrVrs() ), with no duplicated elements, ordered in
     increasing sense and Inf<Index>()-terminated). */
   {
    // by default there are no extra variables, hence nothing has to be done
    }

/*--------------------------------------------------------------------------*/

   virtual void GetICaps( FRow ICps , cIndex_Set nms = 0 ,
			  cIndex strt = 0 , Index stp = Inf<Index>() ) = 0;

/**< Reads the Individual Capacities of the flow variables in the MMCF
   instance. The format of the returned ICps[] depends on nms[], strt and stp
   exactly as in GetCosts() [see above]. */

/*--------------------------------------------------------------------------*/

   virtual void GetMCaps( FRow MCps , cIndex_Set nms = 0 ,
			  cIndex strt = 0 , Index stp = Inf<Index>() ) = 0;

/**< Reads the Mutual Capacities of the flow variables in the MMCF instance. 

   The mutual capacities written in MCps[] are all and only those of arcs
   whose name is
   - comprised between strt (included) and min( stp , m ) (excluded);
   - contained in nms[] if nms != 0 (in this case, it has to be a vector of
     indices in [ 0 .. m ), with no duplicated elements, ordered in increasing
     sense and Inf<Index>()-terminated). */

/* @}-----------------------------------------------------------------------*/
/*------------- METHODS FOR ADDING / REMOVING / CHANGING DATA --------------*/
/*--------------------------------------------------------------------------*/
/** @name Adding/removing/changing data
   @{ */

   virtual void ChgCosts( cCRow NwCsts , cIndex_Set nms = 0 ,
			  cIndex strt = 0 , Index stp = Inf<Index>() ) = 0;

/**< Changes the costs of the flow variables in the MMCF instance.

   Flow variables are numbered from 0 to m * K - 1 according to the format of
   FF in GetSol() [see above]. If nms == 0, the new cost for the flow
   variable strt + i will be taken from NwCsts[ i ] for all 0 <= i <
   min( m * K , stp ) - strt. If nms != 0, it has to point to a vector of
   indices in [ 0 .. m * K ) (with no duplicated elements, ordered in
   increasing sense and Inf<Index>()-terminated); for all j such that strt <=
   nms[ j ] < min( stp , m * K ), the new cost for the flow variable nms[ j ]
   will be taken from NwCsts[ j ]. */

/*--------------------------------------------------------------------------*/

   virtual void ChgXtrCsts( cCRow NwXtrCs , cIndex_Set nms = 0 ,
			    cIndex strt = 0 , Index stp = Inf<Index>() )

/**< Changes the costs of the "extra" variables in the MMCF instance; these
   variables are numbered from 0 to NrXtrVrs() - 1 according to the format
   of XV in GetSol() [see above]. If nms == 0, the new cost for the extra
   variable strt + i will be taken from NwXtrCs[ i ] for all 0 <= i <
   min( NrXtrVrs() , stp ) - strt. If nms != 0, it has to point to a vector
   of indices in [ 0 .. NrXtrVrs() ) (with no duplicated elements, ordered in
   increasing sense and Inf<Index>()-terminated); for all j such that strt <=
   nms[ j ] < min( stp , NrXtrVrs() ), the new cost for the extra variable
   nms[ j ] will be taken from NwXtrCs[ j ]. */
   {
    // by default there are no extra variables, hence nothing has to be done
    }

/*--------------------------------------------------------------------------*/

   virtual void ChgICaps( cFRow NwICps , cIndex_Set nms = 0 ,
			  cIndex strt = 0 , Index stp = Inf<Index>() ) = 0;

/**< Changes the Individual Capacities of the flow variables in the MMCF
   instance. The format of the returned NwICaps[] depends on nms[], strt and
   stp exactly as in ChgCosts() [see above]. */

/*--------------------------------------------------------------------------*/

   virtual void ChgMCaps( cFRow NwMCps , cIndex_Set nms = 0 ,
			  cIndex strt = 0 , Index stp = Inf<Index>() ) = 0;

/**< Changes the Mutual Capacities of the flow variables in the MMCF instance. 
   If nms == 0, the new mutual capacity for arc i will be taken from 
   NwMCps[ i ] for all 0 <= i < min( m , stp ) - strt. If nms != 0, it has
   to point to a vector of indices in [ 0 .. m ) (with no duplicated elements,
   ordered in increasing sense and Inf<Index>()-terminated); for all j such
   that strt <= nms[ j ] < min( stp , m ), the mutual capacity for the arc
   nms[ j ] will be taken from NwMCps[ j ]. */

/*--------------------------------------------------------------------------*/

   virtual void ChgIntVar( cIndex k = Inf<Index>() , bool IntVld = true ,
                           cIndex_Set nms = 0 , Index strt = 0 ,
                           Index stp = Inf<Index>() ) = 0;

/**< Changes the status of (a subset of) the - flow and extra - variables of
   the problem in terms of the integrality constraints imposed upon them.

   ChgIntVar( k , true/false , ... ) for k < K says that some of the variables
   of the commodity k are/aren't integer-valued. The variables to which the
   change is applied are the flow variables of commodity k whose index is:
   - in the set nms[] (which contains indices in [0, m ), ordered in
     increasing sense and Inf<Index>()-terminated);
   - comprised between strt (included) and min( stp , m ) (excluded).
   nms == 0 means "all in the interval [strt, stp)".

   The status of all other variables remains unchanged.

   ChgIntVar( K , ... ) can be used to set the integrality of the "extra"
   (non-flow) variables; of course, in this case the indices must be in the
   range [0, NrXtrVrs() ).

   ChgIntVar( k , ... ) with k > K (e.g. k == Inf<Index>()) applies the
   change to all flow variables, irrespective of the commodity to which they
   belong; in this case, the indices in nms[], strt and stp must be in the
   [ 0 , m * K ) range. */

/*--------------------------------------------------------------------------*/

#if( CHGARCS_MMCF )

   virtual void CloseArcs( cIndex_Set whch ) = 0;

   virtual void OpenArcs( cIndex_Set whch ) = 0;

 #if( CHGARCS_MMCF > 1 )

   virtual void CloseArcs( cIndex k , cIndex_Set whch ) = 0;

   virtual void OpenArcs( cIndex k , cIndex_Set whch ) = 0;

 #endif
#endif

/**< Respectively "Close" and "Open" the arcs indicated in whch, that must
   point to a vector of indices in [ 0 .. m - 1 ] (with no duplicated
   elements, ordered in increasing sense and Inf<Index>()-terminated). The
   first forms close/open arcs for *all* the commodities simultaneously, while
   the second forms do it only for the commodity 'k'.

   Closing an arc for a given commodity 'k' is conceptually equivalent to
   setting its Individual Capacity to zero [see ChgICaps() above]; closing an
   arc for all the commodities is conceptually equivalent to closing the arc
   for all the commodities individually and/or to setting its Mutual Capacity
   to zero [see ChgMCaps() above]. Opening an arc is conceptually equivalent
   to restoring the capacities to their initial value.

   However, the solver may implement the above operations more efficiently;
   furthermore, the original capacities are kept by the solver itself when
   an arc is closed, so that they do not need to be stored somewhere else for
   when it is re-opened.

   Opening an arc that has not previusly been Closed is an error: however, it
   is admitted to call CloseArcs( { i } ) even if CloseArcs( k , { i } ) has
   already been called. The "individual" closure has "precedence" over the
   "global" one, that is if OpenArcs( { i } ) is called afterwards then the
   arc remain closed on all the commodities k such that CloseArcs( k , { i } )
   had been called previously, until OpenArcs( k , { i } ) is explicitly
   called for each of these commodities.

   OpenArcs( k , { i } ) should not be called for an arc that has been
   "globally" closed with CloseArcs( { i } ). */

/*--------------------------------------------------------------------------*/

   virtual void AddExtraVars( cIndex NXV , cCRow XCst = 0 ,
			      cFRow XUb = 0 , cFRow XLb = 0 ,
			      bool IntVar = false , cIndex_Set nms = 0 )

/**< AddExtraVars() adds to the problem NXV extra variables. Costs, upper and
   lower bounds of the new variables are given respectively in the read-only
   vectors pointed by XCst, XUb and XLb; 0 means "all zero" for XCst and
   XLb, and it means "all +infinity" (= no upper bound) for XUb.

   IntVar == true means that some of the new variables are integer-valued. If
   nms != 0, it must point to a vector of indices (ordered in increasing
   sense and Inf<Index>() terminated) containing the indices of those of the
   new extra variables just being added (thus, indices in the [0, NXV) range)
   that are integer-valued; otherwise, all the new variables ar
   integer-valued.

   Note that this method only adds variables: AddExtraConstr() [see below] is
   provided for adding "extra" linear constraints, which may (or may not)
   involve these extra variables. */
   {
    throw(
     MMCFException( "MMCFClass: adding extra variables is not allowed" ) );
    }

/*---------------------------------------------------------------------------*/

   virtual void AddExtraConstr( cIndex NXC , int *IBeg , int *Indx ,
				double *Vals , cFRow XLr = 0 , cFRow XUr = 0 )

/**< AddExtraConstr adds to the problem NXC new linear constraints. IBeg, Indx
   and Vals must contain the description of the new constraints: each
   constraint is represented by the set of indices of variables with nonzero
   coefficient and the corresponding coefficients. The indices and
   coefficients of the i-th constraint, i = 0 .. NXC - 1, must be in Indx and
   Vals, respectively, in the positions between IBeg[ i ] (included) and
   IStp[ i + 1 ] (excluded). The indices corresponding to each constraint,
   in Indx, must be ordered in increasing sense and without duplications.

   The mapping between indices and variables is the following: the variable
   corresponding to arc j for commodity k has the index k * m + j; the i-th
   "extra" variable has the index K * m + i (the order is "commodity-wise",
   with the "extra" variables following).

   Note that IBeg, Indx and Vals are *not* pointers to const: this is done to
   allow the implementation to *overwrite* those vectors, should it find this
   useful. Furthermore, the implementation *need not* to restore the original
   containts of the vectors: no assumption should be done by the caller on the
   contents of IBeg, Indx and Vals after the call. The rationale for this is
   that these data structures may be "big", and that they should be of no use
   after the call, so by allowing the implementation to exploit that memory
   some serious memory problem may be avoided.

   XLr and XUr are, respectively, the lower and upper range of the extra
   constraints: if XLr[ i ] == - F_INF, then the row have sense <=, while if
   XUr[ i ] == F_INF, then the row have sense >=. XLr == 0 means "all
   equal to - F_INF", and XUr == 0 means "all equal to F_INF". It is an
   error if XLr[ i ] == - F_INF *and* XUr[ i ] == F_INF for some i (this is
   not a constraint), and therefore at least one between XLr and XUr must be
   non-0. */
   {
    throw(
     MMCFException( "MMCFClass: adding extra constraints is not allowed" ) );
    }

/* @}-----------------------------------------------------------------------*/
/*------------------------------ DESTRUCTOR --------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Destructor
    @{ */

   virtual ~MMCFClass()

/**< Destructor of the class: it must be virtual. */
   {
    delete MMCFt;
    }

/* @}-----------------------------------------------------------------------*/
/*--------------------- PROTECTED PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

 protected:

/*--------------------------------------------------------------------------*/
/*----------------------- PROTECTED DATA STRUCTURES  -----------------------*/
/*--------------------------------------------------------------------------*/
/** @name Protected fields of the class
    @{ */

  Index NArcs;      ///< Number of Arcs of the graph
  Index NNodes;     ///< Number of Nodes of the graph
  Index NComm;      ///< Number of Commodities
  Index XtrVrs;     ///< Number of "eXtra" Variables
  Index XtrCnst;    ///< Number of "eXtra" Constraints
  Index NSubPr;     ///< into how many SubProblems the MMCF can be divided

  FONumber OptEps;  ///< an OptEps-optimal solution is required
  FONumber FsbEps;  ///< FsbEps is the max allowed violation of constraints

  Index WhchSP;     ///< WhchSP tells which subproblem is one referring to

  FRow  FSol;       ///< pointer to the memory of the Flow Solution
  MFRow MFSol;      ///< as above, but of the MFNumber type
  Index_Set FBse;   ///< if FBse != 0, then it must be "sparse"
  Index WhchFS;     ///< WhchFS tells its "type" of the Flow Solution

  Row   XSol;       ///< pointer to the memory of the "extra" solution
  Index_Set XBse;   ///< if XBse != 0, then it must be "sparse"

  CRow NPot;        /**< pointer to the memory of the Node Potentials (dual
		       costs for the Flow Conservation constraints (1.k)) */
  Index WhchNP;     ///< WhchNP tells its "type" of the Node Potentials

  CRow RCst;        /**< pointer to the memory of the flow reduced costs
		       (dual costs for the bound constraints (2.k))*/
  Index WhchRC;     ///< WhchRC tells its "type" of the flow reduced costs
  
  CRow MCst;        /**< pointer to the memory of the dual costs for the
		       Mutual Capacity constraints (3)*/
  Row XtDV;         /**< pointer to the memory of the dual costs for the
		       "extra" constraints (4)*/
  Row XtRC;         /**< pointer to the memory of the reduced costs for the
		       "extra" variables (dual costs for constraints (5))*/

  ostream *MMCFLog;  ///< the output stream object
  char MMCFLLvl;     ///< the "level of verbosity" of the log

  OPTtimers *MMCFt;  ///< mainly the MMCFSolve() time, probably

/* @}-----------------------------------------------------------------------*/

 };  // end( class MMCFClass )

/*--------------------------------------------------------------------------*/

#if( OPT_USE_NAMESPACES )
 };  // end( namespace MMCFClass_di_unipi_it )
#endif

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif  /* MMCFClas.h included */

/*--------------------------------------------------------------------------*/
/*------------------------End File MMCFClas.h ------------------------------*/
/*--------------------------------------------------------------------------*/

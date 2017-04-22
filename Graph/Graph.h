/*--------------------------------------------------------------------------*/
/*---------------------------- File Graph.h --------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the Graph class, which provides an unified mean for
 * reading descriptions of (Linear) Multicommodity Min Cost Flow problems and
 * storing them in memory, along with a simple interface that can be used by
 * any MMCF solver to access and change the data.
 *
 * It also provides some pre-processing features, aimed at making the MMCF
 * instance easier to solve, such as identification of redundant mutual
 * capacity constraints.
 *
 * \version 2.01
 *
 * \date 11 - 05 - 2012
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
 * Copyright &copy 1994 - 2012 by Antonio Frangioni
 */

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __Graph
 #define __Graph  /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <exception>

#include <limits>

/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE and USING ----------------------------*/
/*--------------------------------------------------------------------------*/

#if( OPT_USE_NAMESPACES )
namespace MMCFGraph_di_unipi_it
{
#endif

/*--------------------------------------------------------------------------*/
/*----------------------------- CLASS Graph --------------------------------*/
/*--------------------------------------------------------------------------*/
/*---------------------------- GENERAL NOTES -------------------------------*/
/*--------------------------------------------------------------------------*/
/** The class Graph is a service class for solvers of (generalized)      
    Multicommodity Min-Cost Flow problems, i.e., problems of the form     

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
                                                                     
  The flow variables X[ k ] and constrains (1.k), (2.k) and (3) are  
  present in all Multicommodity-type problems; the extra variables Y  
  and the extra constraints (4) and (5) are optional.                 
                                                                      
  Flow and extra variables can be declared to be either continuous or  
  integer-valued.                                                      
                                                                      
  The Graph class allows to read a description of a MMCF problem from  
  file or from memory and to make them available in an uniform way to  
  solvers. This class is intended essentially only for initialization; 
  the solvers will copy all the relevant information in their internal 
  data structures, and Graph can be deleted afterwards. However, Graph 
  also offers some nice pre-processing features, intended to make the  
  instance more easily solvable.                                       
                                                                      
  The base class Graph deals with "generic" MMCF problems, i.e. where  
  there are only constrains (1.k), (2.k) and (3); extra constraints    
  and variables can be added, but there is limited support for them.   
  However, all the methods for handling extra constrains and variables 
  are virtual, hence derived classes can be implemented that deal with 
  special cases of MMCF instances with a special structire of the      
  extra constrains and variables.                                      
                                                                      
  In all the comments below, m has to be understood as the number of   
  arcs in the graph G (#A), n as the number of nodes in the graph G    
  (#N) and K as the number of commodities. */

class Graph
{

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*--------------------------------------------------------------------------*/
/*---------------------------- PUBLIC TYPES --------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Public types
    Graph defines three main public types:

    - Index, the type of arc and node indices;

    - FNumber, the type of arc capacities (individual and mutual) and node
      deficits;

    - CNumber, the type of flow costs.

    By re-defining the types in this section one can reduce the memory
    footprint of the object in case "small" data types (e.g., integer ones)
    can be used. However, *it is the user's responsibility to ensure that
    these types are set to reasonable values*.
    @{ */

/*--------------------------------------------------------------------------*/

 typedef unsigned int    Index;          ///< index of a node or arc ( >= 0 )
 typedef Index          *Index_Set;      ///< set (array) of indices
 typedef Index_Set      *Index_Mat;      ///< set of set of indices

 typedef const Index    cIndex;          ///< a read-only index
 typedef cIndex        *cIndex_Set;      ///< read-only index array

/*--------------------------------------------------------------------------*/

 typedef double          FNumber;        ///< type of arc flow
 typedef FNumber        *FRow;           ///< vector of flows
 typedef FRow           *FMat;           ///< matrix of flows

 typedef const FNumber  cFNumber;        ///< a read-only flow
 typedef cFNumber      *cFRow;           ///< read-only flow array

/*--------------------------------------------------------------------------*/

 typedef double          CNumber;        ///< type of arc flow cost
 typedef CNumber        *CRow;           ///< vector of costs
 typedef CRow           *CMat;           ///< matrix of costs

 typedef const CNumber  cCNumber;        ///< a read-only cost
 typedef cCNumber      *cCRow;           ///< read-only cost array

/*--------------------------------------------------------------------------*/

 typedef double          FONumber; 
 /**< type of the objective function: has to hold sums of products of
    FNumber(s) by CNumber(s) */

/*--------------------------------------------------------------------------*/

 typedef bool           *Bool_Vec;        ///< vector of booleans

/*--------------------------------------------------------------------------*/
/** Very small class to simplify extracting the "+ infinity" value for a
    basic type (FNumber, CNumber, Index); just use Inf<type>(). */

   template <typename T>
   class Inf {
    public:
     Inf() {}
     operator T() { return( std::numeric_limits<T>::max() ); }
    };

/*--------------------------------------------------------------------------*/
/** Small class for exceptions. Derives from std::exception implementing the
   virtual method what() -- and since what is virtual, remember to always
   catch it by reference (catch exception &e) if you want the thing to work.
   MMCFGException class are thought to be of the "fatal" type, i.e., problems
   for which no solutions exists apart from aborting the program. Other kinds
   of exceptions can be properly handled by defining derived classes with
   more information. */

 class MMCFGException : public std::exception {
 public:
  MMCFGException( const char *const msg = 0 ) { errmsg = msg; }

  const char* what( void ) const throw () { return( errmsg ); }

 private:
  const char *errmsg;
  };

/*@} -----------------------------------------------------------------------*/
/*--------------------- PUBLIC METHODS OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/
/*----------------------------- CONSTRUCTORS -------------------------------*/
/*--------------------------------------------------------------------------*/

   Graph( const char *const FN , char FT = 's' );

/**< Constructor of the Class: builds an instance of MMCF reading the data
   from the input file `FN'; `FT' must be:

   - 's'  (default) for single-file "Canadian" format; this format is actually
        for Fixed-Charge MMCF, hence there is a Fixed Charge Cost information
        attached to each arc. Thus, m extra variables are defined, their cost
	is set to this Fixed Charge Cost, their upper and lower bounds are set
	set to 0 and 1, respectively, and their type is set to integer;

   - 'c'  for single-file PPRN format: actually, the format read by this
        constructor is a bit more general than the standard definition, since
        it is permitted for an arc not to have single-commodity capacity (-1);

   - 'p'  for a (PSP) formulation  \
   - 'o'  for a (OSP) formulation  | for Jones-Lustig (JL) standard four-files
   - 'd'  for a (ODP) formulation  | format
   - 'u'                           /

   - 'm'  for the four-files mnetgen format, an extension of the original file
        format produced by mnetgen (see mnetgen.C), partially compatible with
        the JL format

   All the four-files formats read the instance from the 4 files `FN'.nod (the
   instance size), `FN'.sup (node supply information - but see 'u'), `FN'.arc
   (arc information) and `FN'.mut (mutual capacity constraints information).
   The 'o' and 'u' formats are essentially the same, but with 'u' the node
   supply informations are searched in a `FN'.od.

   Actually, there are some restrictions on some of the formats:

   - in the 'c' format, all the coefficients corresponding to the same
     "extra" constraint must appear consecutively in the file, and the
     constraints must be well-ordered (that is, all the information about
     the extra constraint 1 must be found first, then all the information
     about the extra constraint 2 and so on); also, the separating lines
     must contain noting (they can be avoided);

   - in the 'd' ('u') format, it is not allowed to put -1 in the origin or
     destination columns (1st and 2nd) of *.sup (*.od) file;

   - in the 'o' format, it is not allowed to put -1 in the origin columns
     (1st) of *.sup file.

   After the end of the constructor, the data read is in a "raw" form, i.e.
   if arc j is not defined for commodity k then CostKJ( k , j ) == 
   Inf<CNumber>() and CapacityKJ( k , j ) == 0, if no single-commodity upper
   bound is defined for arc j then CapacityKJ( k , j ) == Inf<FNumber>(), if
   no mutual capacity constraint is defined for arc j then TotalCapacityJ( j )
   == Inf<FNumber>() and if node i is not defined for commodity k then
   DeficitKJ( k , i ) == Inf<FNumber>().

   A preprocessing phase [see PreProcess() below] is available that, among
   other things, makes sure that all the Inf<FNumber>() arc capacities (that
   may disturb some solver) are replaced with proper, finite values. That
   preprocessing also ensures some obvious things, such as that all arcs
   entering/leaving a non-existent node for commodity k are non-existent. */

/*--------------------------------------------------------------------------*/

   Graph( cIndex n , cIndex m , cIndex comm , const cFRow* Def ,
          cIndex_Set S , cIndex_Set E , cFRow CapTot , const cFRow* Cap ,
          const cCRow* Cost );

/**< Constructor of the class: builds an instance of MMCF reading the data
   from the provided parameters

   - n       is the number of nodes of the network;

   - m       is the number of arcs (or edges) of the network;

   - comm    is the number of commodities;

   - Def     is the (K x n)-matrix of the node deficits: Def[ k ][ i ] is
             the deficit of node i relative to commodity k; Def can be 0,
             meaning that all deficits are 0, and likewise each Def[ k ] can
             be 0, meaning that all the deficits for commodity k are 0;

   - S       is the m-vector of the starting nodes of the arcs;

   - E       is the m-vector of the ending nodes of the arcs - hence, the
             i-th arc, 0 <= i < m, is S[ i ] --> E[ i ];

   - CapTot  is the m-vector of the arc mutual capacities; CapTot can be 0,
             meaning that all mutual capacities are infinite;

   - Cap     is the (K x m)-matrix of the arc capacities: Cap[ k ][ i ] is
             the upper capacity of arc i relative to commodity k; Cap can be
             0, meaning that all arc capacities are infinite, and likewise
             each Cap[ k ] can be 0, meaning that all arc capacities for
             commodity k are infinite;

   - Cost    is the (K x m)-matrix of the arc costs: Cost[ k ][ i ] is the
             cost of arc i relative to commodity k; arcs with Cost[ k ][ i ]
             == Inf<CNumber>() are intended as non-existent, thus allowing
	     graphs to be different for each commodity: this is stronger then
	     setting U[ k ][ i ] == 0 [see PreProcess() below]; Cost can be
	     0,
             meaning that all costs are 0, and likewise each Cost[ k ] can
             be 0, meaning that all the costs for commodity k are 0;

   Arc i for commodity k is non-existent if Cost[ k ][ i ] == Inf<CNumber>(),
   and has no single-commodity capacity if U[ k ][ i ] == Inf<FNumber>(). Arc
   i has no mutual capacity if CapTot[ i ] == Inf<FNumber>(). Node i for
   commodity k is non-existent if Def[ k ][ i ] == Inf<FNumber>(). All arcs
   entering/leaving a non-existent node for commodity k should be
   non-existent; this is guaranteed by PreProcess() [see below].

   This constructor builds a problem with no extra variables and constraints,
   and with all the (flow) variables continuous. This can be changed later
   with SetIntVar(), SetExtraVar() and SetExtraConstr() [see below].

   Memory saving trick: sometimes, all arcs have the same capacity or cost,
   independently from the commodity. It is clearly possible to construct only
   one m-vector containing the capacities/costs, and copy its pointer k times
   in Cap[]/Cost[]. The same holds for deficits. Note that, internally, the
   Graph object will initially allocate all memory anyway and copy K times
   the information to separate vectors; however, the redundancy will be
   eliminated in PreProcess() [see below], if it "survives" (preprocessing
   can alter the values of the data, thus making initially identical vectors
   different). */

/*--------------------------------------------------------------------------*/

   inline Graph( void ) {};

/**< "base" constructor of the class, that in fact does nothing: is tought as
   a mean for derived class to provide their own initialization of the data
   structures, in case they are incompatible with the ones given in the
   "standard" constructor. An help to implement your own constructor is
   given by the protected method CmnIntlz(), see below. */

/*--------------------------------------------------------------------------*/
/*-------------- METHODS FOR READING THE DATA OF THE PROBLEM ---------------*/
/*--------------------------------------------------------------------------*/

  inline Index NrComm( void );

/**< Returns the number of commodities of the network. */

/*--------------------------------------------------------------------------*/

  inline Index NrArcs( void );

/**< Returns the number of arcs of the network. */

/*--------------------------------------------------------------------------*/

  inline Index NrNodes( void );

/**< Returns the number of nodes of the network. */

/*--------------------------------------------------------------------------*/

   inline Index NrExtraVars( void );

/**< Returns the number of extra (non-flow) variables in the problem. */

/*--------------------------------------------------------------------------*/

   inline Index NrExtraConst( void );

   inline Index NrExtraNonZ( cIndex FrstC = 0 , cIndex LstC = Inf<Index>() );

/**< NrExtraConst() returns the total number of "extra" constraints.

   NrExtraNonZ( i , j ) returns the total number of nonzeroes in the
   representation the "extra" constraints i, i + 1, ... j as linear two-sided
   inequalities. The "names" (indices) of "extra" constraints go from 0 to
   NrExtraConst() - 1. It is illegal to call NrExtraNonZ() if there are no
   "extra" constraints, i.e., NrExtraConst() == 0. */

/*--------------------------------------------------------------------------*/

   enum MCFType { kMCF ,
		  kSPT
                  };

   inline MCFType ProblemType( cIndex k );

/**< Returns the type of the k-th subproblem: it can be a generic Min Cost Flow
   problem or a Shortest Path Tree problem. The type is initially set to kMCF,
   and turned to kSPT (if it is the case) only if PreProcessing() or
   MakeSingleSourced() [see below] are invoked. */

/*--------------------------------------------------------------------------*/

   inline bool Directed( void );

/**< Returns true if the graph is directed, false otherwise. */

/*--------------------------------------------------------------------------*/

   inline bool NamesStartFrom1( void );

/**< Returns true if node "names" (in the arc description) go from 1 to n,
   false if they go from 0 to n - 1. */

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

   inline cFRow TotCapacities( void );

   inline FNumber TotalCapacityJ( cIndex j );

/**<The first form returns a read-only pointer to a m-vector containing the
   mutual capacities of the arcs; the second form returns the mutual capacity
   of arc j. */

/*--------------------------------------------------------------------------*/

   inline cIndex_Set Actives( void );

   inline Index NActives( void );

/**< Actives() returns a read-only pointer to a m-vector (ordered in
   increasing sense and Inf<Index>()-terminated) of the indices of arcs that
   have an associated mutual capacity constraint: it returns 0 if every arc
   has a constraint.

   NActives() returns the number of such mutual capacity constraints:
   Actives() == 0 => NActives() == NrArcs().

   This information is computed in PreProcess() [see below], and therefore
   Actives() will always return 0 if called before PreProcess(). Note that
   ProProcess() finds a finite value for mutual capacities, i.e., after
   ProProcess() TotalCapacityJ( j ) returns something < Inf<FNumber>() even
   if the index j is not contained in the vector returned by Actives(). */

/*--------------------------------------------------------------------------*/

   inline cCRow CostsK( cIndex k );

   inline CNumber CostKJ( cIndex k , cIndex j );

/**< The first form returns a read-only pointer to a m-vector containing the
   costs of the arcs relative to commodity k; the second form returns the
   cost of arc j relative to commodity k.

   If k == NrComm(), the costs of "extra" (non-flow) variables are returned;
   it is illegal to call these methods with k == NrComm() if there are no
   extra variables, i.e. NrExtraVars() == 0. */

/*--------------------------------------------------------------------------*/

   inline cFRow CapacitiesK( cIndex k );

   inline FNumber CapacityKJ( cIndex k , cIndex j );

/**< The first form returns a read-only pointer to a m-vector containing the
   single-commodity capacities relative to commodity k; the second form
   returns the single-commodity capacity of arc j for commodity k.

   For k == NrComm() and k == NrComm() + 1, the lower and upper bounds on the
   "extra" (non-flow) variables are returned, respectively: it is illegal to
   call these methods with k > NrComm() if there are no extra variables, i.e.
   NrExtraVars() == 0. */

/*--------------------------------------------------------------------------*/

   inline cIndex_Set ActivesK( cIndex k );

   inline Index NActivesK( cIndex k );

/**< ActivesK( k ) returns a read-only pointer to the vector (ordered in

   increasing sense and Inf<Index>()-terminated) of the indices of arcs that
   have an associated individual capacity constrint for the commodity k; it
   returns 0 if every arc has its constraint.

   NActivesK( k ) returns the number of such individual capacity constraints
   for commodity k: ActivesK( k ) == 0 => NActivesK( k ) == NrArcs().
   NActivesK( NComm ) returns the total number of such constraints, i.e. the
   sum over all k of NActivesK( k ).

   This information is computed in PreProcess() [see below], and therefore
   Actives( k ) will always return 0 if called before PreProcess(). Note
   that ProProcess() finds a finite value for mutual capacities, i.e. after
   ProProcess() CapacityKJ( k , j ) returns something < Inf<FNumber>() even
   if the index j is not contained in the vector returned by Actives( k ). */

/*--------------------------------------------------------------------------*/

   inline cFRow DeficitsK( cIndex k );

   inline FNumber DeficitKJ( cIndex k , cIndex j );

/**< The first form returns a read-only pointer to a n-vector containing the
   node deficits relative to commodity k; the second form returns the deficit
   of node j relative to commodity k.

   For k == NrComm() and k == NrComm() + 1, the lower and upper ranges of the
   "extra" constraints are returned, respectively: it is illegal to call
   these methods with k >= NrComm() if there are no extra constraints, i.e.
   NrExtraConst() == 0. */

/*--------------------------------------------------------------------------*/

   inline cIndex_Set StartN( void );

   inline Index StartNJ( cIndex j );

/**< The first form returns a read-only pointer to a m-vector containing the
   starting nodes of the arcs; the second form returns the starting node of
   arc j. */

/*--------------------------------------------------------------------------*/

   inline cIndex_Set EndN( void );

   inline Index EndNJ( cIndex j );

/**< The first form returns a read-only pointer to a m-vector containing the
   ending nodes of the arcs; the second form returns the ending node of arc j.
   */

/*--------------------------------------------------------------------------*/

   inline Index NIntVar( cIndex k = Inf<Index>() );

   inline cIndex_Set WIntVar( cIndex k );

/**< Some of the variables of the MMCF-like problem may be constrained to be
   integer-valued. NIntVar( k ) returns the number of the variables for
   commodity k that are constrained to be integer-valued. If 0 < NIntVar( k )
   < NrArcs(), WIntVar( k ) returns a read-only pointer to the vector of
   indices (ordered in increasing sense and Inf<Index>()-terminated) of the
   flow variables for commodity k that are integer-valued. If
   NIntVar( k ) == 0 or NIntVar( k ) == NrArcs(), then WIntVar( k ) returns 0;
   none/all the variables for commodity k are integer-valued.

   For k == NrComm(), these methods provide the same information for the
   "extra" (non-flow) variables; of course, the indices in the vector
   returned by WIntVar( NrComm() ) must be in the range [0, NrExtraVars()).

   NIntVar( k ) for k > NrComm() returns the *total* number of variables of
   the problem which are constrained to be integer-valued, counting both flow
   and extra variables. */

/*--------------------------------------------------------------------------*/

   virtual void ExtraConstr( int *IBeg , int *Indx , double *Vals ,
                             Index FrstC = 0 , Index LstC = Inf<Index>() );

/**< Writes in IBeg, Indx and Vals the description of those "extra" linear
   constraints of the problem whose "names" (indices) are comprised between
   FrstC and LstC (see NrExtraNonZ() above).

   The description is constraint-wise: each constraint is represented by the
   set of indices of variables with nonzero coefficient and the corresponding
   coefficients. The indices and coefficients corresponding to the i-th
   constraint that is returned, i = 0, ..., LstC - FrstC, are written in
   Indx and Vals, respectively, in the positions between IBeg[ i ]
   (included) and IBeg[ i + 1 ] (excluded). Thus, IBeg[ LstC - FrstC + 1 ]
   is also written, its content being the index of the first "free"
   position in Indx and Vals (the indices and values of constraint
   LstC - FrstC go up to IBeg[ LstC - FrstC + 1 ] - 1).

   The indices corresponding to each constraint, written in Indx, are
   ordered in increasing sense (and clearly without duplications).

   The mapping between indices and variables is the following: the variable
   corresponding to arc j (0 <= j < m) for commodity k (0 <= k < NrComm())
   has the index k * m + j; the i-th "extra" variable (0 <= i < NrExtraVars())
   has the index K * m + i (the representation is "commodity-wise", with the
   "extra" variables following). */

/*--------------------------------------------------------------------------*/
/*----------------- METHODS THAT ALLOW TO CHANGE THE GRAPH -----------------*/
/*--------------------------------------------------------------------------*/
/*--                                                                      --*/
/*-- These methods can be used to change the data of the instance w.r.t.  --*/
/*-- the original data set up in the constructors: they *must* be called  --*/
/*-- *before PreProcess()*.                                               --*/
/*--                                                                      --*/
/*--------------------------------------------------------------------------*/

   inline void ProblemType( const MCFType NewType , cIndex k );

/**< Sets the type of the k-th subproblem. */

/*--------------------------------------------------------------------------*/

   inline void Directed( bool DG );

/**< Decides if the graph G of the problem has to be considered as a directed
   (DG == true) or an undirected (DG == false) graph, setting the value that
   is returned by Directed( void ) [see above].
   If this method is not called, true (a directed graph) is assumed. */

/*--------------------------------------------------------------------------*/

   inline void NamesStartFrom1( bool NSF1 );

/**< Decides if node "names" (in the arc description) go from 1 to n or from
   0 to n - 1, setting the value that is returned by NamesStartFrom1( void )
   [see above]. If this method is not called, true ("names" starting from 1)
   is assumed.

   Note that no check is done to verify if the data of the instance actually
   corresponds to this setting, nor node "names" are scaled if any of these
   methods are called: this is user's responsibility. */

/*--------------------------------------------------------------------------*/

   virtual void UpDtTotCap( cFRow NewU = 0 );

   virtual void UpDtTotCapJ( cIndex j , cFNumber NewUj = Inf<FNumber>() );

/**< Updates the total upper capacities, either all of them or just the j-th,
   i.e. the one corresponding to the j-th arc. 0 means "all Inf<FNumber>()". */

/*--------------------------------------------------------------------------*/

   virtual void UpDtArcCstK( cIndex k , cCRow NewCk = 0 );

   virtual void UpdateArcCstKJ( cIndex k , cIndex j , cCNumber NewCkj = 0 );

/**< Updates the arc costs, either all of the commodity k or just the one
   relative to j-th arc for commodity k. 0 means "all 0".

   For k == NrComm(), the costs of the "extra" (non-flow) variables are
   changed: it is illegal to call these methods with k >= NrComm() if there
   are no extra variables, i.e. NrExtraVars() == 0. */

/*--------------------------------------------------------------------------*/

   virtual void UpDtArcCapK( cIndex k , cFRow NewUk = 0 );

   virtual void UpDtArcCapKJ( cIndex k , cIndex j , cFNumber NewUkj = 0 );

/**< Update the single-commodity capacities, either all of the commodity k or
   just the one relative to the j-th arc for commodity k. 0 means "all
   Inf<FNumber>()".

   For k == NrComm() and k == NrComm() + 1, the lower and upper bounds on the
   "extra" (non-flow) variables are changed, respectively: it is illegal to
   call these methods with k >= NrComm() if there are no extra variables, i.e.
   NrExtraVars() == 0. */

/*--------------------------------------------------------------------------*/

   virtual void UpDtNdeDfctK( cIndex k , cFRow NewDk = 0 );

   virtual void UpDtNdeDfctKJ( cIndex k , cIndex j , cFNumber NewDkj = 0 );

/**< Update the node deficits, either all of the commodity k or just the one
   relative to the j-th arc for commodity k. 0 means "all 0".

   For k == NrComm() and k == NrComm() + 1, the lower and upper ranges of the
   "extra" constraints are changed, respectively: it is illegal to call these
   methods with k >= NrComm() if there are no extra constraints, i.e.
   NrExtraConst() == 0. */

/*--------------------------------------------------------------------------*/

   virtual void SetIntVar( cIndex k , bool IntVld = true ,
                           cIndex_Set nms = 0 , cIndex strt = 0 ,
                           Index stp = Inf<Index>() );

/**< Gives the Graph object the information about which among the (flow and
   non-flow) variables of the problem are integer-valued.

   SetIntVar( k , true/false , ... ) says that some of the variables of the
   commodity k are/aren't integer-valued. The variables to which the change
   is applied are the flow variables of commodity k whose index is:
   - in the set nms[] (which contains indices in [0, NrArcs() ), ordered in
     increasing sense and Inf<Index>()-terminated);
   - comprised between strt (included) and min( stp , NrArcs() ) (excluded).
   nms == 0 means "all in the interval [strt, stp)".

   The status of all other variables remains unchanged.

   SetIntVar( NrComm() , ... ) can be used to set the integrality of the
   "extra" (non-flow) variables; of course, in this case the indices must
   be in the range [0, NrExtraVars()).

   The original status of the variable when the object is constructed depends
   on the constructor used and on the arguments passed to it [see above]. */

/*--------------------------------------------------------------------------*/

   virtual void SetExtraVars( cIndex NXV );

/**< Tells the Graph object that there are NXV "extra" (non-flow) variables in
   the problem.

   The costs and lower/upper bounds on these new variables can be set with
   the methods UpDtArcCstK[J]( NComm , ... ) and UpDtArcCapK[J]( NComm /
   NComm + 1 , ... ) [see above]; they are set by default respectively to
   0, 0 and Inf<FNumber>() uniformly. Also; by default all extra variables are
   continuous, although this can be changed later with SetIntVar().

   Note that every previous information about extra variables is lost when
   calling this method, that can be called more than once. Also, a call to
   SetExtraVars( 0 ) deallocates all the memory reserved for "extra"
   variables information. */

/*--------------------------------------------------------------------------*/

   virtual void SetExtraConstr( cIndex NXC , int *IBeg , int *Indx ,
				double *Vals );

/**< Gives to the Graph object a description of the "extra" constraints in the
   problem: NXC is the number of such constraints, and IStp, Indx and Vals
   must contain the description.

   The description is constraint-wise: each constraint is represented by the
   set of indices of variables with nonzero coefficient and the corresponding
   coefficients. The indices and coefficients corresponding to the i-th
   constraint, i = 0, ..., NXC - 1, are to be found in Indx and Vals,
   respectively, in the positions between IBeg[ i ] (included) and
   IBeg[ i + 1 ] (excluded). Thus, IBeg[ NXC ] must also be provided that
   contains the index of the first "free" position in Indx and Vals (the
   indices and values of constraint NXC - 1 go up to IBeg[ NXC ] - 1).

   The indices corresponding to each constraint, in Indx, must be ordered
   in increasing sense (and clearly without duplications).

   The mapping between indices and variables is the following: the variable
   corresponding to arc j (0 <= j < m) for commodity k (0 <= k < NrComm())
   has the index k * m + j; the i-th "extra" variable (0 <= i < NrExtraVars())
   has the index K * m + i (the representation is "commodity-wise", with the
   "extra" variables following).

   This is the very same format that is returned by ExtraConstr() [see above].
   Actually, ExtraConstr() (at least, in the implementation of the base Graph
   class) will return exactly the pointers passed to this method, that will
   be retained inside the Graph object. That is, after the call the vectors
   become "property" of the Graph object and should not be changed. This
   method can be called more than once: at each call, the previous pointers
   (if any) are *deallocated* and substituted with the newly provided ones.
   Also, the pointers are deallocated in the destructor of Graph. A call to
   SetExtraConstr( 0 , ... ) deallocates all the memory reserved for extra
   constraints information.

   The lower/upper bounds on the new constraints can be set with the method
   UpDtNdeDfctK[J]( NComm / NComm + 1 , ... ) [see above]: they are set by
   default respectively to 0 and Inf<FNumber>() uniformly.

   Note that no check is done in SetExtraConstr() about the validity of the
   data contained in the provided vectors (the indices being within the
   ranges and properly ordered). However, PreProcess() may use the extra
   constraints for its purposes, and/or modify them (e.g. by removing
   all references to variables that are declared non-existent). */

/*--------------------------------------------------------------------------*/
/*------------------------ PREPROCESSING METHODS ---------------------------*/
/*--------------------------------------------------------------------------*/

   virtual void PreProcess( cFNumber IncUk = 0 , cFNumber DecUk = 0 ,
                            cFNumber IncUjk = 0 , cFNumber DecUjk = 0 ,
                            cFNumber ChgDfct = 0 , cCNumber DecCsts = 0 );

/**< Performs various pre-processing of the data, trying to make the instance
   more easily solvable. The parameters to be given are the following:

   IncUk , DecUk   => (>= 0) upper bounds on the increase and decrease of the
                      mutual capacities: may be Inf<FNumber>() if unknown;

   IncUjk , DecUjk => (>= 0) same as above for single-commodity capacities;

   ChgDfct         => (>= 0) upper bound on the maximum change, in absolute
                      value, of the node deficits: it must be a finite number,
                      since it is used to generate "loose" but finite
                      individual capacities for arcs that have none;

   DecCsts         => (>= 0) upper bound on the decrease of arc Costs: must
                      be < Inf<CNumber>().

   Giving tight bounds (0 is the best, obviously) may cause the preprocessor
   to find more redundant coupling constraints, to squeeze down individual
   arc capacities, to remove more unused arcs and in general to do a better
   preprocessing; for instance, IncUjk == 0 allows PreProcess() to declare
   un-existent (set the cost to Inf<CNumber>()) any arc with 0 individual
   capacity.

   For all k such that, after the pre-processing, the graph has only a source
   and no (existing) arcs have a "real" capacity, the type of the subproblem
   is set to kSPT: all other problem types (see ProblemType( k ) below) are
   left unchanged.

   Important note: in order for PreProcess() to work, it has to be able to
   guess at least an upper bound on the maximum quantity of each commodity
   in the graph. In order to do that, *all arcs* with potentially *negative
   costs* (ChgCsts is used to estimate that) must have a *finite capacity*.

   PreProcess() will also look for redundancy in the data structures (e.g.
   identical costs/deficits/individual capacities for some commodities) and
   eliminate them, thus possibly saving some memory.

   It can be called *only once*. */
   
/*--------------------------------------------------------------------------*/

   virtual void MakeSingleSourced( bool ToAll = false );

/**< This method scans the graph and makes all the commodities single-sourced
   by adding, if necessary, one node to the graph and connecting it to all
   sources with arcs having individual capacities equal to the capacity of

   each source, no mutual capacity and zero cost. The sources then receive 0
   imbalance and the super-source the sum of all their original imbalances.
   The "name" of the source will be the number of nodes in the graph (plus 1
   if NamesStartFrom1() [see above]).

   This is obviously done only if at least one commodity has more than one
   source, unless ToALL == true: in this case, the super-source is constructed
   anyway and it is linked to *all* nodes i with arcs

     / at 0 cost and capacity - B[ i ][ k ]     if B[ i ][ k ] < 0
     |
     \ at Inf<CNumber>() cost and capacity 0    if B[ i ][ k ] >= 0.

   Hence, in this case exactly n new arcs will be constructed, but the
   setting is "resistent" to changes in the imbalances vector: if ToAll if
   false instead, possibly less then n arcs will be constructed.

   After the call, ProblemType( k ) will return kSPT for all k such that
   the graph (before the call) had at least one source: problems with no
   sources (circulation problems) cannot be converted to SPTs, and their
   type is left unchanged (being set to kMCF by the constructor).

   This method can be called only once, and *before* PreProcess(). */

/*--------------------------------------------------------------------------*/
/*---------------------- MISCELLANEOUS METHODS -----------------------------*/
/*--------------------------------------------------------------------------*/

   virtual void OutMPSFile( const char *const FN , bool FxdClmns = false ,
			    const char *const PN = 0 );

/**< Outputs the current MMCF problem as a Linear Programming problem in the
   standard MPS format: if FxdClmns == true, then the fixed-column MPS format
   is used, else (the default) the "modern" MPS format is used, where fields
   have not a fixed column positions but are just blank-delimited.
   `FN' is taken as the [path]name of the output file, and `PN' (if provided)
   is taken as the name of the problem, to be written in the NAME field of the
   MPS file. */

/*--------------------------------------------------------------------------*/

   virtual FONumber UpperBound( void );

/**< Calculates a (very coarse) Upper Bound to the value of the optimal
   solution of the MMCF currently represented in the Graph object. This is
   useful e.g. to discover if the MMCF has no solutions when using a
   decomposition-based algorithm; in fact, in this case one should wait for
   the Lagrangean to go to +INF to declare the problem unfeasible, but any
   finite UB to the value of the optimal solution can be used as +INF. Note,
   however, that the upper bound returned by this method *can be tight*,
   i.e., it can be actually the value of one feasible solution (in
   particular, it surely is in case of a feasibility problem, i.e., all
   costs equal to zero, which is actually feasible; in this case the method
   will return zero).

   If "extra" (non-flow) variables are defined, they are taken into account
   when calculating the upper bound.

   It is better to call UpperBound() *after* PreProcess() [see below]: a
   tighter bound is returned. */

/*--------------------------------------------------------------------------*/
/*------------------------------ DESTRUCTOR --------------------------------*/
/*--------------------------------------------------------------------------*/

  virtual ~Graph();  /// Destructor

/*--------------------------------------------------------------------------*/
/*-------------------- PROTECTED PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/
/*--                                                                      --*/
/*--  The standard user should not need to care about the following part: --*/
/*--  programmers who need to extend the code (i.e. by deriving a new     --*/
/*--  class) may make use of the following methods and data structures.   --*/
/*--                                                                      --*/
/*--  IT IS OBVIOUSLY DANGEROUS TO MODIFY THE DATA STRUCTURES, while it   --*/
/*--  safe to read them.                                                  --*/
/*--                                                                      --*/
/*--------------------------------------------------------------------------*/

 protected:

/*--------------------------------------------------------------------------*/
/*-------------------------- PROTECTED METHODS -----------------------------*/
/*--------------------------------------------------------------------------*/
 
   inline void CmnIntlz( void );

/**< Called at the end of any constructor, does some initializations that are
   common to them all: it is "protected" for allowing derived classes that
   use the "void" constructor to call it. */

/*--------------------------------------------------------------------------*/
/*--------------------- PROTECTED DATA STRUCTURES --------------------------*/
/*--------------------------------------------------------------------------*/

  Index NNodes;        ///< Number of nodes
  Index NArcs;         ///< Number of arcs
  Index NComm;         ///< Number of commodities
  Index NCnst;         ///< Number of arcs with mutual capacity constraints
  Index NXtrV;         ///< Number of "extra" variables
  Index NXtrC;         ///< Number of "extra" constraints
  Index_Set Startn;    ///< Topology of the graph: starting nodes
  Index_Set Endn;      ///< Topology of the graph: ending nodes
  CMat C;              ///< Matrix of the arc costs
  FMat U;              ///< Matrix of the arc upper capacities
  FMat B;              ///< Matrix of the node deficits
  MCFType *PT;         ///< type of flow subproblem
  FRow UTot;           ///< Vector of mutual capacities
  Index_Set Active;    ///< Set of the arcs for which a mutual capacity
                       ///< constraint is defined
  Index_Mat ActiveK;   ///< Like Active for individual capacities
  Index_Set NamesK;    ///< The dual multipliers relative to commodity K
                       ///< start with NamesK[ k ] and end with
                       ///< NamesK[ k + 1 ]
  Index StrtNme;       ///< The "name" of the first node
  bool DrctdPrb;       ///< true if the problem is directed
  Index_Set NInt;      ///< Number of integer-valued variables
  Index_Mat WIsInt;    ///< Which of the variables are integer-valued
  int *IdxBeg;         ///< Description of "extra" constraints: start
  int *CoefIdx;        ///< Description of "extra" constraints: indices
  double *CoefVal;     ///< Description of "extra" constraints: values
  
  Bool_Vec CIsCpy;     ///< true for each row of C[] that is a copy of another
  Bool_Vec UIsCpy;     ///< true for each row of U[] that is a copy of another
  Bool_Vec BIsCpy;     ///< true for each row of B[] that is a copy of another
  Bool_Vec DIsCpy;     ///< true for each row of D[] that is a copy of another

/*--------------------------------------------------------------------------*/

 };  // end( class Graph )

/*--------------------------------------------------------------------------*/
/*------------------- inline methods implementation ------------------------*/
/*--------------------------------------------------------------------------*/
/*-------------- METHODS FOR READING THE DATA OF THE PROBLEM ---------------*/
/*--------------------------------------------------------------------------*/

inline Graph::Index Graph::NrComm( void )
{
 return( NComm );
 }

/*--------------------------------------------------------------------------*/

inline Graph::Index Graph::NrArcs( void )
{
 return( NArcs );
 }

/*--------------------------------------------------------------------------*/

inline Graph::Index Graph::NrNodes( void )
{
 return( NNodes );
 }

/*--------------------------------------------------------------------------*/

inline Graph::Index Graph::NrExtraVars( void )
{
 return( NXtrV );
 }

/*--------------------------------------------------------------------------*/

inline Graph::Index Graph::NrExtraConst( void )
{
 return( NXtrC );
 }

/*- - - - - - - - - - - - - - -- - - - - - - - - - - - - - - - - - - - - -*/

inline Graph::Index Graph::NrExtraNonZ( cIndex FrstC , cIndex LstC )
{
 if( NXtrC )
  return( IdxBeg[ ( LstC >= NXtrC ? NXtrC : LstC ) ] -
	  IdxBeg[ ( FrstC >= 0 ? FrstC : 0 ) ] );
 else
  return( 0 );
 }

/*--------------------------------------------------------------------------*/

inline Graph::MCFType Graph::ProblemType( cIndex k )
{
 return( PT[ k ] );
 }

/*--------------------------------------------------------------------------*/

inline bool Graph::Directed( void )
{
 return( DrctdPrb );
 }

/*--------------------------------------------------------------------------*/

inline bool Graph::NamesStartFrom1( void )
{
 if( StrtNme )
  return( true );
 else
  return( false );
 }

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

inline Graph::cFRow Graph::TotCapacities( void )
{
 return( UTot );
 }

/*- - - - - - - - - - - - - - - - - - -- - - - - - - - - - - - - - - - - - -*/

inline Graph::FNumber Graph::TotalCapacityJ( cIndex j )
{
 return( UTot[ j ] );
 }

/*--------------------------------------------------------------------------*/

inline Graph::cIndex_Set Graph::Actives( void )
{
 return( Active );
 }

/*- - - - - - - - - - - - - - -- - - - - - - - - - - - - - - - - - - - - -*/

inline Graph::Index Graph::NActives( void )
{
 return( NCnst );
 }

/*--------------------------------------------------------------------------*/

inline Graph::cCRow Graph::CostsK( cIndex k )
{
 return( C[ k ] );
 }

/*- - - - - - - - - - - - - - - - - - -- - - - - - - - - - - - - - - - - - -*/

inline Graph::CNumber Graph::CostKJ( cIndex k , cIndex j )
{
 return( C[ k ][ j ] );
 }

/*--------------------------------------------------------------------------*/

inline Graph::cFRow Graph::CapacitiesK( cIndex k )
{
 return( U[ k ] );
 }

/*- - - - - - - - - - - - - - - - - - -- - - - - - - - - - - - - - - - - - -*/

inline Graph::FNumber Graph::CapacityKJ( cIndex k , cIndex j )
{
 return( U[ k ][ j ] );
 }

/*--------------------------------------------------------------------------*/

inline Graph::cIndex_Set Graph::ActivesK( cIndex k )
{
 return( ActiveK[ k ] );
 }

/*- - - - - - - - - - - - - - -- - - - - - - - - - - - - - - - - - - - - -*/

inline Graph::Index Graph::NActivesK( cIndex k )
{
 if( k >= NComm )
  return( NamesK[ NComm ] - *NamesK );
 else
  return( NamesK[ k + 1 ] - NamesK[ k ] );
 }

/*--------------------------------------------------------------------------*/

inline Graph::cFRow Graph::DeficitsK( cIndex k )
{
 return( B[ k ] );
 }

/*- - - - - - - - - - - - - - - - - - -- - - - - - - - - - - - - - - - - - -*/

inline Graph::FNumber Graph::DeficitKJ( cIndex k , cIndex j )
{
 return( B[ k ][ j ] );
 }

/*--------------------------------------------------------------------------*/

inline Graph::cIndex_Set Graph::StartN( void )
{
 return( Startn );
 }

/*- - - - - - - - - - - - - - - - - - -- - - - - - - - - - - - - - - - - - -*/

inline Graph::Index Graph::StartNJ( cIndex j )
{
 return( Startn[ j ] );
 }

/*--------------------------------------------------------------------------*/

inline Graph::cIndex_Set Graph::EndN( void )
{
 return( Endn );
 }

/*- - - - - - - - - - - - - - - - - - -- - - - - - - - - - - - - - - - - - -*/

inline Graph::Index Graph::EndNJ( cIndex j )
{
 return( Endn[ j ] );
 }

/*--------------------------------------------------------------------------*/

inline Graph::Index Graph::NIntVar( cIndex k )
{
 if( k <= NComm )
  return( NInt[ k ] );
 else {
  register Index cnt = 0;
  for( register Index i = 0 ; i <= NComm ; )
   cnt += NInt[ i++ ];

  return( cnt );
  }
 }

/*- - - - - - - - - - - - - - - - - - -- - - - - - - - - - - - - - - - - - -*/

inline Graph::cIndex_Set Graph::WIntVar( cIndex k )
{
 return( WIsInt[ k ] );
 }

/*--------------------------------------------------------------------------*/
/*----------------- METHODS THAT ALLOW TO CHANGE THE GRAPH -----------------*/
/*--------------------------------------------------------------------------*/

inline void Graph::ProblemType( const Graph::MCFType NewType , cIndex k )
{
 PT[ k ] = NewType;
 }

/*--------------------------------------------------------------------------*/

inline void Graph::Directed( bool DG )
{
 DrctdPrb = DG;
 }

/*--------------------------------------------------------------------------*/

inline void Graph::NamesStartFrom1( bool NSF1 )
{
 if( NSF1 )
  StrtNme = 1;
 else
  StrtNme = 0;
 }

/*--------------------------------------------------------------------------*/

#if( OPT_USE_NAMESPACES )
 };  // end( namespace MMCFGraph_di_unipi_it )
#endif

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif  /* Graph.h included */

/*--------------------------------------------------------------------------*/
/*-------------------------- End File Graph.h ------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*--------------------------- File MMCFFlwB.h ------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the class MMCFFlwBase, which solves relaxations of the
 * Multicommodity Min Cost Flow (MMCF) problem where the Mutual Capacity
 * constraints are relaxed, that is, k independent single-commodity Min Cost
 * Flow or Shortest Path subproblems. This is done by using an array of
 * k appropriate single-commodity Min Cost Flow (MCF) solvers derived from
 * the base class MCFClass [see MCFClass.h]. "extra" variables or constraints
 * are not supported. The class implements the MMCFClass interface, see
 * MMCFClas.h. 
 *
 * \version 3.01
 *
 * \date 11 - 05 - 2012
 *
 * \author Antonio Frangioni \n
 *         Operations Research Group \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * Copyright(C) 1992 - 2012 Antonio Frangioni
 */
/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __MMCFFlwB
 #define __MMCFFlwB  /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "MMCFClas.h"

#include "Graph.h"

#include "MCFClass.h"

/*--------------------------------------------------------------------------*/
/*-------------------------------- MACROS ----------------------------------*/
/*--------------------------------------------------------------------------*/

#define FlowBase_HAVE_SPT 1

/**< FlowBase_HAVE_SPT > 0 means that a Shortest Path Tree (SPT) solver is
   available; this allows "super-disaggregate" solutions whereby a subproblem
   is a specific O-D pair within a single (SPT) (possibly with multiple
   destinations for that single origin).

   FlowBase_HAVE_SPT == 1 means that the (SPT) solver is the SPTree class. */

/*--------------------------------------------------------------------------*/
/*------------------------ NAMESPACE and USINGS ----------------------------*/
/*--------------------------------------------------------------------------*/

#if( OPT_USE_NAMESPACES )
namespace MMCFClass_di_unipi_it
{
 using namespace MMCFGraph_di_unipi_it;
#endif

/*--------------------------------------------------------------------------*/
/*----------------------------- CONSTANTS ----------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*----------------------------- CLASSES ------------------------------------*/
/*--------------------------------------------------------------------------*/

class MMCFFlwBase : public MMCFClass
{

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*--------------------------------------------------------------------------*/
/*--------------------------- PUBLIC METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/
/*---------------------------- CONSTRUCTORS --------------------------------*/
/*--------------------------------------------------------------------------*/
/** Default constructor (with no arguments) of the class: gives some default
   values to the data structures of the class, that must be changed later in
   the constructors of the derived classes or by calling SetGraph(). */

   MMCFFlwBase( void ) : MMCFClass()
   {
    MCFs = 0;
    OptVal = 0;
    Cost = 0;
    Ftmp = 0;

    UBnd = Inf<FONumber>();

    Stts = kError;
    FAgg = 0;

    AFSPD = 0;

    #if( FlowBase_HAVE_SPT )
     DFSPD = 0;
     DFSPDn = 0;
     DFSPDb = 0;
    #endif

    #if( CHGARCS_MMCF > 1 )
     ClsdAK = 0;
    #endif
    }

/*--------------------------------------------------------------------------*/

   MMCFFlwBase( Graph *Gh );

/**< Constructor of the class, given the Graph object. */

/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

   void SetGraph( Graph *Gh );

/**< Sets the MMCF solver to work on the given graph. This is the alternative
   to the graph constructor when the empty constructor has been used. Of
   course, it must be called before any other method of the class. */

/*--------------------------------------------------------------------------*/

   void ReOptimize( bool RO = false , cIndex whch = Inf<Index>() );

/**< Tells whether or not the flow subproblems have to be Re-Optimized, i.e. if
   the latest optimal solution have to be (tentatively) used to warm-start the
   solvers. If RO = true this is done: it is usually a good idea (in fact, it
   is assumed if ReOptimize() is never called).

   This can be set for each commodity independently by passing a whch < k, or
   for all commodities together by passing whch == Inf<Index>(). */

/*--------------------------------------------------------------------------*/

   void SetNrSubP( Index FSP = 0 );

/**< Attempts to set the number of subproblems (as reported by NrSubP()) to
   FSP; the "special" setting FSP == 0 corresponds to NrSubP() == 1 (a
   non-decomposable problem), which is assumed if this method is never called.

   If FSP <= NComm, then the subproblems are (possibly groups of) the natural
   flow problems; if FSP > NComm, each one flow subproblem with SPT structure
   and multiple destinations can be disaggregated to a number of subproblems
   with the same origins and less destinations. This requires
   FlowBase_HAVE_SPT > 0. */

/*--------------------------------------------------------------------------*/

   virtual void SetOptEps( const double OE = 0 );

   virtual void SetFsbEps( const double FE = 0 );

/*--------------------------------------------------------------------------*/
/*-------------------- METHODS FOR SOLVING THE PROBLEM ---------------------*/
/*--------------------------------------------------------------------------*/

   MMCFStatus SolveMMCF( void );

/**< Solves the (k) subproblem(s). If a (MCF) solver reports an error (a
   return status > 0, depending on the specific (MCF) solver used and the
   error occurred), the calculation is immediately interrupted, and HpINF is
   returned by GetOptVal(). */

/*--------------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/

   FONumber GetPVal( void );

/*--------------------------------------------------------------------------*/

   inline FONumber GetDVal( void );

/**< GetDVal() returns the same value as GetPVal(), i.e., just the optimal
   value of the current problem. */

/*--------------------------------------------------------------------------*/

   inline FONumber GetUpprBnd( bool &HvSol );

/**< GetUpprBnd() returns the "infeasibility" Upper Bound calculated by Graph:
   if GetLwrBnd() >= GetUpprBnd(), then the problem is unfeasible. However,
   the Upper Bound computed by Graph become meaningless as soon as any
   Chg***() method is called: since the class does not have a "complete"
   picture of the problem, no attempt is done to re-compute a meaningful
   value. */

/*--------------------------------------------------------------------------*/

   inline FONumber GetLwrBnd( bool &HvSol );

/**< GetLwrBnd() returns the best Lower Bound available, i.e., the same value
   as GetDVal(). */

/*--------------------------------------------------------------------------*/

   bool GetPSol( void );

/*--------------------------------------------------------------------------*/

   bool PSolIsFNumber( void )
   {
    return( true );
    }

/*--------------------------------------------------------------------------*/

   bool GetDSol( void );

   inline bool GetLBSol( void );

/**< GetLBSol() does exactly the same as GetDSol() [see GetLwrBnd() above].

   Important note: the potentials may strongly depend on some decisions taken
   at the (MCF)s level. For instance, the potentials might be different in
   case arcs that are "non-existent" in a (MCF) are just dealt with as such,
   or rather "simulated" setting their capacity to 0. Some of the codes under
   the MCFClass interface might therefore provide diferent potentials by
   changing the value of the macro DYNAMIC. Since Reduced Costs are tied to
   node potentials, the same problem can occur there.

   Note that, obviously, the dual costs of Mutual Capacity constraints cannot
   be asked to this solver. */

/*--------------------------------------------------------------------------*/

   FONumber CostOf( void );

/**< This currently only works for primal solutions. */

/*--------------------------------------------------------------------------*/
/*-------------- METHODS FOR READING THE DATA OF THE PROBLEM ---------------*/
/*--------------------------------------------------------------------------*/

   void GetCosts( CRow Csts , cIndex_Set nms = 0 , cIndex strt = 0 ,
		  Index stp = Inf<Index>() );

/*--------------------------------------------------------------------------*/

   void GetICaps( FRow ICps , cIndex_Set nms = 0 , cIndex strt = 0 ,
		  Index stp = Inf<Index>() );

/*--------------------------------------------------------------------------*/
/** Mutual Capacities cannot be asked to this class. */

   void GetMCaps( FRow MCps , cIndex_Set nms = 0 , cIndex strt = 0 ,
		  Index stp = Inf<Index>() )
   {
    throw( MMCFException( "MMCFFlwBase::GetMCaps() should not be called" ) );
    }

/*--------------------------------------------------------------------------*/
/*------------- METHODS FOR ADDING / REMOVING / CHANGING DATA --------------*/
/*--------------------------------------------------------------------------*/

   void ChgCosts( cCRow NwCsts , cIndex_Set nms = 0 , cIndex strt = 0 ,
		  Index stp = Inf<Index>() );

/*--------------------------------------------------------------------------*/

   void ChgICaps( cFRow NwICps , cIndex_Set nms = 0 , cIndex strt = 0 ,
		  Index stp = Inf<Index>() );

/*--------------------------------------------------------------------------*/
/** Mutual Capacities do not exist, hence they cannot be changed. */

   void ChgMCaps( cFRow NwMCps , cIndex_Set nms = 0 , cIndex strt = 0 ,
		  Index stp = Inf<Index>() )
   {
    throw( MMCFException( "MMCFFlwBase::ChgMCaps() should not be called" ) );
    }

/*--------------------------------------------------------------------------*/

   void ChgIntVar( cIndex k = Inf<Index>() , bool IntVld = true ,
		   cIndex_Set nms = 0 , Index strt = 0 ,
		   Index stp = Inf<Index>() );

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
/*------------------------------ DESTRUCTOR --------------------------------*/
/*--------------------------------------------------------------------------*/

   virtual ~MMCFFlwBase();

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
/*-------------------------- PROTECTED TYPES -------------------------------*/
/*--------------------------------------------------------------------------*/

   typedef MCFClass *PMCF;      ///< pointer to a MCFClass (or derived from)
   typedef PMCF     *PPMCF;     ///< array of such poinetrs
   typedef bool     *Bool_Vec;  ///< vector of booleans

/*--------------------------------------------------------------------------*/
/*-------------------------- PROTECTED METHODS -----------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*----------------------- PROTECTED DATA STRUCTURES  -----------------------*/
/*--------------------------------------------------------------------------*/

  PPMCF  MCFs;        ///< The array of k flow solvers
  FONumber *OptVal;   ///< The Optimal Value(s)
  CRow Cost;          ///< Temporary for cost-based calculations
  FRow Ftmp;          ///< Temporary for flow-based calculations

  FONumber UBnd;      /**< The Upper Bound over the optimal value of the
			   (MMCF) (HpINF if no such bound is known)*/
  MMCFStatus Stts;    ///< Status of the optimization
  FRow FAgg;          ///< Aggregate Flow solution

  Index_Set AFSPD;    /**< In each "aggregated flow" subproblem, the
			   involved commodities are
			   AFSPD[ ws - 1 ] <= k < AFSPD[ ws ]*/
  #if( FlowBase_HAVE_SPT )
   Index_Set DFSPD;   /**< DFSPD[ ws - 1 ] is the commodity to which all
			   destinations of subproblem ws belong, for
			   ws = 1, ..., NSubPr*/
   Index_Set DFSPDn;  /**< DFSPDn[ ws - 1 ] is the number of destinations of
			   subproblem ws (== 0 ==> it is a MCF)*/
   Index_Set *DFSPDb; /**< DFSPDb[ ws - 1 ] is the set of destinations of
			   subproblem ws ( == 0 ==> it is a MCF)*/
  #endif

  #if( CHGARCS_MMCF > 1 )
   Bool_Vec ClsdAK;   ///< which arcs have been closed by CloseArc( k , ... )
  #endif

/*--------------------------------------------------------------------------*/
/*--------------------- PRIVATE PART OF THE CLASS --------------------------*/
/*--------------------------------------------------------------------------*/
/*--                                                                      --*/
/*-- Nobody should ever look at this part: everything that is under this  --*/
/*-- advice may be changed without notice in any new release of the code. --*/
/*--                                                                      --*/
/*--------------------------------------------------------------------------*/

 private:

/*--------------------------------------------------------------------------*/
/*-------------------------- PRIVATE METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/

   inline void SolveFlwK( cIndex k );

/*--------------------------------------------------------------------------*/

   inline void SolveSubP( cIndex wh );

/*--------------------------------------------------------------------------*/
/*----------------------- PRIVATE DATA STRUCTURES  -------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/

 };  // end( class MMCFFlwBase )

/*--------------------------------------------------------------------------*/
/*------------------- inline methods implementation ------------------------*/
/*--------------------------------------------------------------------------*/

inline MMCFClass::FONumber MMCFFlwBase::GetDVal( void )
{
 return( MMCFFlwBase::GetPVal() );
 }

/*--------------------------------------------------------------------------*/

inline MMCFClass::FONumber MMCFFlwBase::GetUpprBnd( bool &HvSol )
{
 HvSol = false;
 return( UBnd );
 }

/*--------------------------------------------------------------------------*/

inline MMCFClass::FONumber MMCFFlwBase::GetLwrBnd( bool &HvSol )
{
 HvSol = true;
 return( MMCFFlwBase::GetPVal() );
 }

/*--------------------------------------------------------------------------*/

inline bool MMCFFlwBase::GetLBSol( void )
{
 return( GetDSol() );
 }

/*--------------------------------------------------------------------------*/

#if( OPT_USE_NAMESPACES )
 };  // end( namespace MMCFClass_di_unipi_it )
#endif

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif  /* MMCFFlwB.h included */

/*--------------------------------------------------------------------------*/
/*------------------------End File MMCFFlwB.h ------------------------------*/
/*--------------------------------------------------------------------------*/

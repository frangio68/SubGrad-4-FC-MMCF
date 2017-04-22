/*--------------------------------------------------------------------------*/
/*------------------------- File KnpsFiOrcl.h ------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Definition of the KnpsFiOrcl class which implements the FiOracle interface
 * and solves the knapsack Lagrangian relaxation of the (Fixed-Charge)
 * Multicommodity Min Cost Flow Problem ((FC-)MMCF).
 *
 * \version 1.21
 *
 * \date 3 - 04 - 2013
 *
 * \author Antonio Frangioni \n
 *         Operations Research Group \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Enrico Gorgone \n
 *         Operations Research Group \n
 *         Dipartimento di Elettronica Informatica e Sistemistica \n
 *         Universita' della Calabria \n
 *
 * \author Vitor Barbosa \n
 *         Algoritmi Research Unit \n
 *         Universidade do Minho \n
 *
 * \author Filipe Alvelos \n
 *         Grupo de Optimizacao e Investigacao Operacional \n
 *         Departamento de producao e Sistemas \n
 *         Universidade do Minho \n
 *
 * Copyright &copy 2000 - 2012 by Antonio Frangioni, Enrico Gorgone
 */
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __KnpsFiOrcl
 #define __KnpsFiOrcl /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "FiOracle.h"

#include "MMCFClas.h"
#include "Graph.h"

/*--------------------------------------------------------------------------*/
/*-------------------------------- MACROS ----------------------------------*/
/*--------------------------------------------------------------------------*/

#define LOG_FI 0
/**< If LOG_FI > 0, the FlwFiOrcl class may produce (depending on an input
   parameter) a log of (more or less) important results. */

#define LowB_FI 1
/** If LOG_FI > 0, the FlwFiOrcl class using the global Lower Bound.*/

#define SPF_SUB 0
/**< indicate how to manage the subgradients in the method Gi():

    0: do not sparsify the subgradients;
    1: sparsify the sungradients. */

/*--------------------------------------------------------------------------*/
/*------------------------------ NAMESPACE ---------------------------------*/
/*--------------------------------------------------------------------------*/

#if( OPT_USE_NAMESPACES )
namespace NDO_di_unipi_it
{
 using namespace MMCFGraph_di_unipi_it;
#endif

/*--------------------------------------------------------------------------*/
/*----------------------------- CLASS KnpsFiOrcl ---------------------------*/
/*--------------------------------------------------------------------------*/
/// Instantiation of FiOracle for Knapsack Lagrangian relaxation of FC-MMFC
/** This class instantiates the FiOracle interface [see FiOracle.h] to solve
 * the knapsack Lagrangian relaxation of the (Fixed-Charge) Multicommodity Min
 * Cost Flow Problem ((FC-)MMCF). The problem formulation is below:
 * \f[
 *  \phi(\upsilon)=
 *   \min \sum_{k \in K} \sum_{(i,j) \in A}
 *        \left( c_{ij}^k + \upsilon_i^k - \upsilon_j^k \right) x_{ij}^k +
 *        \sum_{(i,j) \in A} f_{ij}y_{ij} -
 *        \sum_{n \in N} \sum_{k \in K}
 *        \left( \upsilon_i^k - b_i^k \right)
 * \f]
 * \f[
 *   \sum_{k \in K} x_{ij}^k \leq u_{ij} y_{ij} \quad,\quad  (i,j) \in A
 * \f]
 * \f[
 *   0 \leq x_{ij}^k \leq u_{ij}^k y_{ij} \quad,\quad (i,j) \in A, k \in K
 * \f]
 * \f[
 *   y_{ij} \in \{0,1\} \quad,\quad  (i,j) \in A
 * \f]
 * The \f$\upsilon\f$ vector is the Lambda values described in the FiOracle
 * interface. They are dual variables associated to the relaxed flow
 * conservation constraints. */ 

class KnpsFiOrcl : public FiOracle , public MMCFClass
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
/** @name Constructor
   @{*/

   KnpsFiOrcl( Graph *g , istream *iStrm = 0 );

/**< Constructor of the class. The parameter `iStrm', if provided, is taken as
   a pointer to a istream from which the algorithmic parameters for the knapsack
   relaxation are sequentially read in the following order. Each parameter
   must be placed at the beginning of a separate line, max 255 characters
   long, with all the rest of the line up to the first newline character
   (apart from a separating whitespace) being available for comments.
   Any line whose first character is '#' and any blank line is ignored.
   If 0 is passed, the file ends before reaching a given parameter, or
   some parameter is in the wrong format, each non-specified parameter is
   given a default value, shown in [] below.

   `iStrm' is passed to the constructor of FiOracle [see FiOracle.h], which
   reads the general algorithmic parameters out of it; since the constructor
   KnpsFiOrcl is executed after the one of FiOracle, the following parameters
   specific for the SubGrad have to be found in the stream after those of the
   base class:

      -# agg [0]   1 if the function is considered as one function,
                   0 if decomposed

      -# Eps  [1e-6] precision of Fi()

      -# DoAggregation [0] 1 if the variables generation is adopted, 0
                          otherwise.                                       */

/*@}*-----------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

void SetFiLog( ostream *outs = 0 , const char lvl = 0 );

/**< lvl controls the "level of verbosity" of the code. The first four bits
   of lvl have the following meaning:

   0  =>  no log at all (also assumed if log = 0);

   >0  =>  "basic" log: only the errors are reported; */

/*--------------------------------------------------------------------------*/

   void SetMaxName( cIndex MxNme = 0 );

/*--------------------------------------------------------------------------*/

   void SetAggregate( bool agg = true );

/**< Set aggregation option. Tells if the FiOracle works with Fi aggregated
    or disaggregated. */

/*@}------------------------------------------------------------------------*/
/*-------------- METHODS FOR READING THE DATA OF THE PROBLEM ---------------*/
/*--------------------------------------------------------------------------*/

   Index GetNrFi( void ) const;

/*--------------------------------------------------------------------------*/
/*----------------------- METHODS FOR CHANGING DATA ------------------------*/
/*--------------------------------------------------------------------------*/
   
   void SetLambda( cLMRow Lmbd = 0 );

/*--------------------------------------------------------------------------*/

   bool SetPrecision( HpNum Eps );

/*--------------------------------------------------------------------------*/
/*-------------------- METHODS FOR SOLVING THE PROBLEM ---------------------*/
/*--------------------------------------------------------------------------*/

   HpNum Fi( cIndex wFi = Inf<Index>() );

/*--------------------------------------------------------------------------*/

   MMCFStatus SolveMMCF( void ) {

	return kOK;
    }

/*--------------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/

   bool NewGi( cIndex wFi = Inf<Index>() );

/**< The subgradients are calculated using the following formula:
    \f[ gi(\upsilon_i^k) = \sum_{(i,j) \in A} \left( x_{ij}^k -
    \sum_{(j,i)\in A}x_{ji}^k  \right) - b_i^k \f]. */

/*--------------------------------------------------------------------------*/

   Index GetGi( SgRow SubG , cIndex_Set &SGBse , cIndex Name = Inf<Index>() ,
		cIndex strt = 0 , Index stp = Inf<Index>() );

/*--------------------------------------------------------------------------*/

   HpNum GetVal( cIndex Name = Inf<Index>() );

/*--------------------------------------------------------------------------*/

   void SetGiName( cIndex Name );

/*--------------------------------------------------------------------------*/

   FONumber GetPVal( void ) { }
   bool GetPSol( void ) { }

/*--------------------------------------------------------------------------*/

   FONumber GetDVal( void ) { }
   bool GetDSol( void ) { }

/*--------------------------------------------------------------------------*/

   FONumber CostOf( void ) { }

/*--------------------------------------------------------------------------*/
/*-------------------- METHODS FOR READING OTHER RESULTS -------------------*/
/*--------------------------------------------------------------------------*/
/** @name Reading other results
   @{ */

   HpNum GetLowerBound( cIndex wFi = Inf<Index>() );
   void SetLowerBound( HpNum lowB = Inf<HpNum>() );

/*--------------------------------------------------------------------------*/

   FiStatus GetFiStatus( Index wFi = Inf<Index>() );

/* @}-----------------------------------------------------------------------*/
/*-------------- METHODS FOR READING THE DATA OF THE PROBLEM ---------------*/
/*--------------------------------------------------------------------------*/
/** @name Reading the data of the problem
    @{ */

   void GetCosts( CRow Csts , cIndex_Set nms = 0 , cIndex strt = 0 ,
		   Index stp = Inf<Index>() ) { }

   void GetICaps( FRow ICps , cIndex_Set nms = 0 , cIndex strt = 0 ,
		   Index stp = Inf<Index>() ) { }

   void GetMCaps( FRow MCps , cIndex_Set nms = 0 , cIndex strt = 0 ,
		   Index stp = Inf<Index>() ) { }

/*--------------------------------------------------------------------------*/

   HpNum GetGlobalLipschitz( cIndex wFi = Inf<Index>() );

/*@}------------------------------------------------------------------------*/
/*------------- METHODS FOR ADDING / REMOVING / CHANGING DATA --------------*/
/*--------------------------------------------------------------------------*/
/** @name Adding / removing / changing data
   @{ */

   void Deleted( cIndex i = Inf<Index>() );

/*--------------------------------------------------------------------------*/

   void Aggregate( cHpRow Mlt , cIndex_Set NmSt , cIndex Dm ,cIndex NwNm );

/*--------------------------------------------------------------------------*/

   void ChgCosts( cCRow NwCsts , cIndex_Set nms = 0 , cIndex strt = 0 ,
		   Index stp = Inf<Index>() ) { }


   void ChgICaps( cFRow NwICps , cIndex_Set nms = 0 , cIndex strt = 0 ,
		   Index stp = Inf<Index>() ) { }


   void ChgMCaps( cFRow NwMCps , cIndex_Set nms = 0 , cIndex strt = 0 ,
		   Index stp = Inf<Index>() ) { }


   void ChgIntVar( cIndex k = Inf<Index>() , bool IntVld = true ,
           cIndex_Set nms = 0 , Index strt = 0 , Index stp = Inf<Index>() ) { }

/*--------------------------------------------------------------------------*/

#if( CHGARCS_MMCF )

   void CloseArcs( cIndex_Set whch ) { }

   void OpenArcs( cIndex_Set whch ) { }

 #if( CHGARCS_MMCF > 1 )

   virtual void CloseArcs( cIndex k , cIndex_Set whch ) { }

   virtual void OpenArcs( cIndex k , cIndex_Set whch ) { }

 #endif
#endif

/*@} -----------------------------------------------------------------------*/
/*------------------------------ DESTRUCTOR --------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Destructor
   @{ */

   virtual ~KnpsFiOrcl();

/*@}------------------------------------------------------------------------*/
/*-------------------- PROTECTED PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/
/*--                                                                      --*/
/*--  The standard user should not care about the following part: users   --*/
/*--  who need to extend the code by deriving a new class may use these   --*/
/*--  methods and data structures. It is *dangerous* to *modify* the      --*/
/*--  data structures, while it safe to read them.                        --*/
/*--                                                                      --*/
/*--------------------------------------------------------------------------*/

 protected:

/*--------------------------------------------------------------------------*/
/*-------------------------- PROTECTED METHODS -----------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*----------------------- PROTECTED DATA STRUCTURES  -----------------------*/
/*--------------------------------------------------------------------------*/
/** @name Algorithmic parameters
    @{*/

  HpNum Eps;        ///<Precision for some calculations.
  bool Aggrgtd;     ///<Agreggation flag.

/*@}------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Instance data
    @{ */

  Index_Set Startn; ///<Topology of the graph: starting ...
  Index_Set Endn;   ///<... and ending node of each arc

  CRow Costs;   ///< Array of unit costs.
                    /**< Each entry corresponds to the
                     cost of one unit of flow of one commodity
                    traversing one arc. For arc j and commodity k, this
                    cost is OrigCosts[ k * m + j ], where m stands for
                    the total number of arcs.*/
  CRow XtrCosts;   ///< Array of design costs.
                       /**<Each entry is associated with the fixed cost
                       of the corresponding arc being used.*/
  CRow Capacities; ///< Array of individual capacities.
                       /**<Each entry OrigCapacities[ k * m + j ] corresponds to
                       the maximum flow of a given commodity, k,
			           allowed in a given arc, j.*/

  FRow TotCap;     ///<Array of mutual capacity for the arcs.
  FRow Deficits;   ///<deficits
  
  HpNum LowerBound;     ///< A (coarse) lower bound on the value of Fi.

/*@}------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Solutions data
    @{ */
   

  FRow FSolution;    ///< Value of the "actual" flow solution

  Index SolWFi;      ///<Which subproblem belongs last gi created.
  Index_Set OldWFi;  ///<Which subproblem belongs the solution.
  FRow *OldFSols;    ///<History of subgradients (Flow Solution).

  bool DoAggregation; ///< if true, do the aggregation

  HpRow FiLmbd;      ///<Fi[ k ]( Lambda )

/*@}------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Variables management related
    @{ */

  Index_Set SGBse1;     ///< format of Subgradient
  Bool_Vec SlvP;        ///< true if the Lagrangian problem has been solved
  Bool_Vec LsHasChgd;   /**< true if Lambda has changed since the last call
   			                 to NewGi() */

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
/*------------------------ PRIVATE DATA STRUCTURES -------------------------*/
/*--------------------------------------------------------------------------*/

  CRow NwCsts;
  Index_Set *IndNwCsts;
  int *QSStck;   ///< the stack to simulate recoursive calls in QS

  HpNum LpshCnst;   // the Lipschitz constant

/*--------------------------------------------------------------------------*/
/*-------------------------- PRIVATE METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/

   HpNum solve( cIndex wFi = Inf<Index>() );
   ///<Solve wFi arc knapsack problem.

/*--------------------------------------------------------------------------*/

   void quicksort( cIndex arc );

/**< Sort the a vector and corresponding positions stored in the vector b.
    Used to sort the costs of each commodity to find the most atractive
    (negative) costs to use in the knapsack algorithm. */

/*--------------------------------------------------------------------------*/

   void CopyGi( cIndex Name, cIndex wFi , cFRow FlwSol );

/**< Copy the primal solution in the position called Name into the structure
    used for maintaining the previous primal solutions. */

/*--------------------------------------------------------------------------*/

   void FindGlobalLipschitz( void );

/*--------------------------------------------------------------------------*/
/*----------------------- PRIVATE DATA STRUCTURES  -------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/

 };  // end( class KnpsFiOrcl )

/*--------------------------------------------------------------------------*/

#if( OPT_USE_NAMESPACES )
 };  //end namespace NDO_di_unipi_it
#endif

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif  /* KnpsFiOrcl.h included */

/*--------------------------------------------------------------------------*/
/*------------------------- End File KnpsFiOrcl.h --------------------------*/
/*--------------------------------------------------------------------------*/

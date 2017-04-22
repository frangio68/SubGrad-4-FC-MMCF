/*--------------------------------------------------------------------------*/
/*---------------------------- File FlwFiOrcl.h ----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Definition of the FlwFiOrcl class which implements the FiOracle interface
 * and solves the flow Lagrangian relaxation of the (Fixed-Charge)
 * Multicommodity Min Cost Flow Problem ((FC-)MMCF), possibly with "easy"
 * components for the design variables. It derives from MMCFFlwBase (and
 * hence from MMCFClass) and uses its methods to solve the flow relaxation.
 *
 * \version 1.06
 *
 * \date 10 - 03 - 2014
 *
 * \author Antonio Frangioni \n
 *         Operations Research Group \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n \n
 *
 * \author Enrico Gorgone \n
 *         Operations Research Group \n
 *         Diparimento di Elettronica Informatica e Sistemistica \n
 *         Universita' della Calabria \n
 *
 * Copyright(C) 2001 - 2012 by Antonio Frangioni, Enrico Gorgone
 */
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __FlwFiOrcl
 #define __FlwFiOrcl  /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- MACROS ----------------------------------*/
/*--------------------------------------------------------------------------*/

/*------------------------------- LOG_CUT ----------------------------------*/

#define LOG_FI 0
/**< If LOG_FI > 0, the FlwFiOrcl class may produce (depending on an input
   parameter) a log of (more or less) important results. */

#define LowB_FI 1
///< If LowB_FI > 0, the FlwFiOrcl class uses the global Lower Bound.

#define SEP_TYPE 1
/**< Define the type of the separation approach:

    0:  find the maximum violation maxviol of the constraints already involved and
        add the ones whose violation is greater that maxviol
    1:  add the constraints whose violations is greater the epsilon */

#define SPF_SUB 0
/**< indicate how to manage the subgradients in the funciton Gi():

    0: do not sparsify the subgradients;
    1: sparsify the sungradients. */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "FiOracle.h"

#include "MMCFFlwB.h"
#include "Graph.h"

#include "OPTUtils.h"
#include <vector>

/*--------------------------------------------------------------------------*/
/*------------------------------ NAMESPACE ---------------------------------*/
/*--------------------------------------------------------------------------*/

#if( OPT_USE_NAMESPACES )
namespace NDO_di_unipi_it
{
 using namespace MMCFClass_di_unipi_it;
#endif

/*--------------------------------------------------------------------------*/
/*----------------------------- CLASSES ------------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/** This class instantiates the FiOracle interface [see FiOracle.h] to solve
 * the flow Lagrangian relaxation of the (Fixed-Charge) Multicommodity Min
 * Cost Flow Problem ((FC-)MMCF), possibly with "easy" components for the
 * design variables. It derives from MMCFFlwBase (and hence from MMCFClass)
 * and uses its methods to solve the flow relaxation.
 *
 * The problem formulation is:
 * \f[
 *   \phi( \alpha , \beta ) =
 *   \min \left\{
 *        \sum_{k \in K} \sum_{(i, j) \in A}
 *        \left( c_{ij}^k + \alpha_{ij} + \beta_{ij}^k \right)x_{ij}^k +
 *        \left( f_{ij} - \alpha_{ij} u_{ij} -
 *               \sum_{k \in K} \beta_{ij}^k u_{ij}^k
 *               \right) y_{ij}
 *        \right\}
 * \f]
 * \f[
 *   \sum_{(i, j) \in A} x_{ij}^k - \sum_{(j,i) \in A} x_{ij}^k = b_i^k
 *   \quad,\quad i \in N , k \in K
 * \f]
 * \f[
 *   0 \leq x_{ij}^k \leq u_{ij}^k \quad,\quad (i, j) \in A , k \in K
 * \f]
 * \f[
 *   y_{ij} \in \{0, 1\} \quad,\quad (i, j) \in A
 * \f]
 * The pair \f$( \alpha, \beta )\f$ is the Lambda vector described in the
 * FiOracle interface. The \f$\alpha\f$'s are dual variables associated to
 * the relaxed mutual capacity constraints. These variables are always used.
 * The \f$\beta\f$'s are dual variables associated to the relaxed individual
 * capacity constraints. These variables only are used if
 * SetICapRelax( true ) is called. */

class FlwFiOrcl : public FiOracle , public MMCFFlwBase
{

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*--------------------------------------------------------------------------*/
/*---------------------------- PUBLIC TYPES --------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Public Types
	    @{ */

/** Public enum for handling the FlwFiOrcl-specific parameters in
    [see below]. */

    enum FlwParam { kRlx  = 0 , kYiE , kSp1 , kSp2 , kSp3 , kAgg };

/*@} -----------------------------------------------------------------------*/
/*--------------------- PUBLIC METHODS OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/
/*---------------------------- CONSTRUCTOR ---------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructor
    @{*/

   FlwFiOrcl( Graph *g , istream *iStrm = NULL );

/**< Constructor of the class.

   \param g an instance of the class Graph that defines the problem.

   \param iStrm parameters file used to read some parameters not predicted
          in the FiOracle interface. In particular, the following parameters
	  are read:

	  -# SPar1 [0] how often the check for adding variables is performed
	        (0 = never, all variables are there from the beginning)

	  -# SPar2 [0] a variable is removed if it has had value zero for
	         SPar2 consecutive iterations, if 0 = no variables is removed

	  -# SPar3 [0] how often the check form removing variable is performed

	  -# SPar4 [w] relaxation type: 0 weak , b strong formulation and s the
	             individual capacity constraints are the ones and only the ones
	             involved

	  -# YiE   [0] 0: easy components are to be used 1: otherwise

	  -# Aggrgtd [1] 0 if disaggregated master problem 1: aggregated master
	                    problem

	  -# KOld   [0] pre-set number of old solutions that are kept in memory

      -# EpsFi  [1e-6] precision for the design part
      -# EpsCon [1e-6] precision for the constraints

      Note that the separation parameters SPar1-3 are irrelevant if the
      weak relaxation has been required.                                    */

/*@} -----------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
    @{*/

   void SetFiLog( ostream *outs = 0 , const char lvl = 0 );

   /**< lvl controls the "level of verbosity" of the code. The first four bits
      of lvl have the following meaning:

      0  =>  no log at all (also assumed if log = 0);

      >0  =>  "basic" log: only the errors are reported; */


/*--------------------------------------------------------------------------*/

   void SetMaxName( cIndex MxNme );

/*--------------------------------------------------------------------------*/

   void SetAggregate( bool aggrgt = true );

/**< Set aggregation option. Used to tell the FiOracle if it should work with
   Fi aggregated (true) or disaggregated (false) option. */

/*--------------------------------------------------------------------------*/

   void SetRelax( const char sp4 = 'w' , bool YiE = false ,
		  cIndex sp1 = 0 , cIndex sp2 = 0 , cIndex sp3 = 0 );

/**< Set relaxation of Individual Capacities or off.
     Used to tell to the FiOracle if the relaxation of individual
     capacities constraints is considered. In the case of strong
     formulation, set the parameters PAV [0] and PRV [0] for adding and
     removing of the betas. [see FlwFiOrcl( Graph*, istream*=0 ) ]. */

/*--------------------------------------------------------------------------*/

   void SetInitialSet( Index_Set Bse );

/**< Set the initial dictionary. Bse is the base of the current (beta) solution
 * Thus, the dimension is not greater than
 *
 * - MaxNumVar - NArcs, if there are mutual forcing capacity constraints
 * - MaxNumVar, if there are only individual forcing capacity constraints
 */

/*--------------------------------------------------------------------------*/

   void SetEasy ( bool YIsE = false );

/**< Set Easy option.
     Used to tell to the FiOracle that there exist "easy" components of Fi().
     By default YIsEasy has the value to false  and all components are handled
     as "difficult".   */

/*@} -----------------------------------------------------------------------*/
/*-------------- METHODS FOR READING THE DATA OF THE PROBLEM ---------------*/
/*--------------------------------------------------------------------------*/

   Index GetNumVar( void ) const;

/*--------------------------------------------------------------------------*/

   Index GetMaxNumVar( void ) const;

/*--------------------------------------------------------------------------*/

   Index GetNrFi( void ) const;

/*--------------------------------------------------------------------------*/

   bool GetUC( cIndex i );

/*--------------------------------------------------------------------------*/

   HpNum GetGlobalLipschitz( cIndex wFi = Inf<Index>() );

/*--------------------------------------------------------------------------*/

   Index GetBNC( cIndex wFi );

/*--------------------------------------------------------------------------*/
   Index GetBNR( cIndex wFi );

/*--------------------------------------------------------------------------*/

   Index GetBNZ( cIndex wFi );

/*--------------------------------------------------------------------------*/

   void GetBDesc( cIndex wFi , int *Bbeg , int *Bind , double *Bval ,
		  double *lhs , double *rhs , double *cst ,
		  double *lbd , double *ubd );

/*--------------------------------------------------------------------------*/

   Index GetANZ( cIndex wFi , cIndex strt = 0 , Index stp = Inf<Index>() );

/*--------------------------------------------------------------------------*/

   void GetADesc( cIndex wFi , int *Abeg , int *Aind , double *Aval ,
		  cIndex strt = 0 , Index stp = Inf<Index>() );

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
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/

   bool GetPSol( void );

/*--------------------------------------------------------------------------*/

   bool NewGi( cIndex wFi = Inf<Index>() );

/**< The first time NewGi() is called, it computes the "true" subgradients.
   The second time it computes (epsilon-)subgradients associated with a
   primal solution that is obtained by running the heuristic with the
   commodities ordered by decreasing demand. Every subsequent call, it
   computes (epsilon-)subgradients associated with primal solutions that are
   obtained by running the heuristic with the commodities randomly ordered.
   Thus, it can always be called, but it can generate repeated epsilon
   subgradients.

    The "true" subgradients are calculated using the following formulas:
    - If the disaggregated version is used:
    \f[
       gi^k(\alpha_{ij})=gi^k(\beta_{ij}^k)=x_{ij}^k ,\forall\:
                                                      subproblem\: k
    \f]
    \f[ gi^0(\alpha_{i,j})=-u_{ij}y_{ij}\f]
    \f[ gi^0(\beta_{i,j}^k)=b_{ij}^ky_{ij}\f]
    - If the aggregated version is used:
    \f[ gi(\alpha_{i,j})=\left(sum_{k \in K}x_{ij}^k-u_{ij}y_{ij}\right) \f]
    \f[ gi(\beta_{i,j}^k)=x_{ij}^k-b_{ij}^ky_{ij}\f]
   */

/*--------------------------------------------------------------------------*/

   Index GetGi( SgRow SubG , cIndex_Set &SGBse , cIndex Name = Inf<Index>() ,
		cIndex strt = 0 , Index stp = Inf<Index>() );

/*--------------------------------------------------------------------------*/

   HpNum GetVal( cIndex Name = Inf<Index>() );

/*--------------------------------------------------------------------------*/

   void SetGiName( cIndex Name );

/*--------------------------------------------------------------------------*/
/*-------------------- METHODS FOR READING OTHER RESULTS -------------------*/
/*--------------------------------------------------------------------------*/
/** @name Reading other results
    @{ */

  HpNum GetLowerBound( cIndex wFi = Inf<Index>() );

  void SetLowerBound( HpNum lowB = Inf<HpNum>() );

/*--------------------------------------------------------------------------*/

  FiStatus GetFiStatus( Index wFi = Inf<Index>() );

/*--------------------------------------------------------------------------*/

  double AddTime( void );

/**< Returns the user + system time (in seconds) spent during the
     execution for adding of new beta variables. */

/*--------------------------------------------------------------------------*/

  double RemTime( void );

/**< As AddTime() [see above], it returns the time spent for removing
   the beta variables. */


/*--------------------------------------------------------------------------*/

  #if LOG_FI
   Index GetNumVarAvg( void );
  #endif

/*--------------------------------------------------------------------------*/

   bool checkSolution( void );

/*@} -----------------------------------------------------------------------*/
/*-------------- METHODS FOR READING THE DATA OF THE PROBLEM ---------------*/
/*--------------------------------------------------------------------------*/

   inline void GetPar( const int wp , char &value );

   inline void GetPar( const int wp , bool &value );

/**< Read the current value of the "bool" parameters of the FlwFiOracle.
     The enum FlwParam is used (in the obvious way) for selecting the
     parameter to be get.*/

   inline void GetPar( const int wp , Index &value );

/**< Read the current value of the "int" parameters of the FlwFiOracle
     [see above].  The enum FlwParam is used (in the obvious way)
     for selecting the parameter to be get.*/

/*@}------------------------------------------------------------------------*/
/*------------- METHODS FOR ADDING / REMOVING / CHANGING DATA --------------*/
/*--------------------------------------------------------------------------*/

   void Deleted( cIndex i = Inf<Index>() );

/*--------------------------------------------------------------------------*/

   void Aggregate( cHpRow Mlt , cIndex_Set NmSt , cIndex Dm , cIndex NwNm );

/*--------------------------------------------------------------------------*/
/*------------------------------ DESTRUCTOR --------------------------------*/
/*--------------------------------------------------------------------------*/

   virtual ~FlwFiOrcl();

///< Destructor of the class: it must be virtual.

/*--------------------------------------------------------------------------*/
/*-------------------- PROTECTED PART OF THE CLASS -------------------------*/
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

   HpNum EpsFi;     /**< Precision for the design problem.
		       The design solution value is considered (y=1) only if
		       less than or equal to -Eps.
		       This value is set to EFnal (read from
		       parameter input file) by the NDOSolver.*/

   HpNum EpsCon;    /**< Precision for the constraints problem*/

   bool Aggrgtd;    /**< Agreggation flag.
		       Variable that defines if FlwFiOracle outputs
		       aggregated or disaggregated information*/

   bool YIsEasy;    /**< true if the ys are easy components. */

/*@}------------------------------------------------------------------------*/
/** @name Lagrangian data
    @{ */

   Index_Set InvULambda;    /**< Inverse vocabulary of ULambda. This vector
			       contains the name of the variable beta
			       corresponding to the position in Lambda of the
			       index (+ NArcs) in the vector. It terminates
			       with Inf<Index>() */

   Index_Set ZLambdaCount;  /**< Counter of iterations with value zero for
			       each variable. */

   Index_Set SGBse1;       // format of Subgradient

/*@}------------------------------------------------------------------------*/
/** @name Instance data
    @{ */

   Index MaxNumVar;   /**< The maximum number of variables.
			 By default the value of this variable is equal to
			 NrArcs() + NrArcs() * NrComm(), it means the
			 individual capacities constrainsts are relaxed.
			 If SetICapRelax( false ) is called then the value
			 of the variable is NrArcs(), it means only the
			 mutual capacity constraints are relaxed. */

   CRow OrigCosts;   /**< Array of original unit costs.
			Each entry corresponds to the cost of one unit of
			flow of one commodity traversing one arc. For arc j
			and commodity k, this cost is OrigCosts[ k * m + j ],
			where m stands for the total number of arcs.*/

   CRow OrigXtrCosts; /**< Array of original design costs.
			 Each entry is associated with the fixed cost of the
			 corresponding arc; if OrigXtrCosts == 0, then no
			 extra costs are defined for the instance.*/

   CRow OrigCapacities;  /**< Array of original individual capacities.
                          Each entry OrigCapacities[ k * m + j ] corresponds
			  to the maximum flow of a given commodity, k,
			  allowed in a given arc, j.*/

   FRow OrigTotCap;      ///< Array of mutual capacity for the arcs.
   FRow OrigDeficits;    ///< Original deficits.

   HpNum LowerBound;     ///< A (coarse) lower bound on the value of Fi.

/*@}------------------------------------------------------------------------*/
/** @name Solutions data
    @{ */

   Index_Set OldWFi;     ///< Witch subproblem belongs the solution.

   FRow *OldFSols;       ///< History of subgradients (Flow Solution).
   Mat OldXSols;         ///< History of subgradients (Extra Solution).

   SIndex KOld;           ///< number of allocated Old Solution

   FRow FSolution;
   Index SolWFi;         ///< Witch subproblem belongs last gi created.

   Row XSolution;
   Bool_Vec SolvedP;

/*@}------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Variables management related
    @{ */

   Index SPar1;          /**< Variables addition interval.
                            Original constraints are checked every
                            SPar1 and a variable is added for each
                            one currently violated (0 = no adding of
			                variables is allowed). */
   Index SPar2;          /**< Maximum n. of iterations each variable can be 0.
                            If one variable is ParRemVariables iterations
			                with value zero it's removed. */
   Index SPar3;          /**< Remotion of original constraints is checked
                            every SPar3 iterations. */
   char SPar4;           /**< For specifying which kind of constraints have to
                            be relaxed.*/

   Bool_Vec LsHasChgd;   /**< true if Lambda has changed since the last call
			    to NewGi() */

   OPTtimers *Addt;      ///< Timer for variables addition
   OPTtimers *Remt;      ///< Timer for variables deletion

/*@}------------------------------------------------------------------------*/
/*--------------------- PRIVATE PART OF THE CLASS --------------------------*/
/*--------------------------------------------------------------------------*/

 private:

/*--------------------------------------------------------------------------*/
/*-------------------------- PRIVATE METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/

  void XSolA( Row Xsl );

  /* Compute the aggregate extra variables solution. */

  void FSolA( FRow Flw );

  /* Compute the aggregate flow solution. */

/*--------------------------------------------------------------------------*/

  bool checkAddition( void );

// Check if the variables not used are violated by the current solution.

/*--------------------------------------------------------------------------*/

  bool checkRemotion( void );

// Check if the variables are useless.

/*--------------------------------------------------------------------------*/

 void CopyGi( cIndex Name, cIndex wFi , cRow XtrSol , cFRow FlwSol );

/**< Copy the primal solution in the position called Name into the structure
 used for maintaining the previous primal solutions. */

/*--------------------------------------------------------------------------*/

 void FindGlobalLipschitz( void );

/*--------------------------------------------------------------------------*/

  void SolveLagrangian( cIndex SubPName = Inf<Index>() );

/*--------------------------------------------------------------------------*/
/*----------------------- PRIVATE DATA STRUCTURES  -------------------------*/
/*--------------------------------------------------------------------------*/

  Index LstChgItr;  // iteration at which the latest change in the
                    // number of variables (add/remove) happened

  Index LastGiName;

  Index mutual;     // number of mutual constraints

  HpNum L2Cnst;     // the squared Lipschitz constant

  HpRow CoefObj;    // coefficient of the Y part

  vector<Index> *LmbCmp; // for each entry of Lambda tells to which component
                         // is belonging

/*--------------------------------------------------------------------------*/

 };  // end( class FlwFiOrcl )

/*--------------------------------------------------------------------------*/
/*------------------- inline methods implementation ------------------------*/
/*--------------------------------------------------------------------------*/

inline void FlwFiOrcl::GetPar( const int wp , char &value )
{
 switch( wp ) {
  case( kRlx ):  value = SPar4; break;
  default: throw(
    NDOException( "FlwFiOrcl::GetPar( bool ): unknown parameter" ) );
  }
 }  // end( Bundle::GetPar( bool ) )

/*--------------------------------------------------------------------------*/

inline void FlwFiOrcl::GetPar( const int wp , bool &value )
{
 switch( wp ) {
  case( kYiE ):  value = YIsEasy; break;
  case( kAgg ):  value = Aggrgtd; break;
  default: throw(
    NDOException( "FlwFiOrcl::GetPar( bool ): unknown parameter" ) );
  }
 }  // end( Bundle::GetPar( bool ) )

/*--------------------------------------------------------------------------*/

inline void FlwFiOrcl::GetPar( const int wp , Index &value )
{
 switch( wp ) {
  case( kSp1 ):  value = SPar1; break;
  case( kSp2 ):  value = SPar2; break;
  case( kSp3 ):  value = SPar3; break;
  default: throw(
    NDOException( "FlwFiOrcl::GetPar( Index ): unknown parameter" ) );
  }
 }  // end( Bundle::GetPar( Index ) )

/*--------------------------------------------------------------------------*/

#if( OPT_USE_NAMESPACES )
 };  //end( namespace( NDO_di_unipi_it ) )
#endif

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif  /* FlwFiOrcl.h included */

/*--------------------------------------------------------------------------*/
/*------------------------- End File FlwFiOrcl.h ---------------------------*/
/*--------------------------------------------------------------------------*/

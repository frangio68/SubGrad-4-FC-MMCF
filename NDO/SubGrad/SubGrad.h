/*--------------------------------------------------------------------------*/
/*-------------------------- File SubGrad.h --------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 *
 * Definition of the <tt>SubGrad</tt> class, which implements the NDOSolver
 * interface for NonDifferentiable Optimization Solvers, as described in
 * NDOSlver.h, using a (deflected, projected, incremental) subgradient-type
 * algorithm.
 *
 * The user is assumed to be familiar with the algorithm: refer to
 *
 *  A. Frangioni, E. Gorgone, B. Gendron.
 *  <em>On the Computational Efficiency of Subgradient Methods: A Case Study
 *       in Combinatorial Optimization </em> .
 *  Technical Report CIRRELT- 2015-41, Interuniversity Research
 *  Centre on Enterprise Networks, Logistics and Transportation (CIRRELT),
 *  Montreal, Canada, 2015
 *
 * \version 1.00
 *
 * \date 27 - 10 - 2015
 *
 * \author Antonio Frangioni \n
 *         Operations Research Group \n
 *         Department of Informatics \n
 *         University of Pisa \n\n
 *
 * \author Enrico Gorgone \n
 *         Graphs and Mathematical Optimization (GOM) Group \n
 *         Department of Informatics \n
 *         Free University of Brussels \n\n
 *
 * The class requires that the function to be minimized be available under
 * the form of a FiOracle object, as described in FiOracle.h.
 *
 * The class is parametric over the type of both stepsize and direction used:
 * it just relies over the objects of classes Stepsize and Deflection.
 *
 * Copyright &copy 2001 - 2015 by Antonio Frangioni, Enrico Gorgone
 */

/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __SubGrad
 #define __SubGrad  /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- MACROS ----------------------------------*/
/*--------------------------------------------------------------------------*/
/** @defgroup SubGrad_MACROS Compile-time switches in SubGrad.h
    @{ */

#define LOG_SG 1

/**< If LOG_SG > 0, the SubGrad class produces a log of its activities on the
   ostream object and at the "level of verbosity" set with the method
   SetNDOLog() [see below]. */

/* @} end( group( SubGrad_MACROS ) ) */
/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "NDOSlver.h" // NDOSolver is the base class

class SubGrad;  //  The forward declaration is done to allow to use a pointer
                //  to SubGrad inside Stepsize and Deflection.

#include "Deflection.h"  // These classes are defined as friend classes, in
#include "Stepsize.h"    // order to allow them to use directly SubgGrad's
                         // data structures [see below].
/* Note: a few methods of Stepsize and Deflection are implemented at the end of
   SubGrad.C. These are the methods that allow to read protected data of
   SubGrad, which may be useful to the objects. The methods are implemented
   here for the base classes, so that any derived class can use them to
   access to these data. This is crucial, as derived classes from friend
   classes are not friend, and therefore they cannot access data of SubGrad
   unless this capability is explicitly provided by the base classes, who are
   friends. */

#include "CQKnPClass.h"  // This class is used to deal with the constrained
                         // optimization problem if the feasible region
                         // consists of the only knapsack constraint.

#include <vector>     // ... to control the random shuffle
#include <algorithm>

/*--------------------------------------------------------------------------*/
/*----------------------------- NAMESPACE ----------------------------------*/
/*--------------------------------------------------------------------------*/

#if( OPT_USE_NAMESPACES )
namespace NDO_di_unipi_it
{
#endif

 using namespace CQKnPClass_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*----------------------------- CLASSES ------------------------------------*/
/*--------------------------------------------------------------------------*/
/** @defgroup SubGrad_CLASSES Classes in SubGrad.h
    @{ */

/*--------------------------------------------------------------------------*/
/*--------------------------- CLASS SubGrad  -------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/** The <tt>SubGrad</tt> class implements the NDOSolver interface for
    NonDifferentiable Optimization Solvers [see NDOSlver.h], using a unified
    subgradient-type algorithm as described in:

     A. Frangioni, E. Gorgone, B. Gendron.
     <em>On the Computational Efficiency of Subgradient Methods: A Case
         Study in Combinatorial Optimization</em>.
     Technical Report CIRRELT- 2015-41, Interuniversity Research
     Centre on Enterprise Networks, Logistics and Transportation (CIRRELT),
     Montreal, Canada, 2015

    This is in fact a subgradient method (SM) based on abstract rules for
    the computation of both the <em>stepsize</em> \f$\nu_i\f$ and the
    <em>direction</em> \f$d_i\f$. The algorithm in employs the simple
    recurrence formula:
     \f[
      \breve{\lambda}_{i+1} \gets \bar{\lambda}_i - \nu_i d_i
      \quad , \quad
      \lambda_{i+1} \gets {\rm P}_{\Lambda}( \, \breve{\lambda}_{i+1} \, )
    \f]
    where \f${\rm P}\f$ denotes the <em>orthogonal projection</em> on
    \f$\bar{\Lambda}\f$. The point \f$\bar{\lambda}_i\f$ is not necessarily
    the current iterate. For instance, it could be required that
    \f$\bar{\lambda}_{i}\f$ remains unchanged if an ascent direction occurs.
    It recalls somehow the <em>stability center</em> of the bundle methods.

    The class relies on the objects of the friend classes <tt>Stepsize</tt>
    and <tt>Deflection</tt>, which have to return, respectively, the
    stepsize \f$\nu_i\f$ and the deflection coefficient \f$\alpha_i\f$. The
    latter scalar number defines in turn the direction \f$d_i\f$, i.e..
    \f$ d_i = \alpha_i g_i + (1-\alpha_i)d_{i-1} \f$.
    These abstract classes allow us to derive several variants of the method.
    For the sake of simplicity, the original SM characterized by
    \f$\alpha_i = 1\f$ is performed setting the pointer to the object
    Deflection to NULL.

    The class also includes the <em>incremental</em> variant for when the
    function to be maximized is composed by a sum of different functions.
    */

class SubGrad : public NDOSolver
{

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*--------------------------------------------------------------------------*/
/*---------------------------- FRIEND CLASSES ------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Friend classes
	@{ */

/** The classes Stepsize and Deflection are "friends" of SubGrad. This is
    done because Stepsize and Deflection objects may need some protected
    data to work. An issue, however, is that derived classes from friend
    classes are not friend, and therefore actual implementations of
    Stepsize and Deflection cannot access data of SubGrad unless this
    capability is explicitly provided by the base classes, who are friends.
    This is why Stepsize and Deflection define a few methods that allow to
    read protected data of SubGrad: so that any derived class can use them to
    access to these data. These methods are, in fact, implemented at the end
    of SubGrad.C. */

  friend class Stepsize;
  friend class Deflection;
 
/*@} -----------------------------------------------------------------------*/
/*---------------------------- PUBLIC TYPES --------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Public Types
	@{ */

/** Public enum which "extends" the enum NDOSolver::NDOParam for handling
    the SubGrad-specific algorithmic parameters in (the two overloaded
    versions of) SubGrad::SetPar() [see below]. */

   enum SGParam{ kSGPar1 = kLastNDOParam ,
		 kSGPar2 , kSGPar3 , kSGPar4 , kSGPar5 };

/*@} -----------------------------------------------------------------------*/
/*--------------------- PUBLIC METHODS OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*---------------------------- CONSTRUCTOR ---------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructor
    @{ */

   SubGrad( istream *iStrm = 0 );

/**< Constructor of the class. The parameter `iStrm', if provided, is taken
   as a pointer to a istream from which the algorithmic parameters for the
   subgradient algorithm are sequentially read in the following order. Each
   parameter must be placed at the beginning of a separate line, max 255
   characters long, with all the rest of the line up to the first newline
   character (apart from a separating whitespace) being available for
   comments. Any line whose first character is '#' and any blank line is
   ignored. If 0 is passed, the file ends before reaching a given parameter,
   or some parameter is in the wrong format, each non-specified parameter is
   given a default value, shown in [] below.

   `iStrm' is passed to the constructor of NDOSolver [see NDOSolver.h], which
   reads the general algorithmic parameters out of it; since the constructor
   SubGrad is executed after the one of NDOSolver, the following parameters
   specific for the SubGrad have to be found in the stream *after* those of
   the base class:

    -# Index SGPar1 [0] The direction \f$d_i\f$ is assumed to be such a
                    convex combination of the subgradient \f$g_i\f$ and
                    the direction \f$d_{i-1}\f$:
                    \f[
                      d_i = \alpha_i g_i + ( 1 - \alpha_i ) d_{i-1}
                    \f]
                    SGPar1 selects among \f$\{d_i,d_{i-1},g_i\} \f$ the
       vector(s) to be projected over the tangent cone at 
       \f$\bar{\lambda}_i\f$. The field SGPar1 is coded bit-wise, in the
       following way:

        - bit 0:  if 1 the subgradient \f$g_i\f$ is projected, 0 otherwise;

        - bit 1:  if 1 the direction \f$d_{i-1}\f$ is projected, 0 otherwise;

        - bit 2:  if 1 the \f$d_{i}\f$ is projected, 0 otherwise;

       The setting (+7) is redundant. In fact, it is equivalent to (+3)
       because \f$d_i\f$, being a (convex) combination of \f$g_i\f$ and
       \f$d_{i-1}\f$, coincides with its projection.

    -# HpNum SGPar2 [0] To manage the incremental iterations. A sequence
                    of incremental (or inner) iterations NItIncr performed
                    along single-component subgradients could occur before a
       full (or normal, or outer) iteration. The number NItIncr is obtained
       as NItIncr = ceil( <Number of components of \f$f\f$> * SGPar2 ),
       where the number of components of \f$f\f$ includes the 0-th component,
       and then not necessarily coincides with NrFi. Finally, the incremental
       variant only works with no deflection object, i.e setting deflection
       = NULL, because we do not manage with the sum of subgradients of
       different components.

    -# Index SGPar3 [1] The convergence scheme. This field is coded
                        bit-wise in the following way:

                         - bit 0:  1 if the safe rule is used, 0 otherwise

                         - bit 1:  1 for the stepsize-restricted scheme,
                                   0 for the deflection-restricted scheme.

    -# bool SGPar4 [true] SGPar4 enables the use of
                          \f[
                           \hat{\lambda}_{i+1} = \alpha_{i+1}\lambda_i +
                           ( 1 - \alpha_{i+1} ) \hat{\lambda}_i
                          \f]
                          which could have a certain influence on the
                          stopping test [see IsOptimal()].

    -# Index SGPar5 [0] The seed value used in a call to srand. The
                        components are re-shuffled in the incrmental variant,
                        and a random number generator is used.
    */

/*@} -----------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
    @{ */

   void SetStepsize( Stepsize *STP = 0 );

/**< Gives to the <tt>SubGrad</tt> object a pointer to an object of class
   <tt>Stepsize</tt> that will be used to provide \f$\nu_i\f$ during the
   subgradient algorithm.

   The Stepsize object can be changed during the life of a SubGrad object,
   but this change clearly forces the reset of all the information about the
   function accumulated so far. Passing a 0 pointer does exactly this job. */

/*--------------------------------------------------------------------------*/

   void SetDeflection( Deflection *Vol = 0 );

/**< Gives to the <tt>SubGrad</tt> object a pointer to an object of class
   <tt>Deflection</tt> that will be used to provide a deflection coefficient
   \f$ \alpha_i\f$ .

   The Deflection object can be changed during the life of a SubGrad object,
   but this change clearly forces the reset of all the information about the
   function accumulated so far. Passing a 0 pointer does exactly this job.

   In addition, 0 pointer means that the deflection coefficient is kept to 1,
   namely \f$\alpha_i = 1\f$.  */

/*--------------------------------------------------------------------------*/

   void SetQKNP( CQKnPClass *KNP = 0 );

/**< Gives to the SubGrad object a pointer to an object of class CQKnPClass
   that will be used as quadratic knapsack solver during the subgradient
   algorithm.

   CQKnPClass object is employed to project the subgradient/direction when the
   projection problem is a continuous quadratic knapsack one. Moreover, the
   data of the problem are fixed just before the projection, hence no
   information is kept in the Subgradient method.

   The CQKnPClass object can be changed during the life of a SubGrad object,
   but this change clearly implies that the prior object must be discharged.
   Passing a 0 pointer does exactly this job.
 
   The default implementation assigns zero to the pointer of the CQKnPClass object
   directly in the constructor. In fact, CQKnPClass solver comes be helpful 
   whenever in the feasible set the knapsack constraints appear. In contrast to
   what happens in both SetStepsize() and SetDeflection(), the call to the 
   function SetQKNP() is mandatory only if knapsack constraints take place in the
   feasible set. */

/*--------------------------------------------------------------------------*/

   void SetFiOracle( FiOracle *Fi = 0 );

/*--------------------------------------------------------------------------*/

   void SetLambda( cLMRow tLambda = 0 );

/*--------------------------------------------------------------------------*/

   void KeepBestLambda( const bool KBL = true );

/*--------------------------------------------------------------------------*/

   void SetPar( const int wp , const int value );

/**< Extends NDOSolver::SetPar( , cIndex ) for handling the SubGrad-specific
    parameters; the enum SGParam is used (in the obvious way) for selecting
    the parameter to be set. */

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

   void SetPar( const int wp , cHpNum value );

/**< Extends NDOSolver::SetPar( , cHpNum ) for handling the SubGrad-specific
    parameters; the enum SGParam is used (in the obvious way) for selecting
    the parameter to be set. */

/*--------------------------------------------------------------------------*/

   void SetPar( const int wp , const bool value );

/**< Change boolean algorithmic parameters of the SubGrad solver. The enum
   SGParam [see above] is used for selecting the parameter to be set. */

/*--------------------------------------------------------------------------*/

   void SetNDOLog( ostream *outs = 0 , const char lvl = 0 );

/**< lvl controls the "level of verbosity" of the code. The first four bits
    of lvl have the following meaning:

    - 0  =>  no log at all (also assumed if log = 0);

    - 1  =>  "basic" log: only the errors are reported;

    - 2  =>  a detailed step-by-step log of the algorithm is displayed;

    - 4 .. 15 unused, available to derived classes;  */

/*--------------------------------------------------------------------------*/
/*-------------------- METHODS FOR SOLVING THE PROBLEM ---------------------*/
/*--------------------------------------------------------------------------*/
/** @name Solving the problem
    @{ */

   NDOStatus Solve( void );

/**< Tries to minimize the function [see NDOSlver.h]. It implements the
   subgradient algorithm exploiting the protected methods FormD(), SaveDir(),
   FormLambda1(), FiAndGi(), GotoLambda1() [see below].

   Returns   if

   - kOK       optimization has been succesful: a solution that is "optimal"
               (w.r.t. the current parameters settings) has been found;

   - kUnbndd   there has been an error in the FiOracle, i.e. Fi() has returned
               -HpINF, or the function is unbounded below: the latter case can
               be detected only if a lower bound on the min. value of
               \f$f\f$ is available [see FiOracle::GetMinusInfinity()];

   - kUnfsbl   the polyhedral set defined by the constraints is empty: in this
               case, the primal optimal solution is an unbounded *extreme ray*
	           for the dual problem;

   - kStopped  Solve() has been stopped, either by FiOracle::GetFiStatus()
               or because the stepsize has been "too small" during 100 consecutive
               outer iterations [or 100 * NrFi consecutive inner iterations];

   - kStpIter  the max. number of iterations has been exhausted;

   - kStpTime  the max. running time has been reached;

   - kError    There have been some problem in the FiOracle that require to
               stop the optimization.
 
   As for kStopped, "too small" means that \f$ \nu_i \leq 1e-8 * t^*\f$, where
   \f$ t^* \f$ is the optimality related parameter scaling Fi() [see NDOSlver.h].
   There is no reason, in principle, why we couldn't replace \f$1e-8\f$ by
   a parameter, but in order to make the test easier this parameter has been
   fixed to \f$1e-8\f$. We also decided to replace by 100 the parameter saying
   how many outer iterations of consecutive small stepsizes are sufficient
   to stop the algorithm.

   Note that, whatever the exit condition be, the current point is always
   available by calling ReadSol(), and its Fi() value by calling ReadFiVal()
   [see below].   */

/*--------------------------------------------------------------------------*/

   void ReSetAlg( unsigned char RstLvl = 0 );

/**< Resets the internal state of the Subgradient algorithm. Since several
   different things can be reset independently, RstLvl is coded bit-wise:

   - bit 0: if 0 all the algorithmic parameters are reset to the default
     values read by the stream/set by SetPar(), while if 1 they are left
     untouched;

   - bit 1: if 0 the current point is reset to the all-0 vector, while if
     1 it is left untouched;

   - bit 2: if 0 the direction and the subgradient are removed, while if 1
     they are left there;

   - bit 3: if 0 all the constraints are removed, while if 1 all of them
     are left there.

   - bit 4: if 0 the value of Fi() in the stability center is reset to INF
     (i.e., unknown), while if 1 it is left untouched. */

/*@} -----------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/
/** @name Reading the solution
    @{ */

   cLMRow ReadBestSol( cIndex_Set &I , Index &D );

/**< Returns a read-only pointer to the point having the lowest \f$f\f$
   value found so far [see below].*/

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

   HpNum ReadBestFiVal( cIndex wFi = Inf<Index>() );

/**< Returns the best \f$f\f$ value found so far. Independently from which
   "component" of Fi() is chosen, it returns the full function. */

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

   cLMRow ReadSol( cIndex_Set &I , Index &D );

/**< Returns a read-only pointer to the <em>stability center</em>
   \f$\bar{\lambda}_i\f$. If Solve() has returned a kOK and the tStar has
   been properly set, the point returned by ReadSol() - and, a fortiori, the
   one returned by ReadBestSol() - is \f$\epsilon\f$-optimal. */

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

   HpNum ReadFiVal( cIndex wFi = Inf<Index>() );

/**< Independently from which "component" of Fi() is chosen, it returns the
   full function Fi at the <em> stability center </em> \f$\bar{\lambda}_i\f$.
   */

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

   HpNum ReadHatFiVal( void );

/**< Returns \f$\hat{f}\f$, if \f$ \hat{\lambda}_i \f$ is kept in memory.
   Otherwise it returns Inf<HpNum>().   */

 /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

   bool IsOptimal( HpNum eps = 0 );

/**< Returns true if the solution \f$\bar{\lambda}_i\f$ is
   \f$\epsilon\f$-optimal (relative), being \f$\epsilon \f$ set to EpsLin.
   The parameter SGPar4 controls the linearization error to be used in 
   the stopping condition:
 
   -# SGPar4 = true:
      \f[
       t^* \|d_i\| + \max\{ \hat{\epsilon}_i , \epsilon_i \} <= 
       \epsilon * max( 1 , |f_i^{rec}| )
      \f]

   -# SGPar4 = false: The criteria is just
      \f[
       t^* \|d_i\| + \epsilon_i <= \epsilon * max( 1 , |f_i^{rec}| )
      \f]
 
    where \f$ f^{rec}_i = \min \{ \, f_l \,:\, l = 1, \ldots, i \, \}\f$, i.e. 
    the <em> record value </em> on the optimum $f_*$ */

/*--------------------------------------------------------------------------*/

   cHpRow ReadMult( cIndex_Set &I , Index &D , cIndex wFi = Inf<Index>() );

   HpNum ReadLBMult( cIndex wFi = Inf<Index>() );

/*@} -----------------------------------------------------------------------*/
/*-------------- METHODS FOR READING THE DATA OF THE PROBLEM ---------------*/
/*--------------------------------------------------------------------------*/
/** @name Reading the data of the problem
    @{ */

   inline void GetPar( const int wp , int &value );

   inline void GetPar( const int wp , HpNum &value );

   inline void GetPar( const int wp , bool &value );

/*@}------------------------------------------------------------------------*/
/*------------- METHODS FOR ADDING / REMOVING / CHANGING DATA --------------*/
/*--------------------------------------------------------------------------*/
/** @name Adding / removing / changing data
    @{ */

   void AddVariables( Index NNwVrs , cLMRow IVs = 0 );

   void RemoveVariables( cIndex_Set whch = 0 , Index hwmny = 0 );

   void ChgFiV( cIndex wFi = Inf<Index>() );

   void ChgSbG( cIndex strt = 0 , Index stp = Inf<Index>() ,
		   cIndex wFi = Inf<Index>() );

/*@} -----------------------------------------------------------------------*/
/*------------------------------ DESTRUCTOR --------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Destructor
    @{ */

   ~SubGrad();

/*@} -----------------------------------------------------------------------*/
/*---------------------- PROTECTED PART OF THE CLASS -----------------------*/
/*--------------------------------------------------------------------------*/
/*--                                                                      --*/
/*--  The standard user should not care about the following part: users   --*/
/*--  who need to really extend the functionality of the SubGrad class    --*/
/*--  (e.g. by implementing a different subproblem solver) may use these  --*/
/*--  methods and data structures. It is *dangerous* to *modify* the      --*/
/*--  data structures, while it safe to read them.                        --*/
/*--                                                                      --*/
/*--------------------------------------------------------------------------*/

 protected:

/*--------------------------------------------------------------------------*/
/*-------------------------- PROTECTED METHODS -----------------------------*/
/*--------------------------------------------------------------------------*/

   void FiAndGi( cIndex wFi = Inf<Index>() );

/**< Evaluates the function at the new point \f$\lambda_{i+1}\f$, i.e.,
   \f$f(\lambda_{i+1})\f$, and it either computes a subgradient
   \f$g_{i+1} \in \partial f(\lambda_{i+1})\f$, or, if the point is
   infeasible, a constraint. */

/*--------------------------------------------------------------------------*/

   void FormD( void );

/**< The method is used within Solve() [see above], and its job is to compute
   the direction \f$d_i\f$ that appears in the formula of \f$\lambda_{i+1}\f$.
   The direction is actually saved after the variables generation because (i)
   the subgradient keeps in memory just the current direction (ii) by
   generating/removing variables the direction may quickly come to be
   deteriorated. When the variables generation is ended [see GetFiStatus() in
   FiOracle.h], SaveDir() must be called in order to update the direction.

   In addition, the deflection coefficient is computed inside FormD(). As for
   the stepsize, the computation is performed within FormD() only if the scheme
   is <em>deflection-restricted</em>. */

/*--------------------------------------------------------------------------*/

   void SaveDir( void );

/**< The method is combined with FormD() [see above] to carry out the
   direction computation. Moreover, the stepsize is determined into SaveDir()
   when the adopted scheme is <em>stepsize-restricted</em>. */

/*--------------------------------------------------------------------------*/

   void GotoLambda1( void );

/**< Move the current point to \f$\lambda_{i+1}\f$. */

/*--------------------------------------------------------------------------*/

   void FormLambda1( void );

/**< After a (succesfull) call to FormD(), sets the new (unprojected)
   tentative point \f$\breve{\lambda}_{i+1}\f$  as
   \f[
    \breve{\lambda}_{i+1} = \lambda_i - \nu_i d_i \,.
   \f]
   Remark that the point \f$\breve{\lambda}_{i+1}\f$ must be projected before
   calling FiandGi() [see above], i.e., \f$ \lambda_{i+1} = {\rm P}_{\Lambda}
   ( \breve{\lambda}_{i+1} ) \f$. */

/*--------------------------------------------------------------------------*/
/*----------------------- PROTECTED DATA STRUCTURES  -----------------------*/
/*--------------------------------------------------------------------------*/

   Index SGPar1;       ///< projection-strategy parameters
   HpNum SGPar2;       ///< incremental factor
   Index SGPar3;       ///< scheme: stepsize(deflection)-restricted
   bool SGPar4;        ///< control if \f$ \hat{\lambda}_i \f$ is kept
   Index SGPar5;       ///< seed

   Stepsize *stepsize;       ///< pointer to the Stepsize class
   Deflection *deflection;   ///< pointer to the Deflection class
   CQKnPClass *Q2KNP;        ///< pointer to the quadratic knapsack solver

   Index MaxNumVar;     ///< maximum number of variables
   Index MaxNConst;     ///< maximum number of simplex constraints
   bool BoxConst;       ///< true, if there are some box constraints

   LMRow LambdaBar;     ///< the stability center \f$\bar{\lambda}_i\f$
   HpNum FiBar;         ///< full function value at \a LambdaBar

   LMRow Lambda;        /**< the current point \f$ \lambda_i \f$
                             or the trial point \f$ \lambda_{i+1} \f$ */
   HpNum FiLambda;      ///< full function value at \a Lambda

   LMRow LambdaHat;     ///< the point \f$ \hat{\lambda}_i \f$
   HpNum FiHat;         ///< full function value at \a HLmb


   bool LHasChgd;       /**< true if Lambda has changed since the latest call
                             to FiAndGi() */
   bool LHasProj;       /**< true if \a Lambda  has projected in the current
                             iteration */

   bool KpBstL;         ///< if LambdaBest has to be kept
   HpNum FiBest;        ///< the best value of \f$f\f$ found so far
   LMRow LambdaBest;    ///< the best point found so far

   HpNum LowerBound;    ///< Lower Bound over the full function \f$f\f$
   bool TrueLB;         /**< true if LowerBound is a "true" lower bound rather
                               than just the "minus infinity" */

   SgRow Gi;            ///< Gi[ wFi ]( Lambda ), the subgradient
   SgRow dir;           ///< the direction \f$ d_i \f$

   HpNum alpha;         ///< the deflection parameter \f$ \alpha_i \f$
   HpNum step;          ///< the stepsize \f$ \nu_i \f$

   HpNum Sigma;         /**< the linearization error \f$ \sigma_i\f$ of
                             \f$ g_i\f$ at \f$\bar{\lambda}_{i}\f$ */
   HpNum Epsilon;       /**< the linearization error \f$\epsilon_i\f$ of
                             \f$d_i\f$ at \f$\bar{\lambda}_i\f$ */

   HpNum SigmaHat;      /**< the linearization error \f$ \hat{\alpha}_i \f$ of
                             \f$g_i\f$ with respect to \f$\hat{\lambda}_i\f$ */
   HpNum HatEpsilon;    /**< the linearization error \f$ \hat{\epsilon}_i \f$
                             of \f$d_i\f$ with respect to \f$\hat{\lambda}_i\f$ */

   Index_Set SGBase;    ///< the set of indices of Gi[ wFi ]( Lambda )
   Index MaxName;       ///< maximum name to be used in FiOracle->SetGiName()

   HpNum NrmGi;         ///< the (squared) subgradient's norm
   HpNum NrmDir;        ///< the (squared) direction's norm

   HpNum dGk;           ///< scalar product between \f$d_i\f$ and \f$g_i\f$
   HpNum dM1Gk;         ///< scalar product between \f$d_{i-1}\f$ and \f$g_i\f$

   bool dM1GkDone;      /**< true if the scalar product
			   \f$d_{i-1}^{\top} g_i\f$ has been computed */

   Index_Set CnstBeg;   /**< CnstBeg and CnstVol point respectively to a
                             MaxNConst-vector of Index and HpNum, while CnstNxt
                             is a MaxNumVar-vector of SIndex. Each variable is
                             associated to one and only knapsack constraint.
                             CnstVol[ i ] is the volume of the i-th knapsack
                             constraint. CnstBeg[ i ] supplies the starting of
                             the i-th constraint, i.e. the first variable
                             entering the i-th knpasack constraint.
                             CnstNxt[ j ] says the name of the next variable to
                             j entering the same knapsack constraint. If
                             CnstNxt[ j ] == MaxNumVar, then no other variable
                             enters that knapsack constraint, and if
                             CnstNxt[ j ] == -1, the variable does not belong
                             to any knapsack constraint. */

   HpRow CnstVol;       /**< the MaxNConst-vector of HpNum containing the rhs of the
                             knapsack Constraints */
   SIndex_Set CnstNxt;  /**< the MaxNumVar-vector of SIndex saying the next variable
                             appearing in a knapsack constraint [see above] */

   LMRow ub;            ///< upper bounds on the variables
   LMRow lb;            ///< lower bounds on the variables

   bool ZeroComp;       ///< true if Fi() comes with the 0-th component
   Index NItIncr;       /**< number of incremental iterations after an outer
			   iteration */
   bool InnIter;        ///< if true, the current iteration is an inner one

   vector<Index> Seq;   /**< vector containing the randomly shuffled
                           components of the function \f$f\f$ */
   bool DirPos;         /**< indicates where the direction \f$d_i\f$ in
			   located in  the oracle */
   Index_Set MultBse;

   Index CSSCntr;       ///< counter of consecutive SS
   Index CNSCntr;       ///< counter of consecutive NS

   Index CSmallStep;    ///< counter of consecutive \a short step
   bool DoSS;           ///< SS vs NS

   FiOracle::FiStatus fs;  ///< FiOracle status
   bool EmptySet;       ///< true, if the feasible set is empty

/*--------------------------------------------------------------------------*/
/*--------------------- PRIVATE PART OF THE CLASS --------------------------*/
/*--------------------------------------------------------------------------*/

 private:

/*--------------------------------------------------------------------------*/
/*-------------------------- PRIVATE METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/

   inline void Log1( void );

   inline void Log2( void );

   inline void Log3( void );

   inline void Log4( void );

/*--------------------------------------------------------------------------*/

   void EvalConst( void );

   void EvalGi( Index wFi = Inf<Index>() );

/*--------------------------------------------------------------------------*/

   bool ProjFsb( LMRow Lpt , bool & EmptySet );

   void ProjTng( SgRow Gpt , cLMRow Lpt ,
		 cIndex strt = 0 , Index stp = Inf<Index>() );

/*--------------------------------------------------------------------------*/

   void InitStepsize( void );

   void InitDeflection( void );

/*--------------------------------------------------------------------------*/

   void ChgLambdaHat( void );

/*--------------------------------------------------------------------------*/

   void UpdateSigma( void );

/*--------------------------------------------------------------------------*/

   void UpdtLowerBound( void );

/*--------------------------------------------------------------------------*/

   void MemDealloc( void );

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

 };  // end( class SubGrad )

/* @} end( group( SubGrad_CLASSES ) ) */
/*--------------------------------------------------------------------------*/
/*------------------- inline methods implementation ------------------------*/
/*--------------------------------------------------------------------------*/

inline void SubGrad::GetPar( const int wp , int &value )
{
 switch( wp ) {
  case( kSGPar1 ): value = SGPar1; break;
  case( kSGPar3 ): value = SGPar3; break;
  case( kSGPar5 ): value = SGPar5; break;
  default:  NDOSolver::GetPar( wp , value );
  }
 }  // end( SubGrad::GetPar( Index ) )

/*--------------------------------------------------------------------------*/

inline void SubGrad::GetPar( const int wp , HpNum &value )
{
 switch( wp ) {
 case( kSGPar2 ): value = SGPar2; break;
  default:  NDOSolver::GetPar( wp , value );
  }
 }  // end( SubGrad::GetPar( HpNum ) )

/*--------------------------------------------------------------------------*/

inline void SubGrad::GetPar( const int wp , bool &value )
{
 switch( wp ) {
 case( kSGPar4 ): value = SGPar4; break;
 default: throw( NDOException( "SubGrad::GetPar: unknown parameter" ) );
  }
 }  // end( SubGrad::GetPar( HpNum ) )

/*--------------------------------------------------------------------------*/

inline void SubGrad::Log1( void )
{
 #if( LOG_SG )
  if( NDOLLvl > 1 ) {
   *NDOLog << endl << "{" << SCalls << "-" << ParIter << "} ";

   *NDOLog << " FiBar = ";
    if( FiBar == Inf<HpNum>() )
     *NDOLog << " -INF";
    else
      *NDOLog << -FiBar;

    *NDOLog << " ~ FiBest = ";
    if( FiBest == Inf<HpNum>() )
     *NDOLog << " -INF";
    else {
     *NDOLog << -FiBest;

     if( LowerBound > - Inf<HpNum>() )
      *NDOLog << " ~ Gap = " << ( 100 * ( FiBest - LowerBound ) ) /
                                 max( max( LowerBound , -LowerBound ) , 1.0 ) << " %";
     }

   }
 #endif

 }  // end( SubGrad::Log1() )

/*--------------------------------------------------------------------------*/

inline void SubGrad::Log2( void )
{
 #if( LOG_SG )
  if( NDOLLvl > 1 ) {
   *NDOLog << endl << "           ";

   // information relative to the direction  - - - - - - - - - - - - - - - - -
   //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   if( SGPar4 )
    *NDOLog << " ( |D| , Epsilon , HatEpsilon ) = (" << sqrt(NrmDir)
	    << " , " << Epsilon << " , " << HatEpsilon <<  ") ~ GixD = "
	    << dGk;
   else
    *NDOLog << " ( |D| , Epsilon ) = (" << sqrt(NrmDir) << " , "
	    << Epsilon <<  ") ~ GixD = " << dGk ;
   }
 #endif

 }  // end( SubGrad::Log2() )

/*--------------------------------------------------------------------------*/

inline void SubGrad::Log3( void )
{
 #if( LOG_SG )
  if( NDOLLvl > 1 ) {
   *NDOLog << endl << "           ";

   *NDOLog << " Fi1 = ";  // print function value
   if( FiLambda == Inf<HpNum>() )
    *NDOLog << " -INF => STOP.";
   else
    *NDOLog << -FiLambda;
   }
 #endif

 }  // end( SubGrad::Log3() )

/*--------------------------------------------------------------------------*/

inline void SubGrad::Log4( void )
{
 #if( LOG_SG )
  if( NDOLLvl > 1 ) {
   *NDOLog << endl << "           ";
   if( deflection )
    *NDOLog << " ( |Gi| , Sigma ) = (" << sqrt(NrmGi) << " , "
	    << Sigma << ")";
   }
 #endif

 }  // end( SubGrad::Log4() )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#if( OPT_USE_NAMESPACES )
 };  // end( namespace NDO_di_unipi_it )
#endif

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif /* SubGrad.h included */

/*--------------------------------------------------------------------------*/
/*------------------------- End File SubGrad.h -----------------------------*/
/*--------------------------------------------------------------------------*/

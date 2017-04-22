/*--------------------------------------------------------------------------*/
/*--------------------------- File NDOSlver.h ------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Definition of the abstract base class NDOSolver, which defines the
 * interface for convex minimization algorithms.
 *
 * The algorithms embedded under this interface must be capable of minimizing
 * any proper convex (possibly NonDifferentiable, possibly non-exact, possibly
 * subject to constraints) function given only the (somehow poor) knowledge of
 * the function provided by an "oracle" as defined by the class FiOracle [see
 * FiOracle.h].
 *
 * \version 0.70
 *
 * \date 19 - 04 - 2014
 *
 * \author Antonio Frangioni \n
 *         Operations Research Group \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * Copyright &copy 2001 - 2014 by Antonio Frangioni.
 */

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __NDOSolver
 #define __NDOSolver  /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "FiOracle.h"

/*--------------------------------------------------------------------------*/
/*----------------------------- NAMESPACE ----------------------------------*/
/*--------------------------------------------------------------------------*/

#if( OPT_USE_NAMESPACES )
namespace NDO_di_unipi_it
{
 /** @namespace NDO_di_unipi_it
     The namespace NDO_di_unipi_it is defined to hold the NDOSolver and
     FiOracle classes and all the relative stuff. It comprises the
     namespace OPTtypes_di_unipi_it. */
#endif

/*--------------------------------------------------------------------------*/
/*-------------------------- CLASS NDOSolver -------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/** The class NDOSolver provides a general interface between convex
    (NonDifferentiable) optimization solvers and the applications that may
    need to use them. Also, it provides an interface for the interaction
    between the NDO solvers and the FiOracles within more complex
    optimization schemes, such as those suggested in the general notes for
    the class FiOracle [see FiOracle.h].

    The user is assumed to be familiar with the interface of the class
    FiOracle and with the relative algorithmic issues as described in
    FiOracle.h. */

class NDOSolver
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

/** Public enum which is used in (the two overloaded versions of) SetPar()
    [see below]. for selecting the algorithmic parameter to be set.
    kLastNDOParam is provided for allowing derived classes to extend the set
    of parameters they can handle while keeping the same signature for
    SetPar(): by setting kDerivedClassFirstParam = kLastNDOParam, one is
    guaranteed that the two sets of parameters will not collide. */

   enum NDOParam { kMaxItr = 0 , kMaxTme ,
		   ktStar , kEpsLin ,
		   kEInit , kEFnal , kEDcrs , kEStps ,
		   kLastNDOParam
                   };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

/** Public enum describing the possible return values of Solve() [see below].
    */

   enum NDOStatus { kOK = 0 ,
		    kUnbndd , kUnfsbl ,
		    kStopped , kLwPrcsn , kStpIter , kStpTime ,
		    kError
                    };

/*@} -----------------------------------------------------------------------*/
/*--------------------- PUBLIC METHODS OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/
/*---------------------------- CONSTRUCTOR ---------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructor
    @{ */

   inline NDOSolver( istream *iStrm = 0 );

/**< Constructor of the class. The parameter `iStrm', if provided, istaken as
   a pointer to a istream from which the algorithmic parameters for the NDO
   algorithm are sequentially read in the following order. Each parameter
   must be placed at the beginning of a separate line, max 255 carachters
   long, with all the rest of the line up to the first newline carachter
   '\n' (apart from a separating whitespace) being available for comments.
   Any line whose first carachter is '#' and any blank line is ignored.
   If 0 is passed, the file ends before reaching a given parameter, or
   some parameter is in the wrong format, each non-specified parameter is
   given a default value, shown in [] below.

   Only a few parameters are deemed to be general enough to be put in the
   basic NDOSolver class; however, the same stream can be used for passing
   the (possibly many more) algorithmic parameters to the derived classes.
   Since the constructor of NDOSolver is executed before the ones of the
   derived classes, the specific parameters for the derived classes have
   just to be found in the stream after those of the base class.

    Index MaxIter  [0] Maximum number of iterations for each call to Solve()
                   [see below]; if more iterations are required, Solve()
		   should stop returning kStpIter. 0 (the default) means
		   that the check is not done.

    HpNum MaxTime  [0] Maximum running time (in seconds) for each call to
                   Solve() [see below]; if more time required, Solve()
		   should stop returning kStpTime. 0 (the default) means
		   that the check is not done. Note that this check can
		   only be performed if timing of the code is activated
		   [see NDOTime() below].

    HpNum tStar    [1e+2]  Optimality related parameters. Proving that some
    HpNum EpsLin   [1e-6] point Lambda is optimal for a NonDifferentiable
                   Optimization problem involves finding an all-0 subgradient
		   of the function at Lambda. If an all-0 vector is found in
    the epsilon-subdifferential of Lambda, then the point is epsilon-optimal.
    Note that if the minimization problem is subject to constraints, i.e.,
    Fi() has to be minimized only on the points Lambda \in L, the latter
    being a convex set, then the above is referred to a subgradient of the
    "actual function" ( Fi + I_L )( Lambda ), where I_L is the indicator
    function of L (evaluating to 0 inside L and to +INF otherwise). In other
    words, one has to show that there exists a( enspilon-)subgradient of
    Fi() at Lambda that, *after projection on the frontier of L*, is all-0.
    A general stopping condition requires that, if EpsLin is the *relative*
    precision required, a solver can stop if it finds an epsilon-subgradient
    g at Lambda such that
              tStar * || g || + epsilon <= EpsLin * | MaxFi |
    where MaxFi is an estimate of the optimal solution value of the NDO
    problem, tStar is an estimate of the longest step that can be performed
    and || || is a norm-like function. tStar is related to the "scaling" of
    Fi(), and it can be seen as an estimate of the actual decrease that can be
    obtained by moving of an unitary step in the direction of any subgradient.
    Alternatively, the above condition can be seen as a weaker form of
        epsilon <= EpsLin * | MaxFi | / 2
	|| g || <= EpsLin * | MaxFi | / ( 2 * tStar )
    which says that the solver stops when both epsilon and || g || are
    "small", with tStar dictating what "small" means for || g ||. Note that
    each derived class can use different norm-like functions to evaluate one
    or the other of the above conditions. Also, different estimates can be
    used for MaxFi, although using the bast Fi-value found so far is pretty
    common.

    HpNum EInit    [1e-2] The evaluation of function to be minimized may be a
    HpNum EFnal    [1e-6] costly task: in many cases, it requires the solution
    HpNum EDcrs    [.95]  of a - possibly hard - optimization problem. Often,
    Index EStps    [0]    time can be saved if the function is only
                   approximately computed at the beginning of the optimization
		   process; of course, the computation should become more and
    more "exact" as the optimization proceeds. The relative precision required
    to the FiOracle [see SetPrecision() in FiOracle.h] is initially set to
    EInit, and decreased down to EFnal by multiplying it by EDcrs every EStps
    iterations. If the NDO algorithm allows it, EStps can be set to 0 meaning
    that the precision is decreased only if necessary, i.e., when it is
    impossible to proceed otherwise (the algorithm must have some way to
    detect this). The precision is kept fixed if EInit == EFnal. A NDO solver
    which is *not* capable of handling approximate computation of the function
    can ignore the values of these parameters.

   The parameters can be changed (and read) during the normal lifetime of the
   object, see SetPar() (GetPar()) below. */

/*@} -----------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
    @{ */

 /** Passes the FiOracle object to the NDOSolver class.

     This MUST be done PRIOR to ANY CALL to ANY OTHER method of the class!

     This is not done in the constructor in order to allow the user to change
     the function to be minimized during the lifetime of one NDOSolver object.
     Of course, this method must not be called within any other method of the
     class, as it (presumably) causes a complete reset of most internal data
     structures (e.g. the starting point, see below). Note that extensive
     support [see [Add/Remove]Variables() and Chg[FiV/SbG]() below] is given
     for cases where the function changes but the new function is "related" to
     the old one, so that "warm starts" can be attempted; of course, nothing of
     this kind can be expected when the FiOracle is changed with SetFiOracle().

     Passing 0 as the pointer is intended signal the NDOSolver to release as
     much memory as possible and to sit down quietly in its corner, waiting for
     a new FiOracle to become available. After a SetFiOracle( 0 ) call, NO
     OTHER method of the class (apart from the destructor) should be called
     before SetFiOracle() is called again with a non-0 argument. However,
     specific NDO solvers may allow exceptions to this rule. */

   virtual void SetFiOracle( FiOracle *Fi = 0 )
   {
    if( ( Oracle = Fi ) ) {
     NrFi = Fi->GetNrFi();
     NumVar = Fi->GetNumVar();

     SCalls = ParIter = FiEvaltns = GiEvaltns = 0;
     }
    }

/*--------------------------------------------------------------------------*/

   virtual void SetLambda( cLMRow tLambda = 0 ) = 0;

/**< Sets the starting point of the NDO algorithm; if 0 is given, or if
   the method is *not* called, some starting point is chosen by the solver
   (the all-0 vectore being one of the prime candidates). Note that Solve()
   [see below] is expected to repotimize somehow using the results of the
   previous calls, if any, and using the latest "curent point" as the
   starting point is one very typical way of doing it; SetLambda() can be
   used to reset the starting point.

   No calls to SetLambda() are allowed while Solve() is running; the starting
   point can only be changed between two calls to Solve(), and this possibly
   has a nontrivial computational cost.

   Note that tLambda *must* be feasible at least with respect to the
   non-negativity and upper bound constraints on the variables [see GetUC()
   and GetUB() in FiOracle.h]. Since one needs the FiOracle to check it,
   this method must not be called before SetFiOracle() [see above]. */

/*--------------------------------------------------------------------------*/

   virtual void KeepBestLambda( const bool KBL = true ) {};

/**< Some NDO algorithms are not descent, i.e., the point where the algorithm
   stops need not be the point which provides the best Fi-value found so far.
   If KeepBestLambda() is called, the best point found so far is kept,
   otherwise (or if KeepBestLambda( false ) is called afterwards), it is not.
   In any case, the best *Fi-value* found so far is kept [see ReadBestFiVal()
   below]. The default implementation does nothing; this is appropriate for
   a NDO algorithm which *is* of descent. */

/*--------------------------------------------------------------------------*/

   virtual inline void SetPar( const int wp , const int value );

/**< Change "int" algorithmic parameters of the NDO solver. This method is
   virtual because the derived classes may need to react to the changes of
   the parameters, and/or they may (and should) extend this method to allow
   setting the algorithmic-specific parameters.

   The enum NDOParam [see above] is used (in the obvious way) for selecting
   the parameter to be set. Note that the type of `wp' is int rather than
   NDOParam, which is not a problem as enums can always be used where ints
   are required (i.e., SetPar( NDOSolver::kMaxItr , ... ) is a valid call).
   This is done in order to allow derived classes to extend the set of
   parameters they can handle while keeping the same signature for SetPar().
   */

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

   virtual inline void SetPar( const int wp , cHpNum value );

/**< Change "float" algorithmic parameters of the NDO solver. This method is
   virtual because the derived classes may need to react to the changes of
   the parameters, and/or they may (and should) extend this method to allow
   setting the algorithmic-specific parameters.

   The enum NDOParam [see above] is used (in the obvious way) for selecting
   the parameter to be set. Note that the type of `wp' is int rather than
   NDOParam, which is not a problem as enums can always be used where ints
   are required (i.e., SetPar( NDOSolver::kMaxItr , ... ) is a valid call).
   This is done in order to allow derived classes to extend the set of
   parameters they can handle while keeping the same signature for SetPar().
   */

/*--------------------------------------------------------------------------*/
/** The class ouputs "log" information onto the ostream pointed by outs.
    lvl controls the "level of verbosity" of the code; lvl == 0 means that
    nothing at all is printed, and values larger than 0 mean increasing
    amounts of information, the specific effect of each value being derived-
    class-dependent. outs == 0 implies lvl == 0. */

   virtual void SetNDOLog( ostream *outs = 0 , const char lvl = 0 )
   {
    if( ( NDOLog = outs ) )
     NDOLLvl = lvl;
    else
     NDOLLvl = 0;
    }

/*--------------------------------------------------------------------------*/
/** SetNDOTime() allocates an OPTtimers object [see OPTUtils.h] that should
    be used for timing the calls to relevant methods of the class. The time
    can be read with NDOTime() [see below]. By default, or if
    SetNDOTime( false ) is called, no timing is done. Note that, since all the
    relevant methods of the class are pure virtual, NDOSolver can only manage
    the OPTtimers object, but it is due to derived classes to actually
    implement the timing.

    Note that time accumulates over the calls; calling SetNDOTime(), however,
    resets the counters, allowing to time specific groups of calls. */

   virtual void SetNDOTime( const bool TimeIt = true )

   {
    if( TimeIt )
     if( NDOt )
      NDOt->ReSet();
     else
      NDOt = new OPTtimers();
    else {
     delete NDOt;
     NDOt = 0;
     }
    }

/*@} -----------------------------------------------------------------------*/
/*-------------------- METHODS FOR SOLVING THE PROBLEM ---------------------*/
/*--------------------------------------------------------------------------*/
/** @name Solving the problem
    @{ */

   virtual NDOStatus Solve( void ) = 0;

/**< Tries to minimize the function provided by the FiOracle. Note that, in
   principle, a call to Solve() is intended to restart the algorithm from the
   exact state that it reached at the end of the previous call to Solve(), if
   any, unless ReSetAlg() is called [see below]. This allows the user to
   interrupt the optimization process at any time and to resume it seamlessly
   at a later stage. 

   Returns   if

   kOK       Optimization has been succesfull: a solution that is "optimal"
             (w.r.t. the current parameters settings) has been found.

   kUnbndd   Fi() is identically equal to - Infinity (or there has been an
             error in the FiOracle), i.e., Fi() has returned - INF, or it
             is unbounded below; the latter case can be detected only if a
             lower bound on the min. value of Fi is available [see
             FiOracle::GetMinusInfinity()].

   kUnfsbl   The domain of Fi is empty.

   kStopped  Solve() has been stopped for some reason, typically because
             FiOracle::GetFiStatus() returned FiOracle::kFiStop [see
             FiOracle.h]; this is possibly not an error, just maybe a request
             for a pause in the optimization process that might be restored in
             a later moment.

   kLwPrcsn  The NDO algorithm cannot proceed because the function cannot be
             computed with enough precision [see FiOracle::SetPrecision()];
             this usually means that the function has been minimized up to the
             maximum extent that is possible due to the limited precision that
             the FiOracle can provide;

   kStpIter  The max. number of iterations has been exhausted.

   kStpTime  The max. allowed running time has been spent.

   kError    There was a (typically numerical) error of some sort in Solve()
             or in the FiOracle that forced the algorithm to quit.

   Note that, whatever the exit condition be, some "current point" is usually
   available by calling ReadSol(), and its Fi-value by calling ReadFiVal()
   [see below]. */

/*--------------------------------------------------------------------------*/

   virtual void ReSetAlg( unsigned char RstLvl = 0 ) {}

/**< NDO algorithms typically have some parameters (stepsizes, multipliers
   ...) that are updated during the run according to factors such as the
   elapsed number of iterations, the outcome of the optimization process and
   so on. Also, they keep some information (previous subgradients and/or
   directions and/or constraints ...) that they use to drive the subsequent
   optimization.

   When Solve() is called more than once, the "standard interpretation" is
   that the previous optimization had not been really succesfull, so that
   further processing is needed; this situation covers the case when "slight
   changes" in the function require a reoptimization [see
   [Add/Remove]Variables() and Chg[FiV/SbG]() below].

   However, there are times when a "reset" of the internal state of the NDO
   algorithm is required. For instance, different functions defined on the
   variable space may have to be minimized consecutively; of course, this
   requires a complete reset of the internal state of the algorithm (which
   could be obtained with SetFiOracle() [see above], but possibly at a higher
   cost). Alternatively, the "slight" changes in Fi signalled by
   [Add/Remove]Variables() and Chg[FiV/SbG]() could not be slight at all,
   suggesting to restart the optimization basically from scratch. Finally, the
   NDO algorithm may "get stuck" because of some bad decision and not be able
   to terminate the optimization; sometimes, resetting some of the algorithmic
   parameters (stepsizes etc) may help it to recover.

   ReSetAlg() tells the NDO solver to reset the internal state of the NDO
   algorithm. What exactly such a "reset" means is algorithm-specific (the
   standard implementation in the base class does nothing). Indeed, the
   parameter `RstLvl' is provided to allow the caller to specify different
   "levels of reset". The exact meaning of the value of the parameter shall
   be algorithm-specific, but it is intended that larger values mean "more
   conservative" resets (i.e., resetting less things) than smaller values,
   and that 0 means "a complete reset". */

/*@} -----------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/

   virtual cLMRow ReadSol( cIndex_Set &I , Index &D ) = 0;

/**< Returns a read-only pointer to the "current point" of the NDO algorithm,
   which is usually the "best estimate" found so far of an optimal point in
   some algorithmic-specific sense. It may *not* be the point having the
   lowest Fi-value found so far, which is returned by ReadBestSol() [see
   below].

   The "format" of the returned vector depends on what is returned in I and
   D; if I == 0 (=> D == NumVar), then the pointer points to an array whose
   i-th entry is Lambda[ i ], otherwise the point has only (at most) D nonzero
   elements, and its i-th entry represents the value of Lambda[ I[ i ] ]. I[]
   must be ordered in increasing sense and Inf<Index>()-terminated. */

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/** Returns a read-only pointer to the point having the lowest Fi-value
    found so far. Note that it is necessary to call KeepBestLambda() in order
    to be sure that the ReadBestSol() works; if not, it can return 0. If
    the NDO algorithm *is* of descent, the method returns the same point as
    ReadSol() [see above]; the implementation of ReadBestSol() provided by
    the base class works for this case.

    For "format" of the returned vector, also refer to ReadSol(). */

   virtual cLMRow ReadBestSol( cIndex_Set &I , Index &D )
   {
    return( ReadSol( I , D ) );
    }

/*--------------------------------------------------------------------------*/

   virtual HpNum ReadFiVal( cIndex wFi = Inf<Index>() ) = 0;

/**< Returns the Fi-value(s) of the point returned by ReadSol() [see above].
   wFi tells the value of which "component" of Fi is required [see GetNrFi()
   in FiOracle.h]:

   - wFi == 0 requires the value of the linear 0-th component of Fi.

   - 1 <= wFi <= NrFi requires the value of the wFi-th component of Fi.

   - Inf<Index>() > wFi > NrFi requires the value of the full function Fi
     *except* the linear 0-th component.

   - wFi == Inf<Index>() requires the value of the full function Fi. */

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/** Returns the best Fi-value() found so far, i.e., those of the point
    returned by ReadBestSol(). Unlike ReadBestSol(), however, this method is
    supposed to work even if KeepBestLambda() is *not* called (too little
    storage is required here not to do that). If the NDO algorithm *is* of
    descent, the method returns the same values as ReadFiVal() [see above];
    the implementation of ReadBestFiVal() provided by the base class works
    for this case.

    For meaning of wFi, also refer to ReadFiVal(). */

   virtual HpNum ReadBestFiVal( cIndex wFi = Inf<Index>() )
   {
    return( ReadFiVal( wFi ) );
    }

/*--------------------------------------------------------------------------*/
/** This method should return true if the NDOSolver believes that the
    current solution [see ReadSol() above] is eps-optimal (relative). If
    eps == 0, then the optimality tolerances used by the solver to terminate
    are used.

    This method may be called while Solve() is running, e.g. by the FiOracle:
    if IsOptimal(), then the NDOSolver is going to stop at this iteration, so
    the FiOracle can react accordingly. */

   virtual bool IsOptimal( HpNum eps = 0 ) const
   {
    return( ( Result == kOK ) );  // a standard very naive implementation is
                                  // provided, but it should probably be
                                  // redefined in real classes
    }

/*--------------------------------------------------------------------------*/
/** Convex minimization problems have a "dual" viewpoint (which is also
    briefly described in the general notes section of FiOracle.h). This method
    should return "dual" information about the problem, if the NDO algorithm
    is capable of computing it; as this is not always the case, a default
    implementation returning "no such information is available" is provided
    by the base class.

    In general, the dual for the convex minimization problem

        (P)   min{ Fi( Lambda ) }

    is  (D)   inf{ Fi*( z ) : z = 0 }

    where "*" indicates the Fenchel's conjugate operator; in fact, it is
    well-known that

      inf{ Fi( Lambda ) } = - Fi*( 0 ),

    so that the minimization of Fi is equivalent to problem of computing the
    conjugate of Fi in 0. The standard form in which Fi is available is that
    of a "black box" or "oracle" [see FiOracle.h]; from the dual viewpoint, an
    oracle is something that takes a point Lambda in the primal space and
    returns a point z in the dual space (essentially, the space of
    subgradients of Fi) such that Lambda is a *subgradient* of Fi* in z,
    together with the value Fi*( z ). This information is all that can be
    used from a "dual" NDO algorithm in order to solve the dual problem.

    Since Fi* is a convex function, its epigraph is a convex set; each point
    (z, Fi*( z )) of the epigraph can be obtained as a convex combination of
    at most n + 1 other points of the epigraph. From the dual viewpoint, an
    NDO algorihm should produce a set of multipliers Theta[ z ] attached to
    all the subgradients/constraints (items) z generated by the oracle, such
    that

      Sum{z} z * Theta[ z ] = 0  &&  Sum{z} Fi*( z ) * Theta[ z ] = Fi*( 0 ).

    A (read-only) pointer to these Theta[] must be returned by ReadMult(),
    with 0 meaning that they are not available. The "format" of Theta[]
    depends on I: if I == 0, then Theta[] is a "dense" D-vector, i.e.,
    the multiplier of item with "name" i is found in Theta[ i ] for i = 0,
    ..., D - 1, otherwise Theta[ i ] is the multiplier of the item with
    "name" I[ i ] for i = 0, ..., D - 1. I[] must be ordered in increasing
    sense and Inf<Index>()-terminated. The "names" of the items are those
    that are passed to the FiOracle [see SetMaxName() and SetGiName() in
    FiOracle.h].

    In the Lagrangian case, a dual object x[ i ] is attachedto the item z with
    name i; x[ i ] \in X if z is a subgradient, while x[ i ] is an extreme
    ray for X if z is a constraint. Then,

      x = Sum{i \in D[]} Theta[ i ] * x[ i ]

    is an optimal solution for the "convex relaxation" (D) of the original
    problem (OP) [see the general notes section in FiOracle.h], which can be
    constructed with this dual information and the help of the FiOracle.

    When Fi is decomposable, the dual information is naturally "splitted"
    among the components of Fi, since the subgradients are. If 1 <= wFi <=
    NrFi (the number of different components, see SetFiOracle() above) then
    only the dual information relative to that component will be returned,
    i.e., I[] will only contain "names" of items corresponding to component
    wFi; otherwise all the dual information will be returned. In the
    Lagrangian case, a decomposable Fi corresponds to a separable X [see
    FiOracle.h], so that the dual information is divided among the disjoint
    components of X.

    \note For "easy" components of Fi() [see FiOracle::GetBNC() ...], a
          different type of information naturally "takes the place" of the
	  multipliers: the optimal solution x[ wFi ]* of the Lagrangian
	  problem in the current point in terms of the "original variables" of
	  the wFi-th component. Thus, calls to ReadMult( wFi ) with an "easy"
	  wFi have the following different meaning: the vector returned by the
	  method (that may be "dense" or "sparse" as in the standard case)
	  represents x[ wFi ]*, and it is therefore a
	  FiOracle::GetBNC( wFi )-vector. This information, however, can
	  *only* be accessed when calling ReadMult( wFi ) for wFi <= NrFi:
	  in "global" calls (wFi > NrFi) only the multipliers corresponding
	  to non-easy components are returned.

    \note The dual multipliers Theta[] corresponding to subgradients [of any
          given component] should be nonnegative and sum to 1---they are
	  convex combinators---while the multipliers Theta[] corresponding to
	  constraints [for any given component] need only be nonnegative.
	  There is an exception to this rule, however, which happens if the
	  NDO algorithms exploits the information provided by Lower Bounds on
	  the optimal value of Fi and/or its components Fi[ i ]; these bounds
	  have a dual meaning and therefore dual information attached to them,
	  that is returned by the separate method GetLBMult() [see below]. If
	  this happens, the "ordinary" dual multipliers returned by ReadMult()
	  and corresponding to subgradients may sum to a quantity *strictly
	  smaller* than 1.

    If Solve() returns kUnfeasible, the problem is unfeasible; this means that
    it is either *dual unbounded* or *dual empty*. In fact, the dual solution
    obtained as above is then a *feasible ascent extreme ray* for (D), that
    is c( x ) > 0 , A( x ) [<]= 0 and x' + beta * x \in X for each x' \in X
    and beta > 0. Thus, if X is nonempty then (D) is unbounded, otherwise it
    is empty. */

   virtual cHpRow ReadMult( cIndex_Set &I , Index &D ,
			    cIndex wFi = Inf<Index>() )
   {
    I = 0;
    D = 0;
    return( 0 );
    }

/*--------------------------------------------------------------------------*/
/** Some NDO algorithms may exploit the knowledge of Lower Bounds on the
    optimal value of the function [see GetLowerBound() in FiOracle.h] as a
    whole, or of its components. These bound can be thought as all-0
    subgradients, and therefore they have a dual multiplier associated; think
    to the case when the LB *is* the optimal value of the function, so that
    the dual multiplier of the (all-0 subgradient associated to the) LB is 1.
    This method returns the value of the optimal dual multiplier associated to
    the LB; if 1 <= wFi <= NrFi the LB is that on the wFi-th component of Fi,
    otherwise the LB is that on the whole function.

    \note The multipliers returned by ReadMult( wFi ) [see above] should sum
          to 1 - ReadLBMult( wFi ) - ReadLBMult(), except when wFi == 1 so
	  that two multipliers are the same and the sum must be
	  1 - ReadLBMult().

    \note For "easy" components of Fi() [see FiOracle::GetBNC() ...], it
          makes no sense to define lower bounds on the component because
	  they (if any) are already implicit in the available complete
	  description of the function; in fact, FiOracle::GetLowerBound()
	  never produces anything for "easy" components.

    If the NDO algorithm does not expolit the LBs, no dual multipliers are
    associated with them, and this method must return 0. */

   virtual HpNum ReadLBMult( cIndex wFi = Inf<Index>() )
   {
    return( 0 );
    }

/*--------------------------------------------------------------------------*/
/** Returns the number of times that Fi() has been called; if called from
    within Fi(), the present call should be excluded. Note that the class
    offers a protected field for this information, which must obviously be
    properly set by the derived classes. */

   inline Index FiEval( void ) const
   {
    return( FiEvaltns );
    }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/** Returns the number of times that NewGi() has been called; if called from
    within NewGi(), the present call should be excluded. Note that the class
    offers a protected field for this information, which must obviously be
    properly set by the derived classes. */

   inline Index GiEval( void ) const
   {
    return( GiEvaltns );
    }

/*--------------------------------------------------------------------------*/
/** Returns the number of times that Solve() has been called; if called from
    within Solve(), the present call should be included. Note that the class
    offers a protected field for this information, which must obviously be
    properly set by the derived classes. */

   inline Index NrCalls( void ) const
   {
    return( SCalls );
    }

/*--------------------------------------------------------------------------*/
/** Returns the current number of iterations. Note that this counter is *not*
    necessarily reset each time Solve() is called, as multiple calls may
    actually be necessary to "solve the given problem". Thus, ideally this
    counter is only reset when "it is clear that an entirely new problem is
    being solved", such as after a call to ReSetAlg() or when an entirely
    different oracle is passed [see SetFiOracle()]. */

   inline Index NrIter( void ) const
   {
    return( ParIter );
    }

/*--------------------------------------------------------------------------*/
/** If this method is called within any of the methods of the class that are
    "actively timed" (this depends on the subclasses), it returns the user and
    sistem time (in seconds) since the start of that method. If methods that
    are actively timed call other methods that are actively timed, this method
    returns the (...) time since the beginning of the *outer* actively timed
    method. If this method is called outside of any actively timed method, it
    returns the (...) time spent in all the previous executions of all the
    actively timed methods of the class.

    Implementing the proper calls to NDOt->Start() and NDOt->Stop() is due to
    derived classes; these should at least be placed at the beginning and at
    the end, respectively, of Solve() - that is, at least Solve() should be
    "actively timed". */

   void NDOTime( double &t_us , double &t_ss )
   {
    t_us = t_ss = 0;
    if( NDOt )
     NDOt->Read( t_us , t_ss ); 
    }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/** As NDOTime( double & , double & ) [see above], except that returns the
    total (system + user ) time. */

   double NDOTime( void )
   {
    return( NDOt ? NDOt->Read() : 0 );
    }

/*@} -----------------------------------------------------------------------*/
/*-------------- METHODS FOR READING THE DATA OF THE PROBLEM ---------------*/
/*--------------------------------------------------------------------------*/
/** @name Reading the data of the problem
    @{ */

/** Returns the current number of variables of the function; a protected
    field is offered by the base class to keep this information. */

   inline Index GetNumVar( void ) const
   {
    return( NumVar );
    }

/*--------------------------------------------------------------------------*/

   virtual inline void GetPar( const int wp , int &value );

/**< Read the current value of the "int" algorithmic parameters of the NDO
   solver [see SetPar() above]. */

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

   virtual inline void GetPar( const int wp , HpNum &value );

/**< Read the current value of the "float" algorithmic parameters of the NDO
   solver [see SetPar() above]. */

/*@} -----------------------------------------------------------------------*/
/*------------- METHODS FOR ADDING / REMOVING / CHANGING DATA --------------*/
/*--------------------------------------------------------------------------*/
/** @name Adding / removing / changing data
   @{ */

/** Adds `NNwVrs' new variables to the NDO problem. The new variables are
    added at the end of the current variable set, i.e., their "names" will be
    set to NumVar , ... , NumVar + NNwVrs - 1 (where `NumVar' means "the
    number of variables *before* the call"). The total number of variables
    after the call will be < NumVar + NNwVrs when NumVar + NNwVrs >
    FiOracle::GetMaxNumVar() [see FiOracle.h].

    IVs is a NNwVrs-vector containing the initial values for the newly created
    variables; IVs == 0 means that all the initial values are 0.

    After a call to this method, a different function, on a different variable
    space, has to be minimized. Of course, the two functions must not be
    unrelated to one another; more specifically, it is required that the new
    function Fi( Lambda , LambdaNew ) is identical to the old function
    Fi( Lambda ) when restricted to the old space, i.e. that

      Fi( Lambda , 0 ) = Fi( Lambda )  for all Lambda.

    This implies, for instance, that all the prevously collected subgradients
    are still subgradients for the restriction of the new Fi() to the subspace
    where all the new variables are zero. We also require that the every
    subgradient of the old Fi() can be extended to a subgradient of the new
    Fi() by just "filling in" the entries corresponding to the new variables
    with proper values (this is always possible e.g. for polyhedral functions).
    In other words, if the linear function

      L( Lambda ) = G * Lambda + a

    is known to minorize the old Fi(), then there must exist NewG such that

      L( Lambda , LambdaNew ) = G * Lambda + NewG * LambdaNew + a

    minorizes the new Fi(). Even better than that, the "linearization error"
    of L() at any point Lambda

      Fi( Lambda ) - L( Lambda ) >= 0

    is identical to the linearization error of the new L() at [ Lambda , 0 ]

      Fi( Lambda , 0 ) - L( Lambda , 0 ) >= 0.

    Then, the NDOSolver can use GetGi() [see FiOracle.h] to retrieve the new
    entries NewG for the newly created variables of the items "recorded" in
    the FiOracle, if it needs so, while all the already fetched information
    remains valid. Note that, if IVs == 0, the value of the function in
    the new point [ Lambda , 0 ] is also known.

    Therefore, a call to this method assumes that the FiOracle already knows
    about the new set of variables to be created. In particular, the NDOSolver
    can use GetGi() as outlined above, and it can also use GetUB() and GetUC()
    [see FiOracle.h] to retrieve the lower and upper bounds on the variables;
    also, the NDOSolver can assume that, when the method is called,
    Oracle->GetNumVar() [see FiOracle.h] already returns the number of
    variables *after* the addiction of the new NNwVrs ones. Note that the
    likely caller for AddVariables() is *the FiOracle itself*; this is one of
    the reasons why the FiOracle may need a pointer to the NDOSolver that is
    using it [see SetNDOSolver() in FiOracle.h]. If AddVariables() is not
    directly called by the FiOracle, the caller must ensure that the FiOracle
    has been properly updated *before* calling the method. Note that this
    requirement is, for obvious reasons, opposite to what is assumed for
    RemoveVariables() [see below].

    This operation is typically useful in the case of Lagrangian optimization
    where the set of relaxed constraints A( x ) [<]= b is very large, leading
    to a very-large-scale NDO problem. Such a problem could be tackled with a
    row generation scheme, i.e., working with only an "active subset" of the
    full set of constraints and revising it as necessary. In this setting,
    adding variables corresponds to inserting new relaxed constraints in the
    current active set. */

   virtual void AddVariables( Index NNwVrs , cLMRow IVs = 0 )
   {
    throw( NDOException( "AddVariables: adding variables is not supported" )
	   );
    }

/*--------------------------------------------------------------------------*/
/** Removes the variable whose "names" are contained in the vector `whch',
    that should contain `hwnmy' distinct values in the range 0 ... NumVar - 1,
    be ordered in increasing sense and be Inf<Index>()-terminated, i.e.,
    whch[ hwnmy ] must be == Inf<Index>().

    If whch == 0, *all* the variables are eliminated; in this case, hwmny
    is ignored.

    The set of variable "names" is kept contiguous, i.e., it is always the set
    0 ... NumVar - 1 for the value of NumVar *after* the call to the method;
    hence, some of the variables that are not eliminated need to be "renamed".
    This is done as follows: when variable `i' is eliminated, all variables
    with names i + 1 ... NumVar - 1 take names i ... NumVar - 2, respectively
    (i.e., the names are shifted left by one to fill the gap). If multiple
    variables are eliminated, this is repeated for each variable, starting with
    the one with smaller name (the first one in whch) upwards. Note that if the
    *last* variable is eliminated, no renaming has to be done.

    After a call to this method, a different function, on a different variable
    space, has to be minimized. Of course, the two functions must not be
    unrelated to one another; more specifically, it is required that the new
    function Fi( Lambda ) is just the restriction of the old function
    Fi( Lambda , LambdaOld ) to the subspace where all the eliminated
    variables are zero, i.e. that

      Fi( Lambda , 0 ) = Fi( Lambda )  for all Lambda.

    This implies, for instance, that the projection of all the prevously
    collected subgradients on the new space (the elimination of the entries
    corresponding to the removed variables) are still subgradients for the
    new function in the new space. In other words, if the linear function

      L( Lambda , LambdaOld ) = G * Lambda + OldG * LambdaOld + a

    is known to minorize the old Fi(), then

      L( Lambda ) = G * Lambda + a

    minorizes the new Fi(). Even better than that, the "linearization error"
    of L() at any point with LambdaOld = 0

      Fi( Lambda , 0 ) - L( Lambda , 0 ) >= 0

    is identical to the linearization error of L() at Lambda

      Fi( Lambda ) - L( Lambda ) >= 0.

    Also, note that the value of the new function in all the previously
    tested points with LambdaOld = 0 (if any) is known.

    When this method is called, the removed variables must *still be defined*
    in the FiOracle, i.e., the NDOSolver is still allowed to query information
    about the variables being removed from the FiOracle. Of course, after the
    termination of the call to RemovaVariables() the FiOracle must be updated
    to reflect the change in the variables set. Note that the likely caller for
    RemoveVariables() is *the FiOracle itself*; this is one of the reasons why
    the FiOracle may need a pointer to the NDOSolver that is using it [see
    SetNDOSolver() in FiOracle.h]. If ReomveVariables() is not directly
    called by the FiOracle, the caller must ensure that the FiOracle is also
    properly updated *after* that the method returns. Note that this
    requirement is, for obvious reasons, opposite to what is assumed for
    AddVariables() [see above].

    This operation is typically useful in the case of Lagrangian optimization,
    where it corresponds to the deletion of some of the relaxed constraints
    A( x ) [<]= b, hopefully because they have been detected to be redundant.
    */

   virtual void RemoveVariables( cIndex_Set whch = 0 , Index hwmny = 0 )
   {
    throw(
     NDOException( "RemoveVariables: removing variables is not supported" ) );
    }

/*--------------------------------------------------------------------------*/
/** This method signals to the NDO solver that there have been changes in
    the function Fi. These changes are such that the previously obtained
    information about the function is not completely useless, but it needs to
    be properly updated. If 0 <= wFi <= NrFi, the NDO solver is told that only
    the wFi-th component of Fi has changed; if Inf<Index>() > wFi > NrFi then
    all the components except the 0-th have changed, while if wFi ==
    Inf<Index>() then all the components, comprised the 0-th, have changed.

    The class of changes that are signalled by ChgFiV() are those where only
    the Fi-values are affected, but not the first-order information. In the
    Lagrangian case, they correspond to
     = changes in the objective function c( x ) of the Lagrangian problem, or
     = changes in the feasible set X of the Lagrangian problem.

    In both cases, the first-order information (subgradient/constraint)
    corresponding to a particular dual solution/dual extreme ray x that has
    been previously generated can be re-used. If c( x ) changes, the old
    information given by x is still meaningful for the new Fi provided that it
    is properly translated [see GetVal() in FiOracle.h]. For the linear 0-th
    component, this corresponds to a change of the constant `b0'. If X changes,
    x can either be feasible/an extreme ray or not; if it is, then the old item
    associated to x is still valid for Fi with no changes, otherwise it should
    probably be discarded, which is signalled by GetVal() returning - INF.

    In both cases, the changes that are needed to update the information are
    given by just one number for each dual solution x, which can be queried
    by means of the method GetVal() of class FiOracle. Of course, this means
    that the FiOracle already knows about the change. Indeed, the most likely
    caller for ChgFi() is *the FiOracle itself*; this is one of the reasons
    why the FiOracle may need a pointer to the NDOSolver that is using it [see
    SetNDOSolver() in FiOracle.h]. If ChgFiV() is not directly called by the
    FiOracle, the caller must ensure that the FiOracle has been properly 
    informed *before* calling the method. */

   virtual void ChgFiV( cIndex wFi = Inf<Index>() )
   {
    throw( NDOException( "ChgFiV: changing Fi-values is not supported" ) );
    }

/*--------------------------------------------------------------------------*/
/** This method signals to the NDO solver that there have been changes in
    the function Fi. These changes are such that the previously obtained
    information about the function is not completely useless, but it needs to
    be properly updated. If 0 <= wFi <= NrFi, the NDO solver is told that only
    the wFi-th component of Fi has changed; if Inf<Index>() > wFi > NrFi then
    all the components except the 0-th have changed, while if wFi ==
    Inf<Index>() then all the components, comprised the 0-th, have changed.

    ChgSbG() signals that the first-order information relative to the variables
    with "names" comprised between strt and min( stp , NumVar ) - 1 for the
    specified components has changed. The changes are intended not to involve
    the Fi-values, i.e., if Fi-values change together with the first-order
    information then a separate call to ChgFiV() [see above] is required. In
    the Lagrangian case, changes to the first-oredr information correspond to
    changes in the constranints `A[ h ]()'/ right hand side `b' (note that
    these changes *are* typically accompained by changes of the Fi-values).

    A call to this method assumes that the FiOracle already knows about the
    changes. In particular, the NDOSolver can use GetGi() [see FiOracle.h] to
    retrieve the new values for the changed entries of the items of the 
    specified components of Fi, if it needs so. Indeed, the likely caller for
    ChgSbG() is *the FiOracle itself*; this is one of the reasons why the
    FiOracle may need a pointer to the NDOSolver that is using it [see
    SetNDOSolver() in FiOracle.h]. If ChgSbG() is not directly called by the
    FiOracle, the caller must ensure that the FiOracle has been properly
    updated *before* calling the method. */

   virtual void ChgSbG( cIndex strt = 0 , Index stp = Inf<Index>() , 
			cIndex wFi = Inf<Index>() )
   {
    throw( NDOException( "ChgSbG: changing subgradients is not supported" ) );
    }

/*@} -----------------------------------------------------------------------*/
/*------------------------------ DESTRUCTOR --------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Destructor
    @{ */

/** Destructor of the class. Since this is an abstract base class (i.e.,
    derived classes *must* be defined), the destructor must be virtual in
    oder to ensure that "delete p;", where p is a NDOSolver* actually
    pointing to an object of a derived class, works. */

  virtual ~NDOSolver()
  {
   delete NDOt;
   }

/*@} -----------------------------------------------------------------------*/
/*-------------------- PROTECTED PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 protected:

/*--------------------------------------------------------------------------*/
/*----------------------- PROTECTED DATA STRUCTURES  -----------------------*/
/*--------------------------------------------------------------------------*/
/** @name Standard fields
    Although NDOSolver is an abstract base class, it contains some protected
    fields for holding some information that is always going to be there in
    all implementations, so that several (simple) methods of the public
    interface can be given a "standard" implementation that is going to work
    in most cases.
    @{ */

  FiOracle *Oracle;  ///< (pointer to) the oracle for Fi

  Index MaxIter;     ///< maximum number of iterations
  HpNum MaxTime;     ///< maximum time (in seconds) for each call to Solve()

  HpNum tStar;       ///< optimality related parameter: "scaling" of Fi
  HpNum EpsLin;      ///< optimality related parameter: relative precision

  HpNum EInit;       ///< precision-related parameter: initial precision
  HpNum EFnal;       ///< precision-related parameter: final precision
  HpNum EDcrs;       ///< precision-related parameter: rate of decrease
  int EStps;         ///< precision-related parameter: number of steps

  NDOStatus Result;  ///< result of the latest call to Solve()

  Index NumVar;      ///< (current) number of variables
  Index NrFi;        ///< number of components of Fi()

  Index SCalls;      ///< nuber of calls to Solve() (the current included)
  Index ParIter;     ///< nuber of iterations in this run
  Index FiEvaltns;   ///< total number of Fi() calls
  Index GiEvaltns;   ///< total number of Gi() calls


  ostream *NDOLog;   ///< the output stream object for log purposes
  char NDOLLvl;      ///< the "level of verbosity" of the log

  OPTtimers *NDOt;   ///< OPTtimer for timing purposes

/*@} -----------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

 };  // end( class NDOSolver )

/*--------------------------------------------------------------------------*/
/*------------------- inline methods implementation ------------------------*/
/*--------------------------------------------------------------------------*/

inline NDOSolver::NDOSolver( istream *iStrm  )
{
 // initialize algorithmic parameters - - - - - - - - - - - - - - - - - - - -

 DfltdSfInpt( iStrm , MaxIter , Index( 0 ) );
 DfltdSfInpt( iStrm , MaxTime , HpNum( 0 ) );

 DfltdSfInpt( iStrm , tStar  , HpNum( 100 ) );
 DfltdSfInpt( iStrm , EpsLin , HpNum( 1e-6 ) );

 DfltdSfInpt( iStrm , EInit , HpNum( 1e-2 ) );
 DfltdSfInpt( iStrm , EFnal , HpNum( 1e-6 ) );
 DfltdSfInpt( iStrm , EDcrs , HpNum( .95 ) );
 DfltdSfInpt( iStrm , EStps , int( 0 ) );

 // other little initializations- - - - - - - - - - - - - - - - - - - - - - -

 Result = kOK;

 Oracle = 0;
 NumVar = NrFi = 0;

 NDOLog = 0;
 NDOLLvl = 0;

 NDOt = 0;

 }  // end( NDOSolver::NDOSolver )

/*--------------------------------------------------------------------------*/

inline void NDOSolver::SetPar( const int wp , const int value )
{
 switch( wp ) {
  case( kMaxItr ): MaxIter = Index( value ); break;
  case( kEStps ):  EStps = value; break;
  default:         throw(
		    NDOException( "SetPar( Index ): unknown parameter" ) );
  }
 }  // end( NDOSolver::SetPar( Index ) )

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

inline void NDOSolver::SetPar( const int wp , cHpNum value )
{
 switch( wp ) {
  case( kMaxTme ): MaxTime = value; break;
  case( ktStar ):  tStar = value; break;
  case( kEpsLin ): EpsLin = value; break;
  case( kEInit ):  EInit = value; break;
  case( kEFnal ):  EFnal = value; break;
  case( kEDcrs ):  EDcrs = value; break;
  default:         throw(
		    NDOException( "SetPar( HpNum ): unknown parameter" ) );
  }
 }  // end( NDOSolver::SetPar( HpNum ) )

/*--------------------------------------------------------------------------*/

inline void NDOSolver::GetPar( const int wp , int &value )
{
 switch( wp ) {
  case( kMaxItr ): value = int( MaxIter ); break;
  case( kEStps ):  value = EStps; break;
  default:         throw(
		    NDOException( "GetPar( Index ): unknown parameter" ) );
  }
 }  // end( NDOSolver::GetPar( Index ) )

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

inline void NDOSolver::GetPar( const int wp , HpNum &value )
{
 switch( wp ) {
  case( kMaxTme ): value = MaxTime; break;
  case( ktStar ):  value = tStar; break;
  case( kEpsLin ): value = EpsLin; break;
  case( kEInit ):  value = EInit; break;
  case( kEFnal ):  value = EFnal; break;
  case( kEDcrs ):  value = EDcrs; break;
  default:         throw(
		    NDOException( "GetPar( HpNum ): unknown parameter" ) );
  }
 }  // end( NDOSolver::GetPar( HpNum ) )

/*--------------------------------------------------------------------------*/

#if( OPT_USE_NAMESPACES )
 };  // end( namespace NDO_di_unipi_it )
#endif

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif  /* NDOSlver.h included */

/*--------------------------------------------------------------------------*/
/*------------------------- End File NDOSlver.h ----------------------------*/
/*--------------------------------------------------------------------------*/

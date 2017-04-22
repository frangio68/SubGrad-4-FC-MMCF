/*--------------------------------------------------------------------------*/
/*--------------------------- File FiOracle.h ------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * 
 * Definition of the abstract base class FiOracle, which sets the interface
 * for the "black boxes" (oracles) computing [approximate] function values
 * and [epsilon-]subgradients for (constrained) NonDifferentiable
 * Optimization algorithms.
 *
 * \version 0.83
 *
 * \date 02 - 10 - 2014
 *
 * \author Antonio Frangioni \n
 *         Operations Research Group \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * Copyright &copy 2001 - 2012 by Antonio Frangioni.
 */

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef _FiOracle
 #define _FiOracle  /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "OPTtypes.h"

/*--------------------------------------------------------------------------*/
/*------------------------ NAMESPACE and USINGS ----------------------------*/
/*--------------------------------------------------------------------------*/

#if( OPT_USE_NAMESPACES )
namespace NDO_di_unipi_it
{
 /** @namespace NDO_di_unipi_it
     The namespace NDO_di_unipi_it is defined to hold the NDOSolver and
     FiOracle classes and all the relative stuff. It comprises the
     namespace OPTtypes_di_unipi_it. */

 using namespace OPTtypes_di_unipi_it;
#endif

/*--------------------------------------------------------------------------*/
/*--------------------------- CLASS NDOSolver ------------------------------*/
/*--------------------------------------------------------------------------*/
/** Forward declaration of the class NDOSolver, which is defined in
   NDOSlver.h. This class implements the interface for a generic NDO
   algorithm which can use the information provided by FiOracle to minimize
   the function. Since the interaction between the solver and the FiOracle
   can be two-sided, the FiOracle can be provided with a pointer to the
   NDOSolver that is using it [see SetNDOSolver() below], whence the need
   for this forward declaration. */

class NDOSolver;

/*--------------------------------------------------------------------------*/
/*------------------------- CLASS NDOException -----------------------------*/
/*--------------------------------------------------------------------------*/
/** Small class for exceptions. Derives from std::exception implementing the
   virtual method what() -- and since what is virtual, remember to always
   catch it by reference (catch exception &e) if you want the thing to work.
   NDOException class are thought to be of the "fatal" type, i.e., problems
   for which no solutions exists apart from aborting the program. Other kinds
   of exceptions can be properly handled by defining derived classes with
   more information. */

 class NDOException : public exception {
 public:
  NDOException( const char *const msg = 0 ) { errmsg = msg; }

  const char* what( void ) const throw () { return( errmsg ); }

 private:
  const char *errmsg;
  };

/*--------------------------------------------------------------------------*/
/*--------------------------- CLASS FiOracle -------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/** The class FiOracle provides a standard interface between NDO solvers and
    the functions that they have to minimize.

    Given a point Lambda in a finite vector space of proper dimension, the
    "oracle" must be capable of computing the value Fi( Lambda ) of the
    proper convex (possibly nondifferentiable) function Fi. The value is
    allowed to be +Infinity outside of a convex polyhedral set, i.e., the
    effective domain of Fi can be any polyhedron (possibly the whole space).
    We will indicate as Dom(Fi) te set of all points Lambda such that
    Fi( Lambda ) < + INF, i.e., the domain of Fi.

    If Fi is finite in Lambda (i.e., Lambda belongs to Dom(Fi)), the oracle
    must be capable of returning at least one of its subgradients; both
    Fi-values and subgradients may be -- to a certain extent -- subject to
    errors.

    For any point Lambda where Fi evaluates to +Infinity, the oracle must be
    capable of returning a linear constraint which is valid for the domain of
    Fi (i.e., it is satisfied in all points where Fi is finite) and it is
    strictly violated by Lambda. A special kind of constraints, i.e. "box"
    constraints 0 <= Lambda[ i ] <= u[ i ] on some of the variables
    Lambda[ i ], is usually known in advance, and special means are provided
    for informing the NDO solver about these constraints, so that it can
    guarantee feasibility at least w.r.t. those.

    Often, (epsilon-)subgradients and linear constraints, which are both
    carachterized by a real vector in the Lambda-space, will be dealt with
    indifferently, being generally referred to as "items".

    Special support is offered for the case where the function Fi to be
    minimized is "decomposable", i.e., it is given by the sum of k + 1
    independent (convex nondifferentiable) functions

       Fi( Lambda ) = \sum{h = 0 .. k} Fi[ h ]( Lambda )

    where the 0-th component is singled out because it is affine on some
    convex set, i.e.,

                            / b0 + b * Lambda  if Lambda is in Feas0
       Fi[ 0 ]( Lambda ) =  |
                            \ +Infinity        otherwise

    where Feas0 is a convex set; clearly, Dom(Fi) is a subset of Feas0
    (one may alternatively say that Fi is the sum of k + 2 functions, as
    Fi[ 0 ] itself is the sum of an affine function and of the indicator
    function of Feas0, but there is no algorithmic reason for wanting to
    do that). Note that, analogously, each individual component may evaluate
    to +Infinity outside a convex set; in other words, Dom(Fi) is the
    intersection of Dom( Fi[ h ] ) for h = 0 .. k, where Dom( Fi[ 0 ] ) is
    Feas0. Also, note that Feas0 is contained in the hyper-rectangle
    0 <= Lambda[ i ] <= u[ i ] for all variables i which have box constraints
    defined, if any.

    If the oracle is capable of providing subgradients/constraints for each
    "component" Fi[ h ], NDO algorithms can exploit this structure by working
    with this "disaggregate" information rather than with the usual
    "aggregate" information

       Gi( Lambda ) = b + \sum{ h = 1 .. k } Gi[ h ]( Lambda )

    where Gi[ h ]() is the "subgradient function" for the h-th component of
    Fi, and Gi() is the "subgradient function" for the whole Fi. Note that,
    w.l.o.g., we can assume that the special "0-th component" of Fi is a
    linear (affine) function.

    The NDO problem associated with Fi

      (P)   min{ Fi( Lambda ) }

    has a dual problem (in the general case of a decomposable Fi)

      (D)   min  \sum{ h = 1 .. k } Fi*[ h ]( z[ h ] ) + S( s )

            s.t. \sum{ h = 1 .. k } z[ h ] + s + b = 0

    where each Fi*[ h ] is the Fenchel's conjugate function of Fi[ h ],
    and S() is the support function of Feas0, i.e.

       S( s ) = sup{ s * Lambda : Lambda \in Feas0 }

    (this actually being the Fenchel's conjugate function of the indicator
    function of Feas0). Note that if Feas0 is the whole space then S() is
    + INF everythere except for s = 0, so (D) simplifies somewhat.

    Certifying the (epsilon-)optimality of a point Lambda for (P) means
    proving that 0 belongs to the (epsilon-)subdifferential of Fi in Lambda;
    this can be seen as constructing a solution z* of (D) such that
    Fi( Lambda ) = - Fi*[ h ][ z*[ h ] ]. z* is often constructed by taking
    convex combinations of the (epsilon-)subgradients found during the
    optimization process. This is particularly interesting  when Fi is a
    Lagrangian function, that is, an "original problem"

      (OP)   sup{ c( x ) : A( x ) [<]= b , x \in X }

    exists such that Fi is

     Fi( Lambda ) = sup{ c( x ) - Lambda*A( x ) : x \in X } + Lambda * b.

    We will refere to this case as the "Lagrangian case", as the above
    problem, that must be solved by the FiOracle in order to compute the
    value of Fi, is the Lagrangian relaxation of (OP) w.r.t. the constraints
    A( x ) [<]= b. For any (epsilon-)optimal solution x( Lambda ) of the
    Lagrangian relaxation,

       Gi( Lambda ) = b - A( x( Lambda ) )

    is an epsilon-subgradient of Fi in Lambda. In this case, (P) is
    equivalent (if A() and c() are linear functions) to

       (D)   sup{ c( x ) : A( x ) [<]= b , x \in conv( X ) }.

    Thus, in the Lagrangian case any convex combination of (epsilon-)
    subgradients corresponds to a point x \in conv( X ), and in this way an
    optimal solution for the "convex relaxation" (D) of (OP) can be obtained.
    A similar interpretation holds for the information that the FiOracle has
    to provide when a point Lambda outside of the domain of Fi is evaluated,
    i.e., a cutting plane separating Lambda from the domain; these cutting
    planes can be seen to be generated by extreme rays of the unbounded set X,
    so that any x \in conv( X ) is given by the convex combination of feasible
    finite solutions of the Lagrangian relaxation plus the nonnegative
    combination of extreme rays of X.

    The minimization of a (NonDifferentiable) convex function may be only one
    step of a more complicated process, which may require the (approximate)
    minimization of a family of related functions. For instance, if Fi() is a
    Lagrangian function, one may want to solve a set of NDO problems
    corresponding to different restrictions of the same Lagrangian problem,
    e.g. to derive bounds within an implicit enumeration approach to the
    original problem (OP).

    Since all the information about the function is "hidden" in the FiOracle,
    changes in the function can be "triggered" by decisions that happen inside
    the oracle. However, these changes must then be communicated somehow to
    the solver that is using the FiOracle, since typically the solver relies
    on information gathered during the optimization process to drive its
    subsequent decisions, and that information may have become unvalid.
    Pretty often, however, the previous information can be exploited,
    possibily after proper modifications, to "warm-start" the optimization of
    the new function.

    Indeed, part of the interface of the NDOSolver class provides means for
    communicating these changes [see [Add/Remove]Variables(), ChgFiV() and
    ChgSbG() in NDOSlver.h]; thus, even though the interaction between
    NDOSolver and FiOracle is mostly a master-slave one, with the NDOSolver
    acting as the master, there can be times when it is the FiOracle that
    impose changes to the NDOSolver. Note that, in order to call the methods
    of the concerned NDOSolver, the FiOracle must be given a pointer to the
    solver object [see SetNDOSOlver() below]. */

class FiOracle
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

/** Public enum describing the "status" of the FiOracle as returned by
    GetFiStatus() [see below]. */

   enum FiStatus { kFiNorm = 0,
		   kFiError ,
		   kFiChgd ,
		   kFiStop ,
		   kFiCont
                   };

/*@} -----------------------------------------------------------------------*/
/*--------------------- PUBLIC METHODS OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/
/*---------------------------- CONSTRUCTOR ---------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructor
    @{ */

/** Constructor of the class: takes no arguments, since everything that
    concerns the real evaluation of the function must be done in derived
    classes, which will have their parameters. */

   FiOracle( void )
   {
    Slvr = 0;

    NumVar = MaxName = LamBDim = 0;
    Lambda = 0;
    LamBase = 0;
    LHasChgd = true;  // NewGi() has clearly never be called yet

    FiLog = 0;
    FiLLvl = 0;

    Fit = 0;
    }

/*@} -----------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
    @{ */

/** This method is meant to pass to the FiOracle a pointer to the NDOSolver
    object that is using it, and that must be warned if the function Fi
    changes. Passing 0 means that either no solver is currently using this
    FiOracle or that the solver does not want to know about the changes in the
    function. */

   virtual void SetNDOSolver( NDOSolver *NwSlvr = 0 )
   {
    Slvr = NwSlvr;
    }

/*--------------------------------------------------------------------------*/

/** The oracle should ouput any "log" information onto the ostream pointed
    by outs. lvl controls the "level of verbosity" of the code: lvl == 0 means
    that nothing at all is printed, and values larger than 0 mean increasing
    amounts of information, the specific effect of each value being derived-
    class-dependent. outs == 0 implies lvl == 0. */

   virtual void SetFiLog( ostream *outs = 0 , const char lvl = 0 )
   {
    if( ( FiLog = outs ) )
     FiLLvl = lvl;
    else
     FiLLvl = 0;
    }

/*--------------------------------------------------------------------------*/

/** SetFiTime() allocates an OPTtimers object [see OPTtypes.h] that should be
    used for timing the calls to relevant methods of the class. The time can
    be read with FiTime() [see below]. By default, or if SetFiTime( false ) is
    called, no timing is done. Note that, since all the relevant methods ot the
    class are pure virtual, FiOracle can only manage the OPTtimers object, but
    it is due to derived classes to actually implement the timing.

    Note that time accumulates over the calls: calling SetFiTime(), however,
    resets the counters, allowing to time specific groups of calls. */

   virtual void SetFiTime( const bool TimeIt = true )
   {
    if( TimeIt )
     if( Fit )
      Fit->ReSet();
     else
      Fit = new OPTtimers();
    else {
     delete Fit;
     Fit = 0;
     }
    }

/*--------------------------------------------------------------------------*/

/** In some cases, an optimal solution of the dual problem (D) [see the
    general notes] is a mostly welcome thing; typically, this solution is
    given in terms of convex (nonnegative) multipliers which form 0 out of a
    set of subgradients (linear constraints) generated during the run of a NDO
    algorithm [see ReadMult() in NDOSlver.h]. In the Lagrangian case, the
    corresponding convex combination of the corresponding points x \in X gives
    an optimal solution for the "convex relaxation" (D) of the original problem
    (OP); this optimal solution may be a very interesting by-product of the
    optimization process (or even the real target of the whole approach).

    Knowledge of the structure of X is confined in the FiOracle class; thus, it
    is here that (some information about) the solutions x generated during the
    algorithm must be kept if an (approximately) optimal solution of (D) has
    evenctually to be produced. Also, a "naming protocol" is needed to ensure
    the identification between each subgradient/linear constraint and the
    corresponding x \in X/extreme ray of X.

    "Names" are associated with each item in the method SetGiName() [see below]
    for other purposes; these names are also used as labels for the dual
    solutions x. The maximum number of different names that will be used by the
    NDO algorithm corresponds to the maximum number of different dual solutions
    x/extreme rays that the FiOracle may have to store; this method should be
    used to communicate this number to the FiOracle, that can use it to
    properly dimension its internal data structures. The protected field
    `MaxName' of the base class FiOracle is provided for storing this
    information, as in the default implementation of the method.

    Once SetMaxName( n ) has been called, the NDO solver is required to only
    use names 0 .. n - 1 in SetGiName(). SetMaxName() can be called more than
    once to modify the setting, but this may be costly. Also, setting a max
    name smaller than the previous setting obviously discards all the dual
    information corresponding to items with names no longer allowed.

    By default no dual information is kept by the FiOracle, i.e., the maximum
    name is 0; in this case, the NDO solver need not use SetGiName().
    SetMaxName( 0 ) brings back to this situation, typically making the
    FiOracle to discard all the currently available x and deallocate some
    memory.

    Other methods are related with handling of the names and the associated
    dual solutions x, see e.g. Aggregate() or Unused() below. */

   virtual void SetMaxName( cIndex MxNme = 0 )
   {
    MaxName = MxNme;
    }

/*@} -----------------------------------------------------------------------*/
/*-------------- METHODS FOR READING THE DATA OF THE PROBLEM ---------------*/
/*--------------------------------------------------------------------------*/
/** @name Reading the data of the problem
    @{ */

/** Returns the number of variables of the function, i.e., the (maximum)
    size of the vectors to be passed to SetLambda() and SetLamBase() [see
    below].

    A default implementation is given in the base class which just returns
    the content of the protected field "NumVar". */

   virtual Index GetNumVar( void ) const
   {
    return( NumVar );
    }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/** This method provides the only explicit support -- except for the return
    value `kFiChgd' of GetFiStatus() [see below] -- of the class FiOracle for
    the case where the function Fi changes over time. In particular, if the
    number of variables of Fi may change, and a good estimate of the max
    number of variables that can ever appear is available, GetMaxNumVar() can
    be used to return this information to NDO solver, which may then be able
    to manage its memory more efficiently. Some NDO solvers may even refuse
    to handle more variables than the maximum number declared by
    GetMaxNumVar().

    No mechanism, however, is provided in the class FiOracle for changing the
    number of variables of Fi. Methods for adding and removing variables from
    the NDO problem [see [Add/Remove]Variables() in NDOSlver.h] are provided
    by the NDOSolver class, and these methods can be properly invoked by the
    FiOracle using the pointer provided by SetNDOSolver() [see above]. */

   virtual Index GetMaxNumVar( void ) const
   {
    return( GetNumVar() );
    }

/*--------------------------------------------------------------------------*/
/** GetNrFi() returns the number of independent components of Fi(); 1 is the
    minimum number, meaning that the function is not decomposable. Note that
    even non-decomposable functions has a 0-th linear (affine) part: if Fi()
    is a Lagrangian function, for instance, the linear part corresponds to the
    right hand side `b' of the relaxed constraints `A( x ) [<]= b'. In
    general, every function has a linear part, possibly with all-zero
    coefficients (and, therefore, identically zero).

    The return value of GetNrFi() allows to properly select the values of the
    parameter `wFi' in several methods, like Fi() and ***Gi() [see below];
    however, note that all the interface is carefully structured (by giving
    default values to the `wFi' parameters) in such a way that a NDO solver
    which can only deal with non-decomposable functions can easily ignore the
    information about the decomposability of Fi(), and always treat it like a
    unique function.

    A standard implementation is given for non-decomposable oracles. */

   virtual Index GetNrFi( void ) const
   {
    return( 1 );
    }

/*--------------------------------------------------------------------------*/

/**  GetMaxN() returns the number of constraints. This number does not include
           the box constraints. */

   virtual Index MaxNConst( void ) {
    return( 0 );
    }

/*--------------------------------------------------------------------------*/
/** Returns true if the function is continuous. If 1 <= NrFi <= GetNrFi()
    [see above] the return value is about the component NrFi, while if
    NrFi > GetNrFi() then the return value is about the aggregated function
    (NrFi == 0 makes no sense since a linear function is continuous, so an
    oracle may as well decide to ignore it). Note that one would expect that
    the answer for the aggregated function is true if and only it is true for
    all the individual components, although in theory this may not be the
    case.

    A standard implementation is given for oracles of continuous functions. */

   virtual bool IsFiContinuous( cIndex NrFi = Inf<Index>() )
   {
    return( true );
    }

/*--------------------------------------------------------------------------*/
/** Returns true if the function is continuous. If 1 <= NrFi <= GetNrFi()
    [see above] the return value is about the component NrFi, while if
    NrFi > GetNrFi() then the return value is about the aggregated function
    (NrFi == 0 makes no sense since a linear function is convex, so an
    oracle may as well decide to ignore it). Note that one would expect that
    the answer for the aggregated function is true if and only it is true for
    all the individual components, although in theory this may not be the
    case.

    A standard implementation is given for oracles of convex functions. */

   virtual bool IsFiConvex( cIndex NrFi = Inf<Index>() )
   {
    return( true );
    }

/*--------------------------------------------------------------------------*/
/** Returns true if the oracle is able to provide first-order information
    about the function. If 1 <= NrFi <= GetNrFi() [see above] the return value
    is about the component NrFi, while if NrFi > GetNrFi() then the return
    value is about the aggregated function (NrFi == 0 makes no sense since for
    a linear function first-order information must clearly be available, so an
    oracle may as well decide to ignore it). Note that one would expect that
    the answer for the aggregated function is true if and only it is true for
    all the individual components, although in theory this may not be the
    case.

    A standard implementation is given for oracles which are able to provide
    first-order information. */

   virtual bool HasGi( cIndex NrFi = Inf<Index>() )
   {
    return( true );
    }

/*--------------------------------------------------------------------------*/
/** Returns true if the first-order information of the function [see
    HasGi() above] is continuous, i.e., the function is differentiable. If
    1 <= NrFi <= GetNrFi() [see above] the return value is about the
    component NrFi, while if NrFi > GetNrFi() then the return value is about
    the aggregated function (NrFi == 0 makes no sense since the first-order
    information of a linear function is constant, hence continuous, so an
    oracle may as well decide to ignore it). Note that one would expect that
    the answer for the aggregated function is true if and only it is true for
    all the individual components, although in theory this may not be the
    case. Also, note that one would expect this method to return true only if
    HasGi() (for the corresponding component) returns true; however, this
    need not necessarily be the case, as one may know that a function is
    differentiable but still be unable (or unwilling) to compute its gradient.
    Still, the information that the gradient is continuous is important, in
    that the solver can then compute approximate first-order information, e.g.
    by finite differences.

    A standard implementation is given for those oracles which either do not
    provide first-order information, or have it not continuous. */

   virtual bool IsGiContinuous( cIndex NrFi )
   {
    return( false );
    }

/*--------------------------------------------------------------------------*/
/** Returns true if the oracle is able to provide second-order information
    about the function; for this to happen, (the corresponding) HasGi() also
    has to return true, that is, an oracle being able to provide second-order
    information is also necessarily able to provide first-order information.
    If 1 <= NrFi <= GetNrFi() [see above] the return value is about the
    component NrFi, while if NrFi > GetNrFi() then the return value is about
    the aggregated function (NrFi == 0 makes no sense since for a linear
    function second-order information always is the null matrix, so an oracle
    may as well decide to ignore it). Note that one would expect that the
    answer for the aggregated function is true if and only it is true for all
    the individual components, although in theory this may not be the case.

    A standard implementation is given for oracles which are not able to
    provide any form of second-order information. */

   virtual bool HasH( cIndex NrFi )
   {
    return( false );
    }

/*--------------------------------------------------------------------------*/
/** Returns true if the second-order information of the function [see
    HasH() above] is continuous, i.e., the function is twice differentiable.
    If 1 <= NrFi <= GetNrFi() [see above] the return value is about the
    component NrFi, while if NrFi > GetNrFi() then the return value is about
    the aggregated function (NrFi == 0 makes no sense since the second-order
    information of a linear function is constantly equal to the null matrix,
    hence is continuous, so an oracle may as well decide to ignore it). Note
    that one would expect that the answer for the aggregated function is true
    if and only it is true for all the individual components, although in
    theory this may not be the case. Also, note that one would expect this
    method to return true only if HasH() (for the corresponding component)
    returns true; however, this need not necessarily be the case, as one may
    know that a function is twice differentiable but still be unable (or
    unwilling) to compute its Hessian. Still, the information that the
    Hessian is continuous is important, in that the solver can then compute
    approximate second-order information, e.g. by finite differences or via
    some form of quasi-Newton iteration.

    A standard implementation is given for those oracles which either do not
    provide second-order information, or have it not continuous. */

   virtual bool IsHContinuous( cIndex NrFi )
   {
    return( false );
    }

/*--------------------------------------------------------------------------*/
/** Returns the maximum number of dual information stored in the FiOracle(),
    as set by the lastes call to SetMaxName(). A defualt implementation is
    provided which uses the protected field `MaxName'. */

   virtual Index GetMaxName( void ) const
   {
    return( MaxName );
    }

/*--------------------------------------------------------------------------*/
/** The function Fi to be minimized may be unbounded below, i.e., its
    infimum may be - INF. In this case, many NDO algorithms will *never*
    stop, as there is no way of acheiving - INF by a finite number of finite
    improvements. However, there are cases where a "finite - INF" can be
    found, that is, a finite number such that if a smaller Fi-value is found,
    then Fi can be safely declared unbounded below. If Fi is a Lagrangian
    function for a problem (D) that is not known to have a feasible solution,
    a suitable value is any Lower Bound (note that (D) is a max-problem) on
    the value c( x ) of any feasible solution x (s.t. A( x ) [<]= b and
    x \in X).

    GetMinusInfinity() returns such a value, if it can be provided, and - INF
    otherwise. In most cases, such a bound is either available from the very
    beginning or it is not available at all, so calling this method more than
    once is not likely to provide different answers *unless Fi() changes*. In
    fact, changes in Fi() (number of variables etc.) may very well impact on
    the "finite - INF" value, so that this method should be called each
    time this happens. */

   virtual HpNum GetMinusInfinity( void )
   {
    return( - Inf<HpNum>() );
    }

/*--------------------------------------------------------------------------*/

/** The subgradients returned by the FiOracle may happen to be very "sparse",
    i.e., containing very few nonzeroes; sometimes, the maximum number of
    nonzeroes is known in advance. If available, this information may be very
    useful for the NDO solver, as it could save some memory by properly
    dimensioning its data structures.

    GetMaxNZ( wFi ) should return an upper bound on the maximum number of
    nonzero elements in:

    - the subgradient of the "linear part" of Fi if wFi == 0: this should
      always be possible, as that subgradient is constant.

    - the subgradient of the wFi-th component of Fi if 1 <= wFi <= GetNrFi();

    - the aggregated subgradient *excluding* the constant part, i.e. the sum
      of all the subgradients corresponding to wFi = 1 .. GetNrFi(), if
      Inf<Index>() > wFi > GetNrFi();

    - the aggregated subgradient of Fi if wFi == Inf<Index>().

    Of course, returning the maximum number of variables is always legal.

    Important note: by returning any number < GetMaxNumVar(), the FiOracle
    gets an obligation to returning subgradients in a "sparse" form [see
    GetGi() below]. */

   virtual Index GetMaxNZ( cIndex wFi = Inf<Index>() ) const
   {
    return( GetMaxNumVar() );
    }

/*--------------------------------------------------------------------------*/
/** This method has the same meaning as GetMaxNZ() [see above], but it is
    about the sparsity of constraints rather than subgradients, as the two may
    be different. A constraint is returned when Fi [see below] returns + INF.
    Note that, how further specified in the comments to Fi(), even constraints
    are actually associated to the components of Fi.

    GetMaxCNZ( wFi ) should return an upper bound on the maximum number of
    nonzero elements in:

    - a "global" constraint (valid for Feas0 = Dom( Fi[ 0 ]) ) if wFi == 0;

    - a constraint corresponding to the wFi-th component of Fi (valid for
      Dom( Fi[ wFi ) ) if 1 <= wFi <= GetNrFi();

    - the maximum of the returned values for 1, ..., GetNrFi() if
      Inf<Index>() > wFi > GetNrFi(), i.e., the maximum nonber of nonzeroes
      in a constraint valid for any Dom( Fi[ h ] ) for any nonzero h;

    - the maximum of all the above (i.e., any constraint valid for Dom( Fi ))
      if wFi == Inf<Index>().

    Of course, returning the maximum number of variables is always legal.

    Important note: by returning any number < GetMaxNumVar(), the FiOracle
    gets an obligation to returning constraints in a "sparse" form [see
    GetGi() below]. */

   virtual Index GetMaxCNZ( cIndex wFi = Inf<Index>() ) const
   {
    return( GetMaxNumVar() );
    }

/*--------------------------------------------------------------------------*/
/** The variables of the function may be either constrained in sign, i.e.,
    required to be nonnegative, or unconstrained; if Fi is a Lagrangian
    function, for instance, constrained variables correspond to inequality
    constraints while unconstrained variables correspond to equality
    constraints.

    GetUC( i ) returns true if the variable i (0 <= i < GetNumVar()) is
    unconstrained in sign, and it returns false if it is constrained to be
    nonnegative. The status of a variable is not assumed to change over time;
    at least, no support for this is offered by the base FiOracle class, but
    a variable can be deleted and re-created with a different status by using
    the methods [Add/Remove]Variables() of NDOSolver [see NDOSlver.h].

    The default implementation of the method corresponds to all the variables
    being unconstrained. */

   virtual bool GetUC( cIndex i )
   {
    return( true );
    }

/*--------------------------------------------------------------------------*/
/** For some of the variables of the function, an upper bound on the optimal
    value may be known. If this happens for all the variables, and they are
    also all constrained in sign [see GetUC() above], then a compact set to
    which any optimal solution belong is known a priori (this knowledge may be
    *required* by some NDO solvers).

    GetUB( i ) returns the upper bound on variable i; a return value of + INF
    says that variable i has no upper bound.

    Note that a variable which is unconstrained in sign (see GetUC() above)
    can have an upper bound, as well as a variable which is constrained in
    sign can have no upper bound. See the comments on GetUC() above for the
    issue of changes over time of the "status" of a variable.

    The default implementation of the method corresponds to no variable having
    upper bounds. */

   virtual LMNum GetUB( cIndex i )
   {
    return( Inf<LMNum>() );
    }

/*--------------------------------------------------------------------------*/
/** When some of the variables are declared by the FiOracle either to be
    nonnegative [see GetUC() above] or to possess a finite upper bound [see
    GetUB() above], the NDO algorithm should take care to only provide values
    for these variables which satisfy 0 <= Lambda[ i ] <= GetUB( i ) in the
    points it intends to probe the function in [see SetLambda() below].
    However, many FiOracles may resist to small violations to these bounds.
    GetBndEps() returns a "small" number eps such that the FiOracle is
    guaranteed to correctly perform its task as long as

      - eps <= Lambda[ i ] <= ( 1 + eps ) * max( GetUB( i ) , 1 )

    while it may disbehave (e.g. enter in an infinite loop) if the bound
    constraints are more violated than that. The default implementation of the
    method is fine for FiOracles that cannot accept even the tiniest violation
    in the bound constraints they declare.

    The information returned by GetBndEps() is typically useful for another
    purpose, too: it basically says when is that a variable is "zero" for this
    particular FiOracle. Since an error of eps around zero is tolerated by the
    FiOracle, each variable that has a(n absolute) value <= eps can be safely
    assumed to be "zero". Note that the "zero" value for a variable has a
    special meaning for the functions, in particular when adding and deleting
    variables [see AddVariables() and RemoveVariables() in NDOSolver.h]. */

   virtual LMNum GetBndEps( void )
   {
    return( 0 );
    }

/*--------------------------------------------------------------------------*/

/** Returns the Lipschitz constant. If 1 <= NrFi <= GetNrFi() [see above]
    the return value is about the component NrFi, while if NrFi > GetNrFi() then
    the return value is about the aggregated function.  Note that two cases
    are special:

     - Lipschitz constant of aggregated function is less or equal to the sum
     of constants of its components;

     - NrFi == 0 is a special case since norm of b is Lipschitz constant.

    A standard implementation is given for those oracles which either do not
    provide constant Lipschitz, or have it not Lipschitz continuous. */

   virtual HpNum GetGlobalLipschitz( cIndex wFi = Inf<Index>() )

   {
    return( Inf<HpNum>() );
    }

/*--------------------------------------------------------------------------*/
/** @defgroup EasyComp Description of "easy" components of Fi()
    In some cases, the function Fi to be minimized is composed of several
    different components [see GetNrFi() above] which can be subdivided into
    two different classes:

    - "difficult" components, that are known only through a "complicated"
      oracle (the usual situation), and

    - "easy" components, for which a "compact, efficient" descritpion is
      known.

    These methods allow an NDO solver to obtain a compact description of the
    "easy" components of Fi; in particular, the components are considered easy
    if they can be described with a "compact" linear program (any polyhedral
    function can be described a linear program, but here the idea is that the
    linear program must have "few" variables and constraints). That is, the
    h-th component is "easy" if

     Fi[ h ]( Lambda ) = sup{ ( c[ h ] - Lambda * A[ h ] ) x[ h ] :
                              d[ h ] <= B[ h ] x[ h ] <= e[ h ] ,
			      l[ h ] <= x[ h ] <= u[ h ] } 

    i.e., if Fi[ h ] is the Lagrangian function of a "compact" linear program.
    A typical case is when the whole Fi is a Lagrangian function of a complex
    problem, where the Lagrangian problem decomposes into some "easy" (linear)
    and some "hard" (nonlinear/nonconvex) problems. Note that if some
    component has a linear (affine) part, all these are considered to be
    "embedded" in the linear (affine) 0-th component of Fi.

    Some NDO solvers may exploit the knowledge of a complete description of
    the "easy" components of the function; these methods are thought for
    providing the interested solver with the information.

    \note *Important*: "easy" components that reveal themselves cannot be
          considered "standard" components; that is, the "ordinary" Fi(),
	  NewGi(), SetGi(), ... [see below] machinery for computing
	  function values and subgradients is not assumed to work with easy
	  components, and the NDOSolver is *not* allowed to call them for
	  "easy" components. In other words, by declaring a component as
	  "easy" the FiOracle is basically freed by all obbligations regarding
	  it, and it's the NDOSolver's duty to use the provided information to
	  handle it (or refuse to do so and die gracefully). Clearly, whether
	  or not an "easy" component is declared as such is a decision of the
	  FiOrcale, which can then avoid to do so if the NDOSolver it's going
	  to be used with does not support it, or if it can handle it
	  (presumably) more efficiently than a generic linear solver due to
	  its knowledge of the structure of the subproblem.

    \note These methods are intended to be called only with 1 <= wFi <=
          GetNrFi(). Calls with wFi > GetNrFi() make no sense. Calls with
	  wFi == 0 would indeed make sense, as the linear 0-th component of Fi
	  is precisely one very special case of function which has a compact
	  polyhedral description (with only one variable x[ 0 ], an empty
	  B[ 0 ], l[ 0 ] == u[ 0 ] == 1 and A[ 0 ] = b). However, other means
	  are already provided for obtaining this very special descritpion
	  [see GetGi() below].

    \note The existence of a polyhedral description of the components is not
          assumed to vary over time. That is, either a component always has a
	  polyhedral description, or it never has one; in other words,
	  GetBNC( h ) *must always return the same value* throughout all the
	  life of a FiOracle object, or at least as long as the FiOracle is
	  used by one given NDO solver. Changes in the number of variables
	  might in principle happen in some cases (e.g., variables generation
	  within the polyhedral description) but this is *not* supported by
	  the current version of this interface.

    For any component h that has a polyhedral description, the NDO solver
    should in principle call these methods only once. There are exceptions,
    however. The first is if the NDO solver uses variables generation
    techniques, which try to solve a NDO problem by restricting it to a
    (small) "active set" of the variables, minimizing the function in that
    subspace and possibly revising the active set. Such a NDO solver may
    acquire the information about different parts of the matrix A[ h ] in
    different moments, as it only works with a submatrix. The second
    exception is if the NDO solver discovers that one of the components which
    have a polyhedral description is changed (e.g., because the methods
    ChgFiV() and ChgSbG() of class NDOSolver are called, see NDOSlver.h); in
    this case, it may want to call these methods again to acquire the new
    description.

    Since in most cases a polyhedral description of the components is not
    available, these methods are not "pure virtual"; rather, they are given a
    default implementation which corresponds to "no polyhedral description
    available". @{ */

/** Returns the number of variables of the "easy" linear problem which
    describes Fi[ wFi ], that is, the number of columns of the matrices
    B[ wFi ] and A[ wFi ] and the lenght of the vectors x[ wFi ], c[ wFi ],
    l[ wFi ] and u[ wFi ]. If GetBNC( wFi ) returns 0, then no polyhedral
    description of Fi[ wFi ] is available, and no calls to any other of these
    methods with the same wFi are permitted. */

   virtual Index GetBNC( cIndex wFi )
   {
    return( 0 );
    }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/** Returns the number of rows of the matrix B[ wFi ] and the lenght of the
    vectors d[ wFi ] and e[ wFi ]; GetBNR( wFi ) can return 0 if only "box"
    constraints are imposed on the variables x[ wFi ]. */

   virtual Index GetBNR( cIndex wFi )
   {
    return( 0 );
    }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/** Returns the number of nonzeroes in the matrix B[ wFi ]; this is clearly 0
    if GetBNR( wFi ) == 0. */

   virtual Index GetBNZ( cIndex wFi )
   {
    return( 0 );
    }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/** Returns a description of the matrix B[ wFi ] and of the vectors d[ wFi ],
    e[ wFi ], c[ wFi ], l[ wFi ] and u[ wFi ]. Bbeg, Bind and Bval must point
    respectively to a (GetBNC( wFi ) + 1)-vector of ints, a
    GetBNZ( wFi )-vector of ints and a GetBNZ( wFi )-vector of doubles; upon
    return, a description of B[ wFi ] has to be written in these arrays.
    Bbeg[ i ] says where the description of the i-th column of the matrix
    starts in the arrays Bind[] and Bval[], for i = 0 .. GetBNC() - 1; for
    all the indices Bbeg[ i ] <= j < Bbeg[ i + 1 ], Bval[ j ] is the value of
    (B[ h ])[ Bind[ j ] ][ i ], while all the other elements on the i-th
    column are zero. lhs and rhs must point to GetBNR( wFi )-vectors of
    doubles; upon return, the vectors d[ wFi ] and e[ wFi ] have to be written
    there. Finally, cst, lbd and ubd must point to GetBNC()-vectors of
    doubles; upon return, the vectors c[ wFi ], l[ wFi ] and u[ wFi ] have to
    be written there.

    Note that bounds, either of the constraints or on the variables, may be
    infinite. In particular, for any 0 <= i < GetBNR( wFi ), lhs[ wFi ][ i ]
    can be == - Inf<double>() to indicate that the i-th constraint has the
    form B[ wFi ][ i ] x[ wFi ] <= rhs[ wFi ][ i ], while rhs[ wFi ][ i ]
    can be == + Inf<double>() to indicate that the i-th constraint has the
    form lhs[ wFi ][ i ] <= B[ wFi ][ i ] x[ wFi ]; it is assumed that at
    least one of the two bounds is finite, for otherwise there is no i-th
    constraint. Similarly, any 0 <= j < GetBNC( wFi ), lbd[ wFi ][ j ] can
    be == - Inf<double>() and/or ubd[ wFi ][ j ] can be == + Inf<double>(),
    to indicate that x[ wFi ][ j ] has no finite lower/upper bound. In this
    case, both bounds infinite (respectively - and +) can happen, although
    some NDO solvers may require the overall feasible region of the wFi-th
    subproblem to be compact (this makes the function finite everywhere and
    therefore "more regular"). Note that, however, it is always assumed that
    lhs[ wFi ][ i ] <= rhs[ wFi ][ i ] (equality is a definite possibility,
    corresponding to an equality constraint), and lbd[ wFi ][ j ] <=
    ubd[ wFi ][ j ] (equality here being unlikely, for it corresponds to a
    fixed variable that can always be dealt with implicitly), for otherwise
    the component is "trivially empty" (constant - INF, i.e., not proper
    convex).

    The NDO solver is allowed to query only a subset of the information;
    that is, by passing 0 to some of the pointers, the oracle will not write
    the corresponding information. This can be done for all parameters
    individually, except for the first three; that is, if just one among
    Bbeg, Bind and Bval is 0, then all the other ones are to be treated
    as if they were 0, too. It if of course expected that at least one
    among lhs, rhs, cst, lbd and ubd, and/or the three Bbeg, Bind and Bval,
    be nonzero at each call. */

   virtual void GetBDesc( cIndex wFi , int *Bbeg , int *Bind , double *Bval ,
			  double *lhs , double *rhs , double *cst ,
			  double *lbd , double *ubd )
   {
    throw( NDOException( "GetBDesc: this component is not easy" ) );
    }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/** Returns the number of nonzeroes in the matrix A[ wFi ] (whose number of
    columns is GetBNC( wFi ) and whose number of rows is GetNumVar()); more
    precisely, it returns the number of nonzeroes in the submatrix of
    A[ wFi ] containing all the rows with indices between start and min( stp ,
    GetNumVar() ) - 1 (corresponding to the variables with those "names"). */

   virtual Index GetANZ( cIndex wFi , cIndex strt = 0 ,
			 Index stp = Inf<Index>() )
   {
    return( 0 );
    }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/** Return a description of the submatrix of A[ wFi ] containing all the
    rows with indices between start and min( stp , GetNumVar() ) - 1; the
    meaning of Abeg, Aind and Aval is analogous to that of Bbeg, Bind and Bval
    in GetBDesc() [see above]. */

   virtual void GetADesc( cIndex wFi , int *Abeg , int *Aind , double *Aval ,
			  cIndex strt = 0 , Index stp = Inf<Index>() )
   {
    throw( NDOException( "GetADesc: this component is not easy" ) );
    }

/*@} -----------------------------------------------------------------------*/
/** This method allows to read back the pointer to the NDOSolver that has
    been passed to the FiOracle through SetNDOSolver() [see above], if any. */

   virtual NDOSolver *GetNDOSolver( void )
   {
    return( Slvr );
    }

/*@} -----------------------------------------------------------------------*/
/*---------------------- METHODS FOR SETTING LAMBDA ------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Setting Lambda
   @{ */

/** Pass to the FiOracle the point where the function Fi() has to be
    evaluated. The vector Lmbd has to contain the values of all the variables
    in the new point (but it may be in "sparse" format, see SetLamBase()
    below). Passing Lmbd == 0 is equivalent to passing an all-zero vector;
    this is also assumed as default if SetLambda() has not been called yet.

    The FiOracle object is allowed *not* to copy the vector pointed by Lmbd
    in its own memory, but rather to keep a pointer to the original vector
    and keep working with it until a new Lmbd (possibly == 0) is set.
    However, the FiOracle object is *not* allowed to assume that the "old"
    Lmbd is still valid *at the beginning* of the call to SetLambda() where
    it is changed. This also holds, for obvious reasons, for SetLamBase()
    [see below].

    An implementation is given in the base class for those oracles which
    only need to store the pointer in the corresponding protected data
    structure (and remember that it has changed). */

   virtual void SetLambda( cLMRow Lmbd = 0 )
   {
    Lambda = Lmbd;
    LHasChgd = true;
    }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/** The "format" of the vector Lambda set by SetLambda() [see above] depends
    on the argument LamB passed to SetLamBase(). If LmbdB == 0, or
    SetLamBase() is never called, then Lambda is assumed to be in "dense
    format", i.e., the GetNumVar()-vector of the variables. If LmbdB != 0
    instead, then Lambda is in "sparse format", i.e., it contains only the
    *nonzero entries* of the vector of variables. There are LmbdBD
    (<= GetNumVar()) such entries, their indices being contained in LmbdB;
    that is, the "real" Lambda[] is

                   / Lambda[ j ]  if exists j < LmbdBD s.t. LmbdB[ j ] = i
     Lambda[ i ] = |
                   \   0          otherwise.

    Note that while all variables whose indices are not in LmbdB[] surely
    have value 0 the converse may not be true, i.e., some of the Lambda[ j ]
    with j in LmbdB[] may be have value 0. LmbdB[] must be ordered in
    increasing sense and Inf<Index>()-terminated, i.e. LmbdB[ i ] <
    LmbdB[ i + 1 ] for 0 <= i < LmbdBD and LmbdB[ LmbdBD ] == Inf<Index>().

    Once that a "base" is set with SetLamBase(), all the subsequent calls to
    SetLambda() will assume the same format for Lmbd, until the next call to
    SetLamBase(). This means that if a "sparse" Lambda has to be passed to
    the FiOracle, then SetLamBase() has to be called *before* SetLambda(), for
    otherwise the FiOracle will assume that Lambda is "dense" (or sparse with
    a wrong set of nonzeroes) and errors will occur. That is, when SetLamBase()
    is called the FiOracle is *not* allowed to assume that the previously set
    Lambda[] is still valid. However, as for SetLambda(), the FiOracle is
    allowed *not to copy the vector pointed by LmbdB in its own memory*, but
    rather to keep a pointer to the original vector and keep working with it
    until a new LmbdB (possibly == 0) is set. As for Lambda[], the FiOracle
    is *not* allowed to assume that the previously set LmbdB[] is still valid
    at the beginning of the call to SetLamBase() where the base is changed.
    That means that the caller must ensure that it *will not change* the
    content of Lambda[] and LmbdB[] before calling *any other* method of the
    class but these two; however, it is allowed to change that vectors right
    before calling SetLamBase() (or SetLambda() for the sole Lambda[]), i.e.,
    at any time between the last call to another method of the class and a
    call to this.

    An implementation is given in the base class for those oracles which only
    need to store the pointer (and lenght) in the corresponding protected data
    structures. */

   virtual void SetLamBase( cIndex_Set LmbdB = 0 , cIndex LmbdBD = 0 )
   {
    LamBase = LmbdB;
    LamBDim = LamBase ? LmbdBD : NumVar;
    }

/*--------------------------------------------------------------------------*/
/** The computation of the function Fi may be a costly task, e.g. involving
    the solution of a possibly hard optimization problem, as in the Lagrangian
    case. Hence, it may not be always possible, or just smart, to compute the
    value of the function exactly.

    SetPrecision() tells the FiOracle that the value of the function is only
    required with *relative* precision Eps (>= 0): of course, the FiOracle
    can always return values with an higher precision. If SetPrecision() is
    not called, 0 (that is, the best precision that the FiOracle is capable
    of providing anyway) is assumed.

    The return value of SetPrecision() tells the NDO solver whether or not
    the new required precision "changes the life" to the FiOracle. A call to
    SetPrecision( Eps ) is equivalent to the question "Were you, the FiOracle,
    giving me values affected by an absolute error larger than Eps? And, in
    this case, are you capable of providing me new values with a relative
    error smaller than Eps?".

    Thus, the FiOracle should return false if

    - either it was already providing values with precision at least Eps;

    - or it was providing values less accurate than that, but this is just all
      that it can do: its "precision" in computing Fi cannot be taken down to
      the required value.

    A return value of false, therefore, informs the NDO solver that if more
    accurate values than those previously provided are necessary for continuing
    the optimization process, then it should be stopped because such values are
    not available. A return value of true, instead, tells that it is perhaps
    possible to produce more accurate values. In this case, the NDO solver can
    call again Fi() [see below] *without* changing the point to retrieve the
    (hopefully) more accurate values; also, it can ask the FiOracle for more
    (epsilon)-subgradients at the current point, even if it had previously
    declared that it had had enough [see NewGi() below]. However, the NDO
    solver is *not obliged* to verify the claim of the FiOracle, that is it
    can avoid to ask for the values/subgradients in the current point.

    Furthermore, note that a return value of true does not mean that the
    FiOracle *must* be capable of returning more accurate values: it just
    means that it is *possible*. Thus, a FiOracle receiving a call to
    SetPrecision() *should not* try to calculate Fi with the higher precision
    (unless this is very inexpensive) and then answer true/false depending on
    the outcome: if there is a possibility, it should simply answer true and
    wait for the next move from the NDO solver. The reason is that the NDO
    solver may not ask for the values that the FiOracle has computed, hence
    all that work may simply be wasted.

    A final remark about the concept of "precision": there are actually two
    different ways in which a FiOracle can return "approximate" information
    about the function Fi when it is not capable of computing it exactly.

    - The first case is when the FiOracle cannot compute the real value of Fi,
      but it is capable of providing a vector that would be a 0-subgradient if
      the value were the one returned. The typical case it that of a Lagrangian
      function where the Lagrangian subproblem is not solved to optimality, but
      the suboptimal solution is treated as if it were optimal. In this case,
      a higher accuracy means that better solutions of the Lagrangian problem
      must be obtained, actually providing a different (higher) value of the
      function Fi in the same point (and, of course, new subgradients).

    - The second case is when the FiOracle is capable of computing the real
      value of Fi, or a conveniently tight approximation, but it is not
      capable of providing any 0-subgradient corresponding to that value.
      An example is, again in the Lagrangian case, if a relaxation of the
      Lagrangian problem is solved to compute a (tight) upper bound on the
      value of Fi, but no feasible solution of the Lagrangian problem
      (x \in X) can be found that attains such a value, and therefore no
      0-subgradient is available. In this case, a higher accuracy does not
      change the returned value of Fi, but rather obtains epsilon-subgradients
      with smaller epsilons.

    These two cases have different meanings in practice: for instance, in the
    first case the value of Fi is not a valid upper bound on the optimal
    objective function value of (D), while in the second it is. Of course, NDO
    algorithms that need "exact" subgradients do not work in the second case,
    but apart from that here is usually no need to distinguish among the two
    outside the FiOracle. */

   virtual bool SetPrecision( HpNum Eps )
   {
    return( false );
    }

/*@} -----------------------------------------------------------------------*/
/*------------------------ METHODS FOR COMPUTING Fi() -----------------------*/
/*--------------------------------------------------------------------------*/
/** @name Computing Fi()
   @{ */

   virtual HpNum Fi( cIndex wFi = Inf<Index>() ) = 0;

/**< This method must return the value of the function Fi to be minimized in
   the point Lmbd set by SetLambda() and SetLamBase() [see above]. Usually,
   Fi() has to be called *after* SetLambda(), the only exception being if
   Fi( < all-zero vector > ) is required, since Lambda == all-zero is the
   default assumed if SetLambda() is not called.

   wFi tells the value of which of the components of Fi() [see GetNrFi()
   above] has to be returned. Note that all functions, even those that are
   declared non-separable (GetNrFi() == 1), are treated as the sum of a
   generic function Fi[ 1 ]( Lambda ) plus a linear function b * Lambda;
   this is always possible (b can be == 0). The meaning of wFi is:

   - wFi == 0 requires the value of the linear (affine) 0-th component of Fi.

   - 1 <= wFi <= GetNrFi() requires the value of the wFi-th component of Fi,
     which must not be an "easy" one [see GetBNC() etc. above];

   - Inf<Index>() > wFi > GetNrFi() requires the value of the full function Fi
     except that of the "easy" components [see GetBNC() etc. above] and of
     the linear part, i.e. Sum{h = 1 .. GetNrFi(), h not easy} Fi( h ).

   - wFi == Inf<Index>() requires the value of the full function Fi except
     that of the "easy" components (...).

   In case an error occurs in the calculation of Fi(), this can be signalled
   immediately by returning - INF (a proper convex function never evaluates
   to - Infinity); if a "decent" value is available, an alternative is to use
   GetFiStatus() [see below]. Note that an important special case exists when
   it may be reasonable for a FiOracle to return - INF; this is when (one
   component of) Fi() is *not* a proper convex function, that is, it is
   *identically equal to - Infinity*. In the Lagrangian case, this corresponds
   to the fact that (one of) the Lagrangian problem(s) actually is *empty*.
   It is for this reason that if Fi() returns - INF, then the NDO solver
   has the right to declare Fi() unbounded below.

   Fi() may also return + INF, meaning that the point set by SetLambda() is
   outside the feasible region (the domain of the function to be minimized).
   In this case, the information returned by GetGi() [see below] has a
   different meaning. However, note that, when Fi is separable, only some of
   its components may be undefined, while other may still give valid
   subgradients; in other words, the point set by SetLambda() may belong to
   Dom( Fi[ h ] ) for some h, and not for others. Thus, as for subgradients,
   when a constraint is generated the value of wFi in the call to Fi() which
   returned INF allows the NDOSolver to distinguish to which component the
   constraint belongs to.

   Usually, there will be only one call to Fi() (for each component) for each
   call to SetLambda(); however, multiple calls are allowed if the "precision"
   changes [see SetPrecision() above]. */

/*@} -----------------------------------------------------------------------*/
/*------------- METHODS FOR READING SUBGRADIENTS / CONSTRAINTS -------------*/
/*--------------------------------------------------------------------------*/
/** @name Reading subgradients / constraints
   The four methods NewGi(), GetGi(), GetVal() and SetGiName() are the ones
   that allow a NDO solver to obtain first-order information about the
   function Fi at the point Lambda set by SetLambda() [see above].

   The typical use of these methods, in the simplest possible setting, is
   the following:

   - after having set Lambda[] with SetLambda() and computed the value of
     the function with Fi(), the NDO solver asks the FiOracle whether it
     can produce a new item by calling NewGi() (at the first call a return
     value of true is expected, but the method has to be called anyway);

   - the NDO solver then uses GetGi() and GetVal() (in *this order*, see the
     comments to GetVal()] to extract a description of the item from the
     FiOracle;

   - if the item looks "interesting" to the NDO solver, then it calls
     SetGiName() to signal to the FiOracle that some (dual) information
     about it should be retained for later use, otherwise SetGiName() is
     not called;

   - in either case, if the NDO solver wants another item it calls again
     NewGi(), until either it gets a return value of false or it has got
     enough items.

   However, these methods also have other possible uses, as described in
   their interface.
   @{ */

/*--------------------------------------------------------------------------*/
/** This method must be called to ask the FiOracle whether it can produce a
    "new" item corresponding to the point Lambda set by SetLambda() and
    SetLamBase() [see above]. The method can anly be called *after* that
    Fi() (and, therefore, SetLambda()) has been called.

    wFi tells to which of the components of Fi() [see GetNrFi() above] the
    item should correspond. More precisely, if Fi( wFi ) < INF,

    - 1 <= wFi <= GetNrFi() asks for an (epsilon-)subgradient of Fi[ wFi ],
      the wFi-th component of Fi, which must not be an "easy" one [see
      GetBNC() etc. above];

    - Inf<Index>() > wFi > GetNrFi() asks for an aggregated
      (epsilon-)subgradient excluding the constant part and all the easy
      components, i.e., the sum of all the things that would be reported by
      calls with wFi going from 1 to GetNrFi() and avoiding the easy
      components;

    - wFi == Inf<Index>() requires the aggregated (epsilon-)subgradient of
      (the not-easy part of) Fi.

    If Fi( wFi ) == INF, what is required is rather a linear constraint
    separating Lmbd from the domain of Fi[ wFi ]. The method must return true
    if a such new item could be produced, and false otherwise. Note that in
    this case wFi == 0 *does* make sense, as "global" constraints are
    associated with the 0-th component.

    \note *Important*: in general, any oracle should be able to provide on
          e new item for *at least* one of the components of Fi() for each
	  point Lambda. Many NDO solvers actually require that the oracle is
	  capable to provide a new item for *all* the components for each
	  point Lambda.

    \note The NDO solver may, in principle, require a new "aggregated"
          subgradient (e.g. by a call NewGi( Inf<Index>() )) after having
	  required some set, possibly not "covering" all the (not-easy)
	  components of Fi(), of "disaggregated" ones. In this case, the
	  intended behavior is that NewGi() should return true only if a new
	  subgradient can be obtained for *all* components of Fi(), comprised
	  those that already have produced some "disaggregated" one, so that
	  an *entirely new* aggregated subgradient can be produced. The
	  rationale is that if the NDO solver wants to produce an aggregated
	  subgradient out of partly old and partly new information, it can do
	  this by requiring new subgradients to some components individually
	  and then "mix" it with those already obtained for other components.

    Typically, when NewGi( wFi ) returns false, it will keep doing so until
    SetLambda() is called; however, new items could also be produced after a
    call to SetPrecision() [see above] which returns true.

    An implementation is given in the base class for those oracles which can
    only provide a single subgradient for each different value of Lambda[];
    *note that LHasChgd has to be set to false in GetGi() [see below]*. */

   virtual bool NewGi( cIndex wFi = Inf<Index>() )
   {
    return( LHasChgd );
    }

/*--------------------------------------------------------------------------*/

   virtual Index GetGi( SgRow SubG , cIndex_Set &SGBse ,
			cIndex Name = Inf<Index>() ,
			cIndex strt = 0 , Index stp = Inf<Index>() ) = 0;

/**< GetGi() [and GetVal(), see below] can be used to query information
   about the items. The method allows to access both the newly produced item
   corresponding to the last call to NewGi() [see above], if it has not been
   "named" yet [see SetGiName() below], and all "named" items recorded in the
   FiOracle memory, if any [see SetMaxName() above]. This is specified by the
   value of Name. If Name < n, where n is the max name set by the latest call
   to SetMaxName(), then the required information is about the item Name. Note
   that Name must be a valid item name, i.e., a name that has been previously
   used in some call to SetGiName(), and not used in any call to Deleted()
   after the last call of SetGiName() where it was used. If Name == n, then
   the required information is about the constant subgradient of the linear
   0-th component of Fi. Finally, if Name > n then the required information
   is about the newly produced "unnamed" item. Note that, given Name, it is
   immediately known to which component of Fi the associated item belongs.

   GetGi() can be used to retrieve (part of) the linear GetNumVar()-vector
   which carachterizes the item. SubG is a pointer to a vector of SgNum of
   appropriate lenght (see below) where the item have to be written. The
   indices strt and stp tell that what is required is the part of the vector
   corresponding to variables with "names" comprised between strt (included)
   and min{ stp , GetNumVar() } (excluded); the name of a variable is just
   its position in the vector Lmbd. Upon return, the required information has
   to be written in SubG; the "format" of SubG depends on what is returned in
   the read-only pointer SGBse. If SGBse == 0, then SubG[] is in "dense"
   format, i.e., SubG[ i ] is the entry of the subgradient relative to
   variable strt + i for i in 0 .. min( GetNumVar() , stp ) - strt. If
   SGBse != 0, then it must point to the vector of indices of nonzero
   entries of the subgradient (ordered in increasing sense and
   Inf<Index>()-terminated). That is, for i = strt ... stp - 1, the "real"
   SubG is

                   / SubG[ j ]    if exists j < k s.t. SGBse[ j ] = i
     SubG[ i ]  =  |
                   \   0          otherwise.

   where k, the number of nonzeroes of SubG[], is returned by GetGi(). If
   SGBse == 0, then GetGi() should return min{ stp , GetNumVar() } - strt,
   i.e., the lenght of the "dense" SubG[]. Note that, if SGBse != 0, the
   Indices in SGBse[] are in the range [strt, min{ stp , GetNumVar() }).

   According to whether the item is a subgradient or a constraint (which is
   revealed by the return value of Fi()), GetMaxNZ( wFi ) and
   GetMaxCNZ( wFi ) [see above] provide an upper bound on the maximum number
   of nonzeroes in the item; let us denote that by MxNZ. Then, the lenght of
   the vector SubG is required to be at least

      min( MxNZ , min( GetNumVar() , stp ) - strt ) ,

   i.e., the minimum possible amount that is guaranteed to contain all the
   possible nonzeroes. Hence, a non-0 SGBDim *must* always be returned if
   MxNZ < min( GetNumVar() , stp ) - strt ), as in this case the vector SubG
   is *not* long enough to accommodate (the required part of) an item in
   "dense" format. Also, note that the memory for SGBase has to be provided
   by the FiOracle, but the NDO solver is not allowed to keep the pointer
   and must copy the vector, if it so requires, because the content may
   change at the subsequent call to any method of the FiOracle.

   Allowing a NDO solver to query information about items "stored" in the
   FiOracle [see SetGiName() below] has different uses. For instance, the
   NDO solver may use variables generation techniques which try to solve a
   NDO problem by restricting it to a (small) "active set" of the variables,
   minimizing the function in that subspace and possibly revising the active
   set. Thus, at some point the NDO solver may not need some of the
   information related to an item, but this information may still be required
   later. A possibly more important use, however, is tied to the handling of
   changes in the function, which may occur in applications as those alluded
   to in the general notes. The implementation of the methods AddVariables()
   and ChgSbG() of NDOSolver [see NDOSlver.h] typically require calls to
   GetGi() to get information about some components of the items "stored" in
   the FiOracle. Note that these methods of class NDOSolver, in turn, can be
   called *by the FiOracle* itself, e.g. in the method GetFiStatus() [see
   below], so that the FiOracle may well be "prepared" when a call occurs
   simply because the caller was, ultimately, itself. */

/*--------------------------------------------------------------------------*/
/** GetVal() [and GetGi(), see above] can be used to query information
    about the items. The method allows to access both the newly produced item
    corresponding to the last call to NewGi() [see above], if it has not been
    "named" yet [see SetGiName() below], and all "named" items recorded in the
    FiOracle memory, if any [see SetMaxName() above]. This is specified by the
    value of Name. If Name < n, where n is the max name set by the latest call
    to SetMaxName(), then the required information is about the item Name.
    Note that Name must be a valid item name, i.e., a name that has been
    previously used in some call to SetGiName(), and not used in any call to
    Deleted() after the last call of SetGiName() where it was used. If
    Name == n, then the required information is about the constant subgradient
    of the linear 0-th component of Fi. Finally, if Name > n then the required
    information is about the newly produced "unnamed" item. Note that, given
    Name, it is immediately known to which component of Fi the associated item
    belongs.

    The meaning of the value returned by GetVal( Name ) is different if Name
    indicates a subgradient or a linear constraint.

    If Name indicates a constraint, then it is

      SubG * Lambda <= GetVal( Name ),

    i.e., GetVal() returns its right hand side.

    If Name indicates an epsilon-subgradient, then GetVal() returns the value
    of the corresponding linearization of Fi in Lambda, that is, the smallest
    (known) value epsilon >= 0 such that the item Name is an
    epsilon-subgradient of the wFi-th component of Fi. In the Lagrangian case,
    if x is the dual solution corresponding to the epsilon-subgradient Name,
    it is

      epsilon = c[ wFi ]( x ) - Lambda * A[ wFi ]( x ) - Fi[ wFi ]()

    i.e., how much "suboptimal" is the dual solution x for the wFi-th
    Lagrangian problem (x[ i ] is epsilon-optimal). A special case is when
    Name == n, i.e., the unique (sub)gradient of the linear 0-th component
    of Fi is queried; in this case, GetVal() returns the constant `b0'.

    \note When querying the "new" item (Name > n), both GetGi() and GetVal()
          must be called. Although in principle any order of the two calls
	  would be possible, we require that GetGi() is always called (for
	  the "new" item) *before* GetVal(); this may be convenient for the
	  FiOracle, so that, say, it knows without checking that some data
	  structures have been already initialized.

    Allowing a NDO solver to query information about items "stored" in the
    FiOracle [see SetGiName() below] is fundamental for the handling of
    changes in the function, which may occur in applications as those alluded
    to in the general notes. The implementation of the method ChgFiV() of
    NDOSolver [see NDOSlver.h] typically requires calls to GetVal() to get the
    new Fi-values/right hand sides for subgradients/constraints "stored" in
    the FiOracle. Note that, in this case, the special return value - INF
    for GetVal( Name ) tells that the item `Name' is no longer a valid one for
    the new function. Note that that method of class NDOSolver, in turn, can
    be called *by the FiOracle* itself, e.g. in the method GetFiStatus() [see
    below], so that the FiOracle may well be "prepared" when a call occurs
    simply because the caller was, ultimately, itself.

    An implementation is given in the base class for those oracles which
    only return "exact" subgradients and no constraints, and do not store
    past solutions. */

   virtual HpNum GetVal( cIndex Name = Inf<Index>() )
   {
    if( Name < MaxName )
     throw( NDOException( "GetVal: past information is not recorded" ) );

    return( 0 );
    }

/*--------------------------------------------------------------------------*/

   virtual void SetGiName( cIndex Name ) {}

/**< After that a new item has been produced, i.e., a call to NewGi()
   returned true, and that (possibly) GetGi() / GetVal() have been used to
   retrieve information about it, the NDO solver can use SetGiName() to tell
   the FiOracle that a "name" has been assigned to that new item by NDO
   solver. This name can be used later in GetGi() and GetVal() to query
   information about the item - typically, information that was not obtained
   with the first call to GetGi() / GetVal() because the function has changed
   in the meantime. The name can also be used to identify the "dual" solution
   x / extreme ray associated to that item [see SetMaxName() above] in order
   to construct a dual solution x \in conv( X ). If SetGiName() is *not*
   called, the FiOracle is authorized to "dump" the corresponding dual
   information, as this means that the item is not "interesting" for the NDO
   solver (exception to that is the constant subgradient of the linear 0-th
   component of Fi, which is not given a name but still the corresponding -
   unique - dual information may be needed).

   Note that the FiOracle can in principle disregard the "commands" issued by
   SetGiName(), for instance because it knows that no "dual" information will
   ever be required. In fact, this method is given a default implementation
   that does nothing.

   It is possible that SetGiName() is called with a given name more than once
   even if the item has not been explicitly declared as unused [see Deleted()
   below]. When a name of an already existing item is re-used in SetGiName(),
   the dual information corresponding to the old item with that name must
   just be replaced with the dual information corresponding to the new item.

   Note that, in the decomposable Lagrangian case, the choice of wFi in
   NewGi() influences how the FiOracle has to treat the dual information 
   associated to the items. In fact, a decomposable Fi corresponds to a
   feasible set X which is the cartesian product of GetNrFi() disjoint sets
   X[ h ], each being the feasible set of an independent Lagrangian problem,
   plus component-separable objective function

     c( x ) = \sum{ i = 1 .. GetNrFi() } c[ i ]( x[ i ] )

   and constraints

     A( x ) = \sum{ i = 1 .. GetNrFi() } A[ i ]( x[ i ] )

   Hence, Fi( Lmbd ) < INF means that a dual solution x = [ x[ 1 ] ..
   x[ GetNrFi() ] ] is available where each x[ i ] is a (epsilon[ i ]-)optimal
   solution of the i-th Lagrangian subproblem. Hence, if GetNrFi() separate
   calls to NewGi() are issued, then each x[ i ] is used to produce an
   (epsilon[ i ]-)subgradient for its component, so that each x[ i ] will
   receive a different name. If wFi > GetNrFi() instead, x will be considered
   as a unique dual object with an unique name. Analogously, a constraint for
   the h-th component of Fi corresponds to an unbounded ray for the feasible
   set X[ h ] of the h-th Lagrangian problem. Thus, if unbounded rays are
   generated for some components after a call to NewGi() with wFi > GetNrFi(),
   they will have to be considered as an unique unbounded ray for the whole X
   (this is always possible by padding the non-unbounded components with
   zeroes), while the same information will be considered as separate
   unbounded rays if separate calls to NewGi() are issued. Hence, with wFi
   > GetNrFi() it may well happen that the constraint is in fact associated to
   just one of the components, but the caller just does not want to know which
   component this is. This is reasonable, since Fi (as a whole) is +Infinity
   whenever at least one of its components is; in other words, the domain of
   Fi is the intersection of the domains of all its components. Thus,
   although in fact attached to components, the constraints can be thought
   to be "global", as a feasible point Lmbd must satisfy all the constraints
   corresponding to all the components. */

/*--------------------------------------------------------------------------!!

  virtual void GetH( ... )

  < We must decide exactly the signature of GetH(), whether or not we need
  something like "NewH()", whether or not we need "SetHName()", etc.
  In particular, we must decide a reasonable format for the Hessians.
  The best would be an abstract C++ matrix class which can hold different
  types of matrices (diagonal, banded, sparse, dense, ...). One may use
  some existing C++ library for that; however, being this an abstract
  interface it would be better to have it depend on the smallest possible
  number of other libraries (possibly none); if we do use some, it should
  be a "very light" library, i.e., not one such that compiling it is
  very complicated, since many of the solvers would not actually use it.
  In other words, I don't think that lapack++ is the right one.

  However, I don't have a clear opinion on that. Suggestions? */

/*@} -----------------------------------------------------------------------*/
/*-------------------- METHODS FOR READING OTHER RESULTS -------------------*/
/*--------------------------------------------------------------------------*/
/** @name Reading other results
   @{ */

/** In some cases, a Lower Bound on the minimal value of Fi is known; if Fi
    is a Lagrangian function, for instance, the objective function value
    c( x ) of any dual feasible solution x (s.t. A( x ) [<]= b and x \in X)
    gives a valid lower bound. This value can be useful for NDO solver,
    especially if it is tight; for one thing, it ensures that the function is
    not unbounded below.

    GetLowerBound() allows to retrieve this value from the FiOracle, which may
    return - INF to signal that no Lower Bound is known. Note that the value
    may change---typically, increase---over time, e.g. as better feasible
    solutions of the dual problem are generated, so that it may be worth to
    call GetLowerBound() after each call to Fi() to check if a better bound
    has become available. Of course, the value returned by GetLowerBound()
    must be no smaller than that returned by GetMinusInfinity() [see above],
    altough the two can be equal.

    Furthermore, changes in Fi() (number of variables etc.) may possibly render
    any previously obtained Lower Bound invalid, so that this method should be
    called *at least* each time Fi() changes; failure to do that could lead to
    early and incorrect termination of the optimization due to the use of an
    incorrect Lower Bound.

    Lower Bounds may be available separatedly for the wFi-th component of Fi,
    1 <= wFi <= GetNrFi(); these are returned by GetLowerBound( wFi ). Note
    that if wFi is an "easy" component then it makes no sense to provide an
    explicit lower bound, since everything is known for that component and the
    Lower Bound, if any, is already implicit in the description; thus,
    GetLowerBound( wFi ) should always return - INF if wFi is "easy".
    GetLowerBound( 0 ) also makes no sense (i.e., - INF is again a good
    return value) while GetLowerBound( wFi ) for wFi > NrFi must return a
    "global" Lower Bound for the whole Fi.

    \note While "individual" Lower Bounds are not allowed for "easy"
          components, the "global" Lower Bound is to be intended as a bound
	  on the *total* value of Fi(), that is, comprised the value of the
	  "easy" components.

    \note The existence of a "global" Lower Bound does not imply the existence
          of "individual" Lower Bounds; consider the case of the 0-th linear
	  component of Fi, which has no Lower Bound unless the domain of Fi is
	  compact. Vice-versa, "individual" Lower Bounds imply a "global" one
	  only if *every* component has one *comprised* the 0-th (i.e., the
	  0-th component is bounded on the domain, e.g. it is null or the
	  domain is bounded). */

   virtual HpNum GetLowerBound( cIndex wFi = Inf<Index>() )
   {
    return( - Inf<HpNum>() );
    }

/*--------------------------------------------------------------------------*/
/** GetFiStatus() returns the internal status of the FiOracle object. There
    are *two* envisioned uses for this method:

    - The first is as a "hook" for the FiOracle into the main loop of the NDO
      algorithm. That is, the FiOracle can expect the NDO algorithm to call
      GetFiStatus() at least once for each iteration, *before* checking for
      termination. This allows the FiOracle to properly react to things
      happening, e.g. in the case where variables are being dynamically
      generated (see GetMaxNumVar() and GetGi() above). A special setting of
      the wFi parameter is reserved for this use, see below.

    - The second is as a way for the NDOSolver to check whether or not the
      FiOracle actually has computed the latest function value(s) with
      "enough precision" according to the current setting of SetPrecision()
      [see above].

    The two uses are told apart by the value of the parameter wFi, and in
    particular by the fact that wFi > GetNrFi() or not. That is, the second
    usage correspond to the values

    - wFi == 0 requires the status of the full function Fi except that of
      the "easy" components [see GetBNC() etc. above]; this will return
      kFiNorm [see below] if all individual components return kFiNorm,
      kFiStop [see below] if at least one component returns kFiStop and
      none returns kFiError [see below], kFiError if any component returns
      kFiError; note that kFiChgd and kFiCont are not allowed as return
      values.

    - 1 <= wFi <= GetNrFi() requires the status of Fi[ wFi ], which can be
      kFiNorm, kFiStop or kFiError, but not kFiChgd and kFiCont.

    Conversely, the first usage corresponds to any value wFi > GetNrFi().
    This is analogous to wFi == 0, except that kFiChgd and kFiCont can be
    returned. This means that the function is not supposed to change during
    calls to GetFiStatus() with wFi <= GetNrFi(), while it is allowed to
    during calls to GetFiStatus() with wFi > GetNrFi(). A more detailed
    discussion of the meaning of the return status follows:

    kFiNorm    Everything is normal, the function (this particular component)
               has been correctly computed with the precision required by
	       SetPrecision().

    kFiError   There have been some problem in the FiOracle that require to
               stop the optimization (the next results from the oracle may
	       not be correct).

    kFiStop    The meaning of this return value is somwhow different for
               wFi <= GetNrFi() and wFi > GetNrFi(). In the former case, it
	       just means that the oracle is not sure to have computed the
    function value and/or subgradients "accurately enough", as required by
    the parameters set by the latest call to SetPrecision(); whether or not
    this is a problem depends on the NDOSolver, since it may just choose to
    either ignore this, or call again Fi() to "give the FiOracle some more
    time to finish the job". When wFi > GetNrFi() instead, the return value
    instructs the NDOSolver to stop. This may happen  e.g. because the
    FiOracle has a mean for detecting near-optimality of the point Lmbd but
    it cannot provide an almost-zero subgradient to directly prove this to
    the NDO solver. Alternatively, optimization may be stopped because some
    resources (e.g., computation time) have been depleted, and one must live
    with the best solution found. Finally, this may just be a request for a
    "pause" in the optimization process, that may be restarted later; this
    can be useful if the minimization of Fi is -- as it often happens -- just
    a part of a more complex process.

    kFiChgd    Since the last call to GetFiStatus(), something in the function
               Fi has changed; thus, the optimization process has to be
	       restarted. Note that for several types of changes a NDOSolver
    may be able to exploit the information acquired during the (interrupted)
    optimization of the previous function for "warm-starting" the optimization
    of the new function; such changes comprise increase and/or decrease of the
    number of variables [see [Add/Remove]Variables() in NDOSlver.h] and some
    different kinds of changes in the function Fi, such as those handled by
    the methods ChgFiV() and ChgSbG() of class NDOSolver [see NDOSlver.h].
    Note that GetFiStatus() actually is one very good point where to *call*
    all the abovementioned methods of NDOSolver. In turn, those methods of
    NDOSolver may call some other methods of FiOracle to get information about
    the changes, such as GetGi() and GetVal(), or Get[A/B]Desc() [see above].

    kFiCont    As opposed to kFiStop above, the optimization should be
               carried on even though, based on the current data, the NDO
	       solver would have detected optimality; this is useful if the
    FiOracle has chosen to hide some data (e.g., some variables) to the NDO
    solver that is planning to disclose later, but, for some reason, not at
    the current iteration. Note that this may well cause a NDO solver to
    cycle, depending on how its the stopping conditions are implemented. If
    GetFiStatus() returns kFiCont, the NDO solver should not stop before
    having called GetFiStatus() again (e.g., at the subsequent iteration)
    and having got a different return value.

   The implementation provided by the base class does nothing. */

   virtual FiStatus GetFiStatus( Index wFi = Inf<Index>() )
   {
    return( kFiNorm );
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

    Implementing the proper calls to Fit->Start() and Fit->Stop() is due to
    derived classes; these should at least be placed at the beginning and at
    the end, respectively, of Fi(), SetLambda(), GetGi() and GetVal() - that
    is, at least these methods should be "actively timed". */

   void FiTime( double &t_us , double &t_ss )
   {
    t_us = t_ss = 0;
    if( Fit ) Fit->Read( t_us , t_ss ); 
    }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/** As FiTime( double & , double & ) [see above], except that returns the
    total (system + user ) time. */

   double FiTime( void )
   {
    return( Fit ? Fit->Read() : 0 );
    }

/*@} -----------------------------------------------------------------------*/
/*------------- METHODS FOR ADDING / REMOVING / CHANGING DATA --------------*/
/*--------------------------------------------------------------------------*/
/** @name Adding / removing / changing data
   @{ */

   virtual void Deleted( cIndex i = Inf<Index>() ) {}

/**< If so instructed by SetMaxName() [see above], the FiOracle should keep
   "dual" information attached to the subgradients/constraints produced [see
   GetGi() above]. As items become useless for the NDO algorithm, the
   corresponding dual information can be discarded. If called with 0 <= i <
   n, where n is the maximum item name set by the latest call to SetMaxName()
   [see above], Deleted( i ) tells the FiOracle() that the dual information
   corresponding to the item named `i' can be discared; if i > n, then the
   information corresponding to *all items* can be discarded. The FiOracle
   may not need to do nothing in response to a call to Deleted(), as in the
   default implementation.

   Note that there is another way for the FiOracle to discover that a certain
   dual information is outdated, which is simply when its "name" is re-used
   in SetGiName() [see above]; at that point, the dual information
   corresponding to the old item with that name must be replaced with the dual
   information corresponding to the new item. */

/*--------------------------------------------------------------------------*/

   virtual void Aggregate( cHpRow Mlt , cIndex_Set NmSt , cIndex Dm ,
			   cIndex NwNm ) {}

/**< Many NDO algorithms perform operations on the subgradients that they
   obtain from the FiOracle; the most common operation is taking linear or
   convex combinations of some subgradients/constraints.

   In the Lagrangian case, this corresponds to taking linear or convex
   combinations of the corresponding points/extreme rays, generating new
   points/extreme rays. These new dual objects may have to be stored together
   with the "original" ones, if a dual solution is to be obtained in the end.

   This is precisely the meaning of this method; a new dual object must be
   computed by taking a linear combination of Dm dual objects, using as
   multipliers those found in the first Dm positions of Mlt[]. Which dual
   objects have to be used depends on NmSt: if NmSt == 0 then Mlt[ i ]
   is the multiplier of the dual object with name i, for i = 0, ..., Dm - 1,
   otherwise Mlt[ i ] is the multiplier of the dual object with name
   NmSt[ i ], for i = 0, ..., Dm - 1. NmSt[] must be ordered in increasing
   sense and Inf<Index>()-terminated. The new dual object must be given name
   NwNm. Note that *NwNm can be among the names in NmSt[]* (be < Dm if NmSt
   == 0), that is, the NDO solver may want to substitute one existing dual
   object with its linear combination with other dual objects. This method
   can be called several times consecutively, possibly using dual objects
   previously obtained by linear combination as the basis for computing new
   dual objects. Typically, the linear combination is in fact a convex one,
   so that, in the Lagrangian case, the newly obtained dual object is in
   conv( X ). When X is separable in a cartesian product (Fi is decomposable),
   convex combinations usually combine only points belonging to the same
   component.

   Note that the FiOracle may refuse to take care of the dual information
   without warning the NDO solver about it; in fact, an "empty" implementation
   is given by default to all the methods dealing with dual stuff. */

/*@} -----------------------------------------------------------------------*/
/*------------------------------ DESTRUCTOR --------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Destructor
    @{ */

   virtual ~FiOracle()

/**< Destructor of the class. Since this is an abstract base class (i.e.,
   derived classes *must* be defined), the destructor must be virtual in
   oder to ensure that "delete p;", where p is a FiOracle* actually
   pointing to an object of a derived class, works. */
   {
    delete Fit;
    }

/*@} -----------------------------------------------------------------------*/
/*-------------------- PROTECTED PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 protected:

/*--------------------------------------------------------------------------*/
/*----------------------- PROTECTED DATA STRUCTURES  -----------------------*/
/*--------------------------------------------------------------------------*/
/** @name Standard fields
    Although FiOracle is an abstract base class, it contains some protected
    fields for holding some information that is always going to be there in
    all implementations, so that several (simple) methods of the public
    interface can be given a "standard" implementation that is going to work
    in most cases.
    @{ */

  NDOSolver *Slvr;     /**< (pointer to) the NDO solver that is currently
			  using this oracle */

  Index NumVar;        ///< (current) number of variables if Fi()
  Index MaxName;       ///< maximum name to be used in SetGiName()

  cLMRow Lambda;       /**< (pointer to) the point where Fi() has to be
			  evaluated */

  cIndex_Set LamBase;  /**< (pointer to) the set of indices of nonzeroes
			  if Lambda[] is in "sparse" format. */
  Index LamBDim;       ///< length of LamBase[]
  bool LHasChgd;       /**< true if Lambda has changed since the last call
			  to NewGi() */

  ostream *FiLog;      ///< the output stream object for log purposes
  char FiLLvl;         ///< the "level of verbosity" of the log

  OPTtimers *Fit;      ///< OPTtimer for timing purposes

/*@} -----------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

 };  // end( class FiOracle )

/*--------------------------------------------------------------------------*/

#if( OPT_USE_NAMESPACES )
 };  // end( namespace NDO_di_unipi_it )
#endif

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif  /* FiOracle.h included */

/*--------------------------------------------------------------------------*/
/*------------------------- End File FiOracle.h ----------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*--------------------------- File Stepsize.h ------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 *
 * Definition of the abstract <tt>Stepsize</tt> class, which defines the
 * interface for computing the stepsize rule (SR) used by the
 * <tt>SubGrad</tt> class [see SubGrad.h].
 *
 * \version 1.00
 *
 * \date 27 - 10 - 2015
 *
 * \author Antonio Frangioni \n
 *         Operations Research Group \n
 *         Department of Informatics  \n
 *         University of Pisa \n \n
 *
 * \author Enrico Gorgone \n
 *         Graphs and Mathematical Optimization (GOM) Group \n
 *         Department of Informatics \n
 *         Free University of Brussels \n \n
 *
 * Copyright &copy 2001 - 2015 by Antonio Frangioni, Enrico Gorgone.
 */

/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __Stepsize
 #define __Stepsize /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "SubGrad.h"
#include "OPTvect.h"

/*--------------------------------------------------------------------------*/
/*----------------------------- NAMESPACE ----------------------------------*/
/*--------------------------------------------------------------------------*/

#if( OPT_USE_NAMESPACES )
namespace NDO_di_unipi_it
{
#endif

/*--------------------------------------------------------------------------*/
/*-------------------------- CLASS Stepsize --------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/** @defgroup Stepsize_CLASSES Classes in Stepsize.h
    @{ */

/** The class Stepsize provides an interface for the stepsize rule, to be
    used in the SubGrad solver [see SubGrad.h].

    The main functions of the class are NewStep() and GetStepsize(). While
    NewStep() sets up a new step, the value \f$\nu_i\f$ is passed to the
    SubGrad solver by means of a call to GetStepsize().

    The user must extend the class to one or more stepsize rules following
    this interface.  */

class Stepsize
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
	@{ */

/** Constructor of the class. The parameter `iStrm', if provided, is taken
   as a pointer to a istream from which the algorithmic parameters for the
   stepsize formula are sequentially read in the following order. Each must
   be placed at the beginning of a separate line, max 255 characters long,
   with all the rest of the line up to the first newline character '\n'
   (apart from a separating whitespace) being available for comments. Any
   line whose first character is '#' and any blank line is ignored.
   If 0 is passed, the file ends before reaching a given parameter, or
   some parameter is in the wrong format, each non-specified parameter is
   given a default value, shown in [] below.

   Only a few parameters are deemed to be general enough to be put in the
   basic Stepsize class; however, the same stream can be used for passing
   the (possibly many more) algorithmic parameters to the derived classes.
   Since the constructor of Stepsize is executed before the ones of the
   derived classes, the specific parameters for the derived classes have just
   to be found in the stream after those of the base class.

   HpNum LpsFct [1] a general scaling factor (for instance, it is the 
                    Lipschitz factor in the target value SR)  */

   Stepsize( SubGrad *slvr , istream *iStrm = 0 )
   {
    if( ! slvr )
     throw NDOException("Stepsize::Stepsize: no subgradient solver ");
    Solver = slvr;

    if( iStrm )
     DfltdSfInpt( iStrm , LpsFct , HpNum( 1 ) );

    STPLog = 0;
    STPLLvl = 0;

    }  // end( Stepsize::Stepsize() )

/*@} -----------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
   @{ */

/** The class outputs "log" information onto the ostream pointed by outs.
   lvl controls the "level of verbosity" of the code: lvl == 0 means that
   nothing at all is printed, and values larger than 0 mean increasing
   amounts of information, the specific effect of each value being derived-
   class-dependent. outs == 0 implies lvl == 0. */

   virtual void SetSTPLog( ostream *outs = 0 , const char lvl = 0 )
   {
    if( ( STPLog = outs ) )
     STPLLvl = lvl;
    else
     STPLLvl = 0;
    }

/*--------------------------------------------------------------------------*/
/** The method initializes the stepsize rule. */

   virtual void Format( void )
   {
    Beta = MaxBeta = Inf<HpNum>();
    FiLev = - Inf<HpNum>();
    }

/*@}------------------------------------------------------------------------*/
/*-------------------- METHODS FOR STEPSIZE COMPUTATION --------------------*/
/*--------------------------------------------------------------------------*/
/** @name Stepsize computation
    @{ */

/** This method must be called before GetStepsize() [see below]. This is
   indeed the core of every derived class, producing a new step.

   If the derived class implements a target value stepsize rule, the method
   must compute both the scalar \f$ \beta_i \f$ and the level
   \f$ f^{lev}_i \f$.

   Typically, the previous step will be unavailable after the call to
   NewStep() [see GetStepsize()].  */

   virtual void NewStep( void ) = 0;

/*--------------------------------------------------------------------------*/
/** Returns true if the stepsize requires the scalar product
   \f$ d_{i-1}^{\top} g_i\f$. The default implementation don't require this
   value. */

   virtual bool NeedsdkM1Gk( void )
   {
    return( false );
    }

/*@} -----------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/
/** @name Reading the solution
    @{ */

/** This function must be called after NewStep() [see above]. It returns the
    step value \f$\nu_i\f$.
 
    The default implementation offers support to the <em> target value SR </em>.
    For that rule the computation is based on the scalar \f$ \beta_i \f$ and the
    level \f$ f^{lev}_i \f$ [see NewStep()], i.e. the stepsize has the form of:
    \f[
       \nu_i = \beta_i ( \, f_i - f^{lev}_i \, ) / \|d_i\|^2 .
    \f]
 
    If true is passed, is required the stepsize for an inner iteration, otherwise
    for an outer iteration.  */

    virtual HpNum GetStepsize( bool StepIsIncr = false )
    {
     HpNum denomin, step;
     const HpNum STEpZro = 1e-8;
        
     if( StepIsIncr ) {
      HpNum LipCns = GetOracle()->GetGlobalLipschitz();
      denomin = LipCns * LipCns * HpNum( LpsFct * GetNItIcr() );
      }
     else
      denomin = GetDNorm();
        
     if( denomin <= ( STEpZro * max( ABS( ReadFiBar() ) , HpNum(1) ) ) )
      step = 1;
     else {
      HpNum Beta_ = Beta;
      if( Beta_ >= min( 2.0 , MaxBeta ) )
       Beta_ = min( 2.0 - STEpZro , MaxBeta );
      step = Beta_ * ( ReadFiBar() - FiLev ) / denomin;
      }
        
     return( step );
     } // end( GetStepsize() )
    
/*--------------------------------------------------------------------------*/
/** This function must be called after NewStep() [see above].

   The default implementation returns the level \f$ f^{lev}_i \f$ for a
   target value stepsize rule [see below]. It should throw exception for
   other forms of the stepsize rule. */

   virtual HpNum GetLev( void )
   {
    return( FiLev );
    }

/*--------------------------------------------------------------------------*/
/** This function must be called after NewStep() [see above].

   The default implementation returns the scalar \f$ \beta_i \f$ for a
   target value stepsize rule [see below]. It should throw exception for
   other forms of the stepsize rule. */

   virtual HpNum GetBeta( void )
   {
    return( Beta );
    }

/*@} -----------------------------------------------------------------------*/
/*------------- METHODS FOR ADDING / REMOVING / CHANGING DATA --------------*/
/*--------------------------------------------------------------------------*/
/** @name Adding / removing / changing data
    @{ */

/** Changes the <em>maximum value</em> of \f$ \beta_i \f$ for a <em>target
   value stepsize rule</em>.

   Typically, the method implements the safe rule of the stepsize-restricted
   approach. Passing the deflection coefficient \f$ \alpha_i \f$ does exactly
   this job. */

   virtual void SetMaxBeta( const HpNum alpha )
   {
    MaxBeta = alpha;
    }

/*@}------------------------------------------------------------------------*/
/*------------------------------ DESTRUCTOR --------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Destructor
    @{ */

   virtual ~Stepsize() {}

/*@} -----------------------------------------------------------------------*/
/*---------------------- PROTECTED PART OF THE CLASS -----------------------*/
/*--------------------------------------------------------------------------*/

  protected:

/*--------------------------------------------------------------------------*/
/*-------------------------- PROTECTED METHODS -----------------------------*/
/*--------------------------------------------------------------------------*/
/** It tries to provide a new target level \f$ f^{lev}_i \f$. It true is
   returned, the level has been changed. By default, it returns false. */

   virtual bool UpdateTargetLevel( void )
   {
    return false;
    }

/*--------------------------------------------------------------------------*/
/** @name Protected methods to read from Subgrad

    These methods are used to provide access to the information stored into
    the SubGrad object. Stepsize is a "friend" of SubGrad and therefore it
    can read its protected data structures, but classes derived from
    Stepsize are not friend of SubGrad and cannot. This is why these methods
    are defined in the base class (and implemented in SubGrad.C).
    @{ */

/*--------------------------------------------------------------------------*/
/** Returns the deflection coefficient \f$\alpha_i\f$. */

   HpNum GetCoeffDefl( void );

/*--------------------------------------------------------------------------*/
/** Returns the pointer to FiOracle. Thus, stepsize can ask directly the
   FiOracle object for the function information. */

   FiOracle* GetOracle( void );

/*--------------------------------------------------------------------------*/
/** Returns the norm of the subgradient \f$ g_i \f$. */

   cHpNum GetGiNorm( void );

/*--------------------------------------------------------------------------*/
/** Returns the norm of the direction \f$ d_i \f$. */

   cHpNum GetDNorm( void );

/*--------------------------------------------------------------------------*/
/** Returns the scalar product \f$ g_i^{\top} d_i\f$. */

   cHpNum GetdGk( void );

/*--------------------------------------------------------------------------*/
/** Returns the scalar product \f$ g_i^{\top} d_{i-1}\f$. */

  HpNum GetdkM1Gk( void );

/*--------------------------------------------------------------------------*/
/** Returns NItIncr, the parameter for incremental iterations [see SubGrad.h].
 */

   Index GetNItIcr( void );

/*--------------------------------------------------------------------------*/
/** Returns FiLambda [ see SubGrad.h ]. */

   cHpNum ReadFkVal( void );
    
/*--------------------------------------------------------------------------*/
/** Returns FiBar [ see SubGrad.h ]. */
    
   cHpNum ReadFiBar( void );

/*@}------------------------------------------------------------------------*/
/*----------------------- PROTECTED DATA STRUCTURES  -----------------------*/
/*--------------------------------------------------------------------------*/
  
  SubGrad *Solver;  ///< (pointer to) the SubGrad solver
  ostream *STPLog;  ///< the output stream object
  char STPLLvl;     ///< the "level of verbosity"

  HpNum FiLev;      ///< the target level \f$ f^{lev}_i \f$
  HpNum Beta;       ///< beta factor \f$ \beta_i \f$
  HpNum MaxBeta;    ///< maximum value for beta factor

  HpNum LpsFct;     ///< scaling factor

 }; // end( class Stepsize )

/* @} end( group( Stepsize_CLASSES ) ) */
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#if( OPT_USE_NAMESPACES )
 };  // end( namespace NDO_di_unipi_it )
#endif

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif  /* Stepsize.h included */

/*--------------------------------------------------------------------------*/
/*------------------------ End File Stepsize.h -----------------------------*/
/*--------------------------------------------------------------------------*/

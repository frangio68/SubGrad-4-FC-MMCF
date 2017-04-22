/*--------------------------------------------------------------------------*/
/*--------------------------- File Deflection.h ----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Definition of the abstract class <tt>Deflection</tt>, which defines the
 * interface for computing the delflection coefficient within the SubGrad
 * class [see SubGrad.h].
 * 
 * \version 1.00
 *
 * \date 27 - 10 - 2014
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

#ifndef __Deflection
 #define __Deflection /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "SubGrad.h"

/*--------------------------------------------------------------------------*/
/*----------------------------- NAMESPACE ----------------------------------*/
/*--------------------------------------------------------------------------*/

#if( OPT_USE_NAMESPACES )
namespace NDO_di_unipi_it
{
#endif

/*--------------------------------------------------------------------------*/
/*-------------------------- CLASS Deflection -------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/** @defgroup Deflection_CLASSES Classes in Deflection.h
    @{ */

/** The class Deflection provides an interface for the deflection rule (DR),
    to be used in the SubGrad solver [see SubGrad.h].

    The aim of this class is to compute the deflection coefficient
    \f$ \alpha_i \f$, which is employed in the search direction formula as
    follows:
    \f[
     d_i = \alpha_i g_i + (1-\alpha_i)d_{i-1}\;.
    \f]

    The user must extend the class to one or more deflection rules following
    this interface. */

class Deflection
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

/** Constructor of the class. It has no parameters apart the SubGrad object
    using it because there is no "standard" deflection rule (except "no
    deflection", i.e., \f$ \alpha_i = 1 \f$, which however in SubGrad is
    treated by just setting the null Deflection object). */

   Deflection( SubGrad *slvr )
   {
    if( ! slvr )
     throw NDOException( "Deflection::Deflection: no subgradient solver " );
    Solver = slvr;

    VOLLog = 0;
    VOLLLvl = 0;

    }  // end( Deflection::Deflection )

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

   virtual void SetVOLLog( ostream *outs = 0 , const char lvl = 0 )
   {
    if( ( VOLLog = outs ) )
     VOLLLvl = lvl;
    else
     VOLLLvl = 0;
    }

/*--------------------------------------------------------------------------*/
/** The method initializes the DR. It does nothing for the base class. */

   virtual void Format( void ) {}

/*@} -----------------------------------------------------------------------*/
/*-------------------- METHODS FOR DEFLECTION COMPUTATION ------------------*/
/*--------------------------------------------------------------------------*/
/** @name Reading the solution
    @{ */

/** This method must be called before GetDFLCoeff() [see below].
   This is indeed the core of every derived class, producing a new
   deflection coefficient.

   Typically, the previous coefficient will be unavailable after the call to
   NewDEF() [GetDFLCoeff()]. */

   virtual void NewDEF( void ) = 0;

/*@} -----------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/
/** @name Reading the solution
    @{ */

/** Returns the deflection coefficient. This function must be called after
   NewDEF() [see above]. */

   virtual cHpNum GetDFLCoeff( void ) = 0;

/*--------------------------------------------------------------------------*/
/** Returns true if a <em>Serious Step</em> (SS) takes place. Otherwise a
   <em>Null Step (NS)</em> occurs. By default, true is returned.

   Typically, a SS occurs when a \" good \" improvement of the function
   \f$ f \f$ is obtained, usually it should be as good as the expected 
   improvement [see Delta()]. */

   virtual const bool DoSS( void )
   {
    return( true );
    }

/*--------------------------------------------------------------------------*/
/** Returns the expected improvement in the objective function [see DoSS()].
   This value is not useful for all DR. By default Inf<HpNum>() is returned.
   */

   virtual cHpNum Delta( void )
   {
    return( Inf<HpNum>() );
    }

/*@}------------------------------------------------------------------------*/
/*------------------------------ DESTRUCTOR --------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Destructor
    @{ */

   virtual ~Deflection() {}

/*@} -----------------------------------------------------------------------*/
/*---------------------- PROTECTED PART OF THE CLASS -----------------------*/
/*--------------------------------------------------------------------------*/

 protected:

/*--------------------------------------------------------------------------*/
/*-------------------------- PROTECTED METHODS -----------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Protected methods to read from Subgrad

    These methods are used to provide access to the information stored into
    the SubGrad object. Deflection is a "friend" of SubGrad and therefore it
    can read its protected data structures, but classes derived from
    Deflection are not friend of SubGrad and cannot. This is why these methods
    are defined in the base class (and implemented in SubGrad.C).
    @{ */

/*--------------------------------------------------------------------------*/
/** Returns the pointer to FiOracle. Thus, Deflection can ask directly the
   FiOracle object for the function information. */

   FiOracle* GetOracle( void );

/*--------------------------------------------------------------------------*/
/** Returns the stepsize \f$\nu_i\f$. */

   HpNum GetStepsize( void );

/*--------------------------------------------------------------------------*/
/** Returns the norm of the subgradient \f$ g_i \f$. */

   cHpNum GetGiNorm( void );

/*--------------------------------------------------------------------------*/
/** Returns the norm of the direction \f$ d_i \f$. */

   cHpNum GetDNorm( void );

/*--------------------------------------------------------------------------*/
/** Returns the scalar product \f$ g_i^{\top} d_i \f$. */

   cHpNum GetdGk( void );

/*--------------------------------------------------------------------------*/
/** Returns the linearization error \f$ \sigma_i \f$. */

   cHpNum GetSigma( void );

/*--------------------------------------------------------------------------*/
/** Returns the linearization error \f$ \epsilon_i \f$. */

   cHpNum GetEpsilon( void );

/*--------------------------------------------------------------------------*/
/** Returns the full function \f$ f(\lambda_i) \f$. */

   cHpNum ReadFVal( void );

/*@} -----------------------------------------------------------------------*/
/*----------------------- PROTECTED DATA STRUCTURES  -----------------------*/
/*--------------------------------------------------------------------------*/

  SubGrad *Solver;    ///< (pointer to) the SubGrad solver
  ostream *VOLLog;    ///< the output stream object
  char VOLLLvl;       ///< the "level of verbosity"

/*--------------------------------------------------------------------------*/

 };  // end( class Deflection )

/* @} end( group( Deflection_CLASSES ) ) */
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#if( OPT_USE_NAMESPACES )
 };  // end( namespace NDO_di_unipi_it )
#endif

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif  /* Deflection.h included */

/*--------------------------------------------------------------------------*/
/*------------------------ End File Deflection.h ---------------------------*/
/*--------------------------------------------------------------------------*/

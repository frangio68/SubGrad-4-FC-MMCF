/*--------------------------------------------------------------------------*/
/*-------------------------- File ExDualCQKnP.h ----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * 
 * Continuous Quadratic Knapsack Problems (CQKnP) solver based on the
 * standard dual-ascent approach, which extends the DualCQKnP class to
 * non-negative quadratic costs and extended real bounds (and therefore
 * fully conforms to the standard interface for CQKnP solver defined by the
 * abstract base class CQKnpClass).
 *
 * \version 1.06
 *
 * \date 29 - 03 - 2012
 *
 * \author Antonio Frangioni \n
 *         Operations Research Group \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Enrico Gorgone \n
 *         Dipartimento di Elettronica Informatica e Sistemistica \n
 *         Universita' della Calabria \n
 *
 * Copyright &copy 2011 - 2012 by Antonio Frangioni, Enrico Gorgone.
 */

/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __ExDualCQKnP
 #define __ExDualCQKnP /* self-identification: #endif at the end of the file*/

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "DualCQKnP.h"

/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE and USING ----------------------------*/
/*--------------------------------------------------------------------------*/

namespace CQKnPClass_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*---------------------------- CLASSES -------------------------------------*/
/*--------------------------------------------------------------------------*/
/** Continuous Quadratic Knapsack Problems (CQKnP) solver derived from the
    DualCQKnP class (and therefore from CQKnPClass) and extending it, using
    the same standard dual-ascent approach, to non-negative quadratic costs
    and extended real bounds. */

class ExDualCQKnP : public DualCQKnP {

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/
/*--                                                                      --*/
/*--  The following methods and data are the actual interface of the      --*/
/*--  class: the standard user should use these methods and data only.    --*/
/*--                                                                      --*/
/*--------------------------------------------------------------------------*/

 public:

/*--------------------------------------------------------------------------*/
/*--------------------------- PUBLIC METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/
/*---------------------------- CONSTRUCTOR ---------------------------------*/
/*--------------------------------------------------------------------------*/

   ExDualCQKnP( bool sort ) : DualCQKnP( sort ) {};

/*--------------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/

   const double *KNPGetX( void );

/*--------------------------------------------------------------------------*/
/*--------------------- PROTECTED PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/
/*--                                                                      --*/
/*-- Nobody should ever look at this part: everything that is under this  --*/
/*-- advice may be changed without notice in any new release of the code. --*/
/*--                                                                      --*/
/*--------------------------------------------------------------------------*/

 protected:

/*--------------------------------------------------------------------------*/
/*-------------------------- PROTECTED METHODS -----------------------------*/
/*--------------------------------------------------------------------------*/

   bool CheckPFsb( void );

   bool CheckDFsb( void );

/*--------------------------------------------------------------------------*/

  #if DualCQKnP_SANITY_CHECKS
   inline void SanityCheckB( void );  // bounds

   inline void SanityCheckC( void );  // objective function
  #endif

/*--------------------------------------------------------------------------*/

   void SetName( void );

   void PreSort( void );

   void FindDualSol ( void );

 };  // end( class DualCQKnP )

/*--------------------------------------------------------------------------*/

 };  // end( namespace KNPClass_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif  /* ExDualCQKnP.h included */

/*--------------------------------------------------------------------------*/
/*------------------------ End File ExDualCQKnP.h --------------------------*/
/*--------------------------------------------------------------------------*/

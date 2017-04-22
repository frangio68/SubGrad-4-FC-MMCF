/*--------------------------------------------------------------------------*/
/*----------------------------- File Main.C --------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Main() for computing bounds for the Fixed-Charge Multicommodity
 * Capacitated Network Design (FC-MCND) problem under the NDOSolver/FiOracle
 * interface, using the SubGrad solver [see SubGrad.h] and either Flow or
 * Knapsack relaxation [see FlwFiOrcl.h, KnpsFiOrcl.h].
 * 
 * \version 1.00
 *
 * \date 22 - 10 - 2015
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
 * Copyright &copy 2011 - 2015 by Antonio Frangioni , Enrico Gorgone
 */

/*--------------------------------------------------------------------------*/
/*------------------------------- MACROS -----------------------------------*/
/*--------------------------------------------------------------------------*/

#define FileCsv 0
/**< Select the type of the final report:

   - 0 ==> txt file
   - 1 ==> csv file  */

#define ITBreak 0
/**< The user may wish to display the results more times during the execution.
   This job is done changing repeatedly the maximum number of iterations
   (MaxIter). The process can be stopped several times, allowing the reading
   of the results and without inducing any complication

   - 0 ==> MaxIter is defined one time only (by the user)
   - 1 ==> MaxIter in { 100 , 200 , 500 , 1000 , 2000 , 5000 , 10000 }
   - 2 ==> MaxIter in { 1000 , 2000 , 5000 , 10000 , 20000 , 50000 , 100000
                        200000 , 500000 , 1000000 } */

#define WHICH_STEPSOLVER 0
 /**< Select which stepsize rule to be used within the subgradient method (SM):

    - 0 ==> ColorTV
    - 1 ==> FumeroTV
    - 2 ==> Polyak */

#define WHICH_VOLSOLVER 1
 /**< Select which deflction rule to be used within SM:

    - 0 ==> STSubgrad, i.e. the original SM
    - 1 ==> Volume
    - 2 ==> PrimalDual, i.e. the primal-dual SM (PDSM) */

#define WHICH_FIORACLE 1
/**< Select the kind of Lagrangian relaxation:

     - 0 ==> FlwFiOrcl, i.e. flow relaxation
     - 1 ==> KnpsFiOrcl, i.e. knapsack relaxation */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

// NDOSolver header  - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#include "SubGrad.h"

// Deflection header(s), if used  - - -  - - - - - - - - - - - - - - - - - - -
//- - - - - - - - - - - - - - - - - - - - - - - -- - - - - - - - - - - - - - -

#if( WHICH_VOLSOLVER == 1 )
 #include "Volume.h"
#elif( WHICH_VOLSOLVER == 2 )
 #include "PrimalDual.h"
#elif WHICH_VOLSOLVER > 3
 #error "Unable to find the stepsize solver"
#endif

// Stepsize header(s), if used  - -  - - - - - - - - - - - - - - - - - - - - -
//- - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - - - - - - - - -

#if( WHICH_VOLSOLVER != 2 )
 #if( WHICH_STEPSOLVER == 0 )
  #include "ColorTV.h"
 #elif( WHICH_STEPSOLVER == 1 )
  #include "FumeroTV.h"
 #elif( WHICH_STEPSOLVER == 2 )
  #include "Polyak.h"
 #else
  #error "Unable to find the stepsize solver"
 #endif
#endif

// FiOracle header(s)  - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//- - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - - - - - - - - -

#if WHICH_FIORACLE == 0
 #include "FlwFiOrcl.h"
#elif WHICH_FIORACLE == 1
 #include "KnpsFiOrcl.h"
#else
 #error "Unable to find the FiOracle."
#endif

// other includes  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//- - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - - - - - - - - -

#include <fstream>
#include <sstream>

#include "OPTvect.h"

/*--------------------------------------------------------------------------*/
/*-------------------------------- USING -----------------------------------*/
/*--------------------------------------------------------------------------*/

#if( OPT_USE_NAMESPACES )
 using namespace NDO_di_unipi_it;
#endif

/*--------------------------------------------------------------------------*/
/*----------------------------- CONSTANTS ----------------------------------*/
/*--------------------------------------------------------------------------*/

// parameters files  - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#if WHICH_FIORACLE == 0
 const char *const ParF = "ParF.sg";
#else
 const char *const ParF = "ParK.sg";
#endif

// logs files (iteration by iteration) - - - - - - - - - - - - - - - - - - - -
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 const char *const logF = "log.sg";

// output file (final results) - - - - - - - - - - - - - - - - - - - - - - - -
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#if( FileCsv )
 const char *const out = "report.csv";
#else
 const char *const out = "report.txt";
#endif

/*--------------------------------------------------------------------------*/
/*----------------------------- FUNCTIONS ----------------------------------*/
/*--------------------------------------------------------------------------*/

template<class T>
static inline void str2val( const char* const str , T &sthg )
{
 istringstream( str ) >> sthg;
 }

/*--------------------------------------------------------------------------*/
/*-------------------------------- main() ----------------------------------*/
/*--------------------------------------------------------------------------*/

int main( int argc , char **argv )
{
 // read the command-line parameters- - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 char type = 's';
 HpNum optimum = - Inf<HpNum>();

 switch( argc ) {

  case( 4 ): str2val( argv[ 3 ] , optimum ); // pass the optimum of Fi()
  case( 3 ): type = *argv[ 2 ];              // set the type of the input file
  case( 2 ): break;

  default:   cerr << "Usage: " << argv[ 0 ] << " file_name [prb_typ]" << endl;
             return( 1 );
  }

 // enter the try-block - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 try {
  // read the instance - - - - - - - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  Graph *g = new Graph( argv[ 1 ] , type );

  // a call to PreProcess() assures that every arc has a finite capacity
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  g->PreProcess();

  // open the parameters file - - - - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  ifstream ParFile( ParF );
  if( ! ParFile.is_open() )
   cerr << "Warning: cannot open parameters file """ << ParF << """" << endl;

  // the lower bound is set to a percentage of the optimal value according to
  // the following rule - - - - - - - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  double percent;
  DfltdSfInpt( &ParFile , percent , double( 1 ) );

  if( optimum != - Inf<HpNum>() )
   optimum = optimum - percent * ABS(optimum);

  // select the verbosity of all the logs: Stepsize, Deflection , SubGrad and
  // FiOracle - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  int STPlvl, VOLlvl;
  int NDOlvl, Filvl;
  DfltdSfInpt( &ParFile , STPlvl , int( 0 ) );
  DfltdSfInpt( &ParFile , VOLlvl , int( 0 ) );
  DfltdSfInpt( &ParFile , NDOlvl , int( 0 ) );
  DfltdSfInpt( &ParFile , Filvl , int( 0 ) );

  // construct the FiOracle object- - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  #if WHICH_FIORACLE == 0
  FlwFiOrcl *Fi = new FlwFiOrcl( g , &ParFile );
  #else
   KnpsFiOrcl *Fi = new KnpsFiOrcl( g , &ParFile );
  #endif
  Fi->SetLowerBound( optimum );

  // free graph memory (data that will be used is kept in Fi) - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  delete g;

  // construct the NDOSolver object - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  SubGrad *s = new SubGrad( &ParFile );

  // construct the Deflection object  - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  #if( WHICH_VOLSOLVER == 0 )
   Deflection * VOL = 0;
  #elif( WHICH_VOLSOLVER == 1 )
   Volume *VOL = new Volume( s , &ParFile );
  #else
   PrimalDual *VOL = new PrimalDual( s , &ParFile );
  #endif

  // construct the Stepsize object r  - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 #if( WHICH_VOLSOLVER == 2 )
  PrimalDual *STP = VOL;
 #else
  #if( WHICH_STEPSOLVER == 0 )
   ColorTV *STP = new ColorTV( s , &ParFile );
  #elif( WHICH_STEPSOLVER == 1 )
   FumeroTV *STP = new FumeroTV( s , &ParFile );
  #else
   Polyak *STP = new Polyak( s , &ParFile );
  #endif
 #endif

  // close the parameters file  - - - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  ParFile.close();

  // open the log file  - - - - - - - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  ofstream LOGFile( logF , ofstream::out );
  if( ! LOGFile.is_open() )
   cerr << "Warning: cannot open log file """ << logF << """" << endl;

  // set the verbosity for the objects Deflection and Stepsize  - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  STP->SetSTPLog( &LOGFile , STPlvl );
  if( VOL )
   VOL->SetVOLLog( &LOGFile , VOLlvl );

  // pass Stepsize, Volume and CQKnP to the SubGrad    - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  s->SetStepsize( STP );
  s->SetDeflection( VOL );
  s->SetQKNP( 0 );

  // do timing - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  Fi->SetFiTime();
  s->SetNDOTime();

  // pass the FiOracle to the NDOSolver and vice-versa - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  s->SetFiOracle( Fi );
  Fi->SetNDOSolver( s );

  // set the verbosity of the objects FiOracle and SubGrad - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  Fi->SetFiLog( &LOGFile , Filvl );
  s->SetNDOLog( &LOGFile , NDOlvl );

  // open output file  - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  ofstream Out( out , ofstream::app );
  if( ! Out.is_open() )
   cerr << "Warning: cannot open log file """ << out << """" << endl;

  // solve the problem  - - - - - - - - - - - -  - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  s->KeepBestLambda();
  NDOSolver::NDOStatus result;

  for( Index i = 100 ; i <= 100000 ; ) {
   // set several breakpoints by changing MaxIter  - - - - - - - - - - - - -

   #if( ITBreak == 2 )
    s->SetPar( NDOSolver::kMaxItr , int( i * 10 ) );
   #elif( ITBreak == 1 )
    s->SetPar( NDOSolver::kMaxItr , int( i ) );
   #endif
   result = s->Solve();

   #if( FileCsv )
    Out << s->NDOTime() << "\t" << Fi->FiTime() << "\t" << s->FiEval() << "\t";
    #if WHICH_FIORACLE == 0
     Out << Fi->AddTime()<< "\t";
    #endif
    Out << s->NrIter() << "\t";

    Out.precision( 16 );
    switch( result ) {
     case( NDOSolver::kOK ) :
	  Out << -s->ReadBestFiVal() << "\t"<< -s->ReadHatFiVal()
	      << "\t Status: OK";
      break;
     case( NDOSolver::kStpIter ) :
	  Out << -s->ReadBestFiVal() << "\t" << -s->ReadHatFiVal()
	      << "\t Status: kStpIter";
      break;
     case( NDOSolver::kStpTime ) :
	  Out << -s->ReadBestFiVal() << "\t" << -s->ReadHatFiVal()
	      << "\t Status: kStpTime";
      break;
     case( NDOSolver::kStopped ) :
      Out << -s->ReadBestFiVal() << "\t" << -s->ReadHatFiVal()
	  << "\t Status: kStopped";
      break;
     default:
      Out << "\t" << "\t Status: Error";
     }

    Out << endl;

   #else

    Out << "NrIter: " << s->NrIter() << "\t" << "FiEval: " << s->FiEval()
	<< "\t" << "TimeNDO: " << s->NDOTime() << "\t" << "TimeFi: "
	<< Fi->FiTime() << "\t";

    Out.precision( 16 );

    switch( result ) {
     case( NDOSolver::kOK ) :
	 Out << "BestValue: " << -s->ReadBestFiVal()  << endl;
     break;
    case( NDOSolver::kStpIter ) :
	 Out << "BestValue: " << -s->ReadBestFiVal() << endl;
     break;
    case( NDOSolver::kStpTime ) :
	 Out << "BestValue: " << -s->ReadBestFiVal() << endl;
     break;
    case( NDOSolver::kStopped ) :
   	 Out << "BestValue: " << -s->ReadBestFiVal() << endl;
        break;
    default:
	 Out << "\t Status: Error" << endl;
    }

  #endif

  #if( ITBreak == 1 )
   switch( i ) {
    case( 100 ): i = 200; break;
    case( 200 ): i = 500; break;
    case( 500 ): i = 1000; break;
    case( 1000 ): i = 2000; break;
    case( 2000 ): i = 5000; break;
    case( 5000 ): i = 10000; break;
    default: i = 1e6; break;
    }
  #elif( ITBreak == 2 )
   switch( i ) {
    case( 100 ): i = 200; break;
    case( 200 ): i = 500; break;
    case( 500 ): i = 1000; break;
    case( 1000 ): i = 2000; break;
    case( 2000 ): i = 5000; break;
    case( 5000 ): i = 10000; break;
    case( 10000 ): i = 20000; break;
    case( 20000 ): i = 50000; break;
    case( 50000 ): i = 100000; break;
    default: i = 1e6; break;
    }
  #else
   break;
  #endif

  }

 // end resolution- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 // deallocate everything- - - - - - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 s->SetStepsize();
 s->SetDeflection();
 delete STP;

 #if( ( WHICH_VOLSOLVER != 2 ) && VOL )
  delete VOL;
 #endif

 delete s;
 delete Fi;

 // close the log and output file - - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 LOGFile.close();
 Out.close();

 return( 0 );

 } // end try part- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 // managing exceptions - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 catch( exception &e ) {
  cerr << e.what() << endl;
  return( 1 );
  }
 catch(...) {
  cerr << "Error: unknown exception thrown" << endl;
  return( 1 );
  } 

 // the end - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 return( 0 );

 }  // end( Main )

/*--------------------------------------------------------------------------*/
/*-------------------------- End File Main.C -------------------------------*/
/*--------------------------------------------------------------------------*/

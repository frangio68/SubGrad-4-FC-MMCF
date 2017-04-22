/*--------------------------------------------------------------------------*/
/*--------------------------- File Main.C ----------------------------------*/
/*--------------------------------------------------------------------------*/
/*--                                                                      --*/
/*--  Simple main() for testing the MMCFCplex class                       --*/
/*--                                                                      --*/
/*--                            VERSION 2.03                              --*/
/*--                           28 - 06 - 2012                             --*/
/*--                                                                      --*/
/*--                  Original Idea and Implementation by:                --*/
/*--                                                                      --*/
/*--                          Antonio Frangioni                           --*/
/*--                                                                      --*/
/*--                      Operations Research Group                       --*/
/*--                     Dipartimento di Informatica                      --*/
/*--                         Universita' di Pisa                          --*/
/*--                                                                      --*/
/*--                            Enrico Gorgone                            --*/
/*--                                                                      --*/
/*--             Dipartimento di Informatica e Sistemistica               --*/
/*--                      Universita' della Calabria                      --*/
/*--                                                                      --*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*------------------------------- MACROS -----------------------------------*/
/*--------------------------------------------------------------------------*/
   
#define PRINTRESULTS 0
/**< Select which results are printed:
    - 0 ==> none
    - 1 ==> flow variables
    - 2 ==> flow and design variables */

#define FileCsv 0
/**< use a csv file to print final report */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "MMCFCple.h"

#include <fstream>
#include <sstream>   // For istringstream

#include "OPTvect.h"

/*--------------------------------------------------------------------------*/
/*-------------------------------- USING -----------------------------------*/
/*--------------------------------------------------------------------------*/

#if( OPT_USE_NAMESPACES )
using namespace MMCFClass_di_unipi_it;
#endif

const char *const ParF = "ParValue.cpx";

#if( FileCsv )
 const char *const logF = "log.csv";
#else
 const char *const logF = "log.txt";
#endif

const char *const logCPX = "log.cpx";

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
 // initialize and read command line parameters - - - - - - - - - - - - - - -
 // - - - - - -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 HpNum lowB = - Inf<HpNum>();
    
 switch( argc ) {

  case( 3 ): str2val( argv[ 2 ] , lowB );
         
  case( 2 ): break;

  default:   cerr << "Usage: " << argv[ 0 ] << " file_name [typ]"
 		     << endl;
             return( 1 );
  }

 // set the Log of NDData  - - - - - - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 ofstream LOGCPX( logCPX , ofstream::out );
 if( ! LOGCPX.is_open() )
  cerr << "Warning: cannot open log file """ << logCPX << """" << endl;


 // open parameters file - - - - - - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 bool relax;          // true if the integrality is relaxed
 char add;            // add (strong) design constrs if present y*/s/n

 int BBlvl;           // level of verbosity of MIP problem
 char type;           // type of the input file
 double epsilon;      // relative tolerance
 int threads;         // the number of threads

 ifstream ParFile( ParF );
 if( ! ParFile.is_open() )
  cerr << "Warning: cannot open parameters file """ << ParF << """" << endl;

 DfltdSfInpt( &ParFile , type , char('s') );

 DfltdSfInpt( &ParFile , add , char('y') );
 DfltdSfInpt( &ParFile , relax , bool(0) );

 DfltdSfInpt( &ParFile , BBlvl , int(2) );

 DfltdSfInpt( &ParFile , epsilon , double(1e-6) );
 DfltdSfInpt( &ParFile , threads , int(1) );

 ofstream LOGMP( logF , ofstream::app );
 if( ! LOGMP.is_open() )
  cerr << "Warning: cannot open log file """ << logF << """" << endl;

 // read and modify the problem - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 Graph *Gh = new Graph( argv[ 1 ] , type );

 Graph::cIndex k = Gh->NrComm();
 Graph::cIndex m = Gh->NrArcs();
 Graph::cIndex n = Gh->NrNodes();

 if( add == 'n' )
  Gh->SetExtraVars( 0 );     // there are no extra variables

 if( Gh->NrExtraVars() ) {   // if there are extra variables, we assume these
  // are m arc design variables: add the m corresponding design constraints
  // 0 <= \sum_h - x_{ij}^h + u_{ij} y_{ij} <= F_INF, and in case also the
  // k * m strong forcing constraints
  // 0 <= x_{ij}^h + u_{ij}^h y_{ij} <= F_INF

  int nxc;
  double *Vals;
  int *IBeg, *Indx;

  switch( add ) {

  case( 'b' ): nxc = m * ( k + 1 );
         Indx = new int[ m * ( k + 1 ) + 2 * m * k ];
         Vals = new double[ m * ( k + 1 ) + 2 * m * k ];
        break;
  case( 'w' ): nxc = m;
        Indx = new int[ m * ( k + 1 ) ];
        Vals = new double[ m * ( k + 1 ) ];
        break;
  case( 's' ): nxc = m * k;
        Indx = new int[ 2 * m * k ];
        Vals = new double[ 2 * m * k ];
        break;
  case( 'n' ): break;
  default:
	  cerr << endl << "no choice" << endl;
   }

  if( add != 'n' ) {

   IBeg = new int[ nxc + 1 ];

   int ptr = 0; // the position of the cursor into both vectors Vals and Indx
   int cnrt;    // the position of the cursor into IBeg (just for strong ones)

   if( add == 'w' || add == 'b' ) { // weak forcing constraints

    for( Graph::Index j = 0 ; j < m ; j++  ) {
     IBeg[ j ] = ptr;
 	 for(  Graph::Index h = 0 ; h < k ; h++ ) {
 	  Vals[ ptr ] = -1;
 	  Indx[ ptr++ ] = h * m + j;
 	  }

 	 Vals[ ptr ] = Gh->TotalCapacityJ( j );  // mutual capacity of the j-th arc
 	 Indx[ ptr++ ] = k * m + j;              // the j-th extra variable
 	 }
    cnrt= m;
    IBeg[ m ] = ptr;
    }
   else
	cnrt = 0;

   if( add == 's' || add == 'b' ) {  // add strong forcing constraints

    for( Graph::Index j = 0 ; j < m ; j++ )
     for( Graph::Index h = 0 ; h < k ; h++ ) {
      IBeg[ cnrt++ ] = ptr;
      Vals[ ptr ] = -1;
      Indx[ ptr++ ] = h * m + j;
      Vals[ ptr ] = Gh->CapacityKJ( h , j );
      Indx[ ptr++ ] = k * m + j;
      }

    IBeg[ cnrt ] = ptr;
    }

   Gh->SetExtraConstr( nxc , IBeg , Indx , Vals );

   } // end( adding forcing constraints )


  if( relax )
   Gh->SetIntVar( k , false );  // all extra variables are continuous
  else
   Gh->SetIntVar( k , true );
  }

 Gh->PreProcess();    

 // allocate the solver - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 MMCFCplex *mmcf = new MMCFCplex( Gh , &ParFile );

 if( argc > 2 )
  mmcf->SetCplexParam( CPX_PARAM_WORKDIR , argv[ 2 ] );

 // set tolerance  - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 mmcf->SetCplexParam( CPX_PARAM_EPOPT , epsilon );
 mmcf->SetCplexParam( CPX_PARAM_EPGAP , epsilon);

 // pass the number of threads - - - - - - - - - - - - - - - - - - - - - - -

 mmcf->SetCplexParam( CPX_PARAM_THREADS , threads );
 mmcf->SetCplexParam( CPX_PARAM_OBJULIM , (-lowB) * ( 1 - epsilon ) );

 // pass Log- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 mmcf->SetMMCFLog( &LOGCPX , BBlvl );

 // free some memory- - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 delete( Gh );

 // set the timers on - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 mmcf->SetMMCFTime();

 // solve the problem - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 MMCFClass::MMCFStatus Status = mmcf->SolveMMCF();

 // get solution - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 #if( PRINTRESULTS > 0 )
  Index dim = mmcf->NrArcs() * mmcf->NrComm(); // dimension
  double* Flow = new double[ dim ]; // flow

  mmcf->SetMFlwSol( Flow );
  mmcf->GetPSol();  // get flows

  clog << "x = ";
  for( Index i = 0 ; i < dim ; i++ ) // print flow
   clog << Flow[ i ] << "\t";

  clog << endl;

  delete[] Flow;

  double* Costs = new double[ dim ]; // costs
  mmcf->GetCosts( Costs ); // get costs

  clog << "c = ";
  for( Index i = 0 ; i < dim ; i++ )  // print costs
   clog << Costs[ i ] << "\t";

  clog << endl;

  delete[] Costs;

  #if( PRINTRESULTS > 1 )
   if( mmcf->NrXtrVrs() ) {
    Row Xtr = new Number[ mmcf->NrXtrVrs() ];
    mmcf->SetXtrSol( Xtr );
    mmcf->GetPSol();

    clog << " y = ";
    for( Index i = 0 ; i < mmcf->NrXtrVrs() ; i++ )  // print extra vars
     clog << Xtr[ i ] << "\t";

    clog << endl;

    delete[] Xtr;

    double* XtrCs = new double[ mmcf->NrXtrVrs() ];
    mmcf->GetXtrCsts( XtrCs );

    clog << "f = ";
    for( Index i = 0 ; i < mmcf->NrXtrVrs() ; i++ )  // print extra costs
     clog << XtrCs[ i ] << "\t";

    clog << endl;

    delete[] XtrCs;
    }
  #endif
 #endif

 // get the results - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 double tu , ts;
 mmcf->TimeMMCF( tu , ts );       // get the running time

 MMCFClass::FONumber OV1 = mmcf->GetPVal();      // get the primal value
 MMCFClass::FONumber OV2 = mmcf->GetDVal();      // get the dual value

 Index NNodes = mmcf->GetNode();

 // clean up- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 delete mmcf ;

 // output the results- - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 #if( FileCsv )

 LOGMP << tu + ts << "\t" << NNodes << "\t";

 switch( Status ) {
  case( MMCFClass::kOK ) :
   LOGMP.precision( 16 );
   LOGMP <<  OV1 << "\t" <<  OV2 << "\t Status: OK" << endl;
   break;
  case( MMCFClass::kStopped ) :
   LOGMP <<  OV1 << "\t" <<  OV2 << "\t Status: Stopped" << endl;
   break;
  case( MMCFClass::kUnfeasible ) :
   LOGMP << "\t" << "\t Status: Unfeas." << endl;
   break;
  case( MMCFClass::kUnbounded ) :
   LOGMP << "\t" << "\t Status: Unbound." << endl;
   break;
  default :
   LOGMP << "\t" << "\t Status: Error" << endl;
  }

 #else

 LOGMP << "Time:" << tu + ts << " ~ Node limit: " << NNodes << "\t" ;

 switch( Status ) {
  case( MMCFClass::kOK ) :
   LOGMP.precision( 8 );
   LOGMP << "Status: OK, Value: ( " << OV1 << " , "<< OV2 << " ) " << endl;
   break;
  case( MMCFClass::kStopped ) :
   LOGMP << "Status: Stopped: ( " << OV1 << " , "<< OV2 << " ) " << endl;
   break;
  case( MMCFClass::kUnfeasible ) :
   LOGMP << "Status: Unfeas." << endl;
   break;
  case( MMCFClass::kUnbounded ) :
   LOGMP << "Status: Unbound." << endl;
   break;
  default :
   LOGMP << "Status: Error" << endl;
  }

 #endif

 LOGMP.close();
 ParFile.close();

 // the end - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 return( 0 );

 }  // end( main )

/*--------------------------------------------------------------------------*/
/*---------------------------- End File Main.C -----------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*--------------------------- File MMCFCple.C ------------------------------*/
/*--------------------------------------------------------------------------*/
/*--                                                                      --*/
/*--   Multicommodity Min Cost Flow (MMCF) Problems solver, based on      --*/
/*--   calls to the Cplex(TM) 6.6 Callable Libraries for solution of      --*/
/*--   generic LP problems.                                               --*/
/*--                                                                      --*/
/*--                            VERSION 3.04                              --*/
/*--                           21 - 01 - 2014                             --*/
/*--                                                                      --*/
/*--                  Original Idea and Implementation by:                --*/
/*--                                                                      --*/
/*--                           Paola Cappanera                            --*/
/*--                          Antonio Frangioni                           --*/
/*--                                                                      --*/
/*--                      Operations Research Group                       --*/
/*--                     Dipartimento di Informatica                      --*/
/*--                         Universita' di Pisa                          --*/
/*--                                                                      --*/
/*-- Copyright(C) 1996 - 2012 Antonio Frangioni                           --*/
/*--                                                                      --*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- IMPLEMENTATION -------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/


/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/
 
#include "MMCFCple.h"

#include "OPTvect.h"
#include "OPTUtils.h"

#include <assert.h>

/*--------------------------------------------------------------------------*/
/*------------------------------- MACROS -----------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*-------------------------------- USING -----------------------------------*/
/*--------------------------------------------------------------------------*/

#if( OPT_USE_NAMESPACES )
using namespace MMCFClass_di_unipi_it;

using namespace std;
#endif

/*--------------------------------------------------------------------------*/
/*------------------------------ CONSTANTS ---------------------------------*/
/*--------------------------------------------------------------------------*/

static const MMCFClass::CNumber C_INF = Inf<MMCFClass::CNumber>();
static const MMCFClass::FNumber F_INF = Inf<MMCFClass::FNumber>();

static const double EpsXNum = 1e-8;

/*--------------------------------------------------------------------------*/
/*-------------------- IMPLEMENTATION OF MMCFCplex -------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- PUBLIC METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/
/*---------------------------- CONSTRUCTOR ---------------------------------*/
/*--------------------------------------------------------------------------*/

MMCFCplex::MMCFCplex( Graph *Gh , istream *iStrm , CPXENVptr extenv )
{
 // select the algorithm used by cplex- - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 char alg;            // LP/MIP algorithm used
 DfltdSfInpt( iStrm , alg , char('b') );

 // initialization- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 Index YsPopulate;   // define the options for the Lagrangian relaxation
 Index YsLazy;
 double timeLimit;    // time limit

 DfltdSfInpt( iStrm , YsPopulate , Index( 0 ) );
 DfltdSfInpt( iStrm , YsLazy , Index( 1 ) ); // don't change it!!!!!!
 DfltdSfInpt( iStrm , timeLimit , double(1000) );

 Populate = bool( YsPopulate );
 Lazy = bool( YsLazy );

 // pre-processing - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 DfltdSfInpt( iStrm , preind , int(1) );
 DfltdSfInpt( iStrm , represolve , int(-1) );
 DfltdSfInpt( iStrm , prenode , int(0) );
 DfltdSfInpt( iStrm , reduce , int(3) );

 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // setup environment, if necessary - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( extenv )
  env = extenv;
 else {
  if( ! EnvICntr++ ) {
   int ts;
   genv = CPXopenCPLEX( &ts );

   assert( genv );
   assert( ! ts );
   }

  env = genv;
  }

 // select the algorithm used by cplex- - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( Gh->NIntVar() )                   // any variable is integer
  if( alg == 'a' || alg == 'b'  ||     // use a B&B to solve it
      alg == 'c' || alg == 'd'   ||
      alg == 'e' || alg == 'f' )
      alg = 'B';

 MMCFCplex::MMCFCplexA algorithm;
 switch ( alg ) {
  case ( 'a' ): algorithm = MMCFCplex::kPPrim;   break;
  case ( 'b' ): algorithm = MMCFCplex::kPDual;   break;
  case ( 'c' ): algorithm = MMCFCplex::kBarrier; break;
  case ( 'd' ): algorithm = MMCFCplex::kBrCrss;  break;
  case ( 'e' ): algorithm = MMCFCplex::kNPrim;   break;
  case ( 'f' ): algorithm = MMCFCplex::kNDual;   break;
  case ( 'A' ):
   algorithm = 	MMCFCplex::kMIP;
   CPXsetintparam( env , CPX_PARAM_STARTALG , 1 );
   break;
  case ( 'B' ):
   algorithm = 	MMCFCplex::kMIP;
   CPXsetintparam( env , CPX_PARAM_STARTALG , 2 );
   break;
  case ( 'C' ):
   algorithm = 	MMCFCplex::kMIP;
   CPXsetintparam( env , CPX_PARAM_STARTALG , 3 );
   break;
  case ( 'D' ):
   algorithm = 	MMCFCplex::kMIP;
   CPXsetintparam( env , CPX_PARAM_STARTALG , 4 );
   break;
  case ( 'E' ):
   algorithm = 	MMCFCplex::kMIP;
   CPXsetintparam( env , CPX_PARAM_STARTALG , 5 );
   break;
  case ( 'F' ):
   algorithm = 	MMCFCplex::kMIP;
   CPXsetintparam( env , CPX_PARAM_STARTALG , 6 );
   break;
  case ( 'G' ):
   algorithm = 	MMCFCplex::kMIP; // let CPLEX choose
   CPXsetintparam( env , CPX_PARAM_STARTALG , 0 );
   break;
  default:
   throw( MMCFException(
	 "MMCFCplex::SolveMMCF: does not select any algorithm" ) );
   }

  SetCplexA( algorithm );

 // set pre-processing  - - - - - - - - - - - - - - - - - - - - - - - - - - -

 SetCplexParam( CPX_PARAM_PREIND , preind );
 SetCplexParam( CPX_PARAM_REPEATPRESOLVE , represolve );
 SetCplexParam( CPX_PARAM_PRESLVND , prenode );
 SetCplexParam( CPX_PARAM_REDUCE , reduce );

 // set Max Time- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 CPXsetdblparam( env , CPX_PARAM_TILIM , timeLimit );

 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // get size of the problem - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 NNodes = Gh->NrNodes();
 NArcs = Gh->NrArcs();
 NComm = Gh->NrComm();

 NActvA = Gh->NActives();

 if( Gh->Actives() ) {
  ActvArcs = new Index[ NActvA + 1 ];
  VectAssign( ActvArcs , Gh->Actives() , NActvA + 1 );
  }
 else
  ActvArcs = 0;

 Drctd = Gh->Directed();

 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // check if there is a commodity for which at least one node doesn't exist -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 StrtRows = new Index[ NComm + 1 ];
 NDict    = new Index_Set[ NComm ];

 StrtRows[ 0 ] = 0;

 for( Index k = 0 , totNds = 0 ; k < NComm ; ) {
  cFRow Dk = Gh->DeficitsK( k );
  Index_Set NDk = NDict[ k ] = new Index[ NNodes ];
  for( Index j = NNodes ; j-- ; )
   *(NDk++) = ( *(Dk++) < F_INF ? totNds++ : Inf<Index>() );

  if( totNds - StrtRows[ k ] == NNodes ) {
   delete[] NDict[ k ];
   NDict[ k ] = 0;
   }

  StrtRows[ ++k ] = totNds;
  }

 if( StrtRows[ NComm ] == NComm * NNodes ) {
  delete[] NDict;
  NDict = 0;
  delete[] StrtRows;
  StrtRows = 0;
  }

 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // check if there is a commodity for which at least one arc doesn't exist- -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 StrtCols = new Index[ NComm + 1 ];
 ADict    = new Index_Set[ NComm ];

 StrtCols[ 0 ] = 0;

 for( Index k = 0 , totArcs = 0 ; k < NComm ; ) {
  cCRow Ck = Gh->CostsK( k );
  Index_Set ADk = ADict[ k ] = new Index[ NArcs ];
  for( Index j = NArcs ; j-- ; )
   *(ADk++) = ( *(Ck++) < C_INF ? totArcs++ : Inf<Index>() );

  if( totArcs - StrtCols[ k ] == NArcs ) {
   delete[] ADict[ k ];
   ADict[ k ] = 0;
   }

  StrtCols[ ++k ] = totArcs;
  }

 if( StrtCols[ NComm ] == NComm * NArcs ) {
  delete[] ADict;
  ADict = 0;
  }

 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // calculate the size of the LP- - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 int numcols = int( StrtCols[ NComm ] );
 if( ! Drctd )
  numcols *= 2;

 // in undirected graphs the flow on each arc (i,j) is splitted on two
 // different arcs: (i,j) and (j,i). However, the actual flow on the arc
 // (and, therefore the flow that is subject to the mutual capacity 
 // constraint) is meant to be the sum of the two flows.

 int matsz = 2 * numcols;  // for each arc there are two nonzeros: 1 and -1

 // for each comm., count the number of existent arcs with mutual capacity- -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 int NActive = 0;

 if( ActvArcs )  // not all arcs have mutual capacity
  for( Index k = 0 ; k < NComm ; ) {
   cIndex_Set tAA = ActvArcs;
   cCRow Ck = Gh->CostsK( k++ );
   for( Index h ; ( h = *(tAA++) ) < Inf<Index>() ; )
    if( Ck[ h ] < C_INF )
     NActive++;
   }
 else            // all arcs have mutual capacity
  NActive = int( StrtCols[ NComm ] );

 if( ! Drctd )
  NActive *= 2;

 matsz += NActive;  // add the nonzeros relatives to mutual capacity

 Index numrowsE = ( StrtRows ? StrtRows[ NComm ] : NNodes * NComm );

 // numrowsE is the number of rows of the node-arc incidence matrix 

 int numrows = int( numrowsE + NActvA );  // total rows

 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // allocating memory for Cplex structures- - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 double* obj = new double[ numcols ];
 
 double* rhs = new double[ numrows ];

 char* sense = new char[ numrows ];

 int* matbeg = new int[ numcols ];
 
 int* matcnt = new int[ numcols ];
 
 int* matind = new int[ matsz ];
 
 double* matval = new double[ matsz ];
 
 double* ub = new double[ numcols ];

 double* lb = new double[ numcols ];

 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // inputting data for Cplex- - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 const int p = ( Gh->NamesStartFrom1() ? -1 : 0 );  // scaling the node names

 Index w = 0;
 Index h = 0;
 int cntNds = 0;
 for( Index k = 0 ; k < NComm ; k++ ) {
  // initialize the rhs of the flow balance constraints - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  cFRow Dk = Gh->DeficitsK( k );
  for( Index j = 0 ; j < NNodes ; j++ ) {
   cFNumber Dkj = *(Dk++);
   if( Dkj < F_INF ) {
    rhs[ cntNds ] = Dkj;
    sense[ cntNds++ ] = 'E';
    }
   }

  // now consider arcs of comm k- - - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  Index pActive = 0;
  cCRow Ck = Gh->CostsK( k );
  Index_Set NDk = NDict ? NDict[ k ] : 0;
  for( Index j = 0 ; j < NArcs ; j++ ) {
   cCNumber Ckj = *(Ck++);
   if( Ckj < C_INF ) {  // arc j exists for commodity k - - - - - - - - - - -
    obj[ h ] = Ckj;
    matbeg[ h ] = w;

    if( NDk )
     matind[ w ] = NDk[ Gh->StartNJ( j ) + p ];
    else
     matind[ w ] = k * NNodes + Gh->StartNJ( j ) + p;

    matval[ w++ ] = -1;

    if( NDk )
     matind[ w ] = NDk[ Gh->EndNJ( j ) + p ];
    else
     matind[ w ] = k * NNodes + Gh->EndNJ( j ) + p;

    matval[ w++ ] = 1;

    if( ActvArcs )   // not all arcs have a mutual capacity
     if( ActvArcs[ pActive ] == j ) {  // arc j has a mutual capacity
      matind[ w ] = numrowsE + (pActive++);
      matval[ w++ ] = 1;
      matcnt[ h ] = 3;   
      }
     else
      matcnt[ h ] = 2;
    else {           // all arcs have a mutual capacity
     matind[ w ] = numrowsE + (pActive++);
     matval[ w++ ] = 1;
     matcnt[ h ] = 3;   
     }

    lb[ h ] = 0;
    cFNumber cap = Gh->CapacityKJ( k , j );
    ub[ h++ ] = ( cap < F_INF ? cap : CPX_INFBOUND );  // writing ub 
    }
   else     // arc j doesn't exist for comm. k- - - - - - - - - - - - - - - -
    if( ActvArcs && ( ActvArcs[ pActive ] == j ) )
     pActive++;      // advance the pointer

   }  // end( for( j ) )
  }  // end( for( each commodity k ) )

 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // doubling data for undirected graphs - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // for undirected problems: if the arc (i,j) is put in the column h, then
 // the arc (j,i) (that doubles (i,j)) is put in the column h + numcols / 2.
 // note that the the flow variable of (j,i) has the same objective function,
 // capacity and upper bound

 if( ! Drctd ) {
  Index shift = numcols / 2;
  for( Index h = 0 ; h < numcols / 2 ; ) {
   matbeg[ shift ] = matbeg[ h ];
   matcnt[ shift ] = matcnt[ h ];
   ub[ shift ] = ub[ h ];
   lb[ shift ] = lb[ h ];
   obj[ shift++ ] = obj[ h++ ];
   }

  shift = matsz / 2;
  for( Index h = numcols / 2 ; h < numcols ; ) {
   Index stcol = matbeg[ h ];
   matval[ stcol++ ] = 1;
   matval[ stcol ] = -1;

   if( matcnt[ h++ ] == 3 )
    matval[ ++stcol ] = 1;
   }
  }  // end( if( ! Drctd ) )

 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // rhs of the mutual capacity constraints- - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( ActvArcs ) {  // not all arcs have a mutual capacity
  Index pActive = 0;
  for( Index j = 0 , h = numrowsE ; j < NArcs ; j++ )
   if( ActvArcs[ pActive ] == j ) {
    pActive++;
    rhs[ h ] = Gh->TotalCapacityJ( j );
    sense[ h++ ] = 'L';

    }
  }
 else            // all arcs have a mutual capacity
  for( Index j = 0 , h = numrowsE ; j < NArcs ; ) {
   rhs[ h ] = Gh->TotalCapacityJ( j++ );
   sense[ h++ ] = 'L';
   }

 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // creating the LP - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 int status;
 lp = CPXcreateprob( env , &status , "MMCF" );
 assert( ! status );

 status = CPXcopylp( env , lp , numcols , numrows , 1 , obj , rhs , sense ,
		     matbeg , matcnt , matind , matval , lb , ub , 0 );
 assert( ! status );

 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // deleting temporary data structures- - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 delete[] lb;
 delete[] ub;
 delete[] matval;
 delete[] matind;
 delete[] matcnt;
 delete[] matbeg;
 delete[] sense;
 delete[] rhs;
 delete[] obj;

 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // setting the integer variables, if any - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 for( Index k = 0 ; k < NComm ; k++ )
  if( Gh->NIntVar( k ) )
   ChgIntVar( k , true , Gh->WIntVar( k ) );

 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // adding "extra" variables and constraints, if any- - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( ( XtrVrs = Gh->NrExtraVars() ) ) { // "extra" variables  - - - - - - - -

  CRow XtCsts_ = new CNumber[ XtrVrs ];
  VectAssign( XtCsts_ , Gh->CostsK( NComm ) , XtrVrs );

  AddExtraVars( XtrVrs , XtCsts_ ,
		Gh->CapacitiesK( NComm + 1 ) , Gh->CapacitiesK( NComm ) ,
		Gh->NIntVar( NComm ) ? true : false ,
		Gh->WIntVar( NComm ) );

  delete[] XtCsts_;

  }

 if( ( XtrCnst = Gh->NrExtraConst() ) ) {  // "extra" constraints - - - - - -
  cIndex nnz = Gh->NrExtraNonZ();

  int* IBeg = new int[ XtrCnst + 1 ];
  int* Indx = new int[ nnz ];
  double* Vals = new double[ nnz ];

  Gh->ExtraConstr( IBeg , Indx , Vals );

  AddExtraConstr( XtrCnst , IBeg , Indx , Vals ,
		  Gh->DeficitsK( NComm ) , Gh->DeficitsK( NComm + 1 ) );
  delete[] Vals;
  delete[] Indx;
  delete[] IBeg;
  }

 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // final initializations - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 //CA = kNDual;            // default solution algorithm
 PFeas = DFeas = false;  // no fasible solution known yet

 }  // end( MMCFCplex )

/*--------------------------------------------------------------------------*/
/*---------------------- OTHER INITIALIZATIONS -----------------------------*/
/*--------------------------------------------------------------------------*/

void MMCFCplex::SetMMCFLog( ostream *outs , const char lvl ) {

 MMCFClass::SetMMCFLog( outs , lvl );

 #if LOG_MMCF
  CPXsetintparam( env , CPX_PARAM_MIPDISPLAY , MMCFLLvl );

  logfile = CPXfopen ( "cplex.log" , "w");
  int status = CPXsetlogfile( env , logfile );
  assert( ! status );
 #endif
 }

/*--------------------------------------------------------------------------*/
/*-------------------- METHODS FOR SOLVING THE PROBLEM ---------------------*/
/*--------------------------------------------------------------------------*/

MMCFClass::MMCFStatus MMCFCplex::SolveMMCF( void )
{
 // start the timer - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( MMCFt )
  MMCFt->Start();

 PFeas = DFeas = false;

 // solving the problem with the chosen algorithm - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 int status = CPXgetprobtype( env , lp );     // get the current problem type

 #if LOG_MMCF
  int status1 = CPXwriteprob( env , lp , "MMCF.lp" , 0 );
  assert( ! status1 );
 #endif

 if( CA == kMIP ) {                           // using a MIP algorithm
  if( status != CPXPROB_MILP )                 // if it is *not* a MIP
   CPXchgprobtype( env , lp , CPXPROB_MILP );

  // ensure that the problem is a MIP; if the problem was not a MIP, this
  // just adds ctype[] information to the extent that all variables are
  // continuous; if the problem were a MIP but has been turned into a
  // CPXPROB_RELAXED problem, this reverts it to its original MIP form
  }
 else                                         // using a LP algorithm
  if( status == CPXPROB_MILP )                 // if it *is* a MIP
   CPXchgprobtype( env , lp , CPXPROB_LP );   // solve the relaxation

 switch( CA ) {
  case( kPPrim ):   status = CPXprimopt( env , lp );
                    break;
  case( kPDual ):   status = CPXdualopt( env , lp );
                    break;  
  case( kBarrier ): status = CPXbaropt( env , lp );
                    break;
  case( kBrCrss ):  status = CPXhybbaropt( env, lp , 'd' );
                    break;
  case( kMIP ):     status = CPXmipopt( env , lp );
                    IDSol = 0; // reset counter for solution pool
                    break;
  case( kNPrim ):   status = CPXhybnetopt( env , lp , 'p' );
                    break;
  default:          status = CPXhybnetopt( env , lp , 'd' );
  }

 // stop the timer- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( MMCFt )
  MMCFt->Stop();

 status = CPXgetstat( env , lp );

 if( CA != kMIP ) { //- - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( status == CPX_STAT_OPTIMAL   ) { // check optimality- - - - - - - - - -
   PFeas = true;   // For the moment does not consider the cases:
   DFeas = true;   // CPX_STAT_OPTIMAL_FACE_UNBOUNDED (Barrier and
   return( kOK );  // Barrier&Crossover)
   }               // CPX_STAT_OPTIMAL_INFEAS

  if( ( status == CPX_STAT_ABORT_IT_LIM ) ||
	  ( status == CPX_STAT_ABORT_DETTIME_LIM ) ||
	  ( status == CPX_STAT_ABORT_TIME_LIM ) ||
	  ( status == CPX_STAT_ABORT_USER ) ||
      ( status == CPX_STAT_ABORT_OBJ_LIM ) ||
      ( status == CPX_STAT_NUM_BEST ) ) { // find feasible solution - - - - -

   // for only both Barrier and Barrier&Crossover you may consider the
   // following two cases:
   // CPX_STAT_ABORT_PRIM_OBJ_LIM
   // CPX_STAT_ABORT_DUAL_OBJ_LIM

   int foo , pf , df;
   CPXsolninfo( env , lp , &foo , &foo , &pf , &df );

   if( pf )
    PFeas = true;

   if( df )
    DFeas = true;

   return( kStopped );
   }

  } // end if ( CA != kMIP )

 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 switch( CA ) {

  case( kMIP ):  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   if( ( status == CPXMIP_OPTIMAL ) || ( status == CPXMIP_OPTIMAL_TOL ) ) {

    PFeas = true;
    if( Populate ) {

     // CPXmipopt (Iphase): it solves the model to optimality but it retains nodes
     // that might be useful
     // CPXpopulate (IIphase): it generates multiple solutions  by using the
     // information computed and stored in the first phase
     //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

     status =  CPXpopulate( env , lp );

     // please, not consider the following exceptions:
     // CPXMIP_POPULATESOL_LIM: max # solutions generated by populate
     // CPXMIP_OPTIMAL_POPULATED: Populate has completed the enumeration
     //                           of all solutions it could enumerate
     // CPXMIP_OPTIMAL_POPULATED_TOL: [see above]

     // but just verify that a solution pool has been found

     }

    if( !( NumSol = CPXgetsolnpoolnumsolns( env , lp ) ) )
     throw( MMCFException( "MMCFCplex::SolveMMCF: Solution Pool is Empty" ) );

   #if LOG_MMCF
     int status3 = CPXsolwritesolnpoolall( env , lp , "MMCFpool.sol" );
     assert( ! status3 );
   #endif

    return( kOK );
    }

   if( status == CPXMIP_INFEASIBLE )
    return( kUnfeasible );
   else
    if( status == CPX_STAT_UNBOUNDED )
     return( kUnbounded );

   if( ( status == CPXMIP_NODE_LIM_FEAS ) ||
	   ( status == CPXMIP_ABORT_FEAS )    ||
	   ( status == CPXMIP_FAIL_FEAS )     ||
	   ( status == CPXMIP_FAIL_FEAS_NO_TREE ) ||
       ( status == CPXMIP_TIME_LIM_FEAS ) ||
       ( status == CPXMIP_MEM_LIM_FEAS )  ||
       ( status == CPXMIP_SOL_LIM ) ) {
     PFeas = true;
     return( kStopped );
     }

   if( ( status == CPXMIP_NODE_LIM_INFEAS ) ||
	   ( status == CPXMIP_ABORT_INFEAS )    ||
	   ( status == CPXMIP_FAIL_INFEAS ) ||
	   ( status == CPXMIP_FAIL_INFEAS_NO_TREE ) ||
       ( status == CPXMIP_TIME_LIM_INFEAS ) ||
       ( status == CPXMIP_MEM_LIM_INFEAS ) )
     return( kStopped );

     if( status == CPXMIP_DETTIME_LIM_FEAS  ) {
      PFeas = true;
      return( kStopped );
      }

     if( status == CPXMIP_DETTIME_LIM_INFEAS )
      return( kStopped );


   break;

  case( kPPrim ):  // Primal algorithms- - - - - - - - - - - - - - - - - - - - -
  case( kNPrim ):  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  case( kBarrier ):
  case( kPDual ):  // Dual algorithms- - - - - - - - - - - - - - - - - - - - - -
  case( kNDual ):  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  case( kBrCrss ):

   if( status == CPX_STAT_UNBOUNDED )
    return( kUnbounded );
   else
    if( status == CPX_STAT_INFEASIBLE )
     return( kUnfeasible );

  }  // end ( switch )

 // returning the right status- - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 return( kError );  // some conflicts have been occurred

 }  // end( MMCFCplex::SolveMMCF )

/*--------------------------------------------------------------------------*/
/*------------------- METHODS FOR READING RESULTS  -------------------------*/
/*--------------------------------------------------------------------------*/

MMCFClass::FONumber MMCFCplex::GetPVal( void )
{
 double objval_p = Inf<FONumber>();

 if( PFeas )
  if( CA == kMIP )
   CPXgetsolnpoolobjval( env , lp , IDSol , &objval_p );
  else
   CPXgetobjval( env , lp , &objval_p );

 return( objval_p );

 } // end( MMCFCplex::GetPVal )

/*--------------------------------------------------------------------------*/

MMCFClass::FONumber MMCFCplex::GetDVal( void )
{
 double objval_p =  - Inf<FONumber>();

 if( CA == kMIP )
  CPXgetbestobjval( env , lp , &objval_p );
 else
  if( DFeas )
   CPXgetobjval( env , lp , &objval_p );

 return( objval_p );
 } // end( MMCFCplex::GetDVal )

/*--------------------------------------------------------------------------*/

MMCFClass::FONumber MMCFCplex::GetUpprBnd( bool &HvSol )
{
 HvSol = false;              // by default, no integer solution is available
 FONumber UBnd = Inf<FONumber>(); // by default, no UB is known

 if( PFeas ) {

  double thrshld = EpsXNum;  // threshold for evaluating the extra solution

  Row XSolution = new Number[ NArcs ];
  SetXtrSol( XSolution );
  GetExtraSol();

  HvSol = true;
  for( Index j = 0 ; j < NArcs ; j++ )
   if( ( XSolution[ j ] >= thrshld ) && ( ( XSolution[ j ] - 1.0 ) <= -thrshld ) ) {
	HvSol = false;
  	break;
    }

  UBnd = GetPVal();

  XSol = 0;
  delete[] XSolution;

  }

 return( UBnd );  // by default, no UB is known
 }

/*--------------------------------------------------------------------------*/

bool MMCFCplex::GetPSol( void )
{
 if( MFSol )  // the flow solution is required- - - - - - - - - - - - - - - -
 {            //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( ( 0 <= WhchFS ) && ( WhchFS < NComm ) ) {  // one given commodity - - -
   const int start = StrtCols[ WhchFS ];
   const int end = StrtCols[ WhchFS + 1 ] - 1;

   if( CA == kMIP ) // CPX_INCUMBENT_ID == -1 but a copy is always added
	                // with index 0
    CPXgetsolnpoolx( env , lp , IDSol , MFSol , start , end );
   else
    CPXgetx( env , lp , MFSol , start , end );

   if( ! Drctd )  // the problem is undirected: the returned solution for
   {              // x[ ij ][ h ] is the *difference* between the variable
                  // for the "directed" arc and that for the "inverse" arc

    double* twinx = new double[ end - start + 1 ];

    if( CA == kMIP )
     CPXgetsolnpoolx( env , lp , IDSol , twinx  , start + StrtCols[ NComm ] ,
    	  end + StrtCols[ NComm ] );
    else
     CPXgetx( env , lp , twinx , start + StrtCols[ NComm ] ,
	      end + StrtCols[ NComm ] );

    VectSubtract( MFSol , twinx , end - start );

    delete[] twinx;
    }

   cIndex_Set ADk = ADict ? ADict[ WhchFS ] : 0;
   if( ADk )
    TranslateK( MFSol , ADk , NArcs , start );

   }    // end( if( flow of a given commodity ) )
  else {  // either the aggregated or the whole flow solution is required - -
   // fetch the whole flow solution - - - - - - - - - - - - - - - - - - - - -

   int totArcs = StrtCols[ NComm ];
   double *x = new double[ Drctd ? totArcs : 2 * totArcs ];

   if( CA == kMIP )
	CPXgetsolnpoolx( env , lp , IDSol , x , 0 , totArcs - 1 );
   else
    CPXgetx( env , lp ,  x , 0 , totArcs - 1 );

   if( ! Drctd ) {  // the problem is undirected: see above
    if( CA == kMIP )
     CPXgetsolnpoolx( env , lp , IDSol , x + totArcs , totArcs ,
    		 2 * totArcs - 1 );
    else
     CPXgetx( env , lp , x + totArcs , totArcs , 2 * totArcs - 1 );

    VectSubtract( x , x + totArcs , totArcs );
    }

   if( WhchFS == NComm ) {  // the aggregate flow is required - - - - - - - -
    VectAssign( MFSol , MFNumber( 0 ) , NArcs );

    double *tx = x;
    MFRow tMFS = MFSol;
    for( Index k = 0 ; k < NComm ; k++ ) {
     cIndex_Set ADk = ADict ? ADict[ k ] : 0;
     if( ADk ) {  // at least one arc doesn't exist for comm. k
      for( Index j = NArcs ; j-- ; tMFS++ )
       if( *(ADk++) < Inf<Index>() )
	*tMFS += *(tx++);
      }
     else         // all arcs exist for comm. k
      for( Index j = NArcs ; j-- ; )
       *(tMFS++) += *(tx++);

     } // end( for( k ) )
    }      // end( if( the aggregate flow is required ) ) - - - - - - - - - -
   else {  // the whole flow solution is required - - - - - - - - - - - - - -
    VectAssign( MFSol , x , totArcs );

    if( ADict )
     TranslateF( MFSol , ADict , StrtCols , NArcs );
    }

    delete[] x;

    }  // end( else( aggregate or full flow solution ) )- - - - - - - - - - -

   if( FBse )  // the solution is required in "sparse" format
    *Sparsify( MFSol , FBse , ( WhchFS <= NComm ? NArcs : NArcs * NComm ) )
     = Inf<Index>();

   }  // end( if( MFSol ) ) - - - - - - - - - - - - - - - - - - - - - - - - -
      //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( XSol )         // "extra" variables are required - - - - - - - - - - - -
  GetExtraSol();

 // reset internal state- - - - - - - - - - - - - - - - - - - - - - - - - - -

 XSol = 0;
 MFSol = 0;
 FBse = 0;

 if( CA == kMIP )            // remark that the counter for the current solution
  if( IDSol++ < NumSol - 1 ) // is increased too
   return( true );

 return( false );

 }  // end( MMCFCplex::GetPSol )

/*-------------------------------------------------------------------------*/

bool MMCFCplex::GetDSol( void )
{
 if( NPot )  // node potentials are required- - - - - - - - - - - - - - - - -
 {           // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( WhchNP < NComm ) {
   int start = StrtRows ? StrtRows[ WhchNP ] : WhchNP * NNodes ;
   int end = StrtRows ? StrtRows[ WhchNP + 1 ] : ( WhchNP + 1 ) * NNodes;

   CPXgetpi( env , lp , NPot , start , end - 1 );

   cIndex_Set NDk = NDict ? NDict[ WhchNP ] : 0;
   if( NDk )
    TranslateK( NPot , NDk , NNodes , start );
   }
  else {
   int end = StrtRows ? StrtRows[ NComm ] : NComm * NNodes;

   CPXgetpi( env , lp , NPot , 0 , end );

   if( NDict )
    TranslateF( NPot , NDict , StrtRows , NNodes );
   }

  NPot = 0;
  }

 if( RCst )  // flow reduced costs are required - - - - - - - - - - - - - - -
 {           // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if( WhchRC < NComm ) {
   int start = StrtCols[ WhchRC ];
   int end = StrtCols[ WhchRC + 1 ];
	  
   CPXgetdj( env , lp , RCst , start , end - 1 );

   if( ! Drctd ) {
    start += StrtCols[ NComm ];
    end += StrtCols[ NComm ];

    double* twindj = new double[ end - start ];

    CPXgetdj( env , lp , twindj , start , end - 1 );

    VectSubtract( RCst , twindj , end - start );

    delete[] twindj;
    }

   cIndex_Set ADk = ADict ? ADict[ WhchRC ] : 0;
   if( ADk )
    TranslateK( RCst , ADk , NArcs , start );
   }
  else {  // flow reduced costs of ALL commodities
   int totArcs = StrtCols[ NComm ];

   CPXgetdj( env , lp , RCst , 0 , totArcs - 1 );

   if( ! Drctd ) {
    double* twindj = new double[ totArcs ];

    CPXgetdj( env , lp , twindj , totArcs , 2 * totArcs - 1 );

    for( Index i = totArcs - 1 ; i-- ; )
     RCst[ i ] -= twindj[ i ];

    delete[] twindj;
    }

   if( ADict )
    TranslateF( RCst , ADict , StrtCols , NArcs , C_INF );
   }

  RCst = 0;
  }

 if( MCst )  // mutual capacities dual variables are required - - - - - - - -
 {           // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  int start = StrtRows ? StrtRows[ NComm ] : NNodes * NComm;

  CPXgetpi( env , lp , MCst , start , start + NActvA - 1 );

  // expand the reduced cost vector such that MCst[ i ] will contain RC of
  // arc i; if arc j is not active, then MCst[ j ] will be set to C_INF

  if( ActvArcs ) {
   Index_Set Active = ActvArcs;
   for( Index j = NArcs , h = NActvA - 1 ; --j ; )
    if( Active[ h ] == j )
     MCst[ j ] = MCst[ h-- ];
    else
     MCst[ j ] = C_INF;
   }

  MCst = 0;
  }

 if( XtRC )  // reduced costs of extra variables are required - - - - - - - -
 {           // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  int totArcs = StrtCols[ NComm ];

  double * XtRC_ = new double[ XtrVrs ];
  CPXgetdj( env , lp , XtRC_ , totArcs , totArcs + XtrVrs - 1 );
  VectAssign( XtRC , XtRC_ , XtrVrs  );
  delete[] XtRC_;
  XtRC = 0;
  }

 if( XtDV )  // dual costs of extra constraints are required- - - - - - - - -
 {           // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  int start = ( StrtRows ? StrtRows[ NComm ] : NNodes * NComm ) + NActvA;

  double * XtDV_ = new double[ XtrCnst ];
  CPXgetpi( env , lp , XtDV_ , start , start + XtrCnst - 1 );
  VectAssign( XtDV , XtDV_ , XtrVrs  );
  delete[] XtDV_;
  XtDV = 0;
  }

 return( false );

 }  // end( MMCFCplex::GetDSol )

/*---------------------------------------------------------------------------*/
/*---------- METHODS FOR READING THE DATA OF THE PROBLEM  -------------------*/
/*---------------------------------------------------------------------------*/

MMCFClass::FONumber MMCFCplex::CostOf( void )
{
 assert( false );  // not implemented yet
 }

/*---------------------------------------------------------------------------*/

void MMCFCplex::GetCosts( CRow Csts , cIndex_Set nms , cIndex strt ,
			  Index stp )
{
 if( stp > NArcs * NComm )
  stp = NArcs * NComm;

 if( strt == stp )
  return;

 double *obj = 0;

 int cstrt = strt;    // strt translated into cplex names
 int cstp = stp - 1;  // stp - 1 translated into cplex names
 if( ADict ) {
  while( ( cstrt <= cstp ) &&
	 ( ArcPosKJ( cstrt / NArcs , cstrt % NArcs ) == Inf<Index>() ) )
   cstrt++;

  while( ( cstrt <= cstp ) &&
	 ( ArcPosKJ( cstp / NArcs , cstp % NArcs ) == Inf<Index>() ) )
   cstp--;
  }

 if( cstrt <= cstp ) {
  obj = new double[ cstp - cstrt + 1 ];
  CPXgetobj( env , lp , obj , int( cstrt ) , int( cstp ) );
  }

 if( nms ) {
  while( *nms < strt )
   nms++;

  if( ADict )
   for( Index h ; ( h = *(nms++) ) < stp ; ) {
    cIndex pkj = ArcPosKJ( h / NArcs , h % NArcs );
    if( pkj < Inf<Index>() )
     *(Csts++) = obj[ pkj - cstrt ];
    else
     *(Csts++) = C_INF;
    }
  else
   for( Index h ; ( h = *(nms++) ) < stp ; )
    *(Csts++) = obj[ h - cstrt ];
  }
 else
  if( ADict )
   for( Index h = strt ; h < stp ; h++ ) {
    cIndex pkj = ArcPosKJ( h / NArcs , h % NArcs );
    if( pkj < Inf<Index>() )
     *(Csts++) = obj[ pkj - cstrt ];
    else
     *(Csts++) = C_INF;
    }
  else
   VectAssign( Csts , obj , stp - strt );

 delete[] obj;

 }  // end( MMCFCplex::GetCosts )

/*---------------------------------------------------------------------------*/
 
void MMCFCplex::GetXtrCsts( CRow XtrCs , cIndex_Set nms , cIndex strt ,
			    Index stp )
{
 const int start = int( StrtCols[ NComm ] + strt );
 if( stp > XtrVrs )
  stp = XtrVrs;
 const int end = int( StrtCols[ NComm ] + stp );


 double * XtrCs_ = new double[ end - start ];
 CPXgetobj( env , lp , XtrCs_ , start , end - 1 );
 VectAssign( XtrCs , XtrCs_ , XtrVrs  );
 delete[] XtrCs_;

 if( nms ) {
  while( *nms < strt )
   nms++;

  CRow tXC = XtrCs;
  for( Index h ; ( h = *(nms++) ) < stp ; )
   *(tXC++) = XtrCs[ h - strt ];
  }
 }  // end( MMCFCplex::GetXtrCsts )

/*---------------------------------------------------------------------------*/

void MMCFCplex::GetICaps( FRow ICps, cIndex_Set nms , cIndex strt , Index stp )
{
 if( stp > NArcs * NComm )
  stp = NArcs * NComm;

 if( strt == stp )
  return;

 double *ub = 0;

 int cstrt = strt;    // strt translated into cplex names
 int cstp = stp - 1;  // stp - 1 translated into cplex names
 if( ADict ) {
  while( ( cstrt <= cstp ) &&
	 ( ArcPosKJ( cstrt / NArcs , cstrt % NArcs ) == Inf<Index>() ) )
   cstrt++;

  while( ( cstrt <= cstp ) &&
	 ( ArcPosKJ( cstp / NArcs , cstp % NArcs ) == Inf<Index>() ) )
   cstp--;
  }

 if( cstrt <= cstp ) {
  ub = new double[ cstp - cstrt + 1 ];
  CPXgetub( env , lp , ub , int( cstrt ) , int( cstp ) );
  }

 if( nms ) {
  while( *nms < strt )
   nms++;

  if( ADict )
   for( Index h ; ( h = *(nms++) ) < stp ; ) {
    cIndex pkj = ArcPosKJ( h / NArcs , h % NArcs );
    if( pkj < Inf<Index>() )
     *(ICps++) = ub[ pkj - cstrt ];
    else
     *(ICps++) = 0;
    }
  else
   for( Index h ; ( h = *(nms++) ) < stp ; )
    *(ICps++) = ub[ h - cstrt ];
  }
 else
  if( ADict )
   for( Index h = strt ; h < stp ; h++ ) {
    cIndex pkj = ArcPosKJ( h / NArcs , h % NArcs );
    if( pkj < Inf<Index>() )
     *(ICps++) = ub[ pkj - cstrt ];
    else
     *(ICps++) = 0;
    }
  else
   VectAssign( ICps , ub , stp - strt );

 delete[] ub;

 }  // end( MMCFCplex::GetICaps )

/*---------------------------------------------------------------------------*/

void MMCFCplex::GetMCaps( FRow MCps , cIndex_Set nms , cIndex strt ,
			  Index stp )
{
 if( stp > NArcs )
  stp = NArcs;

 if( nms )
  while( *nms < strt )
   nms++;

 if( ! NActvA ) {  // there are no mutual capacity constraints at all
  if( nms )
   for( Index h ; ( h = *(nms++) ) < stp ; )
    *(MCps++) = F_INF;
  else
   for( Index h = strt ; h < stp ; h++ )
    *(MCps++) = F_INF;

  return;
  }

 double *cap = new double[ NActvA ];
 const int strtMC = StrtRows ? StrtRows[ NComm ] : NNodes * NComm;
 CPXgetrhs( env , lp , cap , strtMC , strtMC + NActvA - 1 );

 if( nms )
  if( ActvArcs ) {
   Index j = 0;
   for( Index h ; ( h = *(nms++) ) < stp ; ) {
    while( ActvArcs[ j ] < h )
     j++;

    if( ActvArcs[ j ] == h )
     *(MCps++) = cap[ j++ ];
    else
     *(MCps++) = F_INF;
    }
   }
  else
   for( Index h ; ( h = *(nms++) ) < stp ; )
    *(MCps++) = cap[ h ];
 else
  if( ActvArcs ) {
   Index j = 0;
   for( Index h = strt ; h < stp ; h++ ) {
    while( ActvArcs[ j ] < h )
     j++;

    if( ActvArcs[ j ] == h )
     *(MCps++) = cap[ j++ ];
    else
     *(MCps++) = F_INF;
    }
   }
  else
   for( Index h = strt ; h < stp ; h++ )
    *(MCps++) = cap[ h ];

 delete[] cap;

 }  // end( MMCFCplex::GetMCaps )

/*--------------------------------------------------------------------------*/
/*------------- METHODS FOR ADDING / REMOVING / CHANGING DATA --------------*/
/*--------------------------------------------------------------------------*/

void MMCFCplex::ChgCosts( cCRow NwCsts , cIndex_Set nms , cIndex strt ,
			  Index stp )
{
 if( CPXgetprobtype( env , lp ) == CPXPROB_LP )
  CPXchgprobtype( env , lp , CPXPROB_MILP );

 if( nms ){
  while( *nms < strt ) {
   NwCsts++;
   nms++;
   }

  if( ADict )
   for( Index h ; ( h = *(nms++) ) < stp ; NwCsts++ ) {
    cIndex pkj = ArcPosKJ( h / NArcs , h % NArcs );
    if( pkj < Inf<Index>() )
     CPXchgcoef( env , lp , -1 , int( pkj ) , double( *NwCsts ) );
    }
  else
   for( Index h ; ( h = *(nms++) ) < stp ; )
    CPXchgcoef( env , lp , -1 , int( h ) , double( *(NwCsts++) ) );
  }
 else {
  if( stp > NArcs * NComm )
   stp = NArcs * NComm;

  if( ADict )
   for( Index h = strt ; h < stp ; h++ , NwCsts++ ) {
    cIndex pkj = ArcPosKJ( h / NArcs , h % NArcs );
    if( pkj < Inf<Index>() )
     CPXchgcoef( env , lp , -1 , int( pkj ) , double( *NwCsts ) );
    }
  else
   for( Index h = strt ; h < stp ; h++ )
    CPXchgcoef( env , lp , -1 ,	int( h ) , double( *(NwCsts++) ) );
  }
 }  // end( MMCFCplex::ChgCosts )

/*---------------------------------------------------------------------------*/

void MMCFCplex::ChgXtrBnds( cRow XLr , cRow XUr , cIndex_Set nms ,
		cIndex strt , Index stp ) {

 if( stp > NArcs )
  stp = NArcs;

 if( strt == stp )
  return;

 // some exceptions- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( XUr && XLr || ( !XUr && !XLr ) )
  throw( MMCFException(
    "MMCFCplex::ChgXtrBnds(): this should not happen" ) );

 if( !XtrVrs )
  throw( MMCFException(
	  "MMCFCplex::ChgXtrBnds(): this should not happen" ) );

 // either upper bound or lower bound has to be changed - - - - -  - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 char bt;
 int bi;

 if( XLr )
  bt = 'L';
 else
  bt = 'U';

 // the directed graph is not considered - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( ! Drctd )
  throw( MMCFException(
    "MMCFCplex::ChgXtrBnds(): this should not happen" ) );

 Index start = StrtCols[ NComm ];

 if( nms ) {
  while( *nms < strt ) {
   XLr? XLr++ : XUr++;
   nms++;
   }

  for( Index h ; ( h = *(nms++) ) < stp ; ) {
   bi = h + start;
   if( bt =='L' )
    CPXchgbds( env, lp, 1 , &bi , &bt , XLr++ );
   else
    CPXchgbds( env, lp, 1 , &bi , &bt , XUr++ );
   }

  }
 else {
  if( stp > XtrVrs )
   stp = XtrVrs;

  for( Index i = strt ; i < stp  ; i++ ) {
   bi = i + start;
   if( bt == 'L' )
    CPXchgbds( env, lp, 1 , &bi , &bt , XLr++ );
   else
    CPXchgbds( env, lp, 1 , &bi , &bt , XUr++ );
   }
  }

 } // end( MMCFCplex::ChgXtrBnds( ) )

/*---------------------------------------------------------------------------*/

void MMCFCplex::ChgXtrCsts( cCRow NwXtrCs , cIndex_Set nms , cIndex strt ,
			    Index stp )
{ 
 if( CPXgetprobtype( env , lp ) == CPXPROB_LP )
  CPXchgprobtype( env , lp , CPXPROB_MILP );

 Index start = StrtCols[ NComm ];

 if( nms ) {
  while( *nms < strt ) {
   NwXtrCs++;
   nms++;
   }

  for( Index h ; ( h = *(nms++) ) < stp ; )
   CPXchgcoef( env , lp , -1 , int( start + h ) , *(NwXtrCs++) );
  }
 else {
  if( stp > XtrVrs )
   stp = XtrVrs;

  for( Index i = stp ; i-- > strt ; )
   CPXchgcoef( env , lp , -1 , int( start++ ) , *(NwXtrCs++) );
  } 
 }  // end( MMCFCplex::ChgXtrCsts )

/*--------------------------------------------------------------------------*/

void MMCFCplex::ChgICaps( cFRow NwICps , cIndex_Set nms , cIndex strt ,
			  Index stp )
{
 if( CPXgetprobtype( env , lp ) == CPXPROB_LP )
  CPXchgprobtype( env , lp , CPXPROB_MILP );

 if( stp > NArcs * NComm )
  stp = NArcs * NComm;

 // first count how many individual capacities must be changed- - - - - - - -

 int cnt = 0;
 if( nms ) {
  while( *nms < strt ) {
   NwICps++;
   nms++;
   }

  while( nms[ cnt ] < stp )
   cnt++;
  }
 else
  cnt = stp - strt;

 // allocate temporaries- - - - - - - - - - - - - - - - - - - - - - - - - - -

 char *lu = new char[ cnt ];
 VectAssign( lu , 'U' , cnt );

 int *indices = new int[ cnt ];
 double *bd = new double[ cnt ];

 int* ti = indices;
 double* tbd = bd;

 // construct the vectors - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( nms )
  if( ADict )
   for( Index h ; ( h = *(nms++) ) < stp ; NwICps++ ) {
    cIndex pkj = ArcPosKJ( h / NArcs , h % NArcs );
    if( pkj < Inf<Index>() ) {
     *(ti++) = pkj;
     *(tbd++) = *NwICps;
     }
    else
     cnt--;
    }
  else
   for( Index h ; ( h = *(nms++) ) < stp ; ) {
    *(ti++) = h;
    *(tbd++) = *(NwICps++);
    }
 else
  if( ADict )
   for( Index h = strt ; h < stp ; h++ , NwICps++ ) {
    cIndex pkj = ArcPosKJ( h / NArcs , h % NArcs );
    if( pkj < Inf<Index>() ) {
     *(ti++) = pkj;
     *(tbd++) = *NwICps;
     }
    else
     cnt--;
    }
  else
   for( Index h = strt ; h < stp ; h++ ) {
    *(ti++) = h;
    *(tbd++) = *(NwICps++);
    }

 // call the function - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 CPXchgbds( env , lp , cnt , indices , lu , bd );

 // deallocate the temporaries- - - - - - - - - - - - - - - - - - - - - - - -

 delete[] indices;
 delete[] bd;
 delete[] lu;

 }  // end( MMCFCplex::ChgICaps )

/*--------------------------------------------------------------------------*/

void MMCFCplex::ChgMCaps( cFRow NwMCps, cIndex_Set nms , cIndex strt ,
			  Index stp )
{
 if( CPXgetprobtype( env , lp ) == CPXPROB_LP )
  CPXchgprobtype( env , lp , CPXPROB_MILP );

 cIndex_Set pActive = ActvArcs;
 int strtMC = StrtRows ? StrtRows[ NComm ] : NNodes * NComm;

 if( nms ) {
  while( *nms < strt ) {
   NwMCps++;
   nms++;
   }

  if( pActive )  // nms && pActive - - - - - - - - - - - - - - - - - - - - - -
   for( Index h ; ( h = *(nms++) ) < stp ; NwMCps++ ) {
    while( *pActive < h )
     pActive++;

    if( *pActive == h )
     CPXchgcoef( env , lp , strtMC + ( pActive - ActvArcs ) , -1 , *NwMCps );
    }
  else           // nms && ! pActive - - - - - - - - - - - - - - - - - - - - -
   for( Index h ; ( h = *(nms++) ) < stp ; NwMCps++ )
    CPXchgcoef( env , lp , strtMC + h , -1, *NwMCps ) ;
  }
 else {
  if( stp > NArcs )
   stp = NArcs;

  strtMC += strt;

  if( pActive )  // ! nms && pActive - - - - - - - - - - - - - - - - - - - - -
   for( Index i = strt ; i < stp ; i++ , NwMCps++ ) {
    if( *pActive == i ) {
     CPXchgcoef( env , lp , strtMC++ , -1 , *NwMCps );
     pActive++;
     }
    }
  else           // ! nms && ! pActive - - - - - - - - - - - - - - - - - - - -
   for( Index i = strt ; i < stp ; i++ )
    CPXchgcoef( env , lp , strtMC++ , -1, *(NwMCps++) ) ;
  }
 }  // end( MMCFCplex::ChgMCaps )

/*--------------------------------------------------------------------------*/

void MMCFCplex::ChgIntVar( cIndex k , bool IntVld , cIndex_Set nms ,
			   Index strt , Index stp )
{
 // initialize dimensions - - - - - - - - - - - - - - - - - - - - - - - - - -

 cIndex maxnm = k < NComm ? NArcs : ( k == NComm ? XtrVrs : NComm * NArcs );
 if( stp > maxnm )
  stp = maxnm;

 // count how many variables are (potentially) involved - - - - - - - - - - -

 int cnt = 0;
 if( nms ) {
  while( *nms < strt )
   nms++;

  while( nms[ cnt ] < stp )
   cnt++;
  }
 else
  cnt = stp - strt;

 if( ! cnt )
  return;

 // allocate temporaries- - - - - - - - - - - - - - - - - - - - - - - - - - -

 int *indices = new int[ cnt ];
 char *ctype = new char[ cnt ];

 // now construct the temporaries - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 cIndex shftr = k <= NComm ? k * NArcs : 0;  // "shifter" of arc names
 cnt = 0;                                    // variables actually changed

 if( IntVld ) {  // set variables to integer- - - - - - - - - - - - - - - - -
  // collect the ub[] to tell binary variables from general integer ones

  int cstrt = strt + shftr;    // strt translated into cplex names
  int cstp = stp + shftr - 1;  // stp - 1 translated into cplex names
  if( ADict && ( k != NComm ) ) {
   while( ( cstrt <= cstp ) &&
	  ( ArcPosKJ( cstrt / NArcs , cstrt % NArcs ) == Inf<Index>() ) )
    cstrt++;

   while( ( cstrt <= cstp ) &&
	  ( ArcPosKJ( cstp / NArcs , cstp % NArcs ) == Inf<Index>() ) )
    cstp--;

   cstrt = ArcPosKJ( cstrt / NArcs , cstrt % NArcs );
   cstp = ArcPosKJ( cstp / NArcs , cstp % NArcs );
   }

  double *ub = 0;
  if( cstrt <= cstp ) {
   ub = new double[ cstp - cstrt + 1 ];
   CPXgetub( env , lp , ub , int( cstrt ) , int( cstp ) );
   }

  if( ! nms )
   stp += shftr;

  if( ADict )
   if( nms )
    for( Index h ; ( h = *(nms++) ) < stp ; ) {
     h += shftr;
     cIndex pkj = ArcPosKJ( h / NArcs , h % NArcs );
     if( pkj < Inf<Index>() ) {
      ctype[ cnt ] = ( ub[ pkj - cstrt ] == 1 ? CPX_BINARY : CPX_INTEGER );
      indices[ cnt++ ] = pkj;
      }
     }
   else
    for( Index h = strt + shftr ; h < stp ; h++ ) {
     cIndex pkj = ArcPosKJ( h / NArcs , h % NArcs );
     if( pkj < Inf<Index>() ) {
      ctype[ cnt ] = ( ub[ pkj - cstrt ] == 1 ? CPX_BINARY : CPX_INTEGER );
      indices[ cnt++ ] = pkj;
      }
     }
  else
   if( nms )
    for( Index h ; ( h = *(nms++) ) < stp ; ) {
     ctype[ cnt ] = ( ub[ h + shftr - cstrt ] == 1 ?
		      CPX_BINARY : CPX_INTEGER );
     indices[ cnt++ ] = h + shftr;
     }
   else
    for( Index h = strt + shftr ; h < stp ; h++ ) {
     ctype[ cnt ] = ( ub[ h - cstrt ] == 1 ? CPX_BINARY : CPX_INTEGER );
     indices[ cnt++ ] = h;
     }

  delete[] ub;
  }
 else {        // set variables to continuous - - - - - - - - - - - - - - - -
  if( ! nms )
   stp += shftr;

  if( ADict )
   if( nms )
    for( Index h ; ( h = *(nms++) ) < stp ; ) {
     h += shftr;
     cIndex pkj = ArcPosKJ( h / NArcs , h % NArcs );
     if( pkj < Inf<Index>() ) {
      ctype[ cnt ] = CPX_CONTINUOUS;
      indices[ cnt++ ] = pkj;
      }
     }
   else
    for( Index h = strt + shftr ; h < stp ; h++ ) {
     cIndex pkj = ArcPosKJ( h / NArcs , h % NArcs );
     if( pkj < Inf<Index>() ) {
      ctype[ cnt ] = CPX_CONTINUOUS;
      indices[ cnt++ ] = pkj;
      }
     }
  else
   if( nms )
    for( Index h ; ( h = *(nms++) ) < stp ; ) {
     ctype[ cnt ] = CPX_CONTINUOUS;
     indices[ cnt++ ] = h + shftr;
     }
   else
    for( Index h = strt + shftr ; h < stp ; h++ ) {
     ctype[ cnt ] = CPX_CONTINUOUS;
     indices[ cnt++ ] = h;
     }

  }  // end( else( set vars to continuous )

 // double data for undirected graphs - - - - - - - - - - - - - - - - - - - -

 if( ! Drctd ) {
  cIndex shift = cnt / 2;
  VectAdd( indices + shift , indices , int( StrtCols[ NComm ] ) , shift );
  VectAssign( ctype + shift , ctype , shift );
  }

 // call the function and deallocate the temporaries- - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( CPXgetprobtype( env , lp ) == CPXPROB_LP )
  CPXchgprobtype( env , lp , CPXPROB_MILP );

 if( IntVld ) {
  if( CPXgetprobtype( env , lp ) != CPXPROB_MILP )  // if it is *not* a MIP
   CPXchgprobtype( env , lp , CPXPROB_MILP );      // make it so
  }

 CPXchgctype( env , lp , cnt , indices , ctype );

 if( ! IntVld )  // if some variables have just been set to continuous
  if( ! ( CPXgetnumbin( env , lp ) + CPXgetnumint( env , lp ) ) )
   // and these happened to be the only integer variables left in the problem
   CPXchgprobtype( env , lp , CPXPROB_LP );  // set it to LP

 delete[] indices;
 delete[] ctype;

 }  // end( MMCFCplex::ChgIntVar )

/*--------------------------------------------------------------------------*/

#if( CHGARCS_MMCF )

void MMCFCplex::CloseArcs( cIndex_Set whch )
{
 if( XtrVrs != NArcs )
 throw( MMCFException("MMCFCplex::CloseArcs: no extra variables" ) );

 Index count = 0;
 for( Index j = 0 ; whch[ j++ ] < Inf<Index>() ; )
  count++;

 FRow XUr = new FNumber[ count ];

 // set to 0 the upper bounds on extra variables- - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 VectAssign( XUr , FNumber( 0 ) , count );
 ChgXtrBnds( 0 , XUr , whch );

 delete[] XUr;

 } // end ( MMCFCplex::CloseArcs )

/*--------------------------------------------------------------------------*/

void MMCFCplex::OpenArcs( cIndex_Set whch )
{
 if( XtrVrs != NArcs )
  throw( MMCFException( "MMCFCplex::OpenArcs: no extra variables" ) );

 Index count = 0;
 for( Index j = 0 ; whch[ j++ ] < Inf<Index>() ; )
  count++;

 FRow XUr = new FNumber[ count ];

 // set to 0 the upper bounds on extra variables- - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 VectAssign( XUr , FNumber( 1 ) , count );
 ChgXtrBnds( 0 , XUr , whch );

 delete[] XUr;

 } // end ( MMCFCplex::OpenArcs )

/*--------------------------------------------------------------------------*/

#if( CHGARCS_MMCF > 1 )

void MMCFCplex::CloseArcs( cIndex k, cIndex_Set whch )
{
 assert( false );  // not implemented yet
 }

/*--------------------------------------------------------------------------*/

void MMCFCplex::OpenArcs( cIndex k, cIndex_Set whch )
{
 assert( false );  // not implemented yet
 }

#endif
#endif

/*--------------------------------------------------------------------------*/

void MMCFCplex::AddExtraVars( cIndex NXV , cCRow XCst , cRow XUb ,
			      cRow XLb , bool IntVar , cIndex_Set nms )
{
 // construct the vector telling the type of each "extra" variable- - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 char *ctype = 0;

 if( IntVar ) {
  ctype = new char[ NXV ];

  if( nms )  // some variables are integer, some are not- - - - - - - - - - -
   for( Index i = 0 ; i < NXV ; i++ )
    if( *nms == i ) {
     nms++;
     if( ( ( ! XLb ) || ( XLb[ i ] == 0 ) ) && XUb && ( XUb[ i ] == 1 ) )
      ctype[ i ] = CPX_BINARY;
     else
      ctype[ i ] = CPX_INTEGER;
     }
    else
     ctype[ i ] = CPX_CONTINUOUS;
  else       // all variables are integer - - - - - - - - - - - - - - - - - -
   for( Index i = 0 ; i < NXV ; i++ )
    if( ( ( ! XLb ) || ( XLb[ i ] == 0 ) ) && XUb && ( XUb[ i ] == 1 ) )
     ctype[ i ] = CPX_BINARY;
    else
     ctype[ i ] = CPX_INTEGER;
  }

 // call the function and deallocate the temporaries- - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // in the call, remove the "const" from XCst, XLb and XUb by casting: this
 // is necessary because those guys at Cplex did not care to add the "const"
 // to the definition of their functions, even though constness is respected
 
 if( CPXgetprobtype( env , lp ) == CPXPROB_LP )
  CPXchgprobtype( env , lp , CPXPROB_MILP );

 double *XCst_ = new double[ NXV ];
 VectAssign( XCst_ , XCst , NXV );
 CPXnewcols( env , lp , int( NXV ) , XCst_ , (double *)XLb ,
	     (double *)XUb , ctype , 0 );
 delete[] XCst_;

 delete[] ctype;

 }  // end( MMCFCplex::AddExtraVars )

/*--------------------------------------------------------------------------*/

void MMCFCplex::AddExtraConstr( cIndex NXC , int *IBeg , int *Indx ,
				double *Vals , const double *XLr , const double *XUr )
{
 assert( Drctd );  // only works in the directed case

 // allocate space for data structures to be passed to CPLEX- - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 double* rhs = new double[ NXC ];
 char* sense = new char[ NXC ];

 // compute the "sense" of the constraints- - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 assert( XLr || XUr );

 if( XLr && XUr )
  for( Index i = 0 ; i < NXC ; i++ ) {
   const double XLri = XLr[ i ];
   const double XUri = XUr[ i ];
   if( XLri > -F_INF ) {
    rhs[ i ] = XLri;
    sense[ i ] = ( XUri < F_INF ? ( ( XLri == XUri ) ? 'E' : 'R' ) : 'G' );
    }
   else {
    rhs[ i ] = XUri;
    sense[ i ] = 'L';
    }
   }  // end( for( i ) )
  else
   if( XLr ) {
    VectAssign( rhs , XLr , NXC );
    VectAssign( sense , 'G' , NXC );
    }
   else {
    VectAssign( rhs , XUr , NXC );
    VectAssign( sense , 'L' , NXC );
    }

 // translate the names of variables into Cplex's language- - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( ADict ) {
  int *tIx = Indx;
  int mxnrmvr = NArcs * NComm;
  int dltnm = StrtCols[ NComm ] - mxnrmvr;
  for( Index i = IBeg[ NXC ] ; i-- ; tIx++ )
   if( *tIx >= mxnrmvr )
    *tIx += dltnm;
   else
    *tIx = ArcPosKJ( *tIx / NArcs , *tIx % NArcs );
  }

 // call the function - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( CPXgetprobtype( env , lp ) == CPXPROB_LP )
  CPXchgprobtype( env , lp , CPXPROB_MILP );

 int numrows = CPXgetnumrows( env , lp );

 if( Lazy ) {
  // add extra constraints as "lazy" ones
  SetCplexParam( CPX_PARAM_REDUCE , 0 );
  CPXaddlazyconstraints( env , lp , int( NXC ) , int( IBeg[ NXC ] ) , rhs ,
  	sense , IBeg , Indx , Vals , 0 );
  }
 else
   // add extra constraints as ordinary ones
  CPXaddrows( env , lp , 0 , int( NXC ) , int( IBeg[ NXC ] ) , rhs , sense ,
    IBeg , Indx , Vals , 0 , 0 );


 // add range values for R constraints- - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 
 for( Index i = 0 ; i < NXC ; i++ )
  if( sense[ i ] == 'R' )
   CPXchgcoef( env , lp , numrows + i , -2 , XUr[ i ] - XLr[ i ] );

 // delete temporaries- - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 delete[] rhs;
 delete[] sense;

 }  // end( MMCFCplex::AddExtraConstr )

/*--------------------------------------------------------------------------*/

Index MMCFCplex::GetNode( void ) {

 return( CPXgetnodecnt( env, lp ) );
 }


/*--------------------------------------------------------------------------*/
/*------------------------------ DESTRUCTOR --------------------------------*/
/*--------------------------------------------------------------------------*/

MMCFCplex::~MMCFCplex()
{
 // close the LP- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 CPXfreeprob( env , &lp );

 // release the environment - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( env == genv )
  if( ! --EnvICntr ) {
   CPXcloseCPLEX( &genv );
   genv = 0;
   }

 // deallocate other memory - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( NDict ) {
  for( Index k = NComm ; k-- ; )
   delete[] NDict[ k ];

  delete[] NDict;
  }

 if( ADict ) {
  for( Index k = NComm ; k-- ; )
   delete[] ADict[ k ];

  delete[] ADict;
  }

 delete[] StrtCols;
 delete[] StrtRows;

 delete[] ActvArcs;

 }  // end( MMCFCplex:: ~MMCFCplex )

/*--------------------------------------------------------------------------*/
/*------------------------- PRIVATE METHODS --------------------------------*/
/*--------------------------------------------------------------------------*/

template<class T>
inline void MMCFCplex::TranslateK( T *const vect , const cIndex_Set Dict ,
				   Index n , cIndex strt  ,
				   const T Dflt )
{
 while( n-- )
 if( Dict[ n ] < Inf<Index>() )
  vect[ n ] = vect[ Dict[ n ] - strt ];
 else
  vect[ n ] = Dflt;

 }  // end( MMCFCplex::TranslateK )

/*---------------------------------------------------------------------------*/

template<class T>
inline void MMCFCplex::TranslateF( T *const vect , Index_Mat Dict ,
				   cIndex_Set Strt , Index n ,
				   const T Dflt )
{
 Index fk = 0;          // first commodity for which the standard names are
 while( ! Dict[ fk ] )  // different from the cplex ones
  fk++;

 T* tv = vect + NComm * n;
 for( Index k = NComm ; k-- > fk ; )
  if( Dict[ k ] ) {
   cIndex_Set Dk = Dict[ k ] + n;
   for( Index i = n ; i-- ; )
    if( *(--Dk) < Inf<Index>() )
     *(--tv) = vect[ *Dk ];
    else
     *(--tv) = Dflt;
   }
  else {
   MFRow tvk = vect + Strt[ k + 1 ];
   for( Index i = n ; i-- ; )
    *(--tv) = *(--tvk);
   }

 }  // end( MMCFCplex::TranslateF )

/*--------------------------------------------------------------------------*/

void MMCFCplex::GetExtraSol( void )
{
 int totArcs = StrtCols[ NComm ];

 if( CA == kMIP )
  CPXgetsolnpoolx( env , lp , IDSol , XSol , totArcs ,
	 totArcs + XtrVrs - 1 );
 else
  CPXgetx( env , lp , XSol , totArcs , totArcs + XtrVrs - 1 );

 if( XBse ) {  // sparse format required
  for( Index i = 0 , q = 0 ; i < XtrVrs ; i++ )
   if( XSol[ i ] ) {
    XSol[ q ] = XSol[ i ];
    XBse[ q++ ] = i;
    }
      
  XBse = 0;
  }
 }  // end( MMCFCplex::GetExtraSol )

/*-------------------------------------------------------------------------*/
/*------------------initialize static members------------------------------*/
/*-------------------------------------------------------------------------*/

CPXENVptr MMCFCplex::genv = 0;

Index MMCFCplex::EnvICntr = 0;

/*--------------------------------------------------------------------------*/
/*------------------------ End file MMCFCple.C -----------------------------*/
/*--------------------------------------------------------------------------*/

/*

void MMCFUCFLCplex::GetBounds( Row XLr , Row XUr , cIndex_Set nms ,
		cIndex strt , Index stp ) {
 if( stp > NArcs )
  stp = NArcs;

 // some exceptions- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( strt == stp )
   throw( MMCFException(
    "MMCFCplex::GetBounds(): this should not happen" ) );

 if( XUr && XLr || ( !XUr && !XLr ) )
  throw( MMCFException(
    "MMCFCplex::GetBounds(): this should not happen" ) );

 // the directed graph is not considered - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( ! Drctd )
  throw( MMCFException(
     "MMCFCplex::GetBounds(): this should not happen" ) );

 Index start = StrtCols[ NComm ];
 int bi;
 char bt;

 if( XLr )
  bt = 'L';
 else
  bt = 'U';

 if( nms ) {
  while( *nms < strt ) {
   XLr? XLr++ : XUr++;
   nms++;
   }

  if( bt == 'L' )
   for( Index h ; ( h = *(nms++) ) < stp ; ) {
    bi = h + start;
    CPXgetlb( env , lp, XLr++ , bi , bi );
    }
  else
   for( Index h ; ( h = *(nms++) ) < stp ; ) {
    bi = h + start;
    CPXgetub( env , lp, XUr++ , bi , bi );
    }

  }
 else {
  if( stp > XtrVrs )
   stp = XtrVrs;

  if( bt == 'L' )
   for( Index i = strt ; i < stp  ; i++ ) {
    bi = i + start;
    CPXgetlb( env , lp, XLr++ , bi , bi );
    }
  else
   for( Index i = strt ; i < stp  ; i++ ) {
	bi = i + start;
    CPXgetub( env , lp, XUr++ , bi , bi );
    }
  }

 } //end ( MMCFUCFLCplex::GetBounds ) */

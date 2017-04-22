/*--------------------------------------------------------------------------*/
/*---------------------------- File FlwFiOrcl.C ----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the FlwFiOrcl class, which implements the FiOracle
 * interface, and solves multicommodity minimal cost network flow problems,
 * possibly with fixed charges, using the flow Lagrangian Relaxation.
 *
 * \version 1.10
 *
 * \date 19 - 11 - 2014
 *
 * \author Antonio Frangioni \n
 *         Operations Research Group \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Enrico Gorgone \n
 *         Operations Research Group \n
 *         Dipartimento di Elettronica Informatica e Sistemistica\n
 *         Universita' della Calabria \n
 *
 * Copyright(C) 2001 - 2012 by Antonio Frangioni, Enrico Gorgone
 */

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- IMPLEMENTATION -------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "FlwFiOrcl.h"

#include "OPTvect.h"

#include <fstream>
#include <sstream>

#include "NDOSlver.h"
#include <algorithm>    // std::find

/*--------------------------------------------------------------------------*/
/*-------------------------------- MACROS ----------------------------------*/
/*--------------------------------------------------------------------------*/

#define SubVsVar 1
// it's used in Aggregate() to distinguish the way to aggregate - - - - - - -
//     1 : aggregate for variables
//     0 : aggregate for items

#if( LOG_FI > 0 )
 #define FILOG( y , x ) if( FiLLvl >= y ) *FiLog << x
 #define FILOG2( y , c , x ) if( ( FiLLvl >= y) && c ) *FiLog << x
#else
 #define FILOG( y , x )
 #define FILOG2( y , c , x )
#endif

/*--------------------------------------------------------------------------*/
/*-------------------------------- USING -----------------------------------*/
/*--------------------------------------------------------------------------*/

#if( OPT_USE_NAMESPACES )
  using namespace NDO_di_unipi_it;
#endif

/*--------------------------------------------------------------------------*/
/*--------------------- IMPLEMENTATION OF FlwFiOrcl ------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- PUBLIC METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/
/*---------------------------- CONSTRUCTOR ---------------------------------*/
/*--------------------------------------------------------------------------*/

FlwFiOrcl::FlwFiOrcl( Graph *g , istream *iStrm )
           :
           FiOracle() , MMCFFlwBase( g )
{
 // initialization of relaxation  - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // by default is considered the weakly formulation, the strong linking
 // constraints are not relaxed	and all alphas are there from the beginning

 Index YiE;  // define the options for the Lagrangian relaxation
 Index agg;

 DfltdSfInpt( iStrm , SPar1 , Index( 0 ) );
 DfltdSfInpt( iStrm , SPar2 , Index( 0 ) );
 DfltdSfInpt( iStrm , SPar3 , Index( 0 ) );
 DfltdSfInpt( iStrm , SPar4 , char( 'w' ) );

 DfltdSfInpt( iStrm , YiE , Index( 0 ) );
 DfltdSfInpt( iStrm , agg , Index( 0 ) );

 DfltdSfInpt( iStrm , KOld , SIndex( 0 ) );

 DfltdSfInpt( iStrm , EpsFi , HpNum( 1e-6 ) );
 DfltdSfInpt( iStrm , EpsCon , HpNum( 1e-6 ) );

 // allocate memory for instance data - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 // if the "extra" variables in graph are less than Nr arcs,
 // we does not consider them

 if( g->NrExtraVars() == NArcs )    // if fixed-charge costs are defined
  OrigXtrCosts = new CNumber[ NArcs ];  // f_{ij} : original costs
 else
  OrigXtrCosts = 0;

 // mutual capacities - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( SPar4 != 's' ) {
  OrigTotCap = new FNumber[ NArcs ];             // u_{ij} : mutual capacities

  for( Index j = 0 ; j < NArcs ; j++ )
   OrigTotCap[ j ] = g->TotalCapacityJ( j );
  }
 else
   OrigTotCap = 0;

 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 OrigCosts = new CNumber[ NArcs * NComm ];      // c_{ij}^k : original unit costs
 OrigCapacities = new CNumber[ NArcs * NComm ]; // u_{ij}^k : original individual capacities
 OrigDeficits = new FNumber[ NComm * NNodes ];  // b_i^k : original deficits

 // read instance data- - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 // Associated with each arc there must be a *finite* capacity
 // If the instance does not have it, Graph::PreProcess() computes it
 // Thus, PreProcess MUST be CALLED before the FiOracle constructor.

 if( OrigXtrCosts )
  for( Index j = 0 ; j < NArcs ; j++ )
   OrigXtrCosts[ j ] = g->CostKJ( NComm , j );

 for( Index k = 0 ; k < NComm ; k++ ) {
  for( Index j = 0 ; j < NArcs ; j++ )  {
   OrigCosts[ k * NArcs + j ] = g->CostKJ( k , j );
   OrigCapacities[ k * NArcs + j ] = g->CapacityKJ( k , j );
   }
  for(Index n = 0 ; n < NNodes ; n++ )
   OrigDeficits[ k * NNodes + n ] = g->DeficitKJ( k , n );
  }

 // define all quantities relative to the relaxation type - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 InvULambda = 0;
 LmbCmp = 0;
 ZLambdaCount = 0;
 Addt = Remt = 0;

 SGBse1 = 0;
 XSolution = 0;
 CoefObj = 0;

 SetRelax( SPar4 , YiE , SPar1 , SPar2 , SPar3 );

 // define the aggregation type - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 LsHasChgd = 0;
 SolvedP = 0;
 SetAggregate( agg );

 // allocate memory for primal solution - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 FSolution = new FNumber[ NArcs * NComm ];
 SolWFi = Inf<Index>();

 // other initializations - - - - - - - - - - - - - - - - - - - - - - - - - -

 OldWFi = 0;
 OldFSols = OldXSols = 0;
 MaxName = 0;

 LstChgItr = 0;

 // set relative precision for objective function  - - - - - - - - - - - - - -

 SetOptEps( EpsFi / NComm );
 SetFsbEps( EpsFi );

 FindGlobalLipschitz();
 LowerBound = - Inf<HpNum>();
 LastGiName = Inf<Index>();

 } // end ( FlwFiOrcl::FlwFiOrcl )

/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

void FlwFiOrcl::SetFiLog( ostream *outs , const char lvl ) {

 FiOracle::SetFiLog( outs , lvl );

 #if( LOG_FI )
  if( FiLLvl > 1 )
   *FiLog << endl << "FlowFiOrcl: YisEasy = " << ( YIsEasy? "True ": "False " )
    << " ~ Aggrgtd = " << ( Aggrgtd? "True ": "False " ) << endl <<
	"SPar = " << SPar1 << " | " << SPar2 << " | " <<
	SPar3 << " | " << SPar4 << " ~ EpsFi = " << EpsFi << " ~ EpsCon = "
	<< EpsCon << endl;
 #endif

 } // end ( FlwFiOrcl::SetFiLog )

/*--------------------------------------------------------------------------*/

void FlwFiOrcl::SetMaxName( cIndex MxNme ) {

  if( MxNme > MaxName ) {

   if( !MaxName ) {                          // allocs memory to keep
	                                         // MxNme names in memory
    OldWFi = new Index[ MxNme ];
    OldFSols = new FRow[ MxNme ];
    if( OrigXtrCosts && !YIsEasy )
     OldXSols = new Row[ MxNme ];
    else
     OldXSols = 0;

    Index KOld_ = KOld > 0 ? KOld : -KOld;
    for( Index i = 0 ; i < min( MxNme , KOld_ ) ; i++ )
     if( Aggrgtd ) {
      OldFSols[ i ] = new FNumber[ NArcs * NComm ];
      if( OrigXtrCosts && !YIsEasy )
       OldXSols[ i ] = new Number[ NArcs ];
      OldWFi[ i ] = Inf<Index>();
      }
     else
      if( KOld > 0 ) {  // set OldWFi as the 1st component
       OldFSols[ i ] = new FNumber[ NArcs ];
       if( OrigXtrCosts && !YIsEasy )
        OldXSols[ i ] = 0;
       OldWFi[ i ] = 1;
       }
      else {    // set OldWFi as the 1st component
       OldFSols[ i ] = new FNumber[ NArcs * NComm ];
       if( OrigXtrCosts && !YIsEasy )
        OldXSols[ i ] = new Number[ NArcs ];
       OldWFi[ i ] = Inf<Index>();
       }

    for( Index i = min( MxNme , KOld_ ) ; i < MxNme ; i++ ) {
     OldFSols[ i ] = 0;
     if( OrigXtrCosts && !YIsEasy )
      OldXSols[ i ] = 0;
     }

    MaxName = MxNme;

    }
   else {                                     // reallocs memory to keep
	                                          // MxNme names in memory
	Index_Set newWFi = new Index[ MxNme ];
	VectAssign( newWFi , OldWFi , MaxName );
	delete[] OldWFi;
	OldWFi = newWFi;

	FRow *newFSols = new FRow[ MxNme ];
	VectAssign( newFSols , OldFSols , MaxName );
	delete[] OldFSols;
	OldFSols = newFSols;
    for( Index i = MaxName ; i < MxNme ; i++ )
     OldFSols[ i ] = 0;

    if( OrigXtrCosts && !YIsEasy ) {
     Mat newXSols = new Row[ MxNme ];
     VectAssign( newXSols , OldXSols , MaxName );
     delete[] OldXSols;
     OldXSols = newXSols;
     for( Index i = MaxName ; i < MxNme ; i++ )
      OldXSols[ i ] = 0;
     }
    else {
     if( OldXSols ) {
      for( Index i = 0 ; i < MaxName ; i++ )
       delete[] OldXSols[i];
      delete[] OldXSols;
      OldXSols = 0;
      }
     }

    MaxName = MxNme;
    }

   }
  else
   if( MxNme == 0 ) {

	if( OldWFi ) {
	 delete[] OldWFi;
	 OldWFi = 0;
	 }

	if( OldXSols ) {
     for( Index i = 0 ; i < MaxName ; i++ )
      delete[] OldXSols[ i ];
     delete[] OldXSols;
     OldXSols = 0;
     }

	if( OldFSols ) {
     for( Index i = 0 ; i < MaxName ; i++ )
      delete[] OldFSols[ i ];
     delete[] OldFSols;
     OldFSols = 0;
     }

    MaxName = 0;

    }
    else
     throw( NDOException(
      "FlwFiOrcl::SetMaxName: MaxName reduction only is available to 0."
      ) );

 }  // end( FlwFiOrcl::SetMaxName )

/*--------------------------------------------------------------------------*/

void FlwFiOrcl::SetAggregate( bool aggrgt ) {

 if( LsHasChgd )
  delete[] LsHasChgd;

 if( SolvedP )
  delete[] SolvedP;

 //  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( aggrgt ) {
  Aggrgtd = true;                   // flag for aggregate option
  SetNrSubP( 1 );
  LsHasChgd = 0;

  }
 else {
  Aggrgtd = false;
  SetNrSubP( NComm );   // SolvedP[ wFi -1 ] indicate if the value of the primal
                        // solution of wFi component is computed
  if( YIsEasy )         // LsHasChgd[ wFi - 1 ]  indicate if the
	                    // subgradient of wFi component is computed
   LsHasChgd = new bool[ NComm ];
  else
   LsHasChgd = new bool[ OrigXtrCosts ? NComm + 1 : NComm ];
  }

 if( YIsEasy ) {
  SolvedP = new bool[ NComm ];
  VectAssign( SolvedP , false , NComm );
  }
 else {
  SolvedP = new bool[ OrigXtrCosts ? NComm + 1 : NComm ];
  VectAssign( SolvedP , false , OrigXtrCosts ? NComm + 1 : NComm );
  }

 } // end ( FlwFiOrcl::SetAggregate )

/*--------------------------------------------------------------------------*/

void FlwFiOrcl::SetRelax( const char sp4 , bool YiE , cIndex sp1 , cIndex sp2 ,
		        cIndex sp3 )
{ // set the dimension of the function - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 SPar4 = sp4;
 YIsEasy = YiE;
 SPar1 = sp1;
 SPar2 = sp2;
 SPar3 = sp3;

 mutual = (SPar4 == 's') ? 0 : NArcs;

 if( InvULambda )              // delete InvULambda,
  delete[] InvULambda;         // if there exists
 if( ZLambdaCount )            // delete ZLambdaCount, if there exists
  delete[] ZLambdaCount;

 if( LmbCmp ) {
  for( Index i = 0 ; i < NComm ; i++ )
   LmbCmp[ i ].clear();
  delete[] LmbCmp;
  }

 // some exceptions- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( !OrigXtrCosts )
  if( SPar4 == 's' || SPar4 == 'b' || ( SPar1 > 0 ) || YIsEasy )
   throw( NDOException( "FlwFiOrcl::SetICapRelax: there are not fixed-charges" ) );

 // set the dimension of the problem:    - - - - - - - - - - - - - - - - - - -

 if( OrigXtrCosts ) {
  if( SPar4 == 's' || SPar4 == 'b' ) {
   MaxNumVar = mutual + NArcs * NComm;
   if( SPar1 > 0 )
	NumVar = mutual; // only alphas are involved at beginning
   else
	NumVar = MaxNumVar; // both alphas and betas are involved at beginning
   }
  else
   MaxNumVar = NumVar = NArcs; // only alphas are involved from beginning
  }                            // to the end
 else
  MaxNumVar = NumVar = NArcs;  // pure multicommodity

 // allocate dictionary and define   - - - - - - - - - - - - - - - - - - - - -
 // the dictionary for the betas   - - - - - - - - - - - - - - - - - - - - - -

 if( (SPar4 == 's' || SPar4 == 'b') && ( SPar1 > 0 ) ) {

  InvULambda = new Index[ NArcs * NComm ];
  VectAssign( InvULambda , Index( Inf<Index>() ) , ( NArcs * NComm ) );

  LmbCmp = new vector<Index>[ NComm ];

  if( SPar2 > 0 ) {
   ZLambdaCount = new Index[ NArcs * NComm ];        // 0: the variable is not
   VectAssign( ZLambdaCount , Index( 0 ) , NArcs * NComm );  // defined yet
   }
  else
   ZLambdaCount = 0;
  }
 else {
  InvULambda = 0;
  LmbCmp = 0;
  ZLambdaCount = 0;
  }

 if( SGBse1 ) {
  delete[] SGBse1;
  SGBse1 = new Index[ MaxNumVar + 1 ];
  }
 else
  SGBse1 = new Index[ MaxNumVar + 1 ];

 // timer for betas - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( Addt )
  delete[] Addt;
 if( Remt )
  delete[] Remt;

 if( SPar1 > 0 ) {
  Addt  = new OPTtimers();
  if( SPar2 > 0 )
   Remt = new OPTtimers();
  else
   Remt = 0;
  }
 else
  Addt = Remt = 0;

 // allocate memory for "extra" solution- - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( XSolution )
  delete[] XSolution;

 if( OrigXtrCosts && !YiE ) {
  XSolution = new Number[ NArcs ];
  CoefObj = new HpNum[ NArcs ];
  }
 else {
  XSolution = 0;
  CoefObj = 0;
  }

 } // end ( FlwFiOrcl::SetRelax )

/*--------------------------------------------------------------------------*/

void FlwFiOrcl::SetInitialSet( Index_Set Bse )
{
  // !!!! aggiustare anche LmbIndx

 if( SPar4 == 'w'  )
  throw( NDOException( "FlwFiOrcl::SetInitialSet: there are not individual "
		  "capacity constraints" ) );

 Index count = 0;
 if( Bse ) {
  VectAssign( InvULambda , Index( Inf<Index>() ) , ( NArcs * NComm ) );
  for( Index h ; ( h = *(Bse++) ) < Inf<Index>() ; )
   if( h >= mutual ) {
	InvULambda[ count ] = h  - mutual;
	count++;
    }
  }
 else {
  count = NArcs * NComm;
  }

 // update the NumVar    - - - - - - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 NumVar = count + mutual;

 // update ZLambdaCount  - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( SPar2 > 0 ) {
  if( ZLambdaCount )
   delete[] ZLambdaCount;

  ZLambdaCount = new Index[ NArcs * NComm ];
  VectAssign( ZLambdaCount , Index( 0 ) , NArcs * NComm );
  }

 } // end( FlwFiOrcl::SetInitialSet )

/*--------------------------------------------------------------------------*/

void FlwFiOrcl::SetEasy( bool YIsE ) {

 if( !OrigXtrCosts )
  if( YIsE )
   throw( NDOException( "FlwFiOrcl::SetICapRelax: there are not "
		   "fixed-charges" ) );

 // allocate memory for "extra" solution- - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( XSolution )
  delete[] XSolution;

 if( OrigXtrCosts && !YIsE )
  XSolution = new Number[ NArcs ];
 else
  XSolution = 0;

 YIsEasy = YIsE;

 } // end ( FlwFiOrcl::SetEasy )

/*--------------------------------------------------------------------------*/
/*-------------- METHODS FOR READING THE DATA OF THE PROBLEM ---------------*/
/*--------------------------------------------------------------------------*/

Index FlwFiOrcl::GetNumVar( void ) const {

 return( NumVar );

 } // end ( FlwFiOrcl::GetNumVar )

/*--------------------------------------------------------------------------*/

Index FlwFiOrcl::GetMaxNumVar( void ) const {

 return( MaxNumVar );

 } // end ( FlwFiOrcl::GetMaxNumVar )

/*--------------------------------------------------------------------------*/

Index FlwFiOrcl::GetNrFi( void ) const {

 if( YIsEasy )
  return( ( Aggrgtd ? Index( 1 ) : NComm ) +  1 );
 else
  return( Aggrgtd ? Index( 1 ) : NComm  + ( OrigXtrCosts ? 1 : 0 ) );

 } // end ( FlwFiOrcl::GetNrFi )

/*--------------------------------------------------------------------------*/

bool FlwFiOrcl::GetUC( cIndex i ) {

 return( false );

 } // end ( FlwFiOrcl::GetUC )

/*--------------------------------------------------------------------------*/

HpNum FlwFiOrcl::GetGlobalLipschitz( cIndex wFi )
{
 if( wFi == Inf<Index>() )
  return( sqrt( L2Cnst ) );
 else
  return( Inf<HpNum>() );
 } // end ( FlwFiOrcl::GetGlobalLipschitz )

/*--------------------------------------------------------------------------*/

Index FlwFiOrcl::GetBNC( cIndex wFi ) {

 // in the network design the "easy" subproblem is NSubPr+1, whenever
 // the ys are considered as easy components

 if( OrigXtrCosts && ( wFi == NSubPr + 1 ) && YIsEasy )
  return NArcs;
 else
  return (0);

 } // end ( FlwFiOrcl::GetBNC )

/*--------------------------------------------------------------------------*/

Index FlwFiOrcl::GetBNR( cIndex wFi ) {

 // in the network design the "easy" subproblem is NSubPr+1, whenever
 // the ys are considered as easy components

 if( ( ! OrigXtrCosts ) || ( wFi != NSubPr + 1 ) )
  throw( NDOException( "FlwFiOrcl::GetBNR: this component is not easy" ) );
 if ( !YIsEasy )
  throw( NDOException( "FlwFiOrcl::GetBNR: this component is not "
		  "declared easy" ) );

 // on the variables x[ wFi ] are defined only "box" constraints
 return (0);

 }

/*--------------------------------------------------------------------------*/

Index FlwFiOrcl::GetBNZ( cIndex wFi ) {

 // in the network design the "easy" subproblem is NSubPr+1, whenever
 // the ys are considered as easy components

 if( ( ! OrigXtrCosts ) || ( wFi != NSubPr + 1 ) )
  throw( NDOException( "FlwFiOrcl::GetBNZ: this component is not easy" ) );
 if ( !YIsEasy )
  throw( NDOException( "FlwFiOrcl::GetBNZ: this component is not "
		  "declared easy" ) );

 //	Ge1tBNR( wFi ) == 0
 return (0);
 }

/*--------------------------------------------------------------------------*/

void FlwFiOrcl::GetBDesc( cIndex wFi , int *Bbeg , int *Bind , double *Bval ,
   			  double *lhs , double *rhs , double *cst ,
   			  double *lbd , double *ubd ) {

 // in the network design the "easy" subproblem is NSubPr+1, whenever
 // the ys are considered as easy components
 if( ( ! OrigXtrCosts ) || ( wFi != NSubPr + 1 ) )
  throw( NDOException( "FlwFiOrcl::GetBDesc: this component is not easy" ) );
 if ( !YIsEasy )
  throw( NDOException( "FlwFiOrcl::GetBDesc: this component is not "
		  "declared easy" ) );

 // on the variables x[ wFi ] are defined only "box" constraints
 // and consequently B[wFi] = d[wFi] = e[wFi] = []

 if( Bbeg )
  for( Index j= 0; j <= NArcs; j++ )
   Bbeg[ j ] =  0;

 // - - - - - - -  - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 // c[ wFi ]
 if( cst )
  for( Index j=0; j < NArcs; j++)
   cst[j] = - OrigXtrCosts[ j ];

 // l[ wFi ]
 if( lbd )
  for( Index j=0; j < NArcs; j++)
   lbd[j] = 0;

 // u[ wFi ]
 if( ubd )
  for(Index j=0; j < NArcs; j++)
   ubd[j] = 1;

 } // end (FlwFiOrcl::GetBDesc)

/*--------------------------------------------------------------------------*/

Index FlwFiOrcl::GetANZ( cIndex wFi , cIndex strt , Index stp  ) {

 // in the network design the "easy" subproblem is NSubPr+1, whenever
 // the ys are considered as easy components

 if( ( ! OrigXtrCosts ) || ( wFi != NSubPr + 1 ) )
  throw( NDOException( "FlwFiOrcl::GetANZ: this component is not easy" ) );
 if( !YIsEasy )
  throw( NDOException( "FlwFiOrcl::GetANZ: this component is not "
    "declared easy" ) );

 // the call to this function in not correct
 if( stp <= strt )
  throw( NDOException( "FlwFiOrcl::GetANZ: stop <= start" ) );

 if( stp > NumVar )
   stp = NumVar;

 if( SPar4 == 's' || SPar4 == 'b'  )               // strong formulation
  if( SPar1 == 0 )             // static case: no generation of constraints
   if(  strt != 0  || stp < NumVar )
	throw( NDOException( "FlwFiOrcl::GetANZ: no variables additional"
	  " interval" ) );
   else                        // the betas are present at beginning
	return NArcs * NComm + mutual ;
  else                         // dynamic case: generation of constrains
   if( strt != 0 && ( strt < mutual ) )
	throw( NDOException( "FlwFiOrcl::GetANZ: this should not happen" ) );
   else                               // or first call...
    return( stp - strt );             // only the betas are involved
 else {                        // weak formulation
  if( strt != 0  || stp < NumVar )
   throw( NDOException( "FlwFiOrcl::GetANZ: this should not happen" ) );
  return NArcs;                        // only the alphas are involved
  }

 } // end (FlwFiOrcl::GetANZ)

/*--------------------------------------------------------------------------*/

void FlwFiOrcl::GetADesc( cIndex wFi , int *Abeg , int *Aind , double *Aval ,
  			  cIndex strt , Index stp ) {

 // in the network design the "easy" subproblem is NSubPr+1, whenever
 // the ys are considered as easy components
 if( ( ! OrigXtrCosts ) || ( wFi != NSubPr + 1 ) )
  throw( NDOException( "FlwFiOrcl::GetADesc: this component is not easy" ) );
 if( !YIsEasy )
  throw( NDOException( "FlwFiOrcl::GetADesc: this component is not "
	 "declared easy" ) );

 // the call to this function in not correct
 if( stp <= strt )
  throw( NDOException( "FlwFiOrcl::GetADesc: stop <= start" ) );

 if( stp > NumVar )
  stp = NumVar;

 Index Elem;
 if( SPar4 == 's' || SPar4 == 'b' )                // strong formulation
  if( SPar1 == 0  )             // static case: no generation of constraints
   if(  strt != 0  || stp < MaxNumVar   )
	throw( NDOException( "FlwFiOrcl::GetADesc: no variables additional "
	     "interval" ) );
   else {                       // the betas are present at beginning
    for( Index i = 0; i < NArcs; i++ ) {
     Elem = i * NComm;
     if( SPar4 == 'b' ) {
      Elem += i;
      Abeg[ i ] = Elem;
      Aind[ Elem ] = i;
      Aval[ Elem++ ] = - OrigTotCap[ i ];
      }
     else
      Abeg[ i ] = Elem;
     for ( Index k = 0; k < NComm; k++) {
      Aind[ Elem + k ] = k * NArcs + mutual + i;
      Aval[ Elem + k ] = - OrigCapacities[ k * NArcs + i ];
      }
     }
    Abeg[ NArcs ] = NArcs * NComm + mutual;
    }
  else                          // dynamic case: generation of constrains
   if( strt != 0 && ( strt < mutual ) )
	throw( NDOException( "FlwFiOrcl::GetADesc: this should not happen" ) );
   else
	if( strt == 0 && SPar4 != 's') {           // first call...
	 int count = 0;
	 for( Index i = 0; i < NArcs; i++ ) { // examine each column of
	  Abeg[ i ] = count;                  // the matrix A
	  Aind[ count ] = i;
	  Aval[ count++ ] = - OrigTotCap[ i ];
	  if( stp > NArcs )
	   for( Index j = NArcs; j < stp ;  j++ )
		if( ( ( InvULambda[ j - NArcs ] % NArcs ) == i ) &&
			( j >= NArcs && j < stp ) ) {
		 Aind[ count ] = j;
		 Aval[ count++ ] = -OrigCapacities[ InvULambda[ j - NArcs ] ];
		 }
	  }
	 Abeg[ NArcs ] = count;
	 }
	else {                       // only the betas are involved
	 int count = 0;
   	 for( Index i = 0; i < NArcs; i++ ) { // examine each column of
   	  Abeg[ i ] = count;                  // the matrix A
   	  for( Index j = mutual; j < NumVar;  j++ )
   	   if( ( ( InvULambda[ j - mutual ] % NArcs ) == i ) &&
   		   ( j >= strt && j < stp ) ) {
   	    Aind[ count ] = j;
   	    Aval[ count++ ] = -OrigCapacities[ InvULambda[ j - mutual ] ];
   	    }
   	  }
   	 Abeg[ NArcs ] = count;
     }
 else {                        // weak formulation
  if( strt != 0  || stp < NumVar )
   throw( NDOException( "FlwFiOrcl::GetANZ: this should not happen" ) );
  for ( Index i = 0; i < NArcs; i++ ) {  // only the alphas are involved
   Abeg[ i ] = i;
   Aind[ i ] = i;
   Aval[ i ] = - OrigTotCap[ i ];
   }
  Abeg[ NArcs ] = NArcs;
  }

 } // end (FlwFiOrcl::GetADesc)

/*--------------------------------------------------------------------------*/
/*----------------------- METHODS FOR CHANGING DATA ------------------------*/
/*--------------------------------------------------------------------------*/

void FlwFiOrcl::SetLambda( cLMRow Lmbd ) {

 // update Lambda  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 FiOracle::SetLambda( Lmbd );

 if( !Aggrgtd )
  if( YIsEasy ) {
   VectAssign( SolvedP , false , NComm );
   VectAssign( LsHasChgd , true, NComm );
   }
  else {
   VectAssign( SolvedP , false , ( OrigXtrCosts ? NComm + 1 : NComm ) );
   VectAssign( LsHasChgd , true, ( OrigXtrCosts ? NComm + 1 : NComm ) );
   }
 else {
  if( YIsEasy )
   VectAssign( SolvedP , false , NComm );
  else
   VectAssign( SolvedP , false , ( OrigXtrCosts ? NComm + 1 : NComm ) );
  }

 // update ZLambdaCount   - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( LamBase ) {
  if( ZLambdaCount )
   for( Index j = mutual ; j < NumVar ; j++ )
    if( Lambda[ j ] <= EpsFi )
     ZLambdaCount[ InvULambda[ LamBase[ j ] - mutual ] ]++;
    else
     ZLambdaCount[ InvULambda[ LamBase[ j ] - mutual ] ] = 0;
  }
 else {
  if( ZLambdaCount )
   for( Index j = mutual ; j < NumVar ; j++ )
    if( Lambda[ j ] <= EpsFi )
     ZLambdaCount[ InvULambda[ j - mutual ] ]++;
    else
     ZLambdaCount[ InvULambda[ j - mutual ] ] = 0;
  }

 } // end ( FlwFiOrcl::SetLambda )

/*--------------------------------------------------------------------------*/

bool FlwFiOrcl::SetPrecision( HpNum Eps ) {
  // We  are assuming that the FiOracle is capable of computing
  // the solution with any precision

 EpsFi = Eps;
 return false;

 } // end ( FlwFiOrcl::SetPrecision )

/*--------------------------------------------------------------------------*/
/*-------------------- METHODS FOR SOLVING THE PROBLEM ---------------------*/
/*--------------------------------------------------------------------------*/

HpNum FlwFiOrcl::Fi( cIndex wFi ) {

 // some exceptions- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( YIsEasy && ( wFi == NSubPr + 1) )
  throw( NDOException( "FlwFiOrcl::Fi: this component is easy" ) );

 // "easy" problem exists if fixed charges have been introduced

 // cells for containing the function value - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 HpNum zrVal = 0;       // value of the linear ( affine ) 0-th component
 HpNum FlwVal = 0;      // value of min-cost flow problem
 HpNum FxChrgVal = 0;   // value of fixed charge component

 Index_Set FSolBse = 0;
 HpNum out = 0;         // output variable

 // timer on - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( Fit )
  Fit->Start();

 // value of the linear ( affine ) 0-th component of Fi- - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( !OrigXtrCosts )
  if( (wFi == 0) || (wFi == Inf<Index>())  ) {
   if( LamBase )
	for( Index j = 0; LamBase[ j ] < NArcs;  j++ )
	 zrVal -= OrigTotCap[ LamBase[ j ] ] * Lambda[ j ];
   else
	for( Index j = 0; j < NArcs;  j++ )
     zrVal -= OrigTotCap[ j ] * Lambda[ j ];
   }

 // if( OrigXtrCosts ) ---> zrVal == 0

 // value of fixed-charge component- - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( !YIsEasy && OrigXtrCosts && ( wFi >=  GetNrFi() ) ) {

  SolveLagrangian( NComm );

  // write the solution - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  for( Index j = 0; j < NArcs; j++ )
   if( CoefObj[ j ] <= 0 )
	FxChrgVal += CoefObj[ j ];

  } // end( !YIsEasy && OrigXtrCosts && ( wFi >=  GetNrFi() ) )

 // value of the min-cost flow problem - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( ( wFi != 0 ) && !( OrigXtrCosts && !Aggrgtd && ( wFi == NComm + 1 ) ) ) {

  if( Aggrgtd )
   SolveLagrangian( );
  else
   if( ( wFi >= 1 ) && ( wFi <= NComm ) )
    SolveLagrangian( wFi - 1 );
   else
	SolveLagrangian( );

  FlwVal = GetPVal();
  } // end Flow computation  - - - - - - - - - - - - - - - - - - - - - - - - -

 // timer off- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( Fit )
  Fit->Stop();

 // return the value - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


 if( wFi == 0 )                          // requires the value of the linear
  out = -zrVal;                          // 0-th component of Fi
                                         // zrVal == 0 < --- OrigXtrCosts

 if( wFi >= 1 && wFi <= NSubPr )         // requires the value of the wFi-th component
  out = - FxChrgVal - FlwVal;            // of Fi which must not be an "easy" one
                                         // FxChrgVal is not considered <--- YIsEasy
                                         // FxChrgVal == 0 <--- !OrigXtrCosts
                                         // FxChrgVal == 0 <--- !Aggrgtd

 if( wFi == NSubPr + 1 ) {               // requires the value of the (NComm + 1)-th
  if( OrigXtrCosts )                     // component of Fi
   out = - FxChrgVal;
  else
   out = - FlwVal;
  }

 if( wFi > NSubPr + 1  )                  // wFi < Inf<Index>() :
  out = - zrVal - FxChrgVal - FlwVal;     // requires the value of the full function Fi
                                          // except that of the "easy" components and of
                                          // the linear part
                                          // wFi == Inf<Index>() :
                                          // requires the value of the full function
  	                                      // Fi except that of the "easy" components

 return out;

 } // end ( FlwFiOrcl::Fi )


/*--------------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/


bool FlwFiOrcl::GetPSol( void ) {

 if( XSol ) {
  XSolA( XSol );
  XSol = 0;
  }
 else
  if( FSol ) {
   FSolA( FSol );
   FSol = 0;
   }
  else
   throw( NDOException( "FlwFiOrcl::GetPSol: FSol & XSol is 0" ) );

 return( false );

 } // end ( FlwFiOrcl::GetPSol )

/*--------------------------------------------------------------------------*/

bool FlwFiOrcl::NewGi( cIndex wFi ) {

 bool FreeWFi = true;   // indicate if it is possible produce a "new" item
                        // for the point Lambda

 LastGiName = Inf<Index>();

 Index_Set FSolBse = 0;
 // some exceptions- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( YIsEasy && ( wFi == NComm + 1) )
  throw( NDOException( "FlwFiOrcl::NewGi: this component is easy" ) );

 if( wFi > GetNrFi() && wFi < Inf<Index>() )
  throw( NDOException( "FlwFiOrcl::NewGi: this component is not declared" ) );

 if( wFi == 0 )
  return false;

 // asks for an aggregated (epsilon-) subgradient  - - - - - - - - - - - - - -
 //  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( Aggrgtd )
  FreeWFi = LHasChgd;
 else {
  if( YIsEasy ) {
   if( wFi <= NComm )
    FreeWFi = LsHasChgd[ wFi - 1 ];
   else
    for( Index i = 0; i < NComm; i++ )
     if( !LsHasChgd[ i ] ) {
      FreeWFi = false;
      break;
      }
   }
  else {
   if( wFi <= GetNrFi() )
    FreeWFi = LsHasChgd[ wFi - 1 ];
   else
    for( Index i = 0; i < GetNrFi(); i++ )
     if( !LsHasChgd[ i ] ) {
      FreeWFi = false;
      break;
      }
   }
  }

 // get the optimal solution - - - - - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( FreeWFi ) {

  if( Fit )
   Fit->Start();

  if( Aggrgtd ) {  // aggregate case   - - - - - - - - - - - - - - - - - - - -

   // solve the Flow part  - - - - - - - - - - - - - - - - - - - - - - - - - -

   SolveLagrangian( );
   SetFlwSol( FSolution );

   // solve the Y part     - - - - - - - - - - - - - - - - - - - - - - - - - -

   if( !YIsEasy && OrigXtrCosts ) {
	SolveLagrangian( NComm );
	for( Index j = 0; j < NArcs; j++ )
     if( CoefObj[ j ] <= 0 )
      XSolution[ j ] = 1;
     else
  	  XSolution[ j ] = 0;
    }

   } // end aggregate case   - - - - - - - - - - - - - - - - - - - - - - - - -
  else
   if( ( wFi >= 1 ) && ( wFi <= NComm ) ) {  // pick just a component

    if( SolWFi > NComm ) {
     delete[] FSolution;
     FSolution = new FNumber[ NArcs ];
     }

    SolveLagrangian( wFi - 1 );
    SetFlwSol( FSolution , FSolBse , wFi - 1 ); // pass to the object pointer to the memory
    }
   else {
    if( SolWFi <=  NComm ) {  // get the all solution  - - - - - - - - - - - -
     delete[] FSolution;
     FSolution = new FNumber[ NArcs * NComm ];
     }

    // solve the Y part   - - - - - - - - - - - - - - - - - - - - - - - - - - -

    if( !YIsEasy && OrigXtrCosts && ( wFi >=  GetNrFi() ) ) {

     SolveLagrangian( NComm );

     // write the solution  - - - - - - - - - - - - - - - - - - - - - - - - - -
     // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

     for( Index j = 0; j < NArcs; j++ )
      if( CoefObj[ j ] <= 0 )
       XSolution[ j ] = 1;
      else
   	   XSolution[ j ] = 0;

     } // end( !YIsEasy && OrigXtrCosts && ( wFi >=  GetNrFi() ) )

    // solve the flow part  - - - - - - - - - - - - - - - - - - - - - - - - - - -

    if( !OrigXtrCosts || !( wFi == NComm + 1 ) ) {
     SolveLagrangian( );
     SetFlwSol( FSolution ); // pass to the object pointer to the memory
     } // end ( if )

    }

  // get the solution  - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( !OrigXtrCosts || !( wFi == NComm + 1 ) )
   MMCFFlwBase::GetPSol(); // where there is the flow solution

  // tells to witch subproblem belongs last gi created - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  SolWFi = wFi;

  if( Fit )
   Fit->Stop();

  } // end if( FreeWFi )

 return FreeWFi;

 } // end ( FlwFiOrcl::NewGi )

/*--------------------------------------------------------------------------*/

void FlwFiOrcl::SetGiName( cIndex Name )
{   // some exceptions   - - - - - - - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( Name > MaxName )
  throw( NDOException( "FlwFiOrcl::SetGiName:Name is higher than MaxName." ) );

 if( SPar1 ) { // copy the extra solution    - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( OrigXtrCosts && !YIsEasy ) {
   if( OldXSols[ Name ] )
    if( ( SolWFi >=  NComm + 1 ) ||  Aggrgtd  )
	 Swap( OldXSols[ Name] , XSolution );
    else {
	 delete[] OldXSols[ Name ];
	 OldXSols[ Name ] = 0;
     }
   else
    if( ( SolWFi >=  NComm + 1 ) || Aggrgtd  ) {
	 OldXSols[ Name ] = XSolution;
     XSolution = new Number[ NArcs ];
     }
   }

  // copy the current flow solution in OldFSol  and the name of the item   - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( OldFSols[ Name] )  // OldFSols[ Name] has been already allocated
   if( !Aggrgtd && OrigXtrCosts && ( SolWFi == NComm + 1 ) ) {
    delete[] OldFSols[ Name];
    OldFSols[ Name] = 0;
    OldWFi[ Name ] = SolWFi;
    }
   else {
    Swap( OldFSols[ Name] , FSolution );
    Swap( SolWFi , OldWFi[ Name ] );
    }
  else // pass the FSolution's pointer to OldFSols[ Name]
   if( !Aggrgtd && OrigXtrCosts && ( SolWFi == NComm + 1 ) )
    OldWFi[ Name ] = SolWFi;
   else {
    OldFSols[ Name] = FSolution;
    OldWFi[ Name ] = SolWFi;
    if( Aggrgtd || SolWFi > NComm )
     FSolution = new FNumber[ NArcs * NComm ];
    else
     FSolution = new FNumber[ NArcs ];
    }

  LastGiName = Name;
  }

 } // end ( FlwFiOrcl::SetGiName )

/*--------------------------------------------------------------------------*/

Index FlwFiOrcl::GetGi( SgRow SubG , cIndex_Set &SGBse , cIndex Name,
			cIndex strt , Index stp ) {

 // GetMaxNZ( wFi ) return MaxNumVar, assuming that the maximum number of
 //	nonzeroes is not known in advance

 Index k = 0;  // the number of nonzero of the item

 // some exceptions  - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 // in NewGi are present the following exception:
 // YIsEasy && ( wFi == NSubPr + 1)
 // wFi > GetNrFi() && wFi < Inf<Index>()
 // wFi == 0

 if( stp <= strt )
  throw( NDOException( "FlwFiOrcl::GetGi: stop <= start" ) );

 if( Name < MaxName ) {
  if( ! OldWFi[ Name ] )
   throw( NDOException( "FlwFiOrcl::GetGi: Subgradient does not exist" ) );
  }
 else
  if( Name > MaxName ) {
   if(  strt != 0  || stp < NumVar )
	throw( NDOException( "FlwFiOrcl::GetGi: this should not happen" ) );
   }
  else
   if(  strt > 0  && strt < mutual )
	throw( NDOException( "FlwFiOrcl::GetGi: this should not happen" ) );

 if( stp > NumVar )
  stp = NumVar;

 // SetGiName could call GetGi at the current point several times. At the first
 // time that SetGiName saves the item a name will be assign. Thus, recover the
 // item by using this information - - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( Fit )
  Fit->Start();

 Index actualname;
 if( ( LastGiName < Inf<Index>() ) && ( Name > MaxName ) )
  actualname = LastGiName;
 else
  actualname = Name;

 // computation part - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( actualname == MaxName ) {

  // actualname == MaxName. The required information is about the constant
  // subgradient of the linear 0-th component of Fi  - - - - - - - - - - - - -

  if( OrigXtrCosts ) { // a fixed-cost problem has a all-zero RHS
   *SGBse1 = Inf<Index>();
   SGBse = SGBse1;
   }
  else {  // "pure" MMCF: RHS are the arc capacities
   for( Index i = strt ; i < stp ; i++ )
	if( i < NArcs )
	 SubG[ i ] = OrigTotCap[ i ];
	else
	 SubG[ i - strt ] = 0;
   SGBse = 0;
   k = stp - strt;
   }

  }
 else {

  if( actualname > MaxName ) {

   // asks for an aggregated (epsilon-) subgradient excluding the constant
   // part and all the easy components   - - - - - - - - - - - - - - - - - - -
   //  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   if( SolWFi > GetNrFi() || Aggrgtd ) {

    // subgradient part relative to alphas   - - - - - - - - - - - - - - - - -

	if( SPar4 != 's' )  {
	 if( OrigXtrCosts && !YIsEasy )
	  for( Index j = 0; j < NArcs; j++ )
	   SubG[ j ] = OrigTotCap[ j ] * XSolution[ j ];
	 else
	  for( Index j = 0; j < NArcs; j++ )
	   SubG[ j ] = 0;

	 for( Index j = 0; j < ( NArcs * NComm ) ; j++ )
	  SubG[ j % NArcs ] -= FSolution[ j ];
	 }

    // subgradient part relative to betas  - - - - - - - - - - - - - - - - - -

	if( SPar4 == 's' || SPar4 == 'b' ) {
     if( !SPar1 )
      if( YIsEasy )
       for( Index j = mutual ; j < MaxNumVar ; j++ )
    	SubG[ j ] = - FSolution[ j - mutual ];
      else
       for( Index j = mutual ; j < MaxNumVar ; j++ )
    	SubG[ j ] = - FSolution[ j - mutual ] +
    	  OrigCapacities[ j - mutual ] * XSolution[  j % NArcs ];
     else
      if( stp > mutual )  {
       if( YIsEasy )
        for( Index j = mutual ; j < stp ; j++ )
         SubG[ j ] = - FSolution[ InvULambda[ j - mutual ] ];
       else
        for( Index j = mutual ; j < stp ; j++ )
         SubG[ j ] = - FSolution[ InvULambda[ j - mutual ] ] +
          OrigCapacities[ InvULambda[ j - mutual ] ] *
          XSolution[ InvULambda[ j - mutual ] % NArcs ];
       }
     }

    // requires the aggregated (epsilon-)subgradient of (the not-easy part of)
    // Fi  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    if( SolWFi == Inf<Index>() && !OrigXtrCosts )
    for( Index j = 0 ;j < NArcs; j++ )
     SubG[ j ] += OrigTotCap[ j ];

    } // end ( if: LastWFi >  GetNrFi() )

   else {

	// asks for an (epsilon-)subgradient of Fi[ wFi ] which must not be an
	// "easy" one  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	if( SolWFi <= NComm ) {

     // subgradient part relative to alphas    - - - - - - - - - - - - - - - -

	 if( SPar4 != 's' )
	  for( Index j = 0; j < NArcs; j++ )
	   SubG[ j ] = -FSolution[ j ];

	 // subgradient part relative to betas   - - - - - - - - - - - - - - - - -

	 if( SPar4 == 's' || SPar4 == 'b' ) {
	  if( !SPar1 )
	   for( Index j = mutual ; j < MaxNumVar ; j++ )
		if( ( SolWFi - 1 ) == ( ( j - mutual ) / NArcs ) )
		 SubG[ j ] = - FSolution[ ( j - mutual ) % NArcs ];
		else
	     SubG[ j ] = 0;
	  else
	   if( stp > mutual )  {
		for( Index j = mutual ; j < stp ; j++ )
		 if( ( SolWFi - 1 )  ==  ( InvULambda[ j - mutual ] / NArcs ) )
		  SubG[ j ] = - FSolution[ InvULambda[ j - mutual ] % NArcs ];
		 else
          SubG[ j ] = 0;
	    }
      }

	 }
   	else {

	 // asks for an (epsilon-)subgradient of Fi[ NComm + 1 ] which must not be
	 // an "easy" one: LastWFi == NSubPr + 1 - - - - - - - - - - - - - - - - -
	 //  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   	 // subgradient part relative to alphas    - - - - - - - - - - - - - - - -

   	 if( SPar4 != 's' )
   	  for( Index j = 0 ;j < NArcs; j++ )
       SubG[ j ] =  OrigTotCap[ j ] * XSolution[ j ];

   	 // subgradient part relative to betas   - - - - - - - - - - - - - - - - -

   	 if( SPar4 == 's' || SPar4 == 'b' ) {
   	  if( !SPar1 )
   	   for( Index j = mutual ; j < stp ; j++ )
   		SubG[ j ] = OrigCapacities[  j - mutual ] * XSolution[  j % NArcs ];
   	  else
   	   for( Index j = mutual ; j < stp ; j++ )
   		SubG[ j ] = OrigCapacities[ InvULambda[ j - mutual ] ] *
   		    XSolution[ InvULambda[ j - mutual ] % NArcs ];
   	  }

     }

    } // end ( else: LastWFi <  GetNrFi() )

   } // end ( if:  actualname > MaxName )
  else {

   // asks for an aggregated (epsilon-) subgradient excluding the constant
   // part and all the easy components   - - - - - - - - - - - - - - - - - - -
   //  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   if( OldWFi[ actualname ] >  GetNrFi() ||  Aggrgtd ) {

	// subgradient part relative to alphas   - - - - - - - - - - - - - - - - -
	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    if( ( SPar4 != 's' ) && ( strt < NArcs ) )  {
	 if( OrigXtrCosts && !YIsEasy )
	  for( Index j = strt ; j < min( stp , mutual ) ; j++ )
	   SubG[ j - strt ] = OrigTotCap[ j ] * OldXSols[ actualname ][ j ];
	 else
	  for( Index j = strt ; j < min( stp , mutual ) ; j++ )
	   SubG[ j - strt ] = 0;

	 for( Index j = 0 ; j < ( NArcs * NComm ) ; j++ )
	  if( ( ( j % NArcs ) >= strt ) && ( ( j % NArcs ) < stp ) )
	   SubG[ ( j % NArcs ) - strt ] -= OldFSols[ actualname ][ j ];
	 }

    // part relative to beta - - - - - - - - - - - - - - - - - - - - - - - - -
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    if( SPar4 == 's' || SPar4 == 'b' ) {
     if( !SPar1 )
	  if( YIsEasy )
	   for( Index j = max( strt , mutual ) ; j < stp ; j++ )
	    SubG[ j - strt ] = - OldFSols[ actualname ][ j - mutual ];
	  else
	   for( Index j = max( strt , mutual ) ; j < stp ; j++ )
        SubG[ j - strt ] = - OldFSols[ actualname ][ j - mutual ] +
            OrigCapacities[ j - mutual ] * OldXSols[ actualname ][  j % NArcs ];
	 else
	  if( YIsEasy )
       for( Index j = max( strt , mutual ) ; j < stp ; j++ )
        SubG[ j - strt ] = - OldFSols[ actualname ][ InvULambda[ j - mutual ] ];
	  else
       for( Index j = max( strt , mutual ) ; j < stp ; j++ )
        SubG[ j - strt ] = - OldFSols[ actualname ][ InvULambda[ j - mutual ] ] +
		    OrigCapacities[ InvULambda[ j - mutual ] ] *
		    OldXSols[ actualname ][ InvULambda[ j - mutual ] % NArcs ];
     }

    // requires the aggregated (epsilon-)subgradient of (the not-easy part of)
    // Fi  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    if( ( OldWFi[ actualname ] == Inf<Index>() ) && !OrigXtrCosts )
     for( Index j = strt ;j < stp; j++ )
      SubG[ j - strt ] += OrigTotCap[ j ];

    } // end ( OldWFi[ actualname ] >  GetNrFi() )

   else {

	// asks for an (epsilon-)subgradient of Fi[ wFi ] which must not be an
	// "easy" one  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    if( OldWFi[ actualname ] <= NComm ) {

     // subgradient part relative to alphas  - - - - - - - - - - - - - - - - - -
     //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

     if( ( SPar4 != 's' ) && ( strt < NArcs ) )
      for( Index j = strt ; j < min( stp , mutual ) ; j++ )
       SubG[ j - strt ] = -OldFSols[ actualname ][ j ];

   	 // subgradient part relative to betas   - - - - - - - - - - - - - - - - - -
     //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   	 if( SPar4 == 's' || SPar4 == 'b' ) {
   	  if( !SPar1 )
   	   for( Index j = max( strt , mutual ) ; j < stp ; j++ )
   		if( ( OldWFi[ actualname ] - 1 ) == ( ( j - mutual ) / NArcs ) )
   		 SubG[ j - strt ] = - OldFSols[ actualname ][ ( j - mutual ) % NArcs ];
   		else
   	     SubG[ j - strt ] = 0;
   	  else
   	   for( Index j = max( strt , mutual ) ; j < stp ; j++ )
   		if( ( InvULambda[ j - mutual ] / NArcs ) == ( OldWFi[ actualname ] - 1 ) )
   		 SubG[ j - strt ] = - OldFSols[ actualname ][ InvULambda[ j - mutual ] % NArcs ];
   		else
         SubG[ j - strt ] = 0;
      }

     }
	else {

     // asks for an (epsilon-)subgradient of Fi[ NComm + 1 ] which must not be an
	 // "easy" one: OldWFi[ actualname ] == NSubPr + 1 - - - - - - - - - - - - - - - - -
	 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	 // subgradient part relative to alphas    - - - - - - - - - - - - - - - - - -
     //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	 if( ( SPar4 != 's' ) &&  ( strt < NArcs ) )
	  for( Index j = strt ;j < min( stp , mutual ) ; j++ )
	   SubG[ j - strt ] =  OrigTotCap[ j ] * OldXSols[ actualname ][ j ];

	 // part relative to beta  - - - - - - - - - - - - - - - - - - - - - - - - - -
	 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

     if( SPar4 == 's' || SPar4 == 'b' ) {
      if( !SPar1 )
       for( Index j = max( strt , mutual ) ; j < stp ; j++ )
        SubG[ j - strt ] = OrigCapacities[  j - mutual ] * OldXSols[ actualname ][  j % NArcs ];
      else
       for( Index j = max( strt , mutual ) ; j < stp ; j++ )
        SubG[ j - strt ] = OrigCapacities[ InvULambda[ j - mutual ] ] *
            OldXSols[ actualname ][ InvULambda[ j - mutual ] % NArcs ];
      }

	 }

    } // end ( else: OldWFi[ actualname ] <  GetNrFi() )

   }  // end ( else: actualname < MaxName )

  #if SPF_SUB

   // turn SubG from a "dense" NumVar-vector to a "sparse" one  - - - - - - - -

   Index_Set SGBse2 = Sparsify( SubG , SGBse1 , ( stp - strt ) , strt );
   *SGBse2 = Inf<Index>();
   k = SGBse2 - SGBse1;

   SGBse = SGBse1;

  #else
   k = stp - strt;
   SGBse = 0;
  #endif



  } // end (  Name != MaxName ) )

 if( Fit )
  Fit->Stop();

 // impose that a new call to item is not possible - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( Aggrgtd )
  LHasChgd = false;
 else {
  if( YIsEasy ) {
   if( SolWFi <= NComm )
    LsHasChgd[ SolWFi - 1 ] = false;
   else
    for( Index i = 0; i < NComm ; i++ )
   LsHasChgd[ i ] = false;
   }
  else
   if( SolWFi <= GetNrFi() )
    LsHasChgd[ SolWFi - 1 ] = false;
   else
    for( Index i = 0; i < GetNrFi(); i++ )
     LsHasChgd[ i ] = false;
   }

 return( k );

 } // end ( FlwFiOrcl::GetGi )

/*--------------------------------------------------------------------------*/

HpNum FlwFiOrcl::GetVal( cIndex Name )
{
 HpNum LinErr  = 0;

 if( Name < MaxName ) {
  throw( NDOException( "FlwFiOrcl::SetGiName:aggiustare codice...." ) );

  // ********************* ELIMINARE FSolution e XSolution*************

 // compute the costs for each flow variable   - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 CRow NwCsts = new CNumber[ NArcs * NComm ];
 LMRow tempL = new LMNum[ NumVar ];

 if( LamBase ) {
  VectAssign( tempL , Lambda , LamBDim );
  Densify( tempL , LamBase , LamBDim , NumVar );
  }
 else
  VectAssign( tempL , Lambda , NumVar );

 for( Index j = 0; j < NArcs * NComm; j++ ) {
  NwCsts[ j ] = OrigCosts[ j ];
  if( SPar4 != 's' )
   NwCsts[ j ] += tempL[ j % NArcs ];
   }

 if( SPar4 == 's' || SPar4 == 'b' ) {
  if( SPar1 > 0 ) {
   for( Index j = mutual ; j < NumVar ; j++ )
	NwCsts[ InvULambda[ j - mutual ] ] += tempL[ j ];
    }
  else
   for( Index j = mutual ; j < MaxNumVar ; j++ )
    NwCsts[ j - mutual ] += tempL[ j ];
    }

  // computation of the linearization errors  - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( OldWFi[ Name ] >  GetNrFi() ||  Aggrgtd ) { // asks for the linearization
	                           //  error of all Fi  - - - - - - - - - - - - - -
   FRow tFSols = new FNumber[ NArcs * NComm ];
   VectDiff( tFSols , OldFSols[ Name ] , FSolution , NArcs * NComm );
   LinErr = ScalarProduct( NwCsts , tFSols , NArcs * NComm ) ;
   delete[] tFSols;

   if( !YIsEasy && OrigXtrCosts ) { // compute linearization error of the y
	     // component   - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	CRow NwExtCsts = new CNumber[ NArcs ];
    Row tXSol = new Number[ NArcs ];

	// compute extra costs  - - - - - - - - - - - - - - - - - - - - - - - - -
    //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    VectAssign( NwExtCsts , OrigXtrCosts , NArcs );

    if( SPar4 != 's' )
	 for( Index j = 0; j < NArcs; j++ )
	  NwExtCsts[ j ] -= OrigTotCap[ j ] * tempL[ j ];

	if( SPar4 == 's' || SPar4 == 'b' ) {
	 if( SPar1 > 0 )
	  for( Index j = mutual; j < NumVar; j++ )
       NwExtCsts[ InvULambda[ j - mutual ] % NArcs ] -=
	   	   OrigCapacities[ InvULambda[ j - mutual ] ] * tempL[ j ];
	 else
	  for( Index j = mutual; j < MaxNumVar; j++ )
       NwExtCsts[ j % NArcs ] -= OrigCapacities[ j - mutual ] * tempL[ j ];
	 }

	// compute the error  - - - - - - - - - - - - - - - - - - - - - - - - - -
	//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    VectDiff( tXSol , OldXSols[ Name ] , XSolution , NArcs );
    LinErr += ScalarProduct( NwExtCsts , tXSol , NArcs ) ;

    delete[] tXSol;
    delete[] NwExtCsts;
    } // end y component  - - - - -  - - - - - - - - - - - - - - - - - - - - -

   }
  else { // disaggregated case - - - - - - - - - - - - - - - - - - - - - - - -

   if( OldWFi[ Name ] <= NComm ) { // ask for the linearization error of a
	// flow component  - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	FRow tFSols = new FNumber[ NArcs ];
	VectDiff( tFSols , OldFSols[ Name ] , FSolution , NArcs );
	LinErr = ScalarProduct( NwCsts + NArcs * ( OldWFi[ Name ] - 1 ) ,
			tFSols , NArcs ) ;
	delete[] tFSols;
    }
   else {  // ask for the linearization error of the y component    - - - - -
	//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	CRow NwExtCsts = new CNumber[ NArcs ];
	Row tXSol = new Number[ NArcs ];

	// compute extra costs  - - - - - - - - - - - - - - - - - - - - - - - - -
	//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	VectAssign( NwExtCsts , OrigXtrCosts , NArcs );

	if( SPar4 != 's' )
	 for( Index j = 0; j < NArcs; j++ )
	  NwExtCsts[ j ] -= OrigTotCap[ j ] * tempL[ j ];

	if( SPar4 == 's' || SPar4 == 'b' ) {
	 if( SPar1 > 0 )
	  for( Index j = mutual; j < NumVar; j++ )
	   NwExtCsts[ InvULambda[ j - mutual ] % NArcs ] -=
	   	   OrigCapacities[ InvULambda[ j - mutual ] ] * tempL[ j ];
     else
	  for( Index j = mutual; j < MaxNumVar; j++ )
	   NwExtCsts[ j % NArcs ] -= OrigCapacities[ j - mutual ] * tempL[ j ];
	 }

	// compute the error  - - - - - - - - - - - - - - - - - - - - - - - - - -
	//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	VectDiff( tXSol , OldXSols[ Name ] , XSolution , NArcs );
	LinErr = ScalarProduct( NwExtCsts , tXSol , NArcs ) ;
	delete[] tXSol;
	delete[] NwExtCsts;
    }

   } // end disaggregated case - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  delete[] NwCsts;
  delete[] tempL;

  } // end computation of linearization errors of previous ietms - - - - - - -

 return( LinErr );

 } // end ( FlwFiOrcl::GetVal )

/*--------------------------------------------------------------------------*/
/*-------------------- METHODS FOR READING OTHER RESULTS -------------------*/
/*--------------------------------------------------------------------------*/

HpNum FlwFiOrcl::GetLowerBound( cIndex wFi ) {

 #if LowB_FI

  if( wFi == 0 )
   return( - Inf<HpNum>() );

  if( wFi > GetNrFi() )
   return( LowerBound );
  else
   return( - Inf<HpNum>() );
 #else
  return( - Inf<HpNum>() );
 #endif

 } // end ( FlwFiOrcl::GetLowerBound )

/*--------------------------------------------------------------------------*/

void FlwFiOrcl::SetLowerBound( HpNum lowB ) {

 LowerBound = lowB;

 } // end ( FlwFiOrcl::SetLowerBound )

/*--------------------------------------------------------------------------*/

FiOracle::FiStatus FlwFiOrcl::GetFiStatus( Index wFi )
{
 if( wFi <= GetNrFi() )
  return kFiNorm;

 // nothing to do - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( ! ( SPar1 && ( SPar4 == 's' || SPar4 == 'b' ) ) )
  return( kFiNorm );

 //- - - - - - - -  - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( ! Slvr )
  throw( NDOException( "FlwFiOrcl::GetFiStatus: called by whom??" ) );

 FiStatus status = kFiNorm;
 const Index nrIters = Slvr->NrIter();

 bool add = false;     // boolean for adding
 bool rem = false;     // boolean for removing

 if( MaxName == 0 )
  throw( NDOException("To use the addition of variables the old solutions.\
	 must be saved Change input parameters accordingly."));

 if( ( ( ( nrIters != LstChgItr ) &&  ( ! ( nrIters % SPar1 ) ) )
   || Slvr->IsOptimal() ) && nrIters ) {
  Addt->Start();
  add = checkAddition();
  Addt->Stop();
  }

 if( SPar2 && ( nrIters != LstChgItr ) && ( ! ( nrIters % SPar3 ) ) ) {
  Remt->Start();
  rem = checkRemotion();
  Remt->Stop();
  }

 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if ( rem || add ) {
  status = kFiChgd;        // if variables have just been added/removed
  LstChgItr = nrIters;     // don't let this happen again this iteration
  }

 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 return( status );

 }  // end ( FlwFiOrcl::GetFiStatus )

/*--------------------------------------------------------------------------*/

double FlwFiOrcl::AddTime( void )
{
 return( Addt ? Addt->Read() : 0 );
 }

/*--------------------------------------------------------------------------*/

double FlwFiOrcl::RemTime( void )
{
 return( Remt ? Remt->Read() : 0 );
 }

/*--------------------------------------------------------------------------*/

bool FlwFiOrcl::checkSolution( void ) {


 bool checkOk = true; // boolean for adding

 HpNum MxViol = 0;    // minimum allowed violation
 HpNum lhs;           // lhs of constraint := x_{ij}^k - y_{ij}u_{ij}^k

 // allocate memory for the new solution - - - - - - - - - - - - - - - - - - -

 Row NewXSol = new Number[ NArcs ];
 FRow NewFSol = new FNumber[ NArcs * NComm ];

 VectAssign( NewXSol , Number(0) , NArcs );
 VectAssign( NewFSol , FNumber(0) , NArcs * NComm );


 // get the solution - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 SetXtrSol( NewXSol );
 GetPSol();

 SetFlwSol( NewFSol );
 GetPSol();

 // check if all constraints are satisfied  - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 //mutual capacity constraints  - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( SPar4 != 's' ) {
  MxViol = 0;

  for( Index j = NArcs ; j < NumVar ; j++ ) {
   lhs = 0;
   for( Index k = 0; k < NComm; k++ )
    lhs += NewFSol[ k * NArcs + j ];
   lhs -= NewXSol[ j ] * OrigTotCap[ j ];
   if( lhs < MxViol )
    MxViol = lhs;
   }

  if( MxViol > EpsCon ) {
   checkOk = false;
   FILOG( 1 , endl << "Maximal violation for the mutual constraints is: "
		  << MxViol );
   }
  }

 // individual capacity constraints - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( SPar4 == 's' || SPar4 == 'b' ) {

  MxViol = 0;
  for( Index j = mutual ; j < NumVar ; j++ ) {
   lhs = NewFSol[ InvULambda[ j - mutual ] ] - NewXSol[ InvULambda[ j - mutual ] % NArcs ]
         * OrigCapacities[ InvULambda[ j - mutual ] ];
   if( MxViol < lhs )
    MxViol = lhs;
   }
  MxViol *= ( 1.0 + EpsCon );
  if( MxViol < EpsCon )
   MxViol = EpsCon;

  for( Index k = 0; k < NComm; k++ )
   for( Index j = 0; j < NArcs; j++ ) {
    lhs = NewFSol[ k * NArcs + j ] - NewXSol[ j ] * OrigCapacities[ k * NArcs + j ];
    if( lhs > MxViol ) {
     checkOk = false;
     FILOG( 1 , endl << "Individual capacity constraint ( " << k * NArcs + j
    		 <<  " ) " << " is not satisfied" );
     }
    }
  } // end strong constraints  - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - -  - - - - - - - - - - - - - - - - - - - -

 delete[] NewXSol;
 delete[] NewFSol;

 return( checkOk );

 } // end ( checkSol)

/*--------------------------------------------------------------------------*/
/*------------- METHODS FOR ADDING / REMOVING / CHANGING DATA --------------*/
/*--------------------------------------------------------------------------*/

void FlwFiOrcl::Deleted( cIndex i ) {

 if( i < MaxName )	{

  delete[] OldFSols[ i ];
  OldFSols[ i ] = 0;

  if( OrigXtrCosts && !YIsEasy ) {
   delete[] OldXSols[ i ];
   OldXSols[ i ] = 0;
   }

  }
 else

  SetMaxName( 0 );

 } // end ( FlwFiOrcl::Deleted )

/*--------------------------------------------------------------------------*/

void FlwFiOrcl::Aggregate( cHpRow Mlt , cIndex_Set NmSt , cIndex Dm ,
			   cIndex NwNm )
{
 if( !SPar1 ) // don't do anything if dynamic variables generation is not
  return;     // allowed

 // define the data structure for saving the aggregate primal solution
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 bool ReplacedItem = false; // the item has been substituted  - - - - - - - -
 Index EndVect;             // dimension of the flow vector

 Index InvNwNm;

 Index wFi;
 if( NmSt ) { // get the name of the problem  - - - - - - - - - - - - - - - -
  wFi = OldWFi[ NmSt[ 0 ] ];

  // some exceptions  - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  for( Index i = 0 ; i < Dm ; i++ )
   if( wFi != OldWFi[ NmSt[ i ] ] )
  	throw( NDOException( "FlwFiOrcl::Aggregate: different subgradients." ) );

  // define the length of Flow vector  - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( wFi > NComm || Aggrgtd )
   EndVect = NArcs * NComm;
  else
   EndVect = NArcs;

  // allocate the memory for aggregate solution  - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  for( Index i = 0 ; i < Dm ; i++ )
   if( NmSt[ i ] == NwNm ) {
    InvNwNm = i;
    ReplacedItem = true;
  	break;
    }

  }
 else { // full aggregation case  - - - - - - - - - - - - - - - - - - - - - -
  wFi = OldWFi[ 0 ];

  // some exceptions  - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( wFi != Inf<Index>() )
   throw( NDOException( "FlwFiOrcl::Aggregate: this should not happen" ) );

  for( Index i = 0 ; i < Dm ; i++ )
   if( wFi != OldWFi[ i ] )
    throw( NDOException( "FlwFiOrcl::Aggregate: aggregation of different "
     	 "subgradients." ) );

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  EndVect = NArcs * NComm;

  // allocate the memory for aggregate solution  - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( NwNm < Dm )
   ReplacedItem = true;

  } //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 // some exceptions- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( wFi == 0 )
  throw( NDOException( "FlwFiOrcl::NewGi: don't aggregate the zero-th component" ) );


 // allocate the memory for aggregate solution  - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( !ReplacedItem ) {

  if( OldXSols ) {
   delete[] OldXSols[ NwNm ];
   if( ( wFi >=  NSubPr + 1  ) ||  Aggrgtd )
    OldXSols[ NwNm ] = new Number[ NArcs ];
   else
	OldXSols[ NwNm ] = 0;
   }

  delete[] OldFSols[ NwNm ];
  if( !Aggrgtd && OrigXtrCosts && ( wFi == NSubPr + 1 ) )
   OldFSols[ NwNm ] = 0;
  else
   OldFSols[ NwNm ] = new FNumber[ EndVect ];

  if( NwNm > MaxName )
   throw( NDOException( "FlwFiOrcl::Aggregate:Name is higher than MaxName." ) );

  // tell which is the problem   - - - - - - - - - - - - - - - - - - - - - - -

  OldWFi[ NwNm ] = wFi;

  }  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 // aggregate the solution for variables  - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#if SubVsVar

 if( NmSt ) {
  if( !ReplacedItem ) {
   if( OldXSols && OldXSols[ NwNm ] )  // aggregate extra solution
 	for( Index j = 0; j < NArcs; j++ ) {
     OldXSols[ NwNm ][ j ] = 0;
 	 for( Index i = 0;  i < Dm; i++ )
 	  OldXSols[ NwNm ][ j ] += Mlt[ i ] * OldXSols[ NmSt[i] ][ j ];
 	}
   if( OldFSols[ NwNm ] )  // aggregate flow solution
 	for( Index j = 0; j < EndVect ; j++ ) {
     OldFSols[ NwNm ][ j ] = 0;
     for( Index i = 0; i < Dm; i++ )
      OldFSols[ NwNm ][ j ] += Mlt[ i ] * OldFSols[ NmSt[i] ][ j ];
 	}
   }
  else {
   if( OldXSols && OldXSols[ NwNm ] )  // aggregate extra solution
	for( Index j = 0; j < NArcs; j++ ) {
	 OldXSols[ NwNm ][ j ] *= Mlt[ InvNwNm ];
	 for( Index i = 0; i < InvNwNm ; i++ )
	  OldXSols[ NwNm ][ j ] += Mlt[ i ] * OldXSols[ NmSt[i] ][ j ];
	 for( Index i = InvNwNm + 1 ; i < Dm ; i++ )
	  OldXSols[ NwNm ][ j ] += Mlt[ i ] * OldXSols[ NmSt[i] ][ j ];
	 }
   if( OldFSols[ NwNm ] ) // aggregate flow solution
	for( Index j = 0; j < EndVect ; j++ ) {
	 OldFSols[ NwNm ][ j ] *= Mlt[ InvNwNm ];
	 for( Index i = 0; i < InvNwNm ; i++ )
	  OldFSols[ NwNm ][ j ] += Mlt[ i ] * OldFSols[ NmSt[i] ][ j ];
	 for( Index i = InvNwNm + 1 ; i < Dm ; i++ )
	  OldFSols[ NwNm ][ j ] += Mlt[ i ] * OldFSols[ NmSt[i] ][ j ];
    }
   }
  }
 else { // full aggregation case - - - - - - - - - - - - - - - - - - - - - -
  if( !ReplacedItem ) {
   if( OldXSols  && OldXSols[ NwNm ] )  // aggregate extra solution
    for( Index j = 0; j < NArcs; j++ ) {
     OldXSols[ NwNm ][ j ] = 0;
     for( Index i = 0;  i < Dm; i++ )
      OldXSols[ NwNm ][ j ] += Mlt[ i ] * OldXSols[ i ][ j ];
     }

   if( OldFSols[ NwNm ] ) // aggregate flow solution
	for( Index j = 0; j < EndVect ; j++ ) {
	 OldFSols[ NwNm ][ j ] = 0;
     for( Index i = 0; i < Dm; i++ )
      OldFSols[ NwNm ][ j ] += Mlt[ i ] * OldFSols[ i ][ j ];
    }
   }
  else {
   if( OldXSols && OldXSols[ NwNm ] ) // aggregate extra solution
	for( Index j = 0; j < NArcs; j++ ) {
	 OldXSols[ NwNm ][ j ] *= Mlt[ NwNm ];
	 for( Index i = 0;  i < NwNm; i++ )
	  OldXSols[ NwNm ][ j ] += Mlt[ i ] * OldXSols[ i ][ j ];
	 for( Index i = NwNm + 1 ; i < Dm ; i++ )
	  OldXSols[ NwNm ][ j ] += Mlt[ i ] * OldXSols[ i ][ j ];
	}
   if( OldFSols[ NwNm ] ) // aggregate flow solution
	for( Index j = 0; j < EndVect ; j++ ) {
	 OldFSols[ NwNm ][ j ] *= Mlt[ NwNm ];
	 for( Index i = 0;  i < NwNm ; i++ )
	  OldFSols[ NwNm ][ j ] += Mlt[ i ] * OldFSols[ i ][ j ];
	 for( Index i = NwNm + 1 ; i < Dm ; i++ )
	  OldFSols[ NwNm ][ j ] += Mlt[ i ] * OldFSols[ i ][ j ];
	}
   }
  }  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 // aggregate the solution for items      - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#else

 if( NmSt ) {

  if( !ReplacedItem ) {
   if( OldXSols && OldXSols[ NwNm ] ) {
	VectAssign( OldXSols[ NwNm ] , Number(0) , NArcs );
    for( Index i = 0;  i < Dm; i++ )
   	 for( Index j = 0; j < NArcs; j++ )
      OldXSols[ NwNm ][ j ] += Mlt[ i ] * OldXSols[ NmSt[i] ][ j ];
    }

   if( OldFSols[ NwNm ] ) {
	VectAssign( OldFSols[ NwNm ] , FNumber(0) , EndVect );
    for( Index i = 0; i < Dm; i++ )
     for( Index j = 0; j < EndVect ; j++ )
      OldFSols[ NwNm ][ j ] += Mlt[ i ] * OldFSols[ NmSt[i] ][ j ];
    }

   }
  else {
   if( OldXSols && OldXSols[ NwNm ] ) {
    VectScale( OldXSols[ NwNm ] , Mlt[ InvNwNm ] , NArcs );
    for( Index i = 0;  i < Dm; i++ )
	 if( NmSt[ i ] != NwNm )
	  for( Index j = 0; j < NArcs; j++ )
	   OldXSols[ NwNm ][ j ] += Mlt[ i ] * OldXSols[ NmSt[i] ][ j ];
    }

   if( OldFSols[ NwNm ] ) {
	VectScale( OldFSols[ NwNm ] , Mlt[ InvNwNm ] , EndVect );
    for( Index i = 0; i < Dm; i++ )
     if( NmSt[ i ] != NwNm )
      for( Index j = 0; j < EndVect ; j++ )
       OldFSols[ NwNm ][ j ] += Mlt[ i ] * OldFSols[ NmSt[i] ][ j ];
    }
   }

  }
 else { // full aggregation case - - - - - - - - - - - - - - - - - - - - - -

  if( !ReplacedItem ) {
   if( OldXSols && OldXSols[ NwNm ] ) {
	VectAssign( OldXSols[ NwNm ] , Number(0) , NArcs );
    for( Index i = 0;  i < Dm; i++ )
   	 for( Index j = 0; j < NArcs; j++ )
      OldXSols[ NwNm ][ j ] += Mlt[ i ] * OldXSols[ i ][ j ];
    }

   if( OldFSols[ NwNm ] ) {
	VectAssign( OldFSols[ NwNm ] , FNumber(0) , EndVect );
    for( Index i = 0; i < Dm; i++ )
     for( Index j = 0; j < EndVect ; j++ )
      OldFSols[ NwNm ][ j ] += Mlt[ i ] * OldFSols[ i ][ j ];
    }

   }
  else {
   if( OldXSols && OldXSols[ NwNm ] ) {
    VectScale( OldXSols[ NwNm ] , Mlt[ NwNm ] , NArcs );
    for( Index i = 0;  i < Dm; i++ )
	 if( i != NwNm )
	  for( Index j = 0; j < NArcs; j++ )
	   OldXSols[ NwNm ][ j ] += Mlt[ i ] * OldXSols[ i ][ j ];
    }

   if( OldFSols[ NwNm ] ) {
	VectScale( OldFSols[ NwNm ] , Mlt[ NwNm ] , EndVect );
    for( Index i = 0; i < Dm; i++ )
     if( i != NwNm )
      for( Index j = 0; j < EndVect ; j++ )
       OldFSols[ NwNm ][ j ] += Mlt[ i ] * OldFSols[ i ][ j ];
    }
   }
  } // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#endif

 } // end ( FlwFiOrcl::Aggregate )

/*--------------------------------------------------------------------------*/
/*------------------------------ DESTRUCTOR --------------------------------*/
/*--------------------------------------------------------------------------*/

FlwFiOrcl::~FlwFiOrcl() {


 FILOG2( 1 , Addt , endl << "Variables addition time (sec): " << Addt->Read() );
 FILOG2( 1 , Remt , endl << "Variables remotion time (sec): " << Remt->Read() );
 FILOG( 1 , endl );

 // delete graph varaibles  - - - - - - - - - - - - - - - - - - - - - - - - - -

 delete[] OrigCosts;
 delete[] OrigXtrCosts;
 delete[] OrigTotCap;
 delete[] OrigCapacities;
 delete[] OrigDeficits;

 // delete vocabulary - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 delete[] ZLambdaCount;
 delete[] InvULambda;

 if( LmbCmp ) {
  for( Index i = 0 ; i < NComm ; i++ )
   LmbCmp[ i ].clear();
  delete[] LmbCmp;
  }

 // delete auxilary variables - - - - - - - - - - - - - - - - - - - - - - - - -
 delete[] SolvedP;
 delete[] LsHasChgd;

 delete[] SGBse1;

 // delete primal solution  - - - - - - - - - - - - - - - - - - - - - - - - - -

 delete[] FSolution;

 delete[] XSolution;
 delete[] CoefObj;

 SetMaxName( 0 );

 delete Addt;
 delete Remt;

 } // ( end FlwFiOrcl::~FlwFiOrcl() )

/*--------------------------------------------------------------------------*/
/*-------------------------- PRIVATE METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/

void FlwFiOrcl::XSolA( Row Xsl ) {

 // data structure for either multipliers or "easy" component solution
 // "easy" components solution - - - - - - - - - - - - - - - - - - - - - - - -

 cHpRow mult;
 Index D;
 cIndex_Set I;

 // compute the ys- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( YIsEasy ) {
  mult = Slvr->ReadMult( I , D , ( NSubPr + 1 ) ); // read  x[ wFi ]*
  if( I )   // case: "sparse"
   for( Index i = 0; i < D; i++ )
	Xsl[ I[ i ] ] = mult[ i ];
  else      // case: dense
   for( Index i = 0; i < D; i++ )
	Xsl[ i ] = mult[ i ];
  }
 else {
  mult = Slvr->ReadMult( I , D );
  if( I )  // case: "sparse"
   for( Index i = 0; i < D; i++ ) {
	if( OldXSols[ I[ i ] ] )
	 for( Index j = 0; j < NArcs; j++ )
	  Xsl[ j ] += mult[ i ] * OldXSols[ I[ i ] ][ j ];
    }
  else     // case: dense
   for( Index i = 0; i < D; i++ ) {
	if( OldXSols[ i ] )
	 for( Index j = 0; j < NArcs; j++ )
	  Xsl[ j ] += mult[ i ] * OldXSols[ i ][ j ];
    }
  }

 } // end ( FlwFiOrcl::XSolA )


/*--------------------------------------------------------------------------*/

void FlwFiOrcl::FSolA( FRow Flw ) {

 // data structure for either multipliers or "easy" component solution
 // "easy" components solution - - - - - - - - - - - - - - - - - - - - - - - -

 cHpRow mult;
 Index D;
 cIndex_Set I;

 // compute x = Sum{i \in D[]} Theta[ i ] * x[ i ] - - - - - - - - - - - - - -

 mult = Slvr->ReadMult( I , D );

 if( I ) {
  for( Index i = 0; i < D; i++ )
   if( mult[ i ] > 0 && ( OldWFi[ I[ i ] ] != NSubPr + 1 ) ) {
	if( ( OldWFi[ I[ i ] ] > NSubPr + 1 ) || Aggrgtd )
     for( Index j = 0; j < NArcs * NComm; j++ )
      Flw[ j ] += mult[ i ] * OldFSols[ I[ i ] ][ j ];
	else
	 for( Index j = 0; j < NArcs; j++ )
	  Flw[ ( OldWFi[ I[ i ] ] - 1 ) * NArcs + j ]
		+= mult[ i ] * OldFSols[ I[ i ] ][ j ];
	}
  }
 else
  for( Index i = 0; i < D; i++ )
   if( mult[ i ] > 0 && ( OldWFi[ i ] != NSubPr + 1 ) ) {
    if( ( OldWFi[ i ] > NSubPr + 1 ) || Aggrgtd )
     for( Index j = 0; j < NArcs * NComm; j++ )
      Flw[ j ] += mult[ i ] * OldFSols[ i ][ j ];
    else
     for( Index j = 0; j < NArcs; j++ )
      Flw[ ( OldWFi[ i ] - 1 ) * NArcs + j ]
    	+= mult[ i ] * OldFSols[ i ][ j ];
	}

 } // end ( FlwFiOrcl::FSolA )

/*--------------------------------------------------------------------------*/

bool FlwFiOrcl::checkAddition( void ) {

 // take into account that adding is possible
 //	if the individual capacities
 // constrained relaxed are considered and the generation of betas
 // is also allowed

 bool add = false;    // boolean for adding
 HpNum MxViol = 0;    // minimum allowed violation
 HpNum lhs;           // lhs of constraint := x_{ij}^k - y_{ij}u_{ij}^k

 // allocate memory for the new solution - - - - - - - - - - - - - - - - - - -

 Row NewXSol = new Number[ NArcs ];
 FRow NewFSol = new FNumber[ NArcs * NComm ];

 VectAssign( NewXSol , Number(0) , NArcs );
 VectAssign( NewFSol , FNumber(0) , NArcs * NComm );

 XSolA( NewXSol );
 FSolA( NewFSol );

 // add those betas that correspond to the violated constraints - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 Index count = 0;    // current number of betas

 // in order to avoid the check if any newly generated constraint is already
 // in the pool, we set the minimum allowed violation to a value *strictly
 // larger* than the largest violation of any constraint in the pool

 #if SEP_TYPE

  // mark the variable already involved - - - - - - - - - - - - - - - - - - -

  bool *ULambda = new bool[ NArcs * NComm ];
  VectAssign( ULambda , false , NArcs * NComm );

  for( Index j = mutual; j < NumVar; j++ )
   ULambda[ InvULambda[ j - mutual ] ] = true;

  for( Index j = 0 ; j < NArcs * NComm ; j++ ) {
   lhs = NewFSol[ j ] - NewXSol[ j % NArcs ] * OrigCapacities[ j ];
   MxViol = EpsCon * max( OrigCapacities[ j ] , 1.0 );
   if( lhs > MxViol && !( ULambda[ j ] ) ) {

	InvULambda[ (count + NumVar) - mutual ] = j;

	if( !Aggrgtd )
	 LmbCmp[ j / NArcs ].push_back( count + NumVar - mutual );

	if( ZLambdaCount )
     ZLambdaCount[ j ] = 0;

    L2Cnst += ( OrigCapacities[ j ] * OrigCapacities[ j ] );
    count++;
    }
   }

  delete[] ULambda;

 #else

  if( NumVar ) { // compute the Max violated if some constraints are in the pool
	 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   MxViol = 0;
   for( Index j = mutual ; j < NumVar ; j++ ) {
    lhs = NewFSol[ InvULambda[ j - mutual ] ] - NewXSol[ InvULambda[ j - mutual ] % NArcs ]
           * OrigCapacities[ InvULambda[ j - mutual ] ];
    if( MxViol < lhs )
     MxViol = lhs;
    }

   MxViol *= ( 1.0 + EpsCon );
   if( MxViol < EpsCon )
    MxViol = EpsCon;
  }
  else
   MxViol = EpsCon;

  FILOG( 2 , endl <<  "           " << " MaxViol in the pool = " << MxViol );

  // find constraints to be eliminated  - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  for( Index k = 0; k < NComm; k++ )
   for( Index j = 0; j < NArcs; j++ ) {
     lhs = NewFSol[ k * NArcs + j ] - NewXSol[ j ]
	      * OrigCapacities[ k * NArcs + j ];
     if( lhs > MxViol ) {

	  InvULambda[ (count + NumVar) - mutual ] = k * NArcs + j;

	  if( !Aggrgtd )
	   LmbCmp[ k ].push_back( count + NumVar - mutual );

	  if( ZLambdaCount )
	   ZLambdaCount[ k * NArcs + j ] = 0;

	  L2Cnst += OrigCapacities[ k * NArcs + j ] * OrigCapacities[ k * NArcs + j ];
      count++;
      }
     }
 #endif



 // insert the obtained constraints  - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( count > 0 ) {              // the NDOSolver can assume that, when the
  NumVar += count;              // method is called,  Oracle->GetNumVar()
  FILOG( 2 , endl << "           " << " Added " << count << " constraints " );
  Slvr->AddVariables( count );  // already returns the number of variables *after*
  add = true;                   // the addiction of the new ones
  }

 delete[] NewXSol;
 delete[] NewFSol;
 return add;

 } // end ( FlwFiOrcl::checkAddition )

/*--------------------------------------------------------------------------*/

bool FlwFiOrcl::checkRemotion( void ) {

 // take into account that the removal of some betas
 // is possible when the individual capacities
 // constrained relaxed are considered and the generation of betas
 // is also allowed

 bool rem = false;             // boolean for removing
 Index count = 0;              // number of variables beta to remove

 Index_Set removed = new Index[ ( NumVar - mutual ) + 1 ]; // vocabulary of items
                                                       // which we must eliminate

 if( NumVar > mutual ) {

  // ensure that no variable with Lambda[] > 0 is ever removed - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // it is especially important to avoid removing nonzero variables from the
  // "current point" of the NDO algorithm, as this may easily lead to cycling

  Index LBd;
  cIndex_Set LB;
  cLMRow L = Slvr->ReadSol( LB , LBd );

  // the "inactive" constraints are removed to make space for the new ones
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( LB ) {

   while( *LB < mutual ) {
	LB++;
	L++;
	}

   for( Index h ; ( h = *(LB++) ) < Inf<Index>() ; )
	if( ( *(L++) < EpsCon ) &&
	    ( ZLambdaCount[ InvULambda[ h - mutual ] ] >= SPar2 ) )
	 removed[ count++ ] = h ;

   }
  else {

   for( Index i = mutual ; i < NumVar ; i++ )
    if( ( L[ i ] < EpsCon ) &&
	   ( ZLambdaCount[ InvULambda[ i - mutual ] ] >= SPar2 ) )
	 removed[ count++ ] = i ;

   }
                                     // write Inf<Index>() at the end
  removed[ count ] = Inf<Index>();          // removed[]

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( count > 0  ) {

   // remove the variables in the solver *before* changing the internal
   // information of the Oracle

   Slvr->RemoveVariables( removed , count );
   FILOG( 2 , endl << " Removed constraints = " << count );

   // set to zero the elements in ZLambdaCount whose the indexes belong
   // to set of removed variables

   Index EntrLmbd;
   vector<Index>::iterator it;
   for( Index i = 0 ;  removed[ i ] < Inf<Index>();  i++ ) {
    removed[ i ] -= mutual;

    EntrLmbd = InvULambda[ removed[ i ] ];

    if( !Aggrgtd ) {
     it = find ( LmbCmp[ EntrLmbd / NArcs ].begin() ,
    	    LmbCmp[ EntrLmbd / NArcs ].end(), removed[ i ] );
     LmbCmp[ EntrLmbd / NArcs ].erase( it );
     }

	L2Cnst -= ( OrigCapacities[ EntrLmbd ] * OrigCapacities[ EntrLmbd ] );
    ZLambdaCount[ EntrLmbd ] = 0;

    }

   // "compacts" InvULambda deleting the elements whose indexes
   // are in removed[]

   Compact( InvULambda , removed , NumVar - mutual );


   for( Index j = NumVar - count ; j < MaxNumVar ; j++ )
	InvULambda[ j - mutual ] = Inf<Index>();

   // update NumVar
   NumVar -= count;

   rem = true;

   }

  } // end ( if: NumVar > NArcs )

 delete[] removed;

 return rem;

 } // end ( FlwFiOrcl::checkRemotion )

/*--------------------------------------------------------------------------*/

void FlwFiOrcl::CopyGi( cIndex Name , cIndex wFi , cRow XtrSol , cFRow FlwSol )
{

 // copy the extra solution  - - - - - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( OrigXtrCosts && !YIsEasy ) {
  delete[] OldXSols[ Name ];
  if( wFi >=  NSubPr + 1 ||  Aggrgtd  )  {
   OldXSols[ Name ] = new Number [ NArcs ];
   VectAssign( OldXSols[ Name ] , XtrSol , NArcs );
   }
  else
   OldXSols[ Name ] = 0;
  }

 // copy the flow solution - - - - - - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 Index EndVect;
 if( wFi > NSubPr || Aggrgtd )
  EndVect = NArcs * NComm;
 else
  EndVect = NArcs;

 delete[] OldFSols[ Name ];
 if( !( !Aggrgtd && OrigXtrCosts && wFi == NSubPr + 1 ) ) {
  OldFSols[ Name] = new FNumber[ EndVect ];
  VectAssign( OldFSols[ Name ] , FlwSol , EndVect );
  }
 else
  OldFSols[ Name ] = 0;

 } // end ( FlwFiOrcl::SetGiName )

/*--------------------------------------------------------------------------*/

void FlwFiOrcl::FindGlobalLipschitz( void ) {

 L2Cnst = 0;
 if( OrigXtrCosts )
  for( Index j = 0 ; j < NArcs ; j++ ) {
   if( SPar4 == 'w' || SPar4 == 'b' )
    L2Cnst += ( OrigTotCap[ j ] * OrigTotCap[ j ] );
   if( ( SPar4 == 's' || SPar4 == 'b' ) && ( SPar1 == 0 ) )
    for( Index k = 0 ; k < NComm ; k++ )
     L2Cnst += ( OrigCapacities[ k * NArcs + j ] * OrigCapacities[ k * NArcs + j ] );
   }
 else
  for( Index j = 0 ; j < NArcs ; j++ )
   L2Cnst += ( OrigTotCap[ j ] * OrigTotCap[ j ] );

} // end ( FlwFiOrcl::FindGlobalLipschitz )

/*--------------------------------------------------------------------------*/

void FlwFiOrcl::SolveLagrangian( cIndex SubPName ) {

 if( YIsEasy || !OrigXtrCosts )
  if( SubPName == NComm )
   throw( NDOException( "FlwFiOrcl::SolveLagrangian: this should not be happen" ) );

 if( SubPName == NComm && !SolvedP[ NComm ] ) { // solve the Y subproblem

  // compute original costs   - - - - - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  for( Index j = 0; j < NArcs; j++ )
   CoefObj[ j ] = OrigXtrCosts[ j ];

  if( LamBase )
   for( Index i = 0 ; i < LamBDim ; i++ ) {

    // the following if instruction distinguishes the alpha part of lambda from
    // the beta one  - - - - - - - - - - - - - -  - - - - - - - - - - - - - -
    //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    if( LamBase[ i ] < mutual )
     CoefObj[ LamBase[ i ]  ] -= OrigTotCap[ LamBase[ i ] ] * Lambda[ i ];
    else
     if( SPar4 == 's' || SPar4 == 'b' )
      if( SPar1 > 0 ) {
       CoefObj[ InvULambda[ LamBase[ i ] - mutual ] % NArcs ] -=
    		   OrigCapacities[ InvULambda[ LamBase[ i ] - mutual ] ] * Lambda[ i ];
       }
      else
       CoefObj[ LamBase[ i ] % NArcs ] -=
    		   OrigCapacities[ LamBase[ i ] - mutual ] * Lambda[ i ];
    }
  else
   for( Index i = 0 ; i < NumVar ; i++ ) {

    // the following if instruction distinguishes the alpha part of lambda from
    // the beta one  - - - - - - - - - - - - - -  - - - - - - - - - - - - - -
    //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    if( i < mutual )
	 CoefObj[ i ] -= OrigTotCap[ i ] * Lambda[ i ];
	else
	 if( SPar4 == 's' || SPar4 == 'b' )
	  if( SPar1 > 0 ) {
	   CoefObj[ InvULambda[ i - mutual ] % NArcs ] -=
	    		   OrigCapacities[ InvULambda[ i - mutual ] ] * Lambda[ i ];
	   }
	  else
	    CoefObj[ i % NArcs ] -= OrigCapacities[ i - mutual ] * Lambda[ i ];
	}

  SolvedP[ NComm ] = true;
  } // end Y component solve - - - - - - - - - - - - - - - - - - - - - - - - -
 else
  if( SubPName < NComm ) {  // solve a flow component  - - - - - - - - - - - -

   // set the name of the problem  - - - - - - - - - - - - - - - - - - - - - -
   //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   SetSubP( SubPName + 1 );

   if( !SolvedP[ SubPName ] ) {

    // update the costs of the flow in a particular subproblem  - - - - - - - -
    //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    Index_Set nms = new Index[ NArcs + 1 ];
    CRow NwCsts = new CNumber[ NArcs ];
    Index offset = SubPName * NArcs;
    Index Entijk;

    // set the original costs - - - - - - - - - - - - - - - - - - - - - - - - -
    //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    for( Index j = 0 ; j < NArcs ; j++ ) {
     nms[ j ] = offset + j;
     NwCsts[ j ] = OrigCosts[ offset + j ];
     }

    if( LamBase )
  	 for( Index i = 0 ; i < LamBDim ; i++ ) { // the following if instruction
	  // distinguishes the alpha part of lambda from the beta ones  - - - - - -
	  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	  if( LamBase[ i ] < mutual )
	   NwCsts[ LamBase[ i ] ] += Lambda[ i ];
	  else
	   if( SPar4 == 's' || SPar4 == 'b' ) {
	    if( SPar1 > 0 )
	     Entijk = InvULambda[ LamBase[ i ] - mutual ];
	    else
	     Entijk =  LamBase[ i ] - mutual;
	    if( ( Entijk / NArcs ) == SubPName )
	     NwCsts[ Entijk % NArcs ] += Lambda[ i ];
	    }
	  }
    else {

     // add the alpha - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

     for( Index j = 0 ; j < mutual ; j++ )
      NwCsts[ j ] += Lambda[ j ];

     // add the beta  - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

     if( SPar4 == 's' || SPar4 == 'b' ) {
      if( SPar1 > 0 )
       for( vector<Index>::iterator it = LmbCmp[ SubPName ].begin() ;
    		   it != LmbCmp[ SubPName ].end(); ++it )
        NwCsts[ InvULambda[ *it ] % NArcs ] += Lambda[ *it + mutual ];
      else
       for( Index j = 0 ; j < NArcs ; j++ )
        NwCsts[ j ] += Lambda[ offset + mutual + j ];
      }

     }

    nms[ NArcs ] = Inf<Index>();
    ChgCosts( NwCsts , nms );

    delete[] nms;
    delete[] NwCsts;

    // solve the subproblem - - - - - - - - - - - - - - - - - - - - - - - - -
    //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    SolveMMCF();
    SolvedP[ SubPName ] = true;

    }

   }  // end if( SubPName < NComm && !SolvedP[ NComm ] )
  else {

   for( Index i = 1; i < NComm; i++ )
    if( SolvedP[ 0 ] != SolvedP[ i ] )
	 throw( NDOException( "FlwFiOrcl::SolveLagrangian: this should not be happen" ) );

   // set the name of the problem  - - - - - - - - - - - - - - - - - - - - - -
   //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   if( Aggrgtd )
    SetSubP( 1 );
   else
	SetSubP( NComm + 1 );

   if( !SolvedP[ 0 ] ) {

	// update the costs of the flow subproblems  - - - - - - - - - - - - - - -
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	Index Entijk;
	CRow NwCsts = new CNumber[ NArcs * NComm ];
    VectAssign( NwCsts , OrigCosts , NArcs * NComm );

    if( LamBase )
     for( Index i = 0 ; i < LamBDim ; i++ ) {

      // the following if instruction distinguishes the alpha part of lambda from
	  // the beta one  - - - - - - - - - - - - -  - - - - - - - - - - - - - -
      //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      if( LamBase[ i ] < mutual ) {
       Entijk = LamBase[ i ];
       for( Index k = 0 ; k < NComm ; k++ )
        NwCsts[ k * NArcs + Entijk ] += Lambda[ i ];
       }
      else
       if( SPar4 == 's' || SPar4 == 'b' ) {
        if( SPar1 > 0 )
         Entijk = InvULambda[ LamBase[ i ] - mutual ];
        else
         Entijk =  LamBase[ i ] - mutual;
        NwCsts[ Entijk ] += Lambda[ i ];
        }
      }
    else
     for( Index i = 0 ; i < NumVar ; i++ ) {

      // the following if instruction distinguishes the alpha part of lambda from
	  // the beta one  - - - - - - - - - - - - -  - - - - - - - - - - - - - -
	  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      if( i < mutual )
       for( Index k = 0 ; k < NComm ; k++ )
    	NwCsts[ k * NArcs + i ] += Lambda[ i ];
      else
       if( SPar4 == 's' || SPar4 == 'b' ) {
        if( SPar1 > 0 )
  	     Entijk = InvULambda[ i - mutual ];
        else
         Entijk =  i - mutual;
        NwCsts[ Entijk ] += Lambda[ i ];
        }
      }

    ChgCosts( NwCsts );
    delete[] NwCsts;

    // solve the subproblem - - - - - - - - - - - - - - - - - - - - - - - - -
    //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    SolveMMCF();
    for( Index i = 0; i < NComm; i++ )
     SolvedP[ i ] = true;
    }

   } // end  ( SubPName > NComm )

 } // end ( FlwFiOrcl::SolveLagrangian( cIndex SubPr ) )

/*--------------------------------------------------------------------------*/
/*-------------------------- End File FlwFiOrcl.C --------------------------*/
/*--------------------------------------------------------------------------*/

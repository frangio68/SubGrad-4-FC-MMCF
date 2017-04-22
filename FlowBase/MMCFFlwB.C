/*--------------------------------------------------------------------------*/
/*--------------------------- File MMCFFlwB.C ------------------------------*/
/*--------------------------------------------------------------------------*/
/*--                                                                      --*/
/*-- Multicommodity Min Cost Flow (MMCF) Problems solvers: solves (MMCF)s --*/
/*-- where the Mutual Capacity constraints are relaxed, i.e., k           --*/
/*-- independent single-commodity Min Cost Flow or Shortest Path          --*/
/*-- subproblems. "extra" variables or constraints are not supported. The --*/
/*-- class implements the MMCFClass interface, see MMCFClas.h.            --*/
/*--                                                                      --*/
/*--                            VERSION 3.03                              --*/
/*--                           17 - 12 - 2013                             --*/
/*--                                                                      --*/
/*--                 Original Idea and Implementation by:                 --*/
/*--                                                                      --*/
/*--                          Antonio Frangioni                           --*/
/*--                                                                      --*/
/*--                       Operations Research Group                      --*/
/*--                      Dipartimento di Informatica                     --*/
/*--                          Universita' di Pisa                         --*/
/*--                                                                      --*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- IMPLEMENTATION -------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*-------------------------------- MACROS ----------------------------------*/
/*--------------------------------------------------------------------------*/

#define HAVE_MCF 84

#if( HAVE_MCF )
 #define HAVE_SLV 0
#endif

/** These macros select which types of flow subproblems solvers are available
   for the MMCFFlwBase class; furthermore, for each type of flow problem
   select which solver (derived class of MCFClass) is used.

   HAVE_MCF > 0 means that a Min Cost Flow (MCF) solver is available: the
                actual value tells which solver is used. The first four bits
   of the value can be used to set some special initialization for the solver.
   Possible values are:

     16    the RelaxIV solver
     +1    the Auction() initialization is used (if available)

     32    the MCFCplex solver
     +k    a value of k betweek 1 and 3 set the pricing rule used by the
           network simplex; 0 means automatic (the default)

     48    the MCFZIB solver
     +1    use the dual network simplex (by default, the primal is used)
     +2    use the first-element pricing rule (by default, Dantzig rule is
           used)
     +4    use the Multiple Partial Pricing rule (only for primal simplex)

     64    the CS2 solver.

     80    the MCFSimplex solver
     +1    use the dual network simplex (by default, the primal is used)
     +2    use the first-element pricing rule (by default, Dantzig rule is
           used)
     +4    use the Candidate List Pivot Rule (only for primal simplex)

   HAVE_SLV > 0 means that each MCF solver is actually two of them: using the
   class MCFClone, each solver solves each problem twice using two (possibly
   different) solvers. The first solver is chosen by HAVE_MCF, the second by
   the value of HAVE_SLV; the meaning of the values is the same. */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "MMCFFlwB.h"

#include "OPTvect.h"

#if( HAVE_MCF )
 #if( HAVE_MCF < 32 )
  #include "RelaxIV.h"
  #define MASTERMCF RelaxIV
 #elif( HAVE_MCF < 48 )
  #include "MCFCplex.h"
  #define MASTERMCF MCFCplex
 #elif( HAVE_MCF < 64 )
  #include "MCFZIB.h"
  #define MASTERMCF MCFZIB
 #elif( HAVE_MCF < 80 )
  #include "CS2.h"
  #define MASTERMCF CS2
 #elif( HAVE_MCF < 96 )
  #include "MCFSimplex.h"
  #define MASTERMCF MCFSimplex
 #endif

 #if( HAVE_SLV )
  #if( HAVE_SLV < 32 )
   #include "RelaxIV.h"
   #define SLAVEMCF RelaxIV
  #elif( HAVE_SLV < 48 )
   #include "MCFCplex.h"
   #define SLAVEMCF MCFCplex
  #elif( HAVE_SLV < 64 )
   #include "MCFZIB.h"
   #define SLAVEMCF MCFZIB
  #elif( HAVE_SLV < 80 )
   #include "CS2.h"
   #define SLAVEMCF CS2
  #elif( HAVE_MCF < 96 )
   #include "MCFSimplex.h"
   #define SLAVEMCF MCFSimplex
  #endif

  #include "MCFClone.h"

  #define MCFTYPE MCFClone<MASTERMCF,SLAVEMCF>
 #else
  #define MCFTYPE MASTERMCF
 #endif
#endif

#if( FlowBase_HAVE_SPT == 1 )
 #include "SPTree.h"
#endif

/*--------------------------------------------------------------------------*/
/*-------------------------------- MACROS ----------------------------------*/
/*--------------------------------------------------------------------------*/

#if( HAVE_MCF )
 #define MCFALG( x ) ( x & 1 ? false : true )

 #define MCFZIBPRC( x ) ( x & 2 ? MCFZIB::kFrstElA : \
			( x & 4 ? MCFZIB::kMltPrPr : \
				  MCFZIB::kDantzig ) )

 #define MCFSMXPRC( x ) ( x & 2 ? MCFSimplex::kFirstEligibleArc : \
                        ( x & 4 ? MCFSimplex::kCandidateListPivot : \
			          MCFSimplex::kDantzig ) )
#endif

/*--------------------------------------------------------------------------*/
/*-------------------------------- USING -----------------------------------*/
/*--------------------------------------------------------------------------*/

#if( OPT_USE_NAMESPACES )
using namespace MMCFClass_di_unipi_it;
#endif

/*--------------------------------------------------------------------------*/
/*----------------------------- CONSTANTS ----------------------------------*/
/*--------------------------------------------------------------------------*/

static MMCFClass::cFNumber FlwMltCnst = Eps<MMCFClass::FNumber>();

static MMCFClass::cFNumber DfctCnst = 3;

static MMCFClass::cCNumber CstMltCnst = Eps<MMCFClass::CNumber>() * 100;

/** If flows [costs] are of a float type, the corresponding epsilons in the
   objects of class MCFClass will be set to FlwMltCnst * the average arc
   capacity [ CstMltCnst * max( average absolute value of arc cost , 1 ) ].
   The epsilon for deficits will be set to the above multiplied by m / n
   and multiplied by DfctCnst, meaning that there are no nodes that have
   more neighbours than DfctCnst * (m / n) (the average density). */

/*--------------------------------------------------------------------------*/
/*------------------- IMPLEMENTATION OF MMCFFlwBase ------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- PUBLIC METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/
/*---------------------------- CONSTRUCTORS --------------------------------*/
/*--------------------------------------------------------------------------*/

MMCFFlwBase::MMCFFlwBase( Graph *Gh ) : MMCFClass()
{
 SetGraph( Gh );
 }

/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

void MMCFFlwBase::SetGraph( Graph *Gh )
{
 NComm = Gh->NrComm();
 NNodes = Gh->NrNodes();
 NArcs = Gh->NrArcs();

 // initial memory allocation - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 MCFs = new PMCF[ NComm ];  // vector of pointers to (MCF) solver objects
 Cost = new CNumber[ max( NArcs , NNodes ) ];  // cost-related temp.

 // get the "worst case" UB from Graph- - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                  
 UBnd = FONumber( Gh->UpperBound() );

 // allocate a temporary- - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // the max is for the odd case where NArcs == NNodes - 1, i.e., the graph is
 // still connected but it's a tree, and therefore has more nodes than arcs

 Ftmp = new FNumber[ max( NArcs , NNodes ) ];

 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // main cycle: repeat for each commodity - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 bool HavePtrs = true;
 for( Index k = 0 ; k < NComm ; k++ ) {
  // select capacities to be used - - - - - - - - - - - - - - - - - - - - - -

  cFRow Caps = Gh->CapacitiesK( k );
  if( ! Caps )                  // if k has no individual capacities ...
   Caps = Gh->TotCapacities();  // ... use the mutual ones

  // select the proper solver - - - - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  switch( Gh->ProblemType( k ) ) {
   #if( FlowBase_HAVE_SPT )
    case( Graph::kSPT ):  // (SPT) subproblems- - - - - - - - - - - - - - - -
    {
     MCFs[ k ] = new SPTree( NNodes , NArcs , Gh->Directed() );

     Caps = 0;
     break;
     }
   #endif
   #if( HAVE_MCF )
    default:  // if available, (MCF) subproblems are the default type - - - -
    {
     MCFTYPE *MCFk = new MCFTYPE( NNodes , NArcs );

     #if( HAVE_MCF < 32 )
      #if( AUCTION && ( HAVE_MCF & 1 ) )
       MCFk->SetAuction();
      #endif
     #elif( HAVE_MCF < 48 )
      MCFk->SetPar( CPX_PARAM_NETPPRIIND , int( HAVE_MCF & 48 ) );
     #elif( HAVE_MCF < 64 )
      MCFk->SetAlg( MCFALG( HAVE_MCF ) , MCFZIBPRC( HAVE_MCF ) );
     #elif( HAVE_MCF < 80 )
     #elif( HAVE_MCF < 96 )
      MCFk->SetAlg( MCFALG( HAVE_MCF ) , MCFSMXPRC( HAVE_MCF ) );
     #endif

     #if( HAVE_SLV )
      #if( HAVE_SLV < 32 )
       #if( AUCTION && ( HAVE_SLV & 1 ) )
        MCFk->SlvMCF->SetAuction();
       #endif
      #elif( HAVE_SLV < 48 )
       MCFk->SlvMCF->SetPar( CPX_PARAM_NETPPRIIND , int( HAVE_SLV & 48 ) );
      #elif( HAVE_SLV < 64 )
       MCFk->SlvMCF->SetAlg( MCFALG( HAVE_SLV ) , MCFZIBPRC( HAVE_SLV ) );
      #elif( HAVE_SLV < 80 )
      #elif( HAVE_MCF < 96 )
       MCFk->SlvMCF->SetAlg( MCFALG( HAVE_MCF ) , MCFSMXPRC( HAVE_MCF ) );
      #endif
     #endif

     MCFs[ k ] = MCFk;
     break;
     }
   #else
    default:  // bad type - - - - - - - - - - - - - - - - - - - - - - - - - -

     throw(
        MMCFException( "MMCFFlwBase: no solver for this type of problem" ) );
   #endif

   }  // end switch() - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  // remove INFs from the deficit vector (MCF solvers don't like 'em) - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  VectAssign( Ftmp , Gh->DeficitsK( k ) , NNodes );
  for( FRow tFt = Ftmp + NNodes ; tFt-- > Ftmp ; )
   if( *tFt == Inf<FNumber>() )
    *tFt = 0;

  // load the problem in the solver - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  MCFs[ k ]->LoadNet( NNodes , NArcs , NNodes , NArcs , Caps ,
	              Gh->CostsK( k ) , Ftmp , Gh->StartN() , Gh->EndN() );

  // check for pointers availability- - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( ! ( MCFs[ k ]->MCFGetX() ) )
   HavePtrs = false;

  }  // end for( k ) (main cycle) - - - - - - - - - - - - - - - - - - - - - -
     // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 // auxiliary data structures allocation- - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( HavePtrs ) {  // the temporary is no longer necessary
  delete[] Ftmp;
  Ftmp = 0;
  }

 #if( CHGARCS_MMCF > 1 )
  ClsdAK = new bool[ NArcs * NComm ];
  VectAssign( ClsdAK , false , NArcs * NComm );
 #endif

 // final initializations - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 NSubPr = 1;
 OptVal = new FONumber[ NSubPr + 2 ];
 OptVal[ 0 ] = 0;
 VectAssign( OptVal + 1 , - Inf<FONumber>() , NSubPr + 1 );

 AFSPD = new Index[ 2 ];
 AFSPD[ 0 ] = 0;
 AFSPD[ 1 ] = NComm;
 
 #if( FlowBase_HAVE_SPT )
  DFSPDn = 0;
  DFSPDb = 0;
  DFSPD = 0;
 #endif

 }  // end( SetGraph() )

/*--------------------------------------------------------------------------*/

void MMCFFlwBase::ReOptimize( bool RO , cIndex whch )
{
 if( whch == Inf<Index>() )
  for( Index k = NComm ; k-- ; )
   MCFs[ k ]->SetPar( MCFClass::kReopt , int( RO ? MCFClass::kYes :
					           MCFClass::kNo ) );
 else
  MCFs[ whch ]->SetPar( MCFClass::kReopt , int( RO ? MCFClass::kYes :
						     MCFClass::kNo ) );

 }  // end( ReOptimize )

/*--------------------------------------------------------------------------*/

void MMCFFlwBase::SetNrSubP( Index FSP )
{
 if( ! FSP )
  FSP = 1;

 if( FSP == NSubPr )  // nothing to do
  return;             // do nothing

 // compute how many subproblems can be built - - - - - - - - - - - - - - - -

 #if( FlowBase_HAVE_SPT )
  Index MCFn = 0;      // number of flow subproblems that are MCFs
  Index TotDN = 0;     // total number of destinations in the SPTs
  if( FSP > NComm ) {  // compute how many subproblems can be built
   for( Index k = 0 ; k < NComm ; k++ ) {
    SPTree *SPTk = dynamic_cast<SPTree*>( MCFs[ k ] );
    if( SPTk )
     TotDN += SPTk->DestN();
    else
     MCFn++;
    }

   FSP = min( FSP , Index( MCFn + TotDN ) );
   }
 #else
  FSP = min( FSP , NComm );
 #endif

 #if( FlowBase_HAVE_SPT )
  delete[] DFSPDb;
  delete[] DFSPDn;
  delete[] DFSPD;
  DFSPDb = 0;
  DFSPD = 0;
 #endif
 delete[] AFSPD;

 NSubPr = FSP;

 if( NSubPr == 1 ) {       // the problem is not separable
  AFSPD = new Index[ 2 ];
  AFSPD[ 0 ] = 0;
  AFSPD[ 1 ] = NComm;
  }
 else                      // the problem is separable
  if( NSubPr <= NComm ) {  // subproblems == sets of commodities- - - - - - -
   // construct AFSPD[]: the commodities corresponding to subproblem h are
   // those with AFSPD[ h - 1 ] <= k < AFSPD[ h ] (h = 1 ... NSubPr)
   //
   // distribute the commodities evenly among the NSubPr subproblems: each
   // subproblem gets NComm / NSubPr consecutive (!!) commodities, but
   // NComm % NSubPr subproblems get one more because of roundings

   AFSPD = new Index[ NSubPr + 1 ];
   cIndex i = NComm / NSubPr;
   Index j = NComm % NSubPr;
   AFSPD[ 0 ] = 0;
   for( Index h = 0 ; h++ < NSubPr ; ) {
    AFSPD[ h ] = AFSPD[ h - 1 ] + i;
    if( j ) { AFSPD[ h ]++; j--; }
    }
   }
 #if( FlowBase_HAVE_SPT )
  else {  // subproblems == subsets of O/D pairs in SPTs- - - - - - - - - - -
   // the minimum length DstLen is computed such that by dividing each SPT
   // into subproblems with (roughly) DstLen destinations the total number
   // of subproblems does not exceed NSubPr

   Index DstLen = TotDN / ( NSubPr - MCFn );  // BlkLen >= 1
   for( ; ; DstLen++ )   // start a cycle where the number of subproblems
   {                     // corresponding to a lenght of DstLen is computed
    Index BLFSP = MCFn;  // and DstLen is increased if it is not <= NSubPr
    for( Index k = 0 ; k < NComm ; k++ ) {
     SPTree *SPTk = dynamic_cast<SPTree*>( MCFs[ k ] );
     if( SPTk )
      BLFSP += CeilDiv( SPTk->DestN() , DstLen );
     }
      
    if( BLFSP <= NSubPr ) {
     NSubPr = BLFSP;
     break;
     }
    }

   // now, by dividing each each SPT into subproblems with (roughly)
   // DstLen destinations each, we obtain exactly NSubPr subproblems

   DFSPD = new Index[ NSubPr ];
   DFSPDn = new Index[ NSubPr ];
   DFSPDb = new Index_Set[ NSubPr ];

   Index ws = 0;
   for( Index k = 0 ; k < NComm ; k++ ) {
    SPTree *SPTk = dynamic_cast<SPTree*>( MCFs[ k ] );
    if( SPTk ) {
     Index SPn = CeilDiv( SPTk->DestN() , DstLen );
     cIndex_Set DSTn = SPTk->Dests();
     cIndex i = SPTk->DestN() / SPn;
     Index j = SPTk->DestN() % SPn;
     for( ; SPn-- ; ) {
      // divide evenly the SPTk->DestN() destinations of commodity k
      // among the SPn subproblems that it generates
      Index h = i;
      if( j ) { h++; j--; }

      DFSPD[ ws ] = k;
      DFSPDn[ ws ] = h;
      for( Index_Set tD = DFSPDb[ ws++ ] = new Index[ h ] ; h-- ; )
       *(tD++) = *(DSTn++);

      }  // end( for( SPn ) )
     }  // end( if( SPTk ) )
    else {
     DFSPD[ ws ] = k;
     DFSPDn[ ws ] = 0;
     DFSPDb[ ws++ ] = 0; 
     }
    }  // end( for( k ) )
   }  // end( else( subproblems == subsets of O/D pairs ) )- - - - - - - - -
 #endif

 delete[] OptVal;
 OptVal = new FONumber[ NSubPr + 2 ];
 OptVal[ 0 ] = 0;
 VectAssign( OptVal + 1 , - Inf<FONumber>() , NSubPr + 1 );

 }  // end( SetNrSubP() )

/*--------------------------------------------------------------------------*/

void MMCFFlwBase::SetOptEps( const double OE )
{
 OptEps = ( OE == 0 ? CstMltCnst : OE );

 // set the epsilon for costs in the subproblem solver - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 CRow costsK = new CNumber[ NArcs ];
 for( Index k = 0 ; k < NComm ; k++) {
  MCFs[ k ]->MCFCosts( costsK );

  CNumber Cavg = 0;
  cCRow tC = costsK;
  for( Index i = NArcs ; i-- ; tC++ )
   if( *tC < Inf<CNumber>() )
    Cavg += ABS( *tC );

  MCFs[ k ]->SetPar( MCFClass::kEpsCst ,
		     OptEps * max( Cavg / CNumber( NArcs ) , CNumber( 1 ) ) );
  }

 delete[] costsK;
 }

/*--------------------------------------------------------------------------*/

void MMCFFlwBase::SetFsbEps( const double FE )
{
 FsbEps = ( FE == 0 ? FlwMltCnst : FE );

 // set the epsilon for flows in the subproblem solver - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 FRow capsK = new FNumber[ NArcs ];
 for( Index k = 0 ; k < NComm ; k++ ) {
  MCFs[ k ]->MCFUCaps( capsK );

  FNumber avrg = 0;
  cFRow Caps = capsK;
  for( cFRow tCaps = Caps + NArcs ; tCaps > Caps ; ) {
   cFNumber ttC = *(--tCaps);
   if( ttC < Inf<FNumber>() )
    avrg += ttC;
   }

  MCFs[ k ]->SetPar( MCFClass::kEpsFlw , FsbEps * ( avrg / NArcs ) );
  MCFs[ k ]->SetPar( MCFClass::kEpsDfct ,
		     FsbEps * ( avrg / NNodes ) * DfctCnst );
  }

 delete[] capsK;
 }

/*--------------------------------------------------------------------------*/
/*-------------------- METHODS FOR SOLVING THE PROBLEM ---------------------*/
/*--------------------------------------------------------------------------*/

MMCFClass::MMCFStatus MMCFFlwBase::SolveMMCF( void )
{
 if( MMCFt )
  MMCFt->Start();

 Stts = kOK;

 if( WhchSP )  // not the constant subproblem - - - - - - - - - - - - - - - -
  if( ( NSubPr > 1 ) && ( WhchSP <= NSubPr ) )  // a specific subproblem- - -
   SolveSubP( WhchSP );                         //- - - - - - - - - - - - - -
  else          // all problems together- - - - - - - - - - - - - - - - - - -
  {             //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   // solve them- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   for( Index k = 0 ; k < NComm ; )
    SolveFlwK( k++ );

   // now compute the objective function(s) - - - - - - - - - - - - - - - - -

   #if( FlowBase_HAVE_SPT )
    if( DFSPD )  // subproblems can be subsets of O/D pairs
     for( Index wh = 1 ; wh <= NSubPr ; wh++ ) {
      FONumber OV;
      if( DFSPDn[ wh - 1 ] ) {
       SPTree *SPTk = (SPTree*) MCFs[ DFSPD[ wh - 1 ] ];
       OV = SPTk->MCFGetFO( DFSPDn[ wh - 1 ] , DFSPDb[ wh - 1 ] );
       }
      else
       OV = MCFs[ DFSPD[ wh - 1 ] ]->MCFGetFO();

      OptVal[ wh ] = ( OV < Inf<FONumber>() ? OV : Inf<FONumber>() );
      }
    else         // subproblems are sets of commodities
   #endif
     for( Index wh = 1 ; wh <= NSubPr ; wh++ ) {
      OptVal[ wh ] = 0;

      for( Index k = AFSPD[ wh - 1 ] ; k < AFSPD[ wh ] ; k++ ) {
       FONumber OV = MCFs[ k ]->MCFGetFO();
       if( OV < Inf<FONumber>() )
        OptVal[ wh ] += OV;
       else {
        OptVal[ wh ] = Inf<FONumber>();
	break;
	}
       }
      }

   }  // end( else( all problems together )- - - - - - - - - - - - - - - - -
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( MMCFt )
  MMCFt->Stop();

 return( Stts );

 }  // end( MMCFFlwBase::SolveMMCF )

/*--------------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/

MMCFClass::FONumber MMCFFlwBase::GetPVal( void )
{
 FONumber FO;
 if( WhchSP <= NSubPr )
  FO = OptVal[ WhchSP ];
 else {
   if( OptVal[ NSubPr + 1 ] == - Inf<FONumber>() ) {  // sum not computed yet
    OptVal[ NSubPr + 1 ] = 0;

    for( Index wh = 1 ; wh <= NSubPr ; wh++ )
     if( OptVal[ wh ] == - Inf<FONumber>() )
      throw( MMCFException( "GetPVal(): o.f. value not available" ) );
     else
      if( OptVal[ wh ] == Inf<FONumber>() ) {
       OptVal[ NSubPr + 1 ] = Inf<FONumber>();
       break;
       }
      else
       OptVal[ NSubPr + 1 ] += OptVal[ wh ];
    }

   FO = OptVal[ NSubPr + 1 ];
   if( ( WhchSP > NSubPr + 1 ) && ( FO < Inf<FONumber>() ) )
    FO += OptVal[ 0 ];
  }

 return( FO );
 
 }  // end( MMCFFlwBase::GetPVal )

/*--------------------------------------------------------------------------*/

bool MMCFFlwBase::GetPSol( void )
{
 if( ! FSol ) {
  FBse = 0;
  return( false );
  }

 if( ! WhchSP )
  throw( MMCFException( "Error: GetPSol() for 0-th subproblem" ) );

 if( WhchSP <= NSubPr )   // asking for the solution of one subproblem- - - -
 {                        //- - - - - - - - - - - - - - - - - - - - - - - - -
  if( WhchFS == NComm )
   throw( MMCFException( "GetPsol(): aggregated flows not supported yet" ) );

  if( WhchFS < NComm ) {  // ... of one specific commodity- - - - - - - - - -
   bool IsThis = false;

   #if( FlowBase_HAVE_SPT )
   if( DFSPD ) {         // subproblems can be subsets of O/D pairs
     if( DFSPD[ WhchSP - 1 ] == WhchFS )  // and this commodity is the one
      if( DFSPDn[ WhchSP - 1 ] ) {
       SPTree *SPTk = (SPTree*) MCFs[ WhchFS ];
       SPTk->MCFGetX( DFSPDn[ WhchSP - 1 ] , DFSPDb[ WhchSP - 1 ] ,
		      FSol , FBse );
       }
      else
       IsThis = true;
     }
    else                  // subproblems are sets of commodities
   #endif
     if( ( WhchFS >= AFSPD[ WhchSP - 1 ] ) && ( WhchFS < AFSPD[ WhchSP ] ) )
      IsThis = true;

   if( IsThis )           // this commodity is the one
    MCFs[ WhchFS ]->MCFGetX( FSol , FBse );
   else                   // this commodity is not the one
    if( FBse )            // want it in sparse format
     *FBse = Inf<Index>();
    else                  // want it in dense format
     VectAssign( FSol , FNumber( 0 ) , NArcs );
   }
  else {  // the full solution of the subproblem- - - - - - - - - - - - - - -
   #if( FlowBase_HAVE_SPT )
    if( DFSPD ) {         // subproblems can be subsets of O/D pairs
     if( DFSPDn[ WhchSP - 1 ] ) {
      SPTree *SPTk = (SPTree*) MCFs[ DFSPD[ WhchSP - 1 ] ];
      SPTk->MCFGetX( DFSPDn[ WhchSP - 1 ] , DFSPDb[ WhchSP - 1 ] , FSol ,
		     FBse );
      }
     else
      MCFs[ DFSPD[ WhchSP - 1 ] ]->MCFGetX( FSol , FBse );

     if( DFSPD[ WhchSP - 1 ] )  // names must be translated
      if( FBse )
       for( ; *FBse < Inf<Index>() ; FBse++ )
	*FBse += NArcs * DFSPD[ WhchSP - 1 ];
     }
    else                  // subproblems are sets of commodities
   #endif
    {
     if( ( ! FBse ) && ( WhchSP > 1 ) ) {
      VectAssign( FSol , FNumber( 0 ) , NArcs * AFSPD[ WhchSP - 1 ] );
      FSol += NArcs * AFSPD[ WhchSP - 1 ];
      }

     for( Index k = AFSPD[ WhchSP - 1 ] ; k < AFSPD[ WhchSP ] ; k++ ) {
      MCFs[ k ]->MCFGetX( FSol , FBse );

      if( FBse )
       if( k )    // names must be translated
        for( ; *FBse < Inf<Index>() ; FBse++ , FSol++ )
	 *FBse += NArcs * k;
       else
        for( ; *FBse < Inf<Index>() ; FBse++ )
	 FSol++;
      else
       FSol += NArcs;
      }

     if( ( ! FBse ) && ( WhchSP < NSubPr ) )
      VectAssign( FSol, FNumber( 0 ), NArcs * ( NComm - AFSPD[ WhchSP ] ) );

     }  // end( else( subproblems are sets of commodities ) )
   }  // end( else( the full solution of the subproblem ) - - - - - - - - - -
  }
 else  // asking for a full solution- - - - - - - - - - - - - - - - - - - - -
 {     // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( WhchFS == NComm ) {  // ... aggregated- - - - - - - - - - - - - - - - -
   MCFs[ 0 ]->MCFGetX( FSol );

   for( Index k = 0 ; ++k < NComm ; ) {
    cFRow Fk = MCFs[ k ]->MCFGetX();
    if( ! Fk ) {
     MCFs[ k ]->MCFGetX( Ftmp );
     Fk = Ftmp;
     }

    VectSum( FSol , Ftmp , NArcs );
    }

   if( FBse )            // in sparse format
    *( Sparsify( FSol , FBse , NArcs ) ) = Inf<Index>();
   }
  else                   // ... disaggregated - - - - - - - - - - - - - - - -
   if( FBse ) {          // in sparse format
    Index i = 0;
    for( Index k = 0 ; k < NComm ; i += NArcs ) {
     MCFs[ k++ ]->MCFGetX( FSol , FBse );
     while( *FBse < Inf<Index>() ) {
      *(FBse++) += i;
      FSol++;
      }
     }
    }
   else                  // in dense format
    for( Index k = 0 ; k < NComm ; FSol += NArcs )
     MCFs[ k++ ]->MCFGetX( FSol );
 
  }  // end( else( asking for a full solution ) ) - - - - - - - - - - - - - -
     // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 FSol = 0;
 FBse = 0;
 
 return( false );  // only one solution is available

 }  // end( MMCFFlwBase::GetPSol )

/*--------------------------------------------------------------------------*/

bool MMCFFlwBase::GetDSol( void )
{
 if( ! WhchSP )
  throw( MMCFException( "Error: GetDSol() for 0-th subproblem" ) );

 if( NPot )  // get Node Potentials - - - - - - - - - - - - - - - - - - - - -
 {           // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( WhchSP <= NSubPr ) {  // of one subproblem- - - - - - - - - - - - - - -
   if( WhchNP < NComm ) {   // ... of one specific commodity
    bool IsThis = false;

    #if( FlowBase_HAVE_SPT )
     if( DFSPD ) {         // subproblems can be subsets of O/D pairs
      if( DFSPD[ WhchSP - 1 ] == WhchNP )  // and this commodity is the one
       IsThis = true;
      }
     else                  // subproblems are sets of commodities
    #endif
      if( ( WhchNP >= AFSPD[ WhchSP - 1 ] ) && ( WhchNP < AFSPD[ WhchSP ] ) )
       IsThis = true;

    if( IsThis )           // this commodity is the one
     MCFs[ WhchNP ]->MCFGetPi( NPot );
    else
     VectAssign( NPot , CNumber( 0 ) , NNodes );
    }
   else {                  // the full solution of the subproblem
    Index Intk;
    Index Fnlk;

    #if( FlowBase_HAVE_SPT )
     if( DFSPD ) {
      Intk = DFSPD[ WhchSP - 1 ];
      Fnlk = Intk + 1;
      }
     else
    #endif
     {
      Intk = AFSPD[ WhchSP - 1 ];
      Fnlk= AFSPD[ WhchSP ];
      }

    if( Intk ) {
     VectAssign( NPot , CNumber( 0 ) , NNodes * Intk );
     NPot += NNodes * Intk;
     }

    for( Index k = Intk ; k < Fnlk ; NPot += NNodes )
     MCFs[ k++ ]->MCFGetPi( NPot );

    if( Fnlk < NComm )
     VectAssign( NPot , CNumber( 0 ) , NNodes * ( NComm - Fnlk ) );

    }  // end( else( the full solution of the subproblem )
   }
  else  // asking for a full solution - - - - - - - - - - - - - - - - - - - -
   for( Index k = 0 ; k < NComm ; NPot += NNodes )
    MCFs[ k++ ]->MCFGetPi( NPot );

  NPot = 0;

  }  // end( if( get Node Potentials )- - - - - - - - - - - - - - - - - - - -
     // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( RCst )  // get Flow Reduced Costs- - - - - - - - - - - - - - - - - - - -
 {           // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( WhchSP <= NSubPr ) {  // of one subproblem- - - - - - - - - - - - - - -
   if( WhchRC < NComm ) {   // ... of one specific commodity
    bool IsThis = false;

    #if( FlowBase_HAVE_SPT )
     if( DFSPD ) {         // subproblems can be subsets of O/D pairs
      if( DFSPD[ WhchSP - 1 ] == WhchRC )  // and this commodity is the one
       IsThis = true;
      }
     else                  // subproblems are sets of commodities
    #endif
      if( ( WhchRC >= AFSPD[ WhchSP - 1 ] ) && ( WhchRC < AFSPD[ WhchSP ] ) )
       IsThis = true;

    if( IsThis )           // this commodity is the one
     MCFs[ WhchRC ]->MCFGetRC( RCst );
    else
     VectAssign( RCst , CNumber( 0 ) , NArcs );
    }
   else {                  // the full solution of the subproblem
    Index Intk;
    Index Fnlk;

    #if( FlowBase_HAVE_SPT )
     if( DFSPD ) {
      Intk = DFSPD[ WhchSP - 1 ];
      Fnlk = Intk + 1;
      }
     else
    #endif
     {
      Intk = AFSPD[ WhchSP - 1 ];
      Fnlk= AFSPD[ WhchSP ];
      }

    if( Intk ) {
     VectAssign( RCst , CNumber( 0 ) , NArcs * Intk );
     RCst += NArcs * Intk;
     }

    for( Index k = Intk ; k < Fnlk ; RCst += NArcs )
     MCFs[ k++ ]->MCFGetRC( RCst );

    if( Fnlk < NComm )
     VectAssign( RCst , CNumber( 0 ) , NArcs * ( NComm - Fnlk ) );

    }  // end( else( the full solution of the subproblem )
   }
  else  // asking for a full solution - - - - - - - - - - - - - - - - - - - -
   for( Index k = 0 ; k < NComm ; RCst += NArcs )
    MCFs[ k++ ]->MCFGetRC( RCst );

  RCst = 0;

  }  // end( if( get Flow Reduced Costs ) - - - - - - - - - - - - - - - - - -
     // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( MCst ) throw(
  MMCFException( "MMCFFlwBase::GetDSol() called with non-0 MCst" ) );

 if( XtRC ) throw(
  MMCFException( "MMCFFlwBase::GetDSol() called with non-0 XtRC" ) );

 if( XtDV ) throw(
  MMCFException( "MMCFFlwBase::GetDSol() called with non-0 XtDV" ) );

 return( false );  // only one solution is available

 }  // end( MMCFFlwBase::GetDSol )

/*--------------------------------------------------------------------------*/

MMCFClass::FONumber MMCFFlwBase::CostOf( void )
{
 FONumber CO = 0;

 if( FSol ) {
  if( FBse ) throw( MMCFException(
             "MMCFFlwBase::CostOf(): sparse solutions not supported yet" ) );

  if( WhchSP <= NComm ) {
   if( WhchSP ) {
    cCRow Ck = MCFs[ WhchSP - 1 ]->MCFCosts();
    if( ! Ck ) {
     MCFs[ WhchSP - 1 ]->MCFCosts( Cost );
     Ck = Cost;
     }

    CO += ScalarProduct( Ck , FSol + ( WhchSP - 1 ) * NArcs , NArcs );
    }
   }
  else
   if( WhchSP > NComm ) {
    cFRow FS = FSol;
    for( Index k = 0 ; k < NComm ; k++ , FS += NArcs ) {
     cCRow Ck = MCFs[ k ]->MCFCosts();
     if( ! Ck ) {
      MCFs[ k ]->MCFCosts( Cost );
      Ck = Cost;
      }

     CO += ScalarProduct( Ck , FS , NArcs );
     }
    }

  FSol = 0;
  FBse = 0;
  }

 return( CO );

 }  // end( MMCFFlwBase::CostOf )

/*--------------------------------------------------------------------------*/
/*-------------- METHODS FOR READING THE DATA OF THE PROBLEM ---------------*/
/*--------------------------------------------------------------------------*/

void MMCFFlwBase::GetCosts( CRow Csts , cIndex_Set nms , cIndex strt ,
			    Index stp )
{
 if( stp > NComm * NArcs )
  stp = NComm * NArcs;

 if( strt == stp )
  return;

 cIndex lstk = ( stp - 1 ) / NArcs;
 Index k = strt / NArcs;
 Index lk = k * NArcs;

 if( nms ) {
  while( *nms < strt )
   nms++;

  for( ; k < lstk ; k++ , lk += NArcs )
   for( Index h ; ( h = (*nms) - lk ) < NArcs ; nms++ )
    *(Csts++) = MCFs[ k ]->MCFCost( h );

  for( Index h ; ( h = *(nms++) ) < stp ; )
   *(Csts++) = MCFs[ k ]->MCFCost( h - lk );
  }
 else
  if( k == lstk )
   MCFs[ k ]->MCFCosts( Csts , 0 , strt - lk , stp - lk );
  else {
   MCFs[ k++ ]->MCFCosts( Csts , 0 , strt - lk );
   Csts += NArcs - ( strt - lk );
   lk += NArcs;

   for( ; k < lstk ; Csts += NArcs , lk += NArcs )
    MCFs[ k++ ]->MCFCosts( Csts );

   MCFs[ k ]->MCFCosts( Csts , 0 , 0 , stp - lk );
   }

 }  // end( MMCFFlwBase::GetCosts )

/*--------------------------------------------------------------------------*/

void MMCFFlwBase::GetICaps( FRow ICps , cIndex_Set nms , cIndex strt ,
			    Index stp )
{
 if( stp > NComm * NArcs )
  stp = NComm * NArcs;

 if( strt == stp )
  return;

 cIndex lstk = ( stp - 1 ) / NArcs;
 Index k = strt / NArcs;
 Index lk = k * NArcs;

 if( nms ) {
  while( *nms < strt )
   nms++;

  for( ; k < lstk ; k++ , lk += NArcs )
   for( Index h ; ( h = (*nms) - lk ) < NArcs ; nms++ )
    *(ICps++) = MCFs[ k ]->MCFUCap( h );

  for( Index h ; ( h = *(nms++) ) < stp ; )
   *(ICps++) = MCFs[ k ]->MCFUCap( h - lk );
  }
 else
  if( k == lstk )
   MCFs[ k ]->MCFUCaps( ICps , 0 , strt - lk , stp - lk );
  else {
   MCFs[ k++ ]->MCFUCaps( ICps , 0 , strt - lk );
   ICps += NArcs - ( strt - lk );
   lk += NArcs;

   for( ; k < lstk ; ICps += NArcs , lk += NArcs )
    MCFs[ k++ ]->MCFUCaps( ICps );

   MCFs[ k ]->MCFUCaps( ICps , 0 , 0 , stp - lk );
   }

 }  // end( MMCFFlwBase::GetICaps )

/*--------------------------------------------------------------------------*/
/*------------- METHODS FOR ADDING / REMOVING / CHANGING DATA --------------*/
/*--------------------------------------------------------------------------*/

void MMCFFlwBase::ChgCosts( cCRow NwCsts , cIndex_Set nms , cIndex strt ,
			    Index stp )
{
 if( stp > NComm * NArcs )
  stp = NComm * NArcs;

 if( strt == stp )
  return;

 cIndex lstk = ( stp - 1 ) / NArcs;
 Index k = strt / NArcs;
 Index lk = k * NArcs;

 if( nms ) {
  while( *nms < strt ) {
   NwCsts++;
   nms++;
   }

  for( ; k < lstk ; k++ , lk += NArcs )
   for( Index h ; ( h = (*nms) - lk ) < NArcs ; nms++ )
    MCFs[ k ]->ChgCost( h , *(NwCsts++) );

  for( Index h ; ( h = *(nms++) ) < stp ; )
   MCFs[ k ]->ChgCost( h - lk , *(NwCsts++) );
  }
 else
  if( k == lstk )
   MCFs[ k ]->ChgCosts( NwCsts , 0 , strt - lk , stp - lk );
  else {
   MCFs[ k++ ]->ChgCosts( NwCsts , 0 , strt - lk );
   NwCsts += NArcs - ( strt - lk );
   lk += NArcs;

   for( ; k < lstk ; NwCsts += NArcs , lk += NArcs )
    MCFs[ k++ ]->ChgCosts( NwCsts );

   MCFs[ k ]->ChgCosts( NwCsts , 0 , 0 , stp - lk );
   }

 VectAssign( OptVal + 1 , - Inf<FONumber>() , NSubPr + 1 );
 // this should be done much better: only some of the subproblems may
 // actually need to be recomputed!!

 UBnd = Inf<FONumber>();

 }  // end( MMCFFlwBase::ChgCosts )

/*--------------------------------------------------------------------------*/

void MMCFFlwBase::ChgICaps( cFRow NwICps , cIndex_Set nms , cIndex strt ,
			    Index stp )
{
 if( stp > NComm * NArcs )
  stp = NComm * NArcs;

 if( strt == stp )
  return;

 cIndex lstk = ( stp - 1 ) / NArcs;
 Index k = strt / NArcs;
 Index lk = k * NArcs;

 if( nms ) {
  while( *nms < strt ) {
   NwICps++;
   nms++;
   }

  for( ; k < lstk ; k++ , lk += NArcs )
   for( Index h ; ( h = (*nms) - lk ) < NArcs ; nms++ )
    MCFs[ k ]->ChgUCap( h , *(NwICps++) );

  for( Index h ; ( h = *(nms++) ) < stp ; )
   MCFs[ k ]->ChgUCap( h - lk , *(NwICps++) );
  }
 else
  if( k == lstk )
   MCFs[ k ]->ChgUCaps( NwICps , 0 , strt - lk , stp - lk );
  else {
   MCFs[ k++ ]->ChgUCaps( NwICps , 0 , strt - lk );
   NwICps += NArcs - ( strt - lk );
   lk += ( strt - lk );

   for( ; k < lstk ; NwICps += NArcs , lk += NArcs )
    MCFs[ k++ ]->ChgUCaps( NwICps );

   MCFs[ k ]->ChgUCaps( NwICps , 0 , 0 , stp - lk );
   }

 VectAssign( OptVal + 1 , - Inf<FONumber>() , NSubPr + 1 );
 // this should be done much better: only some of the subproblems may
 // actually need to be recomputed!!

 UBnd = Inf<FONumber>();

 }  // end( MMCFFlwBase::ChgICaps )

/*--------------------------------------------------------------------------*/

void MMCFFlwBase::ChgIntVar( cIndex k , bool IntVld , cIndex_Set nms ,
			     Index strt , Index stp  )
{
 if( IntVld && ( k < NComm ) && ( Eps<FNumber>() == FNumber( 0 ) ) )
  throw(
   MMCFException( "MMCFFlwBase::ChgIntVar() for non-integer FNumbers" ) );

 }  // end( MMCFFlwBase::ChgIntVar )

/*--------------------------------------------------------------------------*/

#if( CHGARCS_MMCF == 1 ) /*- - - - - - - - - - - - - - - - - - - - - - - - -*/
                         /*- - - - - - - - - - - - - - - - - - - - - - - - -*/

/* Generic note about [Close/Open]Arcs(): if UBnd is a valid upper bound on
   the largest possible value of any feasible solution, then it remains valid
   when arcs are closed (the set of feasible solutions can only shrink).
   Thus, there is no need to update UBnd when closing arcs. Therefore, there
   is also no need of updating it when *opening* arcs, provided that *UBnd is
   computed without taking into account closed/open arcs*. */

 void MMCFFlwBase::CloseArcs( cIndex_Set whch )
 {
  for( Index k = NComm ; k-- ; ) {
   cIndex_Set tw = whch;
   for( Index h ; ( h = *(tw++) ) < Inf<Index>() ; )
    MCFs[ k ]->CloseArc( h );
   }
  
  VectAssign( OptVal + 1 , - Inf<FONumber>() , NSubPr + 1 );
  // this should be done much better: only some of the subproblems may
  // actually need to be recomputed!!

  }  // end( MMCFFlwBase::CloseArcs )

/*--------------------------------------------------------------------------*/

 void MMCFFlwBase::OpenArcs( cIndex_Set whch )
 {
  for( Index k = NComm ; k-- ; ) {
   cIndex_Set tw = whch;
   for( Index h ; ( h = *(tw++) ) < Inf<Index>() ; )
    MCFs[ k ]->OpenArc( h );
   }
  
  VectAssign( OptVal + 1 , - Inf<FONumber>() , NSubPr + 1 );
  // this should be done much better: only some of the subproblems may
  // actually need to be recomputed!!

  }  // end( MMCFFlwBase::OpenArcs )

#elif( CHGARCS_MMCF > 1 )  /*- - - - - - - - - - - - - - - - - - - - - - - -*/
                           /*- - - - - - - - - - - - - - - - - - - - - - - -*/

 void MMCFFlwBase::CloseArcs( cIndex_Set whch )
 {
  Bool_Vec tCAK = ClsdAK;
  for( Index k = NComm ; k-- ; tCAK += NArcs ) {
   cIndex_Set tw = whch;
   for( Index h ; ( h = *(tw++) ) < Inf<Index>() ; )
    if( ! tCAK[ h ] )
     MCFs[ k ]->CloseArc( h );
   }
  
  VectAssign( OptVal + 1 , - Inf<FONumber>() , NSubPr + 1 );
  // this should be done much better: only some of the subproblems may
  // actually need to be recomputed!!

  }  // end( MMCFFlwBase::CloseArcs )

/*--------------------------------------------------------------------------*/

 void MMCFFlwBase::OpenArcs( cIndex_Set whch )
 {
  Bool_Vec tCAK = ClsdAK;
  for( Index k = NComm ; k-- ; tCAK += NArcs ) {
   cIndex_Set tw = whch;
   for( Index h ; ( h = *(tw++) ) < Inf<Index>() ; )
    if( ! tCAK[ h ] )
     MCFs[ k ]->OpenArc( h );
   }
  
  VectAssign( OptVal + 1 , - Inf<FONumber>() , NSubPr + 1 );
  // this should be done much better: only some of the subproblems may
  // actually need to be recomputed!!

  }  // end( MMCFFlwBase::OpenArcs )

/*--------------------------------------------------------------------------*/

 void MMCFFlwBase::CloseArcs( cIndex k , cIndex_Set whch )
 {
  Bool_Vec tCAK = ClsdAK + k * NArcs;
  cIndex_Set tw = whch;
  for( Index h ; ( h = *(tw++) ) < Inf<Index>() ; )
   if( ! tCAK[ h ] ) {
    MCFs[ k ]->CloseArc( h );
    tCAK[ h ] = true;
    }
  
  VectAssign( OptVal + 1 , - Inf<FONumber>() , NSubPr + 1 );
  // this should be done much better: only some of the subproblems may
  // actually need to be recomputed!!

  }  // end( MMCFFlwBase::CloseArcs( k ) )

/*--------------------------------------------------------------------------*/

 void MMCFFlwBase::OpenArcs( cIndex k , cIndex_Set whch )
 {
  Bool_Vec tCAK = ClsdAK + k * NArcs;
  cIndex_Set tw = whch;
  for( Index h ; ( h = *(tw++) ) < Inf<Index>() ; )
   if( tCAK[ h ] ) {
    MCFs[ k ]->OpenArc( h );
    tCAK[ h ] = false;
    }
  
  VectAssign( OptVal + 1 , - Inf<FONumber>() , NSubPr + 1 );
  // this should be done much better: only some of the subproblems may
  // actually need to be recomputed!!

  }  // end( MMCFFlwBase::OpenArcs( k ) )

#endif  /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
        /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

/*--------------------------------------------------------------------------*/
/*------------------------------ DESTRUCTOR --------------------------------*/
/*--------------------------------------------------------------------------*/

MMCFFlwBase::~MMCFFlwBase()
{
 delete[] OptVal;

 delete[] AFSPD;

 #if( CHGARCS_MMCF > 1 )
  delete[] ClsdAK;
 #endif

 for( Index i = NComm ; i-- ; )
  delete MCFs[ i ];

 delete[] Ftmp;

 delete[] Cost;

 delete[] MCFs;

 }  // end( ~MMCFFlwBase )

/*--------------------------------------------------------------------------*/
/*-------------------------- PROTECTED METHODS -----------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*-------------------------- PRIVATE METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/

inline void MMCFFlwBase::SolveFlwK( cIndex k )
{
 MCFs[ k ]->SolveMCF();

 if( MCFs[ k ]->MCFGetStatus() ) {
  OptVal[ k ] = Inf<FONumber>();
  Stts = kUnfeasible;
  }
 /*!!
 else {
  MCFs[ k ]->CheckPSol();
  MCFs[ k ]->CheckDSol();
  }
  !!*/
 }  // end( SolveFlwK )

/*--------------------------------------------------------------------------*/

inline void MMCFFlwBase::SolveSubP( cIndex wh )
{
 if( OptVal[ wh ] > - Inf<FONumber>() ) {  // subproblem already solved
  if( OptVal[ wh ] == Inf<FONumber>() )    // nothing to do
   Stts = kUnfeasible;
  return;
  }

 #if( FlowBase_HAVE_SPT )
  FONumber OV;
  if( DFSPD )  // subproblems can be subsets of O/D pairs - - - - - - - - - -
  {            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   Index k = DFSPD[ wh - 1 ];
   if( DFSPDn[ wh - 1 ] ) {  // ... and this one is
    // check that this SPT has not been already solved- - - - - - - - - - - -
    // WARNING: we assume that all subproblems corresponding to the same
    // commodity are consecutive

    for( Index whm = wh - 1 ; whm ; whm-- ) {
     if( DFSPD[ whm - 1 ] != k )
      break;

     if( OptVal[ whm ] > - Inf<FONumber>() ) {
      k = Inf<Index>();
      break;
      }
     }

    if( k < Inf<Index>() )
     for( Index whp = wh + 1 ; whp <= NSubPr ; whp++ ) {
      if( DFSPD[ whp - 1 ] != k )
       break;

      if( OptVal[ whp ] > - Inf<FONumber>() ) {
       k = Inf<Index>();
       break;
       }
      }

    }  // end( if( this one is a part of a SPT ) )

   if( k < Inf<Index>() )  // if not, solve it- - - - - - - - - - - - - - - -
    SolveFlwK( k );

   // compute the objective function- - - - - - - - - - - - - - - - - - - - -

   if( DFSPDn[ wh - 1 ] ) {
    SPTree *SPTk = (SPTree*) MCFs[ DFSPD[ wh - 1 ] ];
    OV = SPTk->MCFGetFO( DFSPDn[ wh - 1 ] , DFSPDb[ wh - 1 ] );
    }
   else
    OV = MCFs[ DFSPD[ wh - 1 ] ]->MCFGetFO();

   OptVal[ wh ] = ( OV < Inf<FONumber>() ? OV : Inf<FONumber>() );
   }
  else         // subproblems are sets of commodities - - - - - - - - - - - -
               // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 #endif
  {
   OptVal[ wh ] = 0;

   for( Index k = AFSPD[ wh - 1 ] ; k < AFSPD[ wh ] ; k++ ) {
    SolveFlwK( k );

    if( OptVal[ wh ] < Inf<FONumber>() ) {
     FONumber OV = MCFs[ k ]->MCFGetFO();
     if( OV < Inf<FONumber>() )
      OptVal[ wh ] += OV;
     else
      OptVal[ wh ] = Inf<FONumber>();
     }
    }
   }

 }  // end( SolveSubP )

/*--------------------------------------------------------------------------*/
/*------------------------ End File MMCFFlwB.C -----------------------------*/
/*--------------------------------------------------------------------------*/

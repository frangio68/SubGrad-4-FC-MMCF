/*--------------------------------------------------------------------------*/
/*---------------------------- File KnpsFiOrcl.C ---------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the KnpsFiOrcl class which implements the FiOracle
 * interface and solves the knapsack Lagrangian relaxation of the
 * (Fixed-Charge) Multicommodity Min Cost Flow Problem ((FC-)MMCF).
 *
 * \version 1.12
 *
 * \date 05 - 09 - 2014
 *
 * \author Antonio Frangioni \n
 *         Operations Research Group \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Enrico Gorgone \n
 *         Operations Research Group \n
 *         Dipartimento di Elettronica Informatica e Sistemistica \n
 *         Universita' della Calabria \n
 *
 * \author Vitor Barbosa \n
 *         Algoritmi Research Unit \n
 *         Universidade do Minho \n
 *
 * \author Filipe Alvelos \n
 *         Grupo de Optimizacao e Investigacao Operacional \n
 *         Departamento de producao e Sistemas \n
 *         Universidade do Minho \n
 *
 * Copyright &copy 2000 - 2012 by Antonio Frangioni, Enrico Gorgone
 */
/*--------------------------------------------------------------------------*/
/*--------------------------- IMPLEMENTATION -------------------------------*/
/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "KnpsFiOrcl.h"

#include "OPTvect.h"

#include <stdlib.h>

/*--------------------------------------------------------------------------*/
/*-------------------------------- MACROS ----------------------------------*/
/*--------------------------------------------------------------------------*/

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
/*--------------------- IMPLEMENTATION OF KnpsFiOrcl -----------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- PUBLIC METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/
/*---------------------------- CONSTRUCTOR ---------------------------------*/
/*--------------------------------------------------------------------------*/


KnpsFiOrcl::KnpsFiOrcl( Graph *g , istream *iStrm )
           :
           FiOracle() , MMCFClass()
{
 Index agg;
 DfltdSfInpt( iStrm , agg , Index( 0 ) );
 DfltdSfInpt( iStrm , Eps , HpNum( 1e-6 ) );
 DfltdSfInpt( iStrm , DoAggregation , bool( false ) );

 // allocate memory for instance data  - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 NArcs = g->NrArcs();
 NNodes = g->NrNodes();
 NComm = g->NrComm();

 // if the "extra" variables in graph are less than Nr arcs throw exception

 if( g->NrExtraVars() == NArcs )
  XtrCosts = new CNumber[ NArcs ];  // f_{ij} : original costs
 else
  throw( NDOException( "KnpsFiOrcl::NewGi: this component is not declared" ) );

 TotCap = new FNumber[ NArcs ];              // u_{ij} : mutual capacities
 Deficits = new FNumber[ NComm * NNodes ];  // b_i^k : original deficits
 Costs = new CNumber[ NArcs * NComm ];      // c_{ij}^k : original unit costs
 Capacities = new CNumber[ NArcs * NComm ]; // u_{ij}^k : original individual
                                                  //   capacities
 Startn = new Index[ NArcs ];
 Endn = new Index[ NArcs ];

 // read instance data - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 for( Index j = 0 ; j < NArcs ; j++ ) {
  XtrCosts[ j ] = g->CostKJ( NComm , j );
  TotCap[ j ] = g->TotalCapacityJ( j );
  }
  
 for( Index k=0 ; k < NComm ; k++ ) {
  for( Index j = 0 ; j < NArcs ; j++ ) {
   Costs[ k * NArcs + j ] = g->CostKJ( k , j );
   Capacities[ k * NArcs + j ] = g->CapacityKJ( k , j );
   }
  for( Index n = 0 ; n < NNodes ; n++ )
   Deficits[ k*NNodes + n ] = g->DeficitKJ( k , n );
  }

 if( g->NamesStartFrom1()==true)
  for( Index j = 0 ; j < NArcs ; j++ ) {
   Startn[ j ]= g->StartNJ( j ) - 1;
   Endn[ j ] = g->EndNJ( j ) - 1;
   }
 else
  for( Index j = 0 ; j < NArcs ; j++ ) {
   Startn[ j ] = g->StartNJ( j );
   Endn[ j ] = g->EndNJ( j );
   }

 // define the aggregation type  - - - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 LsHasChgd = 0;
 SetAggregate( agg );

 // allocate memory for primal solution  - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 FSolution = new FNumber[ NArcs * NComm ];
 SolWFi = Inf<Index>();

 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 NumVar = NNodes * NComm;
 SGBse1 = new Index[ NumVar + 1 ];
                                // SolvedP[ wFi -1 ] = false indicates
 SlvP = new bool[ NArcs ];     // that the primal solution of wFi
 FiLmbd = new HpNum[ NArcs ];  // component is not computed yet

 // allocate memory for sorting of costs - - - - - - - - - - - - - - - - - - -

 NwCsts = new CNumber[ NArcs * NComm ];
 IndNwCsts = new Index_Set[ NArcs ];
 for( Index j = 0 ; j < NArcs ; j++ )
  IndNwCsts[ j ] = new Index[ NComm + 1 ];
 QSStck = new int[ NComm ];

 // other initializations  - - - - - - - - - - - - - - - - - - - - - - - - - -

 OldWFi = 0;
 MaxName = 0;

 FindGlobalLipschitz();
 LowerBound = - Inf<HpNum>();

 } // end( KnpsFiOrcl::KnpsFiOrcl )

/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

void KnpsFiOrcl::SetFiLog( ostream *outs , const char lvl )
{
 FiOracle::SetFiLog( outs , lvl );

 #if( LOG_FI )
  if( FiLLvl > 1 )
   *FiLog << endl << " Aggrgtd = " << ( Aggrgtd? "True ": "False " ) << endl
    << " ~ EpsFi = " << Eps  << endl;
 #endif

 } // end ( FlwFiOrcl::SetFiLog )

/*--------------------------------------------------------------------------*/

void KnpsFiOrcl::SetMaxName( cIndex MxNme )
{
 if( MxNme > MaxName ) {

  if( MaxName == 0 ) {  // allocs memory to keap MxNme names in memory - - - -

   OldWFi = new Index[ MxNme ];
   OldFSols = new FRow[ MxNme ];
   for( Index i = 0 ; i < MxNme ; i++ )
    OldFSols[ i ] = 0;

   MaxName = MxNme;
   }
  else { // reallocs memory to keap MxNme names in memory - - - - - - - - - -

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

   MaxName = MxNme;
   }

  } // end ( MxNme > MaxName )
 else {

  if( MxNme == 0 ) {

   if( OldWFi ) {
    delete[] OldWFi;
    OldWFi = 0;
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
   throw( NDOException( "KnpsFiOrcl::SetMaxName: MaxName reduction "
		   "only is available to 0." ) );
  }  

 } // end ( KnpsFiOrcl::SetMaxName )

/*--------------------------------------------------------------------------*/

void KnpsFiOrcl::SetAggregate( bool agg )
{
 if( LsHasChgd )
  delete[] LsHasChgd;

 if( agg == true ) {
  Aggrgtd = true;
  LsHasChgd = 0;
  }
 else { // LsHasChgd[ wFi - 1 ] = true indicates that the subgradient of wFi
  Aggrgtd = false;   //  component is not computed yet
  LsHasChgd = new bool[ NArcs ];
  }

 } // end( KnpsFiOrcl::SetAggregate )

/*--------------------------------------------------------------------------*/
/*-------------- METHODS FOR READING THE DATA OF THE PROBLEM ---------------*/
/*--------------------------------------------------------------------------*/

Index KnpsFiOrcl::GetNrFi( void ) const
{
 return(  Aggrgtd ? Index( 1 ) : NArcs );

 } // end ( KnpsFiOrcl::GetNrFi )

/*--------------------------------------------------------------------------*/
/*----------------------- METHODS FOR CHANGING DATA ------------------------*/
/*--------------------------------------------------------------------------*/

void KnpsFiOrcl::SetLambda( cLMRow Lmbd ) {

 // update Lambda  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 FiOracle::SetLambda( Lmbd );

 if( !Aggrgtd )
  VectAssign( LsHasChgd , true , NArcs );
 VectAssign( SlvP , false , NArcs );

 } // end ( FlwFiOrcl::SetLambda )

/*--------------------------------------------------------------------------*/

bool KnpsFiOrcl::SetPrecision( HpNum EpsA )
{
 Eps = EpsA;
 return false;

 } // end ( KnpsFiOrcl::SetPrecision )

/*--------------------------------------------------------------------------*/
/*-------------------- METHODS FOR SOLVING THE PROBLEM ---------------------*/
/*--------------------------------------------------------------------------*/

HpNum KnpsFiOrcl::Fi( cIndex wFi )
{
 if( Fit )
  Fit->Start();

 HpNum out = 0;

 // value of the linear ( affine ) 0-th component of Fi- - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( (wFi == 0) || ( wFi == Inf<Index>() ) ) {
  if( LamBase )
   for( Index i = 0 ; i < LamBDim ; i++ )
    out += Lambda[ i ] * Deficits[ LamBase[ i ] ];
  else
   for( Index i = 0 ; i < NumVar ; i++ )
	out += Lambda[ i ] * Deficits[ i ];
  }

 // evaluate the other components  - - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( wFi != 0 ) {
  if( Aggrgtd || wFi > NArcs ) {
   for( Index j = 0; j < NArcs ; j++ )
	if( !SlvP[ j ] ) {
	 out += solve( j + 1 );
	 SlvP[ j ] = true;
	 }
	else
	 out += FiLmbd[ j ];
   }
  else {
   if( !SlvP[ wFi - 1 ] ) {
    out += solve( wFi );
    SlvP[ wFi - 1 ] = true;
    }
   else
    out += FiLmbd[ wFi-1 ];
   }
  }

 if( Fit )
  Fit->Stop();

 return( -out );

 } // end ( KnpsFiOrcl::Fi )

/*--------------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/

bool KnpsFiOrcl::NewGi( cIndex wFi )
{
 bool FreewFi = true;  // if FreewFi is true, it is possible produce a "new"
	                   // item for the point Lambda
 SolWFi = wFi;        // tells to witch subproblem belongs last gi created

 // some exceptions- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( wFi > NArcs && wFi < Inf<Index>() )
  throw( NDOException( "KnpsFiOrcl::NewGi: this component is not declared" ) );

 if( wFi == 0 )
  return false;

 // asks for an aggregated (epsilon-) subgradient  - - - - - - - - - - - - - -
 //  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( Fit )
  Fit->Start();

 if( Aggrgtd ) {
  for( Index j = 0; j < NArcs ; j++ )
   if( !SlvP[ j ] )
	throw( NDOException( "KnpsFiOrcl::NewGi: Fi() hasn't been called yet" ) );
  FreewFi = LHasChgd;
  }
 else {
  if( wFi <= NArcs ) {
   FreewFi = LsHasChgd[ wFi - 1 ];
   if( !SlvP[ wFi - 1 ] ) {
	solve( wFi );
    SlvP[ wFi - 1 ] = true;
    }
   }
  else
   for( Index j = 0; j < NArcs ; j++ ) {
    if( !LsHasChgd[ j ] ) {
     FreewFi = false;
     break;
     }
    if( !SlvP[ j ] )
  	 throw( NDOException( "KnpsFiOrcl::NewGi: Fi() hasn't been called yet" ) );
    }
  }

 if( Fit )
  Fit->Stop();

 return FreewFi;
 
 } // end ( KnpsFiOrcl::NewGi )

/*--------------------------------------------------------------------------*/

Index KnpsFiOrcl::GetGi( SgRow SubG , cIndex_Set &SGBse , cIndex Name,
		cIndex strt , Index stp )
{
 // GetMaxNZ( wFi ) return MaxNumVar, assuming that the maximum number of
 //	nonzeroes is not known in advance

 Index k = 0;  // the number of nonzero of the item

 // some exceptions  - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( stp <= strt )
  throw( NDOException( "KnpsFiOrcl::GetGi: stop <= start" ) );

 if(  strt != 0  || stp < NumVar )
  throw( NDOException( "KnpsFiOrcl::GetGi: this should not happen" ) );

 if( Name < MaxName && !OldWFi[ Name ] )
  throw( NDOException( "KnpsFiOrcl::GetGi: Subgradient asked do not exist" ) );

 if( stp > NumVar )
  stp = NumVar;

 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( Fit )
  Fit->Start();

 VectAssign( SubG , SgNum( 0 ) , NumVar );

 if( Name == MaxName ) {

  // Name == MaxName. The required information is about the constant
  // subgradient of the linear 0-th component of Fi  - - - - - - - - - - - - -

  // *** aggiustare come in FlowFi *** potrebbe essere che nella generazione dinamica
  // la componente zero viene chiesta poco alla volta!!!! ????

  for(Index j = 0 ; j < NumVar ; j++ )
   SubG[ j ] =- Deficits[ j ];
  SGBse = 0;
  k = NumVar;
  }

 else {

  if( Name > MaxName ) {

   // asks for an aggregated (epsilon-) subgradient
   //  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   if( SolWFi > NArcs || Aggrgtd ) {

	for( Index j = 0 ; j < NArcs ; j++ )
	 for( Index k = 0 ; k < NComm ; k++ ) { // it is excluded the constant part
      SubG[ Startn[ j ] + k*NNodes ] -= FSolution[ j + NArcs*k ];
      SubG[ Endn[ j ] +  k*NNodes ] += FSolution[ j + NArcs*k ];
	  }

    // asks for the constant subgradient of the linear 0-th component of Fi
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    if( SolWFi == Inf<Index>() )
     for( Index j = 0 ;j < NumVar ; j++ )
      SubG[ j ] -= Deficits[ j ];

    }  // end ( LastWFi >  NArcs )
   else {

	for( Index k = 0 ; k < NComm ; k++ ) {
	 SubG[ Startn[ SolWFi - 1 ] + k*NNodes ] -= FSolution[ (SolWFi - 1) + NArcs*k ];
	 SubG[ Endn[ SolWFi - 1 ] + k*NNodes ] += FSolution[ (SolWFi - 1) + NArcs*k ];
	 }

    } // end if( 1 <= LastWFi <= NArcs )

   } // end ( if:  Name > MaxName )
  else {

   // asks for an aggregated (epsilon-) subgradient
   //  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   if( OldWFi[ Name ] >  NArcs || Aggrgtd ) {

	for( Index j = 0 ; j < NArcs ; j++ )
	 for( Index k = 0 ; k < NComm ; k++ ) { // it is excluded the constant part
      SubG[ Startn[ j ] + k*NNodes ] -= OldFSols[ Name ][ j + NArcs*k ];
	  SubG[ Endn[ j ] +  k*NNodes ] += OldFSols[ Name ][ j + NArcs*k ];
	  }

	// asks for the constant subgradient of the linear 0-th component of Fi
	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	if( SolWFi == Inf<Index>() )
	 for( Index j = 0 ;j < NumVar ; j++ )
	  SubG[ j ] -= Deficits[ j ];

	} // end ( OldWFi[ Name ] >  NArcs )
   else {
	for( Index k = 0 ; k < NComm ; k++ ) {
	 Index arc = OldWFi[ Name ] - 1;
	 SubG[ Startn[ arc ] + k*NNodes ] -= OldFSols[ Name ][ k ];
	 SubG[ Endn[ arc ] + k*NNodes ] += OldFSols[ Name ][ k ];
	 }

    } // end if( 1 <= OldWFi[ Name ] <= NArcs )
   } // end ( else: Name < MaxName )

  #if SPF_SUB
  Index_Set SGBse2 = Sparsify( SubG , SGBse1 , NumVar );
  *SGBse2 = Inf<Index>();
  k = SGBse2 - SGBse1;
  SGBse = SGBse1;
  #else
   k = NumVar;
   SGBse = 0;
  #endif

  } // end (  Name != MaxName ) )

 if( Fit )
  Fit->Stop();

 // tell that a new call to item is not possible  - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( Aggrgtd )
  LHasChgd = false;
 else {
  if( SolWFi <= NArcs )
   LsHasChgd[ SolWFi - 1 ] = false;
  else
   for( Index i = 0; i < NArcs ; i++ )
    LsHasChgd[ i ] = false;
  }

 return( k );

 } // end ( KnpsFiOrcl::GetGi )

/*--------------------------------------------------------------------------*/

HpNum KnpsFiOrcl::GetVal( cIndex Name )
{
 if( Name < MaxName )
  throw( NDOException( "KnpsFiOrcl::GetVal: Not implemented yet" ) );

 return( 0 );
 }  // end ( KnpsFiOrcl::GetVal )

/*--------------------------------------------------------------------------*/

void KnpsFiOrcl::SetGiName( cIndex Name )
{ // copy the primal solution in the position called Name into the structure
  // used for maintaining the previous primal solutions - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( DoAggregation ) // don't do anything if dynamic variables generation is not
  CopyGi( Name , SolWFi , FSolution ); // allowed

 } // end ( KnpsFiOrcl::SetGiName )

/*--------------------------------------------------------------------------*/

HpNum KnpsFiOrcl::GetLowerBound( cIndex wFi )
{
 #if LowB_FI
  if( wFi > GetNrFi() )
   return( LowerBound );
  else
   return( - Inf<HpNum>() );
  // *** da rivedere ?!?!?!? ***
  /*Index j;
  HpNum temp = g->UpperBound();
  if( temp != HpINF ) {
   for( j = 0 ; j < NArcs ; j++ )
    temp += OrigXtrCosts[ j ];
   }
  return( - temp );*/
 #else
  return( - Inf<HpNum>() );
 #endif

 }  // end ( KnpsFiOrcl::GetLowerBound )

/*--------------------------------------------------------------------------*/

void KnpsFiOrcl::SetLowerBound( HpNum lowB ) {

 LowerBound = lowB;

 } // end ( FlwFiOrcl::SetLowerBound )

/*--------------------------------------------------------------------------*/

FiOracle::FiStatus KnpsFiOrcl::GetFiStatus( Index wFi )
{
 return kFiNorm;
 }

/*--------------------------------------------------------------------------*/
/*-------------- METHODS FOR READING THE DATA OF THE PROBLEM ---------------*/
/*--------------------------------------------------------------------------*/

HpNum KnpsFiOrcl::GetGlobalLipschitz( cIndex wFi )
{
 if( wFi == Inf<Index>() )
  return( LpshCnst );
 else
  return( Inf<HpNum>() );

 } // end ( KnpsFiOrcl::GetGlobalLipschitz )

/*--------------------------------------------------------------------------*/
/*------------- METHODS FOR ADDING / REMOVING / CHANGING DATA --------------*/
/*--------------------------------------------------------------------------*/

void KnpsFiOrcl::Deleted( cIndex i )
{
 if( i < MaxName )	{

  delete[] OldFSols[ i ];
  OldFSols[ i ] = 0;

  }
 else
  SetMaxName( 0 );

 } // end ( KnpsFiOrcl::Deleted )

/*--------------------------------------------------------------------------*/

void KnpsFiOrcl::Aggregate( cHpRow Mlt , cIndex_Set NmSt , cIndex Dm ,
		cIndex NwNm )
{
 if( !DoAggregation ) // don't do anything if dynamic variables generation is not
  return;     // allowed

 // define the data structure for saving the aggregate primal solution
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 Index wFi;

 FRow FlwSol = new FNumber[ NArcs * NComm ];
 VectAssign( FlwSol , FNumber(0) , NArcs * NComm );

 if( NmSt ) {

   wFi = OldWFi[ NmSt[ 0 ] ];

  // some exceptions - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  for( Index i = 0 ; i < Dm ; i++ )
   if( wFi != OldWFi[ NmSt[ i ] ] )
    throw( NDOException( "KnpsFiOrcl::Aggregate: aggregation of different "
    		"subgradients" ) );

  if( wFi == 0 )
   throw( NDOException( "KnpsFiOrcl::Aggregate: this should not happen. " ) );

  // compute aggregation for flow variables: - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( wFi > NArcs || Aggrgtd )
   for( Index i = 0; i < Dm; i++ )
    for( Index j = 0; j < NArcs * NComm; j++)
     FlwSol[ j ] += Mlt[ i ] * OldFSols[ NmSt[ i ] ][ j ];
  else
   for( Index i = 0; i < Dm; i++ )
	for( Index k = 0; k < NComm ; k++)
	 FlwSol[ (wFi-1) + k*NArcs ] += Mlt[ i ] * OldFSols[ NmSt[ i ] ][ k ];

  }
 else {

  // some exceptions  - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  //  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  wFi = OldWFi[ 0 ];

  if( wFi != Inf<Index>() )
   throw( NDOException( "KnpsFiOrcl::Aggregate: this should not happen." ) );

  for( Index i = 0 ; i < Dm ; i++ )
   if( wFi != OldWFi[ i ] )
    throw( NDOException( "KnpsFiOrcl::Aggregate: aggregation of different "
    	 "subgradients." ) );

  // compute aggregation for flow variables:   - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  for( Index i = 0; i < Dm; i++ )
   for( Index j = 0; j < NArcs * NComm; j++)
    FlwSol[ j ] += Mlt[ i ] * OldFSols[ i ][ j ];
  }

 // set the name of the new item   - - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 CopyGi( NwNm , wFi , FlwSol );

 delete[] FlwSol;
 } // end ( FlwFiOrcl::Aggregate )

/*--------------------------------------------------------------------------*/
/*------------------------------ DESTRUCTOR --------------------------------*/
/*--------------------------------------------------------------------------*/

KnpsFiOrcl::~KnpsFiOrcl()
{
 // delete graph varaibles  - - - - - - - - - - - - - - - - - - - - - - - - - -

 delete [] Startn;
 delete [] Endn;

 delete[] Costs;
 delete[] XtrCosts;
 delete[] TotCap;
 delete[] Capacities;
 delete[] Deficits;

 // delete primal solution  - - - - - - - - - - - - - - - - - - - - - - - - - -

 delete[] FSolution;
 SetMaxName( 0 );

 // delete memory for sorting - - - - - - - - - - - - - - - - - - - - - - - - -

 delete[] QSStck;
 delete[] NwCsts;
 for( Index j = 0 ; j < NArcs ;  j++ )
  delete[] IndNwCsts[ j ];
 delete[] IndNwCsts;

 // delete other variables   - - - - - - - - - - - - - - - - - - - - - - - - - -

 delete[] LsHasChgd;
 delete[] SGBse1;
 delete[] SlvP;
 delete[] FiLmbd;

} // end ( KnpsFiOrcl::~KnpsFiOrcl )

/*--------------------------------------------------------------------------*/
/*-------------------------- PRIVATE METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/

HpNum KnpsFiOrcl::solve( cIndex wFi )
{

 HpNum out = 0; // optimum value of the knapsack problem
 FNumber rhs;   // capacity of knapsack constraint
 HpRow tFSol = new HpNum[ NComm ];  // temporary flow solution

 if( wFi == Inf<Index>() ) {

  if( LamBase ) {

   LMRow tempL = new LMNum[ NumVar ];
   VectAssign( tempL , Lambda , LamBDim );
   Densify( tempL , LamBase , LamBDim , NumVar );

   for( Index j = 0 ; j < NArcs ; j++ ) {
    Index_Set indCst = IndNwCsts[ j ];
    for( Index k = 0 ; k < NComm ; k++  ) {
     Index kj = k * NArcs + j;
     NwCsts[ kj ] = Costs[ kj ] + tempL[ k * NNodes + Startn[ j ] ]
               - tempL[ k * NNodes + Endn[ j ] ];
     if( NwCsts[ kj ] < - Eps )
      *( indCst++ ) =  kj;
     }
    *indCst = Inf<Index>();
    }
    delete[] tempL;
   }
  else
   for( Index j = 0 ; j < NArcs ; j++ ) {
    Index_Set indCst = IndNwCsts[ j ];
    for( Index k = 0 ; k < NComm ; k++  ) {
     Index kj = k * NArcs + j;
     NwCsts[ kj ] = Costs[ kj ] + Lambda[ k * NNodes + Startn[ j ] ]
             - Lambda[ k * NNodes + Endn[ j ] ];
     if( NwCsts[ kj ] < - Eps )
      *( indCst++ ) =  kj;
     }
    *indCst = Inf<Index>();
    }

  // sort NwCsts - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  for( Index j = 0 ; j < NArcs ; j++ )
   quicksort( j );

  HpNum tout;
  for( Index j = 0 ; j < NArcs ; j++ ) {

   rhs = TotCap[ j ];
   VectAssign( tFSol , HpNum(0) , NComm );
   tout = 0;

   // solve the flow part - - - - - - - - - - - - - - - - - - - - - - - - - - -
   // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   for( Index_Set tI = IndNwCsts[ j ] ; *tI < Inf<Index>() ; tI++  ) { // examine the
	Index comm = *tI / NArcs;   // vector of the ordered new costs and pick ones
	                                 // corresponding to the arc wFi -1
	if( NwCsts[ *tI ] < -Eps &&  rhs > Eps ) {
	 tFSol[ comm ] = min( Capacities[ *tI ] , rhs );
	 rhs -= tFSol[ comm ];
	 tout += tFSol[ comm ] * NwCsts[ *tI ];
	 }
	}

   // solve the fixed charge part - - - - - - - - - - - - - - - - - - - - - - -
   // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   if( ( XtrCosts[ j ] + tout ) + Eps < 0 ) {
	tout += XtrCosts[ j ];
	for( Index k = 0 ; k < NComm ; k++ )
	 FSolution[ j + NArcs*k ] = tFSol[ k ];
	out += tout;
	}
   else {
	for( Index k = 0 ; k < NComm ; k++ )
	 FSolution[ j + NArcs*k ] = 0;
	tout = 0;
	}
   FiLmbd[ j ] = tout;
   }

  } // end if( wFi = Inf<Index>() )
 else {
  if( LamBase ) {

   LMRow tempL = new LMNum[ NumVar ];
   VectAssign( tempL , Lambda , LamBDim );
   Densify( tempL , LamBase , LamBDim , NumVar );

   Index_Set indCst = IndNwCsts[ wFi -1 ];
   for( Index k = 0 ; k < NComm ; k++  ) {
    Index kj = k * NArcs + ( wFi -1 ) ;
    NwCsts[ kj ] = Costs[ kj ] + tempL[ k * NNodes + Startn[ wFi -1 ] ]
                - tempL[ k * NNodes + Endn[ wFi -1 ] ];
    if( NwCsts[ kj ] < - Eps )
     *( indCst++ ) =  kj;
    }
   *indCst = Inf<Index>();
   delete[] tempL;
   }
  else {
   Index_Set indCst = IndNwCsts[ wFi -1 ];
   for( Index k = 0 ; k < NComm ; k++  ) {
    Index kj = k * NArcs + ( wFi -1 );
    NwCsts[ kj ] = Costs[ kj ] + Lambda[ k * NNodes + Startn[ wFi -1 ] ]
              - Lambda[ k*NNodes + Endn[ wFi -1 ] ];
    if( NwCsts[ kj ] < - Eps )
      *( indCst++ ) =  kj;
    }
   *indCst = Inf<Index>();
   }

  // sort NwCsts  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  quicksort( wFi -1 );

  rhs = TotCap[ wFi - 1 ];
  VectAssign( tFSol , HpNum(0) , NComm );

  // solve the flow part  - - - - - - - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  for( Index_Set tI = IndNwCsts[ wFi -1 ] ; *tI < Inf<Index>() ; tI++  ) { // examine the
   Index comm = *tI / NArcs;   // vector of the ordered new costs and pick ones
                               // corresponding to the arc wFi -1
   if( NwCsts[ *tI ] < -Eps &&  rhs > Eps ) {
    tFSol[ comm ] = min( Capacities[ *tI ] , rhs );
    rhs -= tFSol[ comm ];
    out += tFSol[ comm ] * NwCsts[ *tI ];
    }
   }

  // solve the fixed charge part  - - - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( ( XtrCosts[ wFi - 1 ] + out ) + Eps < 0 ) {
   out += XtrCosts[ wFi-1 ];
   for( Index k = 0 ; k < NComm ; k++ )
    FSolution[ (wFi-1) + NArcs*k ] = tFSol[ k ];
   }
  else {
   for( Index k = 0 ; k < NComm ; k++ )
    FSolution[ (wFi-1) + NArcs*k ] = 0;
   out = 0;
   }

  FiLmbd[ wFi-1 ] = out;
  } // end if( 1 <= wFi <= NArcs )

 delete[] tFSol;
 return( out );

 }  // end ( KnpsFiOrcl::solve )

/*--------------------------------------------------------------------------*/

void KnpsFiOrcl::quicksort( cIndex arc )
{
 // a non-recoursive QuickSort implementation - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 Index nSort;

 // find the smallest element - - - - - - - - - - - - - - - - - - - - - - - -

 Index_Set InMin , tI;
 double Min = Inf<double>();
 for( tI = IndNwCsts[ arc ] ; *tI < Inf<Index>() ; tI++ ) {
  const double tCh = NwCsts[ *tI ];
  if( tCh < Min ) {
   Min = tCh;
   InMin = tI;
   }
  }

 nSort = tI - IndNwCsts[ arc ];
 if( nSort <= 1 )   // no elements to sort
  return;

 Swap( *InMin , *IndNwCsts[ arc ] );

 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 int *Top = QSStck;
 int l = 0;
 int u = nSort - (++l);

 double pvt;  // pivot element: at the beginning, it is the smallest one

 int lwr;
 int upr;

 for(;;) {
  lwr = l + 1;
  upr = u;

  // set the pivot element - - - - - - - - - - - - - - - - - - - - - - - - - -
  pvt = NwCsts[ IndNwCsts[ arc ][ l ] ];

  // put the elements smaller than pivot at left of pivot, while the bigger
  // ones at right of pivot

  while( lwr <= upr ) {
   while( lwr <= upr )
    if( NwCsts[ IndNwCsts[ arc ][ lwr ] ] <= pvt )
     lwr++;
    else
     break;

   while( upr >= lwr )
    if( NwCsts[ IndNwCsts[ arc ][ upr ] ] > pvt )
     upr--;
    else
     break;

   if( lwr < upr )
	Swap( IndNwCsts[ arc ][ lwr++ ] , IndNwCsts[ arc ][ upr-- ] );
   }

  Swap( IndNwCsts[ arc ][ l ] , IndNwCsts[ arc ][ upr ] );    // swap the pivot at middle point

  // the following is the "recoursive" part, by using the stack- - - - - - - -

  if( --upr > l )    // left recoursion
   if( upr - 1 > l ) {
    *(Top++) = l;
    *(Top++) = upr;
    }
   else
    if( NwCsts[ IndNwCsts[ arc ][ upr ] ] < NwCsts[ IndNwCsts[ arc ][ l ] ] )
     Swap( IndNwCsts[ arc ][ l ] , IndNwCsts[ arc ][ upr ] );

  upr++;

  if( u > ++upr )    // right recoursion
   if( u > upr + 1 ) {
    *(Top++) = upr;
    *(Top++) = u;
    }
   else
    if( NwCsts[ IndNwCsts[ arc ][ u ] ] < NwCsts[ IndNwCsts[ arc ][ upr ] ] )
     Swap( IndNwCsts[ arc ][ u ] , IndNwCsts[ arc ][ upr ] );

  if( Top == QSStck )
   break;

  u = *(--Top);
  l = *(--Top);

  }  // end for( ever )
 }  // end( KnpsFiOrcl::quicksort )

/*--------------------------------------------------------------------------*/

void KnpsFiOrcl::CopyGi( cIndex Name , cIndex wFi , cFRow FlwSol )
{  // some exceptions  - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( Name >= MaxName )
  throw( NDOException( "KnpsFiOrcl::SetGiName:Name is higher than MaxName." ) );

 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 OldWFi[ Name ] = wFi;

 if( wFi <= 0 )
	 throw( NDOException( "KnpsFiOrcl::SetGiName:Name is higher than MaxName." ) );

 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 delete[] OldFSols[ Name ];
 if( wFi >  NArcs || Aggrgtd ) {
  OldFSols[ Name ] = new FNumber[ NArcs * NComm ];
  VectAssign( OldFSols[ Name ] , FlwSol , NArcs * NComm );
  }
 else {
  OldFSols[ Name ] = new FNumber[ NComm ];
  for(Index k = 0; k < NComm; k++ )
 	OldFSols[ Name ][ k ] = FlwSol[ k * NArcs + (wFi - 1) ];
  }

 } // end( KnpsFiOrcl::CopyGi )

/*--------------------------------------------------------------------------*/

void KnpsFiOrcl::FindGlobalLipschitz( void ) {

 HpNum Lip1, Lip2, LipConTmp;
 LpshCnst = 0;

 for( Index k = 0 ; k < NComm ; k++ )
  for( Index n = 0 ; n < NNodes ; n++ ) {
   Lip1 = - Deficits[ k * NNodes + n ];
   Lip2 = - Deficits[ k * NNodes + n ];
   for( Index j = 0 ; j < NArcs ; j++ )
	if( Endn[ j ] == n )
	 Lip1 += Capacities[ k * NArcs + j ];
    else
	 if( Startn[ j ] == n )
	  Lip2 -= Capacities[ k * NArcs + j ];
   LipConTmp = max( ABS( Lip1 ) , ABS( Lip2 ) );
   LpshCnst += LipConTmp * LipConTmp;
   }

 LpshCnst = sqrt( LpshCnst );
 } // end ( KnpsFiOrcl::FindGlobalLipschitz )

/*--------------------------------------------------------------------------*/
/*-------------------------- End File KnpsFiOrcl.C -------------------------*/
/*--------------------------------------------------------------------------*/

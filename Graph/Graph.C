/*--------------------------------------------------------------------------*/
/*---------------------------- File Graph.C --------------------------------*/
/*--------------------------------------------------------------------------*/
/*---                                                                     --*/
/*--  The Graph class provide an unified mean for reading descriptions of --*/
/*--  (Linear) Multicommodity Min Cost Flow Problems and storing them in  --*/
/*--  memory, along with a simple interface that can be used to access    --*/
/*--  and change the data.                                                --*/
/*--                                                                      --*/
/*--                            VERSION 2.01                              --*/
/*--                           11 - 05 - 2012                             --*/
/*--                                                                      --*/
/*--                  Original Idea and Implementation by:                --*/
/*--                                                                      --*/
/*--                          Antonio Frangioni                           --*/
/*--                      Operations Research Group                       --*/
/*--                     Dipartimento di Informatica                      --*/
/*--                         Universita' di Pisa                          --*/
/*--                                                                      --*/
/*--                           Paola Cappanera                            --*/
/*--                      Operations Research Group                       --*/
/*--                Dipartimento di Sistemi e Informatica                 --*/
/*--                        Universita' di Firenze                        --*/
/*--                                                                      --*/
/*--           Copyright (C) 1994 - 2012 by Antonio Frangioni             --*/
/*--                                                                      --*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- IMPLEMENTATION -------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "Graph.h"

#include <fstream>

#include <string.h>

#include <assert.h>

/*--------------------------------------------------------------------------*/
/*-------------------------------- USING -----------------------------------*/
/*--------------------------------------------------------------------------*/

#if( OPT_USE_NAMESPACES )
 using namespace MMCFGraph_di_unipi_it;
#else
 using namespace std;
#endif

/*--------------------------------------------------------------------------*/
/*-------------------------------- MACROS ----------------------------------*/
/*--                                                                      --*/
/*--      Some small macro definitions, used throughout the code.         --*/
/*--                                                                      --*/
/*--------------------------------------------------------------------------*/

#define GOODN( n ) if( ( n <= 0 ) || ( Index( n ) > NNodes ) )	\
                    throw( MMCFGException( "Invalid actual node name" ) )

#define GOODN2( n ) if( ( ( n <= 0 ) || ( Index( n ) > NNodes ) ) && \
			( n != -1 ) ) \
                     throw( MMCFGException( "Invalid generic node name" ) )

#define GOODP( k ) if( ( k <= 0 ) || ( Index( k ) > NumProd ) ) \
                    throw( MMCFGException( "Invalid actual product name" ) )

#define GOODP2( k ) if( ( ( k <= 0 ) || ( Index( k ) > NumProd ) ) && \
			( k != -1 ) ) \
                     throw( MMCFGException( "Invalid generic product name" ) )

#define GOODL( l ) if( Index( l ) > NCnst ) \
                    throw( MMCFGException( "Invalid link name" ) )

/*--------------------------------------------------------------------------*/
/*------------------------------ CONSTANTS ---------------------------------*/
/*--------------------------------------------------------------------------*/

static const Graph::CNumber C_INF = Graph::Inf<Graph::CNumber>();
static const Graph::FNumber F_INF = Graph::Inf<Graph::FNumber>();

/*--------------------------------------------------------------------------*/
/*--------------------------- LOCAL FUNCTIONS ------------------------------*/
/*--------------------------------------------------------------------------*/

template<class T>
static inline void VectAssign( T *const g , const T x , Graph::cIndex n )
{
 // g[ i ] = x for each i = 0 .. n - 1

 for( T *tg = g + n ; tg > g ; )
  *(--tg) = x;
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

template<class T1, class T2>
static inline void VectAssign( T1 *g1 , const T2 *g2 , Graph::Index n )
{
 // g1 := g2

 for( ; n-- ; )
  *(g1++) = *(g2++);
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

template<class T>
inline void VectSubtract( T *const g , const T x , Graph::cIndex n )
{
 // g[ i ] -= x for each i = 0 .. n - 1; useful for *unsigned* data types
 // where VectSum( g , -x , n ) would not work

 for( T *tg = g + n ; tg > g ; )
  *(--tg) -= x;
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

template<class T>
static inline bool EqualVect( const T *g1 ,  const T *g2 , Graph::Index n )
{
 // returns true <=> the two n-vectors g1 and g2 are element-wise identical

 for( ; n-- ; )
  if( *(g1++) != *(g2++) )
   return( false );

 return( true );
 }

/*--------------------------------------------------------------------------*/
/*----------------------- IMPLEMENTATION OF Graph  -------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- PUBLIC METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/
/*----------------------------- CONSTRUCTORS -------------------------------*/
/*--------------------------------------------------------------------------*/

Graph::Graph( const char *const FN , char FT )
{
 // check parameters- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( ( FT != 'm' ) && ( FT != 'p' ) && ( FT != 'd' ) && ( FT != 'o' ) &&
     ( FT != 'u' ) && ( FT != 's' ) && ( FT != 'c' ) )
  throw( MMCFGException( "Graph::Graph( file ): invalid file type" ) );

 bool FourFiles = ( ( FT != 's' ) && ( FT != 'c' ) );

 // in principle there is no "extra" stuff- - - - - - - - - - - - - - - - - -

 NXtrV = NXtrC = 0;
 IdxBeg = CoefIdx = 0;
 CoefVal = 0;

 // reading general informations- - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 cIndex l = strlen( FN );
 char *Name = new char[ l + 5 ];  // temporary string containing the constant
 strcpy( Name , FN );             // part of the pathname + space for `.XXX'

 if( FourFiles )
  strcpy( Name + l , ".nod" );

 ifstream inFile( Name );
 if( ! inFile.is_open() )
  throw( MMCFGException( "Graph::Graph( file ): can't open file" ) );

 if( FourFiles ) {
  inFile >> NComm;
  inFile >> NNodes;
  inFile >> NArcs;
  inFile >> NCnst;
  }
 else {
  inFile >> NNodes;
  inFile >> NArcs;
  inFile >> NComm;

  if( FT == 'c' ) {   // in the PPRN format, read description of side - - - -
   inFile >> NXtrC;   // constraints and prepare the data structures

   Index NNZ;
   inFile >> NNZ;

   if( NNZ ) {
    IdxBeg  = new int[ NXtrC ];
    CoefIdx = new int[ NNZ ];
    CoefVal = new double[ NNZ ];
    }
   }

  NCnst = NArcs;
  }

 if( NNodes <= 1 )
  throw( MMCFGException( "Graph::Graph( file ): wrong node number" ) );
 if( NArcs <= 0 )
  throw( MMCFGException( "Graph::Graph( file ): wrong arc number" ) );
 if( NComm <= 0 )
  throw( MMCFGException( "Graph::Graph( file ): wrong commodity number" ) );
 if( NCnst > NArcs )
  throw( MMCFGException( "Graph::Graph( file ): wrong constraints number" ) 
	 );

 if( FourFiles )
  inFile.close();

 Index_Set Origins;
 Index_Set Destins;

 Index_Set StartOfK;
 Index_Set TempIdx;

 Index NumProd = NComm;

 // format-dependent part - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( FourFiles )  // preparing to read the supply file
  if( FT == 'u' ) {
   strcpy( Name + l , ".od" );
   FT = 'd';
   }
  else
   strcpy( Name + l , ".sup" );

 // determining the actual number of commodities for (OSP) or (ODS)- - - - - -
 // formulations: in the first case, a commodity is a pair ( product , - - - -
 // origin ), while in the second case it is a triplet ( product ,-  - - - - -
 // origin , destination ) - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( ( FT == 'd' ) || ( FT == 'o' ) ) {
  StartOfK = new Index[ NumProd + 1 ];
  VectAssign( StartOfK , Index( 0 ) , NumProd + 1 );

  inFile.clear();        // ensure failbits are not left dirty
  inFile.open( Name );   // commodities can be told from supplies
  if( ! inFile.is_open() )
   throw( MMCFGException( "Graph::Graph( file ): can't open file" ) );

  int origin;   // the *.sup file is read once here just to count the number
  int dest;     // of commodities: the actual data reading will be done later
  int comm;
  FNumber flow;

  if( FT == 'd' )  // in (ODS) count the different O/D pairs - - - - - - - - -
   while( inFile >> origin ) {
    GOODN( origin );

    inFile >> dest;
    GOODN( dest );

    inFile >> comm;
    GOODP2( comm );

    inFile >> flow;

    if( comm != -1 )
     StartOfK[ comm ]++;
    else
     for( Index i = NumProd ; i ; )
      StartOfK[ i-- ]++;
    }
  else  // in (OSP) count the number of different Origins- - - - - - - - - - -
   while( inFile >> origin ) {
    GOODN( origin );

    inFile >> dest;
    GOODN2( dest );

    inFile >> comm;
    GOODP2( comm );

    inFile >> flow;

    if( dest == -1 ) {
     if( comm != -1 )
      StartOfK[ comm ]++;
     else
      for( Index i = NumProd ; i ; )
       StartOfK[ i-- ]++;
     }
    }

  // really construct StartOfK - - - - - - - - - - - - - - - - - - - - - - - -

  *StartOfK = 0;
  for( Index i = 2 ; i <= NumProd ; i++ )
   StartOfK[ i ] += StartOfK[ i - 1 ];

  NComm = StartOfK[ NumProd ]; // note that NComm can "surprisingly" be
                               // < NumProd if some product does not appear

  Origins = new Index[ NComm ];

  inFile.close();

  }  // end if( (OSP) or (ODP) )

 // allocating and initializing memory- - - - - - - - - - - - - - - - - - - -
 // (note that this part is again common) - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 C = new CRow[ NComm + 1 ];  // allocate costs- - - - - - - - - - - - - - - -

 if( FT == 'c' )
  for( Index i = 0 ; i < NComm ; i++ )
   C[ i ] = new CNumber[ NArcs ];
 else
  for( Index i = 0 ; i < NComm ; i++ ) {
   C[ i ] = new CNumber[ NArcs ];         // arcs are un-existent unless
   //   VectAssign( C[ i ] , C_INF , NArcs );  // otherwise stated
   VectAssign( C[ i ] , C_INF , NArcs );  // otherwise stated
   }

 C[ NComm ] = 0;

 U = new FRow[ NComm + 2 ];  // allocate capacities - - - - - - - - - - - - -

 if( FT == 'c' )
  for( Index i = 0 ; i < NComm ; i++ )
   U[ i ] = new FNumber[ NArcs ];
 else
  for( Index i = 0 ; i < NComm ; i++ ) {
   U[ i ] = new FNumber[ NArcs ];                // arcs are un-existent
   VectAssign( U[ i ] , FNumber( 0 ) , NArcs );  // unless otherwise stated
   }

 U[ NComm ] = U[ NComm + 1 ] = 0;

 B = new FRow[ NComm + 2 ];  // allocate deficits - - - - - - - - - - - - - -

 if( FT == 'c' )
  for( Index i = 0 ; i < NComm ; i++ )
   B[ i ] = new FNumber[ NNodes ];
 else
  for( Index i = 0 ; i < NComm ; i++ ) {
   B[ i ] = new FNumber[ NNodes ];                // nodes all have 0 deficit
   VectAssign( B[ i ] , FNumber( 0 ) , NNodes );  // unless otherwise stated
   }

 B[ NComm ] = B[ NComm + 1 ] = 0;

 // allocate start/end nodes and mutual capacities- - - - - - - - - - - - - -

 Startn = new Index[ NArcs ];
 Endn   = new Index[ NArcs ];

 UTot   = new FNumber[ NArcs ];

 // allocate info on integrality of the variables - - - - - - - - - - - - - -

 NInt = new Index[ NComm + 1 ];
 VectAssign( NInt , Index( 0 ) , NComm + 1 );

 WIsInt = new Index_Set[ NComm + 1 ];
 VectAssign( WIsInt , Index_Set( 0 ) , NComm + 1 );

 // reading supply/demand infos, or everything in one-files format- - - - - -
 // (this part is partly splitted again)- - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( FourFiles ) {
  if( FT == 'm' )
   TempIdx = new Index[ NCnst ];
  else {
   TempIdx = new Index[ max( Index( NumProd + 1 ) , NComm ) ];

   if( ( FT == 'o' ) || ( FT == 'd' ) )
    VectAssign( TempIdx , StartOfK , NumProd + 1 );
   }

  inFile.clear();        // ensure failbits are not left dirty
  inFile.open( Name );  // the right name is already there
  if( ! inFile.is_open() )
   throw( MMCFGException( "Graph::Graph( file ): can't open file" ) );
  }

 switch( FT ) {

 case( 's' ): // Canadian format- - - - - - - - - - - - - - - - - - - - - - -
 {            //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  // allocate the data structures for "extra" things- - - - - - - - - - - - -

  CRow FC = C[ NComm ] = new CNumber[ NXtrV = NArcs ];

  U[ NComm ] = new FNumber[ NArcs ];
  VectAssign( U[ NComm ] , FNumber( 0 ) , NArcs );      // "extra" variables
  U[ NComm + 1 ] = new FNumber[ NArcs ];                // are in the ...
  VectAssign( U[ NComm + 1 ] , FNumber( 1 ) , NArcs );  // ... [0, 1] range
  NInt[ NComm ] = NArcs;                                // ... and integer

  for( Index i = 0 ; i < NArcs ; i++ ) {  // read arc-related info- - - - - -
   inFile >> Endn[ i ];
   GOODN( Endn[ i ] );

   inFile >> Startn[ i ];
   GOODN( Startn[ i ] );
   if( Startn[ i ] == Endn[ i ] )
    throw( MMCFGException( "Graph::Graph( file ): self-loop" ) );

   inFile >> FC[ i ];

   FNumber f;
   inFile >> f;

   UTot[ i ] = ( f >= 0 ? f : F_INF );

   Index h;
   inFile >> h;

   for( ; h-- ; ) {
    Index k;
    inFile >> k;
    GOODP( k );

    inFile >> C[ --k ][ i ];
    inFile >> f;

    U[ k ][ i ] = ( f >= 0 ? f : F_INF );

    }  // end for( h )
   }  // end for( i )

  for( Index k ; inFile >> k ; ) {  // read node-related info - - - - - - - -
   GOODP( k );

   Index i;
   inFile >> i;
   GOODN( i );

   inFile >> B[ --k ][ --i ];
   }

  break;

  }  // end case( s ) - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 case( 'c' ): // PPRN format- - - - - - - - - - - - - - - - - - - - - - - - -
 {            //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  for( Index k = 0 ; k < NComm ; k++ )  // read all costs
   for( Index i = 0 ; i < NArcs ; )
    inFile >> C[ k ][ i++ ];

  for( Index k = 0 ; k < NComm ; k++ )  // read all capacities
   for( Index i = 0 ; i < NArcs ; ) {
    FNumber f;
    inFile >> f;

    U[ k ][ i++ ] = ( f >= 0 ? f : F_INF );
    }

  for( Index k = 0 ; k < NComm ; k++ )  // read all supplies
   for( Index i = 0 ; i < NNodes ; )
    inFile >> B[ k ][ i++ ];

  for( Index i = 0 ; i < NArcs ; ) {      // read all total capacities
   FNumber f;
   inFile >> f;

   UTot[ i++ ] = ( f >= 0 ? f : F_INF );
   }

  for( Index i = 0 ; i < NArcs ; i++ ) {  // read graph topology
   inFile >> Startn[ i ];
   GOODN( Startn[ i ] );

   inFile >> Endn[ i ];
   GOODN( Endn[ i ] );

   if( Startn[ i ] == Endn[ i ] )
    throw( MMCFGException( "Graph::Graph( file ): self-loop" ) );
   }

  if( NXtrC ) {  // if there are "extra" constraints- - - - - - - - - - - - -
   FRow EL = B[ NComm ]     = new FNumber[ NXtrC ];
   FRow EU = B[ NComm + 1 ] = new FNumber[ NXtrC ];
  
   for( Index i = 0 ; i < NXtrC ; ) {  // read "extra" Uppr./Lwr. bounds
    inFile >> EL[ i ];
    inFile >> EU[ i++ ];
    }

   Index currc = 0;
   Index currpos = 0;

   for( Index j ; inFile >> j ; ) {  // read extra constraints description:
    Index k;                         // j = arc name
    inFile >> k;                     // commodity name

    CoefIdx[ currpos ] = (--k) * NArcs + (--j);

    Index h;
    inFile >> h;                     // constraint name
    h--;

    while( h > currc ) {
     IdxBeg[ currc ] = currpos;
     currc++;
     }

    inFile >> CoefVal[ currpos++ ];    // the coefficient

    }  // end( for( ! eof() ) )

   IdxBeg[ currc ] = currpos;

   }  // end( if( extra constraints ) )

  break;

  }  // end case( c ) - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 case( 'm' ): // mnetgen format - - - - - - - - - - - - - - - - - - - - - - -
 {            //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  for( Index node ; inFile >> node ; ) {
   GOODN( node );

   int comm;
   inFile >> comm;
   GOODP2( comm );

   FNumber flow;
   inFile >> flow;

   if( comm == -1 )
    for( Index k = 0 ; k < NComm ; )
     B[ k++ ][ node - 1 ] = - flow;
   else
    B[ comm - 1 ][ node - 1 ] = - flow;
   }

  break;

  }  // end case( m ) - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 case( 'p' ): // JL (PSP) format- - - - - - - - - - - - - - - - - - - - - - -
 {            //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  for( int origin ; inFile >> origin ; ) {
   GOODN2( origin );

   int dest;
   inFile >> dest;
   GOODN2( dest );
   if( ( ( origin == -1 ) && ( dest == -1 ) ) ||
       ( ( origin != -1 ) && ( dest != -1 ) ) )
    throw( MMCFGException( "Graph::Graph( file ): exactly one p/d in PSP" ) );

   int comm;
   inFile >> comm;
   GOODP2( comm );

   FNumber flow;
   inFile >> flow;

   if( comm != -1 ) {
    comm--;

    if( origin < 0 )
     B[ comm ][ dest - 1 ] = flow;
    else
     B[ comm ][ origin - 1 ] = - flow;
    }
   else
    if( origin < 0 )
     for( Index i = NumProd ; i-- ; )
      B[ i ][ dest - 1 ] = flow;
    else
     for( Index i = NumProd ; i-- ; )
      B[ i ][ origin - 1 ] = - flow;

   }  // end( for( ! eof() ) )

  break;

  }  // end case( p ) - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 case( 'o' ): // JL (OSP) format- - - - - - - - - - - - - - - - - - - - - - -
 {            // while reading supplies, commodity "names" are assigned - - -

  for( int origin ; inFile >> origin ; ) {
   int dest;
   inFile >> dest;

   int comm;
   inFile >> comm;

   FNumber flow;
   inFile >> flow;

   if( comm != -1 ) {
    // origin or destination node for the given pair ( product , origin )

    Index i = StartOfK[ --comm ];

    while( ( i < TempIdx[ comm ] ) && ( Origins[ i ] != Index( origin ) ) )
     i++;  // seek the name of the commodity

    if( i == TempIdx[ comm ] ) {  // a "new" commodity
     Origins[ i ] = origin;
     TempIdx[ comm ]++;
     }

    if( dest == -1 )  // it is an origin
     B[ i ][ origin - 1 ] = - flow;
    else
     B[ i ][ dest - 1 ] = flow;
    }
   else {  // comm == -1
    // origin or destination node for all the commodities ( k , origin )
    // for each product k

    Index i;
    for( Index k = NumProd ; k-- ; ) {
     i = StartOfK[ k ];

     while( ( i < TempIdx[ k ] ) && ( Origins[ i ] != Index( origin ) ) )
      i++;  // seek the name of the commodity

     if( i == TempIdx[ k ] ) {  // a "new" commodity
      Origins[ i ] = origin;
      TempIdx[ i ]++;
      }

     if( dest == -1 )  // it is an origin
      B[ i ][ origin - 1 ] = - flow;
     else
      B[ i ][ dest - 1 ] = flow;

     }  // end for( k )
    }  // end else( comm == -1 )
   }  // end while( ! eof() )

  break;

  }  // end case( o ) - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 default: // 'd' == JL (ODS) format - - - - - - - - - - - - - - - - - - - - -
 {        // while reading supplies, commodity "names" are assigned - - - - -

  Destins = new Index[ NComm ];

  for( int origin ; inFile >> origin ; ) {
   int dest;
   inFile >> dest;

   int comm;
   inFile >> comm;

   FNumber flow;
   inFile >> flow;

   if( comm != -1 ) {
    comm--;

    Origins[ TempIdx[ comm ] ] = origin;
    Destins[ TempIdx[ comm ] ] = dest;

    B[ TempIdx[ comm ] ][ dest - 1 ] = flow;
    B[ TempIdx[ comm ]++ ][ origin - 1 ] = - flow;
    }
   else
    for( Index i = NumProd ; i-- ; ) {
     Origins[ TempIdx[ i ] ] = origin;
     Destins[ TempIdx[ i ] ] = dest;

     B[ TempIdx[ i ] ][ dest - 1 ] = flow;
     B[ TempIdx[ i ]++ ][ origin - 1 ] = - flow;
     }
    }  // end while( ! eof() )
   }   // end default()- - - - - - - - - - - - - - - - - - - - - - - - - - - -
  }    // end switch( FT ) - - - - - - - - - - - - - - - - - - - - - - - - - -
       //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 inFile.close();

 if( FourFiles ) {
  // continue only in the multi-file formats - - - - - - - - - - - - - - - - -

  VectAssign( UTot , F_INF , NArcs );

  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // now the though part: reading arc infos - - - - - - - - - - - - - - - - -
  // (again, this part is splitted) - - - - - - - - - - - - - - - - - - - - -

  strcpy( Name + l , ".arc" );

  inFile.clear();        // ensure failbits are not left dirty
  inFile.open( Name );
  if( ! inFile.is_open() )
   throw( MMCFGException( "Graph::Graph( file ): can't open file" ) );

  if( FT == 'm' )   // mnetgen format - - - - - - - - - - - - - - - - - - - -
  {                 // it is dealt with separatedly, since it's simpler: the
   Index who;       // name of the arc (who) is explicitely given
   
   while( inFile >> who ) {
    if( ( who <= 0 ) || ( who > NArcs ) )
     throw( MMCFGException( "Graph::Graph( file ): invalid arc name" ) );
    who--;

    Index from;
    inFile >> from;
    GOODN( from );

    Index to;
    inFile >> to;
    GOODN( to );
    if( from == to )
     throw( MMCFGException( "Graph::Graph( file ): self-loop" ) );

    int comm;
    inFile >> comm;
    GOODP2( comm );

    CNumber cost;
    inFile >> cost;

    FNumber cap;
    inFile >> cap;
    if( cap < 0 )
     cap = F_INF;

    Index ptr;
    inFile >> ptr;
    GOODL( ptr );

    if( ptr )
     TempIdx[ ptr - 1 ] = who;

    Startn[ who ] = from;
    Endn[ who ] = to;

    if( comm == -1 )
     for( Index k = 0 ; k < NComm ; ) {
      C[ k ][ who ] = cost;
      U[ k++ ][ who ] = cap;
      }
    else {
     C[ --comm ][ who ] = cost;
     U[ comm ][ who ] = cap;
     }
    }   // end while()
   }    // end if( mnetgen )
  else {  // the three JL formats - - - - - - - - - - - - - - - - - - - - - -
          //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Index unbndld = 0;  // counter of unbundled links so far

   for( Index from ; inFile >> from ; ) {  // first, the usual data reading
    GOODN( from );

    Index to;
    inFile >> to;
    GOODN( to );
    if( from == to )
     throw( MMCFGException( "Graph::Graph( file ): self-loop" ) );

    int comm;
    inFile >> comm;
    GOODP2( comm );

    CNumber cost;
    inFile >> cost;

    int cap;
    inFile >> cap;

    int origin;
    inFile >> origin;
    GOODN2( origin );

    int dest;
    inFile >> dest;
    GOODN2( dest );

    Index ptr;
    inFile >> ptr;
    GOODL( ptr );

    // now, the main part: from the triplet (origin, destination, product)
    // plus the file type (p, d, o) a list of applicable commodities is
    // constructed and put into TempIdx: then the arc will be replicated over
    // all the commodities of the list
    // in all the three fields, a "-1" takes the place of a wildcard

    Index TmpCommCntr = 1;

    switch( FT ) {

    case( 'p' ):  // the simplest, only 2 subcases- - - - - - - - - - - - - -

     if( comm != -1 )      // a specific commodity (prod)
      *TempIdx = comm - 1;
     else                  // all commodities (!?)
      for( Index i = TmpCommCntr = NComm ; i-- ; )
       TempIdx[ i ] = i;

     break;

    case( 'o' ):  // tougher, 4 subcases- - - - - - - - - - - - - - - - - - -

     if( comm != -1 ) {
      Index i = StartOfK[ comm - 1 ];

      if( origin != -1 ) {          // a specific commodity ( prod , origin )
       while( ( i < StartOfK[ comm ] ) &&
	      ( Origins[ i ] != Index( origin ) ) )
        i++;

       assert( i != StartOfK[ comm ] );
       *TempIdx = i;
       }
      else                          // all commodities with a given product
       for( TmpCommCntr = 0 ; i < StartOfK[ comm ] ; )
        TempIdx[ TmpCommCntr++ ] = i++;
      }
     else
      if( origin != -1 )            // all commodities with a given origin
       for( Index i = TmpCommCntr = 0 ; i < NumProd ; ) {
        Index k = StartOfK[ i++ ];

        while( ( k < StartOfK[ i ] ) &&
	       ( Origins[ k ] != Index( origin ) ) )
         k++;

        if( k < StartOfK[ comm ] )
         TempIdx[ TmpCommCntr++ ] = k;
        }
      else                         // all commodities
       for( Index i = TmpCommCntr = NComm ; i-- ; )
        TempIdx[ i ] = i;

     break;

    default:     // == 'd', the toughest: 8 subcases- - - - - - - - - - - - -

     if( comm != -1 ) {
      Index i = StartOfK[ comm - 1 ];

      if( dest != -1 ) {
       if( origin != -1 ) {      // a specific commodity (prod, origin, dest)
        while( ( i < StartOfK[ comm ] ) && ( Destins[ i ] != Index( dest ) )
               && ( Origins[ i ] != Index( origin ) ) )
         i++;

        assert( i != StartOfK[ comm ] );
        *TempIdx = i;
        }
       else                      // all commodities with a given (prod, dest)
        for( TmpCommCntr = 0; i < StartOfK[ comm ] ; i++ )
         if( Destins[ i ] == Index( dest ) )
          TempIdx[ TmpCommCntr++ ] = i;
       }
      else {                     // dest == -1
       if( origin != -1 ) {      // all commodities with a given (prod, orig)
        for( TmpCommCntr = 0; i < StartOfK[ comm ] ; i++ )
         if( Origins[ i ] == Index( origin ) )
          TempIdx[ TmpCommCntr++ ] = i;
        }
       else                      // all commodities with a given (prod)
        for( TmpCommCntr = 0; i < StartOfK[ comm ] ; i++ )
         TempIdx[ TmpCommCntr++ ] = i;
       }
      }
     else                       // comm == -1
      if( dest != -1 ) {
       if( origin != -1 ) {     // all commodities with a given (orig, dest)
        for( Index i = TmpCommCntr = 0 ; i < NComm ; i++ )
         if( ( Destins[ i ] == Index( dest ) ) &&
	     ( Origins[ i ] == Index( origin ) ) )
          TempIdx[ TmpCommCntr++ ] = i;
        }
       else                    // all commodities with a given (dest)
        for( Index i = TmpCommCntr = 0 ; i < NComm ; i++ )
         if( Destins[ i ] == Index( dest ) )
          TempIdx[ TmpCommCntr++ ] = i;
       }
     else                      // dest == -1
      if( origin != -1 ) {     // all commodities with a given (origin)
       for( Index i = TmpCommCntr = 0 ; i < NComm ; i++ )
        if( Origins[ i ] == Index( origin ) )
         TempIdx[ TmpCommCntr++ ] = i;
        }
      else                     // all commodities
       for( Index i = TmpCommCntr = NComm ; i-- ; )
        TempIdx[ i ] = i;

     }  // end switch( FT )

    // now "filling" the proper arc for each commodity

    for( Index i = TmpCommCntr ; i-- ; ) {
     comm = TempIdx[ i ];
     Index who;

     if( ptr ) {      // if ptr != 0 it's easy
      who = ptr - 1;  // ( ptr - 1 ) is already the correct name

      Startn[ who ] = from;
      Endn[ who ] = to;
      }
     else {           // otherwise find the "name" of arc (from, to)
      Index k = 0;    // and put it into who

      while( ( k < unbndld ) &&
             ( ( Startn[ NCnst + k ] != from ) ||
               ( Endn[ NCnst + k ] != to ) ||
               ( C[ comm ][ NCnst + k ] < C_INF ) ) )
       k++;

      // search for an arc (from, to) already defined among the unbundled ones
      // and whose "instance" relative to commodity comm has not already been
      // taken: this is not the only way of accomodating unbundled arcs, (in
      // case of multiple instances of an unbundled arc (i, j)), but it is
      // easy to see that all the resulting problems, however you distribute
      // the instances to arcs, are equivalent

      who = NCnst + k;
      assert( who < NArcs );

      if( k == unbndld ) {  // if no such arc exists ...
       unbndld++;           // ... a new one is created
       Startn[ who ] = from;
       Endn[ who ] = to;
       }
      }  // end else( ! ptr )

     C[ comm ][ who ] = cost;
     U[ comm ][ who ] = ( cap >= 0 ? cap : F_INF );

     } // end for( all comm. )
    }  // end while( ! eof() )

   if( NCnst + unbndld < NArcs )
    NArcs = NCnst + unbndld;

   }   // end else( JL formats )

  inFile.close();

  // reading mutual capacities- - - - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  strcpy( Name + l , ".mut" );
  inFile.clear();        // ensure failbits are not left dirty
  inFile.open( Name );
  if( ! inFile.is_open() )
   throw( MMCFGException( "Graph::Graph( file ): can't open file" ) );

  for( Index i = 0 ; i < NCnst ; ) {
   Index j;
   inFile >> j;

   FNumber f;
   inFile >> f;

   if( FT == 'm' )
    j = TempIdx[ i++ ];
   else
    j = i++;

   UTot[ j ] = ( f >= 0 ? f : F_INF );
   }

  // temporary deallocation and final things- - - - - - - - - - - - - - - - -

  if( ( FT == 'd' ) || ( FT == 'o' ) ) {
   delete[] StartOfK;
   delete[] Origins;

   if( FT == 'd' )
    delete[] Destins;
   }

  delete[] TempIdx;

  }  // end if( FourFiles )

 delete[] Name;


 // common initializations- - - - - - - - - - - - - - - - - - - - - - - - - -

 CmnIntlz();

 }  // end( Graph( char* , char ) )

/*--------------------------------------------------------------------------*/

Graph::Graph( cIndex n , cIndex m , cIndex comm , const cFRow* Def ,
              cIndex_Set S , cIndex_Set E , cFRow CapTot , const cFRow* Cap ,
              const cCRow* Cost )
{
 // first initializations - - - - - - - - - - - - - - - - - - - - - - - - - -

 NNodes = n;
 NArcs = NCnst = m;
 NComm = comm;

 NXtrV = NXtrC = 0;
 IdxBeg = CoefIdx = 0;
 CoefVal = 0;

 C = new CRow[ NComm + 1 ];  // allocate costs- - - - - - - - - - - - - - - -

 for( Index k = 0 ; k < NComm ; k++ ) {
  C[ k ] = new CNumber[ NArcs ];

  if( Cost && Cost[ k ] )
   VectAssign( C[ k ] , Cost[ k ] , NArcs );
  else
   VectAssign( C[ k ] , CNumber( 0 ) , NArcs );
  }

 C[ NComm ] = 0;

 U = new FRow[ NComm + 2 ];  // allocate capacities - - - - - - - - - - - - -

 for( Index k = 0 ; k < NComm ; k++ ) {
  U[ k ] = new FNumber[ NArcs ];

  if( Cap && Cap[ k ] )
   VectAssign( U[ k ] , Cap[ k ] , NArcs );
  else
   VectAssign( U[ k ] , F_INF , NArcs );
  }

 U[ NComm ] = U[ NComm + 1 ] = 0;

 B = new FRow[ NComm + 2 ];  // allocate deficits - - - - - - - - - - - - - -

 for( Index k = 0 ; k < NComm ; k++ ) {
  B[ k ] = new FNumber[ NNodes ];

  if( Def && Def[ k ] )
   VectAssign( B[ k ] , Def[ k ] , NNodes );
  else
   VectAssign( B[ k ] , FNumber( 0 ) , NNodes );
  }

 B[ NComm ] = B[ NComm + 1 ] = 0;

 // allocate start/end nodes and mutual capacities- - - - - - - - - - - - - -

 Startn = new Index[ NArcs ];
 VectAssign( Startn , S , NArcs );

 Endn = new Index[ NArcs ];
 VectAssign( Endn , E , NArcs );

 UTot = new FNumber[ NArcs ];

 if( CapTot )
  VectAssign( UTot , CapTot , NArcs );
 else
  VectAssign( UTot , F_INF , NArcs );

 // allocate info on integrality of the variables - - - - - - - - - - - - - -

 NInt = new Index[ NComm + 1 ];
 VectAssign( NInt , Index( 0 ) , NComm + 1 );

 WIsInt = new Index_Set[ NComm + 1 ];
 VectAssign( WIsInt , Index_Set( 0 ) , NComm + 1 );

 // common initializations- - - - - - - - - - - - - - - - - - - - - - - - - -

 CmnIntlz();

 }  // end( Graph( mem ) )

/*--------------------------------------------------------------------------*/
/*-------------- METHODS FOR READING THE DATA OF THE PROBLEM ---------------*/
/*--------------------------------------------------------------------------*/

void Graph::ExtraConstr( int *IBeg , int *Indx , double *Vals ,
			 Index FrstC , Index LstC )
{
 if( NXtrC ) {
  FrstC = max( FrstC , Index( 0 ) );
  LstC = min( LstC , Index( NXtrC - 1 ) );

  VectAssign( IBeg , IdxBeg + FrstC , LstC - FrstC + 2 );
  if( IdxBeg[ FrstC ] )
   VectSubtract( IBeg , IdxBeg[ FrstC ] , LstC - FrstC + 2 );

  cIndex nnz = IdxBeg[ LstC + 1 ] - IdxBeg[ FrstC ];

  VectAssign( Indx , CoefIdx + IdxBeg[ FrstC ] , nnz );
  VectAssign( Vals , CoefVal + IdxBeg[ FrstC ] , nnz );
  }
 }

/*--------------------------------------------------------------------------*/
/*----------------- METHODS THAT ALLOW TO CHANGE THE GRAPH -----------------*/
/*--------------------------------------------------------------------------*/

void Graph::UpDtTotCap( cFRow NewU )
{
 if( NewU )
  VectAssign( UTot , NewU , NArcs );
 else
  VectAssign( UTot , F_INF , NArcs );
 }

/*- - - - - - - - - - - - - - - - - - -- - - - - - - - - - - - - - - - - - -*/

void Graph::UpDtTotCapJ( cIndex j , cFNumber NewUj )
{
 UTot[ j ] = NewUj;
 }

/*--------------------------------------------------------------------------*/

void Graph::UpDtArcCstK( cIndex k , cCRow NewCk )
{
 if( NewCk )
  VectAssign( C[ k ] , NewCk , ( k < NComm ? NArcs : NXtrV ) );
 else
  VectAssign( C[ k ] , CNumber( 0 ) , ( k < NComm ? NArcs : NXtrV ) );
 }

/*- - - - - - - - - - - - - - - - - - -- - - - - - - - - - - - - - - - - - -*/

void Graph::UpdateArcCstKJ( cIndex k , cIndex j , cCNumber NewCkj )
{
 C[ k ][ j ] = NewCkj;
 }

/*--------------------------------------------------------------------------*/

void Graph::UpDtArcCapK( cIndex k , cFRow NewUk )
{
 if( NewUk )
  VectAssign( U[ k ] , NewUk , ( k < NComm ? NArcs : NXtrV ) );
 else
  VectAssign( U[ k ] , F_INF , ( k < NComm ? NArcs : NXtrV ) );
 }

/*- - - - - - - - - - - - - - - - - - -- - - - - - - - - - - - - - - - - - -*/

void Graph::UpDtArcCapKJ( cIndex k , cIndex j , cFNumber NewUkj )
{
 U[ k ][ j ] = NewUkj;
 }

/*--------------------------------------------------------------------------*/

void Graph::UpDtNdeDfctK( cIndex k , cFRow NewDk )
{
 if( NewDk )
  VectAssign( B[ k ] , NewDk , ( k < NComm ? NNodes : NXtrC ) );
 else
  VectAssign( B[ k ] , FNumber( 0 ) , ( k < NComm ? NNodes : NXtrC ) );
 }

/*- - - - - - - - - - - - - - - - - - -- - - - - - - - - - - - - - - - - - -*/

void Graph::UpDtNdeDfctKJ( cIndex k , cIndex j , cFNumber NewDkj )
{
 B[ k ][ j ] = NewDkj;
 }

/*--------------------------------------------------------------------------*/

void Graph::SetIntVar( cIndex k , bool IntVld , cIndex_Set nms ,
                       cIndex strt , Index stp )
{
 cIndex mxind = ( k >= NComm ? NXtrV : NArcs );
 if( stp > mxind )
  stp = mxind;

 if( nms )
  while( *nms < strt )
   nms++;

 if( IntVld ) {  // set variables to "integer"- - - - - - - - - - - - - - - -
  if( NInt[ k ] == mxind )  // all variables are already integer
   return;                  // making some more is impossible

  if( ( ! nms ) && ( ! strt ) && ( stp == mxind ) )  // set all variables
   NInt[ k ] = mxind;
  else {                            // set only a subset of the variables
   // first turn WIsInt[ k ] into a 0/1 vector

   Index_Set WIIk = WIsInt[ k ];
   if( ! WIIk ) {
    WIsInt[ k ] = WIIk = new Index[ mxind ];
    VectAssign( WIIk , Index( 0 ) , mxind );
    }
   else {
    cIndex_Set tWIIk = WIIk + NInt[ k ] - 1;
    for( Index i = mxind ; i-- ; )
     if( *tWIIk == i ) {
      WIIk[ i ] = Index( 1 );
      tWIIk--;
      }
     else
      WIIk[ i ] = Index( 0 );
    }

   // now set new 1s where necessary

   if( nms )
    while( *nms < stp )
     WIIk[ *(nms++) ] = Index( 1 );
   else
    for( Index i = strt ; i < stp ; )
     WIIk[ i++ ] = Index( 1 );

   // finally, transform back WIsInt[ k ] into a vector of indices

   Index ncnt = 0;
   for( Index i = 0 ; i < mxind ; i++ )
    if( WIIk[ i ] )
     WIIk[ ncnt++ ] = i;

   NInt[ k ] = ncnt;
   }

  if( NInt[ k ] == mxind ) {
   delete[] WIsInt[ k ];
   WIsInt[ k ] = 0;
   }
  else
   WIsInt[ k ][ NInt[ k ] ] = Inf<Index>();
  }
 else {          // set variables to "continuous" - - - - - - - - - - - - - -
  if( ! NInt[ k ] )  // all variables are already continuous
   return;           // making some more is impossible

  if( ( ! nms ) && ( ! strt ) && ( stp == mxind ) )  // set all variables
   NInt[ k ] = 0;
  else {                            // set only a subset of the variables
   // first turn WIsInt[ k ] into a 0/1 vector

   Index_Set WIIk = WIsInt[ k ];
   if( ! WIIk ) {
    WIsInt[ k ] = WIIk = new Index[ mxind ];
    VectAssign( WIIk , Index( 1 ) , mxind );
    }
   else {
    cIndex_Set tWIIk = WIIk + NInt[ k ] - 1;
    for( Index i = mxind ; i-- ; )
     if( *tWIIk == i ) {
      WIIk[ i ] = Index( 1 );
      tWIIk--;
      }
     else
      WIIk[ i ] = Index( 0 );
    }

   // now set new 0s where necessary

   if( nms )
    while( *nms < stp )
     WIIk[ *(nms++) ] = Index( 0 );
   else
    for( Index i = strt ; i < stp ; )
     WIIk[ i++ ] = Index( 0 );

   // finally, transform back WIsInt[ k ] into a vector of indices

   Index ncnt = 0;
   for( Index i = 0 ; i < mxind ; i++ )
    if( WIIk[ i ] )
     WIIk[ ncnt++ ] = i;

   NInt[ k ] = ncnt;
   }

  if( ! NInt[ k ] ) {
   delete[] WIsInt[ k ];
   WIsInt[ k ] = 0;
   }
  else
   WIsInt[ k ][ NInt[ k ] ] = Inf<Index>();
  }
 }  // end( SetIntVar() )

/*--------------------------------------------------------------------------*/

void Graph::SetExtraVars( cIndex NXV )
{
 if( NXtrV != NXV ) {
  delete[] U[ NComm + 1 ];
  delete[] U[ NComm ];

  delete[] C[ NComm ];

  delete[] WIsInt[ NComm ];
  WIsInt[ NComm ] = 0;

  NInt[ NComm ] = 0;

  NXtrV = NXV;
  if( NXtrV ) {
   C[ NComm ] = new CNumber[ NXtrV ];
   VectAssign( C[ NComm ] , CNumber( 0 ) , NXtrV );

   U[ NComm ] = new FNumber[ NXtrV ];
   VectAssign( U[ NComm ] , FNumber( 0 ) , NXtrV );

   U[ NComm + 1 ] = new FNumber[ NXtrV ];
   VectAssign( U[ NComm + 1 ] , F_INF , NXtrV );
   }
  else {
   C[ NComm ] = 0;
   U[ NComm ] = U[ NComm + 1 ] = 0;
   }
  }
 }  // end( Graph::SetExtraVars() )

/*--------------------------------------------------------------------------*/

void Graph::SetExtraConstr( cIndex NXC , int *IBeg , int *Indx ,
			    double *Vals )
{
 delete[] IdxBeg;
 delete[] CoefIdx;
 delete[] CoefVal;

 delete[] B[ NComm + 1 ];
 delete[] B[ NComm ];

 if( ( NXtrC = NXC ) ) {
  IdxBeg = IBeg;
  CoefIdx = Indx;
  CoefVal = Vals;

  B[ NComm ] = new FNumber[ NXtrC ];
  VectAssign( B[ NComm ] , FNumber( 0 ) , NXtrC );

  B[ NComm + 1 ] = new FNumber[ NXtrC ];
  VectAssign( B[ NComm + 1 ] , F_INF , NXtrC );
  }
 else {
  IdxBeg = CoefIdx = 0;
  CoefVal = 0;

  B[ NComm ] = B[ NComm + 1 ] = 0;
  }
 }

/*--------------------------------------------------------------------------*/
/*------------------------ PREPROCESSING METHODS ---------------------------*/
/*--------------------------------------------------------------------------*/

void Graph::PreProcess( cFNumber IncUk , cFNumber DecUk , cFNumber IncUjk ,
                        cFNumber DecUjk , cFNumber ChgDfct ,
                        cCNumber DecCsts )
{
 if( ChgDfct >= F_INF )
  throw( MMCFGException( "Graph::PreProcess: infinite ChgDfct" ) );
 if( DecCsts > C_INF )
  throw( MMCFGException( "Graph::PreProcess: infinite DecCsts" ) );
 
 // allocate (temporary) data structures- - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 Active = new Index[ NArcs ];

 FRow tmpv = new FNumber[ NComm ];
 Index_Set srck = new Index[ NComm ];

 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // squeeze rhss, declare arcs as "non-existent", etc.- - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 // ensure that all arcs entering/leaving a non-existent node do not exist- -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 for( Index k = 0 ; k < NComm ; k++ )
  for( Index i = 0 ; i < NArcs ; i++ )
   if( ( B[ k ][ Startn[ i ] - StrtNme ] == F_INF ) ||
       ( B[ k ][ Endn[ i ] - StrtNme ] == F_INF ) )
    C[ k ][ i ] = C_INF;

 // ensure that all non-existent arcs have zero capacity- - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 for( Index k = 0 ; k < NComm ; k++ )
  for( Index i = 0 ; i < NArcs ; i++ )
   if( C[ k ][ i ] == C_INF )
    U[ k ][ i ] = 0;

 // a *very* rough estimate of the max. flow across any arc is computed for
 // each commodity, and it is stored in tmpv[ k ] - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 FNumber maxU = 0;
 for( Index k = NComm ; k-- ; ) {
  // first, count the flow out of the sources

  Index srcs = 0;    // meanwhile, the sources are counted
  FNumber maxUk = 0;
  cFRow Bk = B[ k ];
  for( cFRow tBk = Bk + NNodes ; tBk-- > Bk ; )
   if( *tBk < 0 ) {
    srcs++;
    maxUk -= *tBk;
    }

  // now the contribution of arcs with potentially negative costs

  cFRow Uk = U[ k ];
  cFRow tUT = UTot + NArcs;
  CRow Ck = C[ k ] + NArcs;
  for( cFRow tUk = Uk + NArcs ; tUk-- > Uk ; ) {
   cFNumber tMF = min( *tUk , *(--tUT) );

   if( *(--Ck) < DecCsts ) {
    if( tMF >= F_INF )
     throw( MMCFGException(
	     "Graph::PreProcess: negative cost, infinite capacity" ) );
     maxUk += tMF;
     }
   }

  srck[ k ] = srcs;
  maxUk += ( ( NNodes + 1 ) / 2 ) * ChgDfct;  // count potential changes in
                                              // the deficits
  maxU += ( tmpv[ k ] = maxUk );
  }

 // detection of redundant mutual capacity constraints is attempted, and- - -
 // all the mutual capacity upper bounds are turned to finite values- - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 for( Index i = NCnst = 0 ; i < NArcs ; i++ ) {
  if( ( ! IncUk ) && ( ! UTot[ i ] ) ) {   // if mutual capacities can not
   for( Index k = NComm ; k-- ; ) {        // increase, and UTot[] == 0 ...
    C[ k ][ i ] = C_INF;          // ... this arc does not exist
    U[ k ][ i ] = 0;
    }

   continue;
   }

  if( DecUk == F_INF ) {     // all mutual capacity constraints exist
   if( UTot[ i ] == F_INF )  // but those that are declared non-so
    UTot[ i ] = maxU;                 // ensure that UTot is "finite" anyway
   else
    Active[ NCnst++ ] = i;

   continue;
   }

  // compute is an upper bound on the max quantity of flow (of any commodity)
  // on arc i: if capacities can increase indefinitely, the only bound is
  // the total quantity of flow in the graph

  FNumber Ui = 0;
  if( IncUjk < F_INF )
   for( Index k = NComm ; k-- ; )
    if( U[ k ][ i ] == F_INF )
     Ui += tmpv[ k ];
    else
     Ui += min( tmpv[ k ] , U[ k ][ i ] + IncUjk );
  else
   Ui = maxU;

  // note: when e.g. the mutual capacity and the sum of all the individual
  // capacities of an arc are identical, the arc is marked as "inactive"; this
  // is an arbitrary choice, since one could as well keep it and eliminate all
  // the individual capacities

  if( UTot[ i ] >= Ui - DecUk )
   UTot[ i ] = Ui;
  else
   Active[ NCnst++ ] = i;

  }  // end for( i )

 if( NCnst < NArcs )
  Active[ NCnst ] = Inf<Index>();

 // now a squeeze of single-commodity capacities is attempted, and SPTs are -
 // definitively recognized - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // meanwhile, construct the "active" individual capacity constraints

 for( Index k = 0 ; k < NComm ; k++ ) {
  cIndex_Set tA = Active;

  delete[] ActiveK[ k ];
  ActiveK[ k ] = new Index[ NArcs ];
  Index_Set tAK = ActiveK[ k ];
  Index cnt = 0;  // active individual capacity constraints

  for( Index i = 0 ; i < NArcs ; i++ ) {
   bool Ai = ( *tA == i );    // true if arc i is "active"
   if( Ai )
    tA++;

   if( C[ k ][ i ] == C_INF )  // a non-existent arc
    continue;

   if( ( ! IncUjk ) && ( ! U[ k ][ i ] ) ) {
    // an arc that can be declared non-existent by its capacity
    // (that will never increase)
    C[ k ][ i ] = C_INF;
    continue;
    }

   if( DecUjk < F_INF ) {
    // if individual capacities cannot decrease forever, then the
    // individual capacity constraint of some existing arc can be
    // declared redundant

    if( U[ k ][ i ] >= tmpv[ k ] + DecUjk ) {
     // the constraint is redundant because there will never be that much
     // flow in the graph

     U[ k ][ i ] = min( tmpv[ k ] , UTot[ i ] );  // give it a "nice"
     continue;                                    // finite value anyway
     }

    if( ( IncUk < F_INF ) &&
	( Ai && ( U[ k ][ i ] >= UTot[ i ] + IncUk + DecUjk ) ) ) {
     // if mutual capacities cannot increase forever, some individual
     // capacities may be declared redundant by the mutual capacity
     // note that the mutual capacity of an arc can be used to declare
     // that the individual capacity is redundant only if the arc is
     // "active", as "inactive" arcs precisely mean that no mutual
     // capacity constraint is imposed on them (i.e., the value of
     // UTot[ i ] is not really meaningful and can be ignored)

     U[ k ][ i ] = UTot[ i ];  // give it a "nice" finite value anyway
     continue;
     }
    }

   *(tAK++) = i;
   cnt++;

   }  // end for( i )

  if( ( ! cnt ) && ( srck[ k ] == 1 ) )
   PT[ k ] = kSPT;

  NamesK[ k + 1 ] = NamesK[ k ] + cnt;

  if( cnt >= NArcs ) {  // all individual capacity constraints are active
   delete[] ActiveK[ k ];
   ActiveK[ k ] = 0;
   }
  else {                // some are active, some are not
   tAK = new Index[ cnt + 1 ];
   VectAssign( tAK , ActiveK[ k ] , cnt );
   tAK[ cnt ] = Inf<Index>();
   delete[] ActiveK[ k ];
   ActiveK[ k ] = tAK;
   }
  }   // end for( k )

 if( NCnst >= NArcs ) {
  delete[] Active;
  Active = 0;
  }

 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // find and eliminate redundancies in the data structures- - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 // examine B[] - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 BIsCpy = new bool[ NComm ];
 VectAssign( BIsCpy , bool( false ) , NComm );

 bool cpy = false;
 for( Index k = 1 ; k < NComm ; k++ )
  for( Index i = 0 ; i < k ; i++ )
   if( ( ! BIsCpy[ i ] ) && EqualVect( B[ k ] , B[ i ] , NNodes ) ) {
    BIsCpy[ k ] = cpy = true;
    delete[] B[ k ];
    B[ k ] = B[ i ];
    break;
    }

 if( ! cpy ) {
  delete[] BIsCpy;
  BIsCpy = 0;
  }

 // examine U[] and UTot- - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 UIsCpy = new bool[ NComm ];
 VectAssign( UIsCpy , bool( false ) , NComm );

 cpy = false;
 if( EqualVect( U[ 0 ] , UTot , NArcs ) ) {
  UIsCpy[ 0 ] = cpy = true;
  delete[] U[ 0 ];
  U[ 0 ] = UTot;
  }

 for( Index k = 1 ; k < NComm ; k++ )
  for( Index i = 0 ; i < k ; i++ )
   if( ( ! UIsCpy[ i ] ) && EqualVect( U[ k ] , U[ i ] , NArcs ) ) {
    UIsCpy[ k ] = cpy = true;
    delete[] U[ k ];
    U[ k ] = U[ i ];
    break;
    }

 if( ! cpy ) {
  delete[] UIsCpy;
  UIsCpy = 0;
  }

 // examine C[] - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 CIsCpy = new bool[ NComm ];
 VectAssign( CIsCpy , bool( false ) , NComm );

 cpy = false;
 for( Index k = 1 ; k < NComm ; k++ )
  for( Index i = 0 ; i < k ; i++ )
   if( ( ! CIsCpy[ i ] ) && EqualVect( C[ k ] , C[ i ] , NArcs ) ) {
    CIsCpy[ k ] = cpy = true;
    delete[] C[ k ];
    C[ k ] = C[ i ];
    break;
    }

 if( ! cpy ) {
  delete[] CIsCpy;
  CIsCpy = 0;
  }

 // cleanup - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 delete[] srck;
 delete[] tmpv;

 }  // end( Graph::PreProcess )

/*--------------------------------------------------------------------------*/

void Graph::MakeSingleSourced( bool ToAll )
{
 Index_Set NewArcs = new Index[ NNodes ];
 Index_Set CmmStts = new Index[ NComm ];

 // count the sources and set NNewArcs == Inf<Index>() if there is at least one
 // commodity with more than two sources- - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 Index NNewArcs = 0;
 Index k = NComm;

 for( ; k-- ; ) {
  CmmStts[ k ] = NNodes;

  for( Index i = NNodes ; i-- ; )
   if( B[ k ][ i ] < 0 )
    if( CmmStts[ k ] == NNodes )
     CmmStts[ k ] = i;
    else {
     CmmStts[ k ] = NNewArcs = Inf<Index>();
     break;
     }

  if( CmmStts[ k ] != NNodes )  // not a circulation subproblem
   PT[ k ] = kSPT;
  }

 // now CmmStts[ k ] contains Inf<Index>() if the commodity has more than 1
 // source, NNodes if it has no sources and its only source name otherwise
 // phase two: construct the new arcs, if necessary - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( ToAll )
  for( Index i = NNewArcs = NNodes ; i-- ; )
   NewArcs[ i ] = i + NArcs;
 else {
  if( ! NNewArcs ) {
   delete[] CmmStts;
   delete[] NewArcs;
   return;
   }

  for( k = NNodes ; k-- ; )
   NewArcs[ k ] = 0;

  for( k = NComm ; k-- ; )
   if( CmmStts[ k ] == Inf<Index>() )
    for( Index i = NNodes ; i-- ; )
     if( B[ k ][ i ] < 0 )
      NewArcs[ i ]++;

  // now NewArcs[ i ] counts how many times i is a source for a commodity,
  // excluded those that have *only* i as a source

  for( NNewArcs = k = 0 ; k < NNodes ; k++ )
   if( NewArcs[ k ] )
    NewArcs[ k ] = NArcs + NNewArcs++;
  }

 // now NewArcs[ i ] contains the name of the arc super-source -> i, if any
 // resize the data structures- - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 for( k = 0 ; k < NComm ; ) {
  CRow NewCk = new CNumber[ NArcs + NNewArcs ];

  VectAssign( NewCk , C[ k ] , NArcs );
  delete[] C[ k ];
  C[ k++ ] = NewCk;
  }

 for( k = 0 ; k < NComm ; ) {
  FRow NewUk = new FNumber[ NArcs + NNewArcs ];

  VectAssign( NewUk , U[ k ] , NArcs );
  delete[] U[ k ];
  U[ k++ ] = NewUk;
  }

 for( k = 0 ; k < NComm ; ) {
  FRow NewBk = new FNumber[ NNodes + 1 ];
  NewBk[ NNodes ] = 0;

  VectAssign( NewBk , B[ k ] , NNodes );
  delete[] B[ k ];
  B[ k++ ] = NewBk;
  }

 Index_Set tNdp = new Index[ NArcs + NNewArcs ];
 VectAssign( tNdp , Startn , NArcs );
 delete[] Startn;
 Startn = tNdp;

 tNdp = new Index[ NArcs + NNewArcs ];
 VectAssign( tNdp , Endn , NArcs );
 delete[] Endn;
 Endn = tNdp;

 FRow tUt = new FNumber[ NArcs + NNewArcs ];
 VectAssign( tUt , UTot , NArcs );
 delete[] UTot;
 UTot = tUt;

 // construct the new arcs- - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 for( k = 0 ; k < NNodes ; k++ )
  if( NewArcs[ k ] ) {
   Startn[ NewArcs[ k ] ] = NNodes + StrtNme;
   Endn[ NewArcs[ k ] ] = k + StrtNme;
   UTot[ NewArcs[ k ] ] = F_INF;
   }

 for( k = 0 ; k < NComm ; k++ ) {
  FRow Uk = U[ k ];
  CRow Ck = C[ k ];

  if( CmmStts[ k ] < NNodes ) {  // single-sourced commodity
   Uk += NArcs;
   Ck += NArcs;

   for( Index i = NNewArcs ; i-- ; )
   {
    *(Uk++) = 0;
    *(Ck++) = C_INF;
    }
   }
  else {  // multiple-sourced or no-sourced commodity
   FRow Bk = B[ k ];

   for( Index i = NNodes ; i-- ; ) {
    Index h = NewArcs[ i ];

    if( h )
     if( Bk[ i ] < 0 ) {
      Bk[ NNodes ] += Bk[ i ];
      Uk[ h ] = - Bk[ i ];
      Ck[ h ] = 0;
      Bk[ i ] = 0;
      }
     else {
      Uk[ h ] = 0;
      Ck[ h ] = C_INF;
      }
    }  // end for( i )
   }
  }  // end for( k )

 delete[] CmmStts;
 delete[] NewArcs;

 NArcs += NNewArcs;
 NNodes++;

 }  // end( MakeSingleSourced );

/*--------------------------------------------------------------------------*/
/*---------------------- MISCELLANEOUS METHODS -----------------------------*/
/*--------------------------------------------------------------------------*/

void Graph::OutMPSFile( const char *const FN , bool FxdClmns ,
                        const char *const PN )
{
 if( ! DrctdPrb )
  throw( MMCFGException(
	  "Graph::OutMPSFile: undirected graphs not supported yet" ) );
 
 // opening file - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 fstream of( FN , ios::out );

 // writing problem name - - - - - - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 of << "NAME      ";
 if( PN )
  of << PN;

 // writing ROWS section: flow balance constraints - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // note: the flow balance constraint for node i of commodity k is named
 // "f<k>_<i>" independently from the fact that there may be un-existent
 // nodes; thus, there may be "holes" in the constraint names

 of << endl << "ROWS" << endl << " N  obj" << endl;

 for( Index k = 0 ; k < NComm ; k++ )
  for( Index i = 0 ; i < NNodes ; i++ )
   if( B[ k ][ i ] < F_INF )
    of << " E  f" << k << "_" << i + StrtNme << endl;

 // writing ROWS section: mutual capacity constraints- - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // note: the mutual capacity constraint for arc i is named "m<i>"
 // independently from the fact that there may be "inactive" arcs; thus,
 // there may be "holes" in the constraint names

 if( Active ) {
  cIndex_Set tA = Active;
  for( Index i ; ( i = *(tA++) ) < Inf<Index>() ; )
   of << " L  m" << i << endl;
  }
 else
  for( Index i = 0 ; i < NArcs ; i++ )
   of << " L  m" << i << endl;

 // writing COLUMNS section- - - - - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // note: the flow on arc i of commodity k is named "x<k>_<i>" independently
 // from the fact that there may be un-existent arcs; thus, there may be
 // "holes" in the column names

 of << "COLUMNS" << endl;

 if( FxdClmns )  // fixed-columns format - - - - - - - - - - - - - - - - - -
 {               //- - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  of.setf( ios::fixed , ios::floatfield );
  of.setf( ios::left , ios::adjustfield );
  of.setf( ios::showpoint );
  of.precision( 9 );

  for( Index k = 0 ; k < NComm ; k++ ) {
   cIndex_Set tA = Active;
   for( Index i = 0 ; i < NArcs ; i++ ) {
    if( C[ k ][ i ] < C_INF ) {
     of.width( 4 );
     of << "    x" << k << "_" << i;

     of << " obj     ";
     if( C[ k ][ i ] >= 0 )
      of << " ";
 
     of.width( 9 );
     of << C[ k ][ i ] << "     f";
     of.width( 4 );
     of << k << "_" << Startn[ i ] << "       -1." << endl << "    x" << k
	<< "_" << i << " f" << k << "_" << Endn[ i ] << "        1.";

     if( ( ! tA ) || ( *tA == i ) ) {
      of.width( 9 );
      of << "     m" << i << "        1.";
      }

     of << endl;

     }  // end( if( C[ k ][ i ] < C_INF ) )

    if( tA && ( *tA == i ) )
     tA++;

    }  // end( for( i ) )
   }  // end( for( k ) )

  // writing RHS section: flow balance constraints - - - - - - - - - - - - -

  bool first = true;

  of.precision( 8 );
  of << "RHS" << endl;

  for( Index k = 0 ; k < NComm ; k++ )
   for( Index i = 0 ; i < NNodes ; i++ )
    if( B[ k ][ i ] < F_INF ) {
     if( first )
      of << "    rhs ";

     of << "     f";
     of.width( 4 );
     of << k << "_" << i + StrtNme << " ";
     if( B[ k ][ i ] >= 0 )
      of << " ";
     of.width( 9 );
     of << B[ k ][ i ];

     if( ! first )
      of << endl;

     first = first ? false : true;
     }

  // writing RHS section: mutual capacity constraints - - - - - - - - - - -

  of.precision( 9 );
  of.width( 9 );

  cIndex_Set tA = Active;
  for( Index i = 0 ; i < NArcs ; i++ ) {
   if( tA )
    if( *tA == i )
     tA++;
    else
     continue;

   if( first )
    of << "    rhs ";

   of << "     m" << i << UTot[ i ];

   if( ! first )
    of << endl;

   first = first ? false : true;
   }

  if( ! first )
   of << endl;

  // writing BOUNDS section - - - - - - - - - - - - - - - - - - - - - - - - -

  of << "BOUNDS" << endl;

  for( Index k = 0 ; k < NComm ; k++ ) {
   cIndex_Set tAK = ActiveK[ k ];
   for( Index i = 0 ; i < NArcs ; i++ ) {
    if( tAK )
     if( *tAK == i )
      tAK++;
     else
      continue;

    if( C[ k ][ i ] == C_INF )
     continue;

    of << " UP BOUND     x";
    of.width( 4 );
    of << k << "_" << i << " ";
    of.width( 9 );
    of << U[ k ][ i ] << endl;
    }
   }
  }
 else  // "modern" format, non-fixed-columns- - - - - - - - - - - - - - - - -
 {     // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  for( Index k = 0 ; k < NComm ; k++ ) {
   cIndex_Set tA = Active;
   for( Index i = 0 ; i < NArcs ; i++ ) {
    if( C[ k ][ i ] < C_INF ) {
     of << " x" << k << "_" << i << "\tobj\t" << C[ k ][ i ] << "\tf" << k
	<< "_" << Startn[ i ] << "\t-1" << endl << " x" << k << "_" << i
	<< "\tf" << k << "_" << Endn[ i ] << "\t1";

     if( ( ! tA ) || ( *tA == i ) )
      of << "\tm" << i << "\t1";

     of << endl;
     }

    if( tA && ( *tA == i ) )
     tA++;
    }
   }

  // writing RHS section: flow balance constraints - - - - - - - - - - - - -

  bool first = true;
  of << "RHS" << endl;

  for( Index k = 0 ; k < NComm ; k++ )
   for( Index i = 0 ; i < NNodes ; i++ )
    if( B[ k ][ i ] && ( B[ k ][ i ] < F_INF ) ) {
     if( first )
      of << " rhs";

     of << "\tf" << k << "_" << i + StrtNme << "\t" << B[ k ][ i ];

     if( ! first )
      of << endl;

     first = first ? false : true;
     }

  // writing RHS section: mutual capacity constraints - - - - - - - - - - -

  if( Active ) {
   cIndex_Set tA = Active;
   for( Index i ; ( i = *(tA++) ) < Inf<Index>() ; ) {
    if( first )
     of << " rhs";

    of << "\tm" << i << "\t" << UTot[ i ];

    if( ! first )
     of << endl;

    first = first ? false : true;
    }
   }
  else
   for( Index i = 0 ; i < NArcs ; i++ ) {
    if( first )
     of << " rhs";

    of << "\tm" << i << "\t" << UTot[ i ];

    if( ! first )
     of << endl;

    first = first ? false : true;
    }

  if( ! first )
   of << endl;

  // writing BOUNDS section - - - - - - - - - - - - - - - - - - - - - - - - -

  of << "BOUNDS" << endl;

  for( Index k = 0 ; k < NComm ; k++ )
   if( ActiveK[ k ] ) {
    cIndex_Set tAK = ActiveK[ k ];
    for( Index i ; ( i = *(tAK++) ) < Inf<Index>() ; )
     if( C[ k ][ i ] < C_INF )
      of << " UP bound\tx" << k << "_" << i << "\t" << U[ k ][ i ] << endl;
     }
   else
    for( Index i = 0 ; i < NArcs ; i++ )
     if( C[ k ][ i ] < C_INF )
      of << " UP bound\tx" << k << "_" << i << "\t" << U[ k ][ i ] << endl;

  }  // end else( FxdClmns )- - - - - - - - - - - - - - - - - - - - - - - - -
     // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 // end - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 of << "ENDATA" << endl << endl;

 }  // end( OutMPSFile )

/*--------------------------------------------------------------------------*/

Graph::FONumber Graph::UpperBound( void )
{
 FONumber UB = 0;

 for( Index k = NComm ; k-- ; ) {
  // calculate the total amount of the k-th commodity - - - - - - - - - - -

  FNumber Fk = 0;
  cFRow dfcts = B[ k ];
  Index i = NNodes;

  for( ; i-- ; )
   if( dfcts[ i ] > 0 )
    Fk += dfcts[ i ];

  // now the contribution of commodity k to the UB - - - - - - - - - - - - -

  cFRow Uk = U[ k ];
  cCRow Ck = C[ k ];

  for( dfcts = UTot , i = NArcs ; i-- ; )
   if( Ck[ i ] > 0 )
    UB += Ck[ i ] * min( Fk , min( dfcts[ i ] , Uk[ i ] ) );

  }  // end for( k )

 if( NXtrV ) {  // add the contribution of the "extra" variables - - - - - -
  cCRow Ce = C[ NComm ];
  cFRow LBe = U[ NComm ];
  cFRow UBe = U[ NComm + 1 ];
  for( Index i = NXtrV ; i-- ; LBe++ , UBe++ ) {
   cCNumber Cei = *(Ce++);
   if( Cei > 0 )
    UB += Cei * (*UBe);
   else
    UB += Cei * (*LBe);
   }
  }

 return( UB );

 }  // end( UpperBound )

/*--------------------------------------------------------------------------*/
/*------------------------------ DESTRUCTOR --------------------------------*/
/*--------------------------------------------------------------------------*/

Graph::~Graph()
{
 delete[] IdxBeg;
 delete[] CoefIdx;
 delete[] CoefVal;

 for( Index k = NComm ; k-- ; )
  delete[] ActiveK[ k ];

 delete[] ActiveK;
 delete[] NamesK;

 delete[] Active;

 delete[] WIsInt;
 delete[] NInt;

 delete[] UTot;

 delete[] Endn;
 delete[] Startn;

 // deleting B[]- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( BIsCpy ) {
  for( Index k = NComm ; k-- ; )
   if( ! BIsCpy[ k ] )
    delete[] B[ k ];

  delete[] BIsCpy;
  }
 else
  for( Index k = NComm ; k-- ; )
   delete[] B[ k ];

 delete[] B;


 // deleting U[]- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 delete[] U[ NComm + 1 ];
 delete[] U[ NComm ];

 if( UIsCpy ) {
  for( Index k = NComm ; k-- ; )
   if( ! UIsCpy[ k ] )
    delete[] U[ k ];

  delete[] UIsCpy;
  }
 else
  for( Index k = NComm ; k-- ; )
   delete[] U[ k ];

 delete[] U;

 // deleting C[]- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 delete[] C[ NComm ];

 if( CIsCpy ) {
  for( Index k = NComm ; k-- ; )
   if( ! CIsCpy[ k ] )
    delete[] C[ k ];

  delete[] CIsCpy;
  }
 else
  for( Index k = NComm ; k-- ; )
   delete[] C[ k ];

 delete[] C;

 delete[] PT;

 }  // end( ~Graph )

/*--------------------------------------------------------------------------*/
/*-------------------------- PROTECTED METHODS -----------------------------*/
/*--------------------------------------------------------------------------*/

inline void Graph::CmnIntlz( void )
{
 // some initializations that are common to all the constructors- - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 StrtNme = 1;
 Active = 0;
 DrctdPrb = true;

 PT = new MCFType[ NComm ];

 Index k = NComm;
 for( ; k-- ; )
  PT[ k ] = kMCF;

 CIsCpy = UIsCpy = DIsCpy = BIsCpy = 0;

 // find arcs that might have individual capacity constraints - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 ActiveK = new Index_Set[ NComm ];
 NamesK  = new Index[ NComm + 1 ];

 for( *NamesK = NCnst , k = 0 ; k < NComm ; k++ ) {
  // first, count how many active constraints there are - - - - - - - - - - -

  Index i = 0;
  Index cnt = 0;
  cCRow Ck = C[ k ];
  cFRow Uk = U[ k ];
  for( ; i++ < NArcs ; Uk++ )
   if( ( *(Ck++) < C_INF ) && ( *Uk < F_INF ) )
    cnt++;

  // second, (if necessary) construct the actual vector of indices - - - - - -

  NamesK[ k + 1 ] = NamesK[ k ] + cnt;

  if( cnt < NArcs ) {
   ActiveK[ k ] = new Index[ cnt + 1 ];

   Index_Set tAKk = ActiveK[ k ];
   for( Ck = C[ k ] , Uk = U[ k ] , i = 0 ; i < NArcs ; i++ , Uk++ )
    if( ( *(Ck++) < C_INF ) && ( *Uk < F_INF ) )
     *(tAKk++) = i;

   *tAKk = Inf<Index>();
   }
  else
   ActiveK[ k ] = 0;

  }  // end( for( k ) )

 }  // end( CmnIntlz )

/*--------------------------------------------------------------------------*/
/*-------------------------- End File Graph.C ------------------------------*/
/*--------------------------------------------------------------------------*/

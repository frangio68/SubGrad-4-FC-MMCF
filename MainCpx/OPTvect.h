/*--------------------------------------------------------------------------*/
/*----------------------- File OPTvect.h -----------------------------------*/
/*--------------------------------------------------------------------------*/
/*--                                                                      --*/
/*--  A bunch of little useful template inline functions.                 --*/
/*--                                                                      --*/
/*--  Scalar functions:                                                   --*/
/*--                                                                      --*/
/*--   ABS(), sgn(), min(), max(), Swap(), CeilDiv()                      --*/
/*--                                                                      --*/
/*--  Number-returning vector functions:                                  --*/
/*--                                                                      --*/
/*--   Norm(), OneNorm(), INFNorm(), SumV()                               --*/
/*--   MaxVecV(), MinVecV(), MaxVecI(), MinVecI()                         --*/
/*--   ScalarProduct(), ScalarProduct[B[B]]()                             --*/
/*--                                                                      --*/
/*--  Vector operations:                                                  --*/
/*--                                                                      --*/
/*--   Vect[M]Assign[B[B]](), VectSum[B[B]](), VectSubtract[B[B]](),      --*/
/*--   Vect[I]Scale[B[B]](), VectAdd(), VectDiff(), VectMult(),           --*/
/*--   VectDivide(), VectXcg[B[B]]()                                      --*/
/*--                                                                      --*/
/*--  Sparse/dense vector transformations:                                --*/
/*--                                                                      --*/
/*--   Sparsify(), SparsifyT(), SparsifyAT(), Densify(), Compact()        --*/
/*--                                                                      --*/
/*--  Array manipulation:                                                 --*/
/*--                                                                      --*/
/*--   Merge(), ShiftVect(), RotateVect(), ShiftRVect(), RotateRVect()    --*/
/*--                                                                      --*/
/*--  Searching:                                                          --*/
/*--                                                                      --*/
/*--   Match(), EqualVect()                                               --*/
/*--   BinSearch(), BinSearch1(), BinSearch2()                            --*/
/*--   HeapIns(), HeapDel()                                               --*/
/*--                                                                      --*/
/*--                             Version 2.61                             --*/
/*--                            08 - 08 - 2012                            --*/
/*--                                                                      --*/
/*--                 Original Idea and Implementation by:                 --*/
/*--                                                                      --*/
/*--                          Antonio Frangioni                           --*/
/*--                                                                      --*/
/*--                       Operations Research Group                      --*/
/*--                      Dipartimento di Informatica                     --*/
/*--                          Universita' di Pisa                         --*/
/*--                                                                      --*/
/*--             Copyright (C) 1992 - 2012 by Antonio Frangioni           --*/
/*--                                                                      --*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __OPTvect
 #define __OPTvect  /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "OPTtypes.h"

/*--------------------------------------------------------------------------*/
/*--------------------------- NAMESPACE ------------------------------------*/
/*--------------------------------------------------------------------------*/

#if( OPT_USE_NAMESPACES )
namespace OPTtypes_di_unipi_it
{
#endif

/*--------------------------------------------------------------------------*/
/*----------------------------- FUNCTIONS ----------------------------------*/
/*--------------------------------------------------------------------------*/
/*--                                                                      --*/
/*--  Scalar functions, like ABS(), min(), max(), sgn() ...               --*/
/*--                                                                      --*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

template<class T>
inline T ABS( const T x )
{
 return( x >= T( 0 ) ? x : -x );
 }

/*--------------------------------------------------------------------------*/

template<class T>
inline T sgn( const T x )
{
 return( x ? ( x > T( 0 ) ? T( 1 ) : T( -1 ) ) : T( 0 ) );
 }

/*--------------------------------------------------------------------------*/
/* max() and min() are already defined on gcc 3.x */

#ifndef __GNUC__

template<class T>
inline T min( const T x , const T y )
{ 
 return( x <= y ? x : y );
 }

template<class T>
inline T max( const T x , const T y )
{
 return ( x >= y ? x : y );
 }

#else
#if __GNUC__ < 3

template<class T>
inline T min( const T x , const T y )
{ 
 return( x <= y ? x : y );
 }

template<class T>
inline T max( const T x , const T y )
{
 return ( x >= y ? x : y );
 }

#endif
#endif

/*--------------------------------------------------------------------------*/

template<class T>
inline void Swap( T &v1 , T &v2 )
{
 const T temp = v1;

 v1 = v2;
 v2 = temp;
 }

/*--------------------------------------------------------------------------*/

template<class T1, class T2>
inline T1 CeilDiv( const T1 x , const T2 y )
{
 // returns the ceiling of (smallest integer number not smaller than) x / y,
 // which have both to be integer types because the `%' operation is used

 T1 temp = x / y;
 if( x % y )
  temp++;

 return( temp );
 }

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*--                                                                      --*/
/*-- Vector functions: norms, scalar products, assignments of vectors,    --*/
/*-- multiplication by a scalar and so on. The types involved in these    --*/
/*-- operations must support the elementary arithmetical operations '+',  --*/
/*-- '-', '*' and '/', as well as the assigment '='.                      --*/
/*--                                                                      --*/
/*-- Vectors are usually just a set of n (a parameter of the function)    --*/
/*-- consecutive elements of a certain type; however, in many cases it    --*/
/*-- is useful to "restrict" them to a subset of their entries. Hence,    --*/
/*-- (read-only) vectors of indices (cIndex_Set) play a special role; for --*/
/*-- a vector g, with g{B} (where B is a vector of indices) we indicate   --*/
/*-- the "restricted" vector [ g[ B[ i ] ]. A typical reason for dealing  --*/
/*-- with "restricted" vectors is that "sparse" vectors, those having few --*/
/*-- nonzeroes w.r.t. their lenght, are very common.                      --*/
/*--                                                                      --*/
/*-- Vector of indices are meant to be "infinity-terminated", i.e., an    --*/
/*-- Inf<Index>() must be found immediately after the last "valid" Index, and    --*/
/*-- sometimes they are required to be ordered, typically in increasing   --*/
/*-- sense.                                                               --*/
/*--                                                                      --*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*--       Number-returning functions (norms, scalar products ...)        --*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

template<class T>
inline T Norm( const T *g , Index n )
{
 // returns Sum{ i = 0 .. n - 1 } g[ i ]^2 = g * g

 T t = 0;
 for( ; n-- ; ) {
  const T tmp = *(g++);
  t += tmp * tmp;
  }

 return( t );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

template<class T>
inline T Norm( const T *g , cIndex_Set B )
{
 // Norm( g{B} )

 T t = 0;
 for( Index h ; ( h = *(B++) ) < Inf<Index>() ; ) {
  const T tmp = g[ h ];
  t += tmp * tmp;
  }

 return( t );
 }

/*--------------------------------------------------------------------------*/

template<class T>
inline T OneNorm( const T *g , Index n )
{
 // returns Sum{ i = 0 .. n - 1 } ABS( g[ i ] )

 T t = 0;
 for( ; n-- ; )
  t += ABS( *(g++) );

 return( t );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

template<class T>
inline T OneNorm( const T *g , cIndex_Set B )
{
 // OneNorm( g{B} )

 T t = 0;
 for( Index h ; ( h = *(B++) ) < Inf<Index>() ; )
  t += ABS( g[ h ] );

 return( t );
 }

/*--------------------------------------------------------------------------*/

template<class T>
inline T INFNorm( const T *g , Index n )
{
 // returns Max{ i = 0 .. n - 1 } ABS( g[ i ] )

 T t = 0;
 for( ; n-- ; ) {
  const T tmp = *(g++);
  if( t < tmp )
   t = tmp;
  }

 return( t );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

template<class T>
inline T INFNorm( const T *g , cIndex_Set B )
{
 // INFNorm( g{B} )

 T t = 0;
 for( Index h ; ( h = *(B++) ) < Inf<Index>() ; ) {
  const T tmp = g[ h ];
  if( t < tmp )
   t = tmp;
  }

 return( t );
 }

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

template<class T>
inline T SumV( const T *g , Index n )
{
 // returns Sum{ i = 0 .. n - 1 } g[ i ]

 T t = 0;
 for( ; n-- ; )
  t += *(g++);

 return( t );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

template<class T>
inline T SumV( const T *g , cIndex_Set B )
{
 // SumV( g{B} )

 T t = 0;
 for( Index h ; ( h = *(B++) ) < Inf<Index>() ; )
  t += g[ h ];

 return( t );
 }

/*--------------------------------------------------------------------------*/

template<class T>
inline T MaxVecV( const T *g , Index n )
{
 // returns Max{ i = 0 .. n - 1 } g[ i ]; n *must* be > 0

 T max = *g;
 for( ; --n ; )
  if( *(++g) > max )
   max = *g;

 return( max );
 }

/*--------------------------------------------------------------------------*/

template<class T>
inline T MinVecV( const T *g , Index n )
{
 // returns Min{ i = 0 .. n - 1 } g[ i ]; n *must* be > 0

 T min = *g;
 for( ; --n ; )
  if( *(++g) < min )
   min = *g;

 return( min );
 }

/*--------------------------------------------------------------------------*/

template<class T>
inline Index MaxVecI( const T *g , cIndex n )
{
 // returns the *index* of the maximum element of the n-vector v

 T max = *g;
 Index maxi = 0;
 for( Index i = maxi ; ++i < n ; )
  if( *(++g) > max )
  {
   maxi = i;
   max = *g;
   }

 return( maxi );
 }

/*--------------------------------------------------------------------------*/

template<class T>
inline Index MinVecI( const T *g , cIndex n )
{
 // returns the *index* of the minimum element of the n-vector v

 T min = *g;
 Index mini = 0;
 for( Index i = mini ; ++i < n ; )
  if( *(++g) < min )
  {
   mini = i;
   min = *g;
   }

 return( mini );
 }

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

template<class T1, class T2>
inline T1 ScalarProduct( const T1 *g1 , const T2 *g2 , Index n )
{
 // returns Sum{ i = 0 .. n - 1 } g1[ i ] * g2[ i ]

 T1 t = 0;
 for( ; n-- ; )
  t += (*(g1++)) * (*(g2++));

 return( t );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

template<class T1, class T2>
inline T1 ScalarProduct( const T1 *g1 , const T2 *g2 , cIndex_Set B )
{
 // ScalarProduct( g1 , g2{B} )

 T1 t = 0;
 for( Index h ; ( h = *(B++) ) < Inf<Index>() ; )
  t += (*(g1++)) * g2[ h ];

 return( t );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

template<class T1, class T2>
inline T1 ScalarProduct( const T1 *g1 , cIndex_Set B1 ,
                         const T2 *g2 , cIndex_Set B2 )
{
 // ScalarProduct( g1{B} , g2{B} ), B = intersection of B1 and B2

 T1 t = 0;
 Index h = *B1;
 Index k = *B2;
 if( ( h < Inf<Index>() ) && ( k < Inf<Index>() ) )
  for(;;)
   if( k < h )
   {
    if( ( k = *(++B2) ) == Inf<Index>() )
     break;
    g2++;
    }
   else
    if( h < k )
    {
     if( ( h = *(++B1) ) == Inf<Index>() )
      break;
     g1++;
     }
    else
    {
     t += (*(g1++)) * (*(g2++));
     if( ( k = *(++B2) ) == Inf<Index>() )
      break;
     if( ( h = *(++B1) ) == Inf<Index>() )
      break;
     }

 return( t );

 }  // end( ScalarProduct( g1 , B1 , g2 , B2 )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

template<class T1, class T2>
inline T1 ScalarProductB( const T1 *g1 , const T2 *g2 , cIndex_Set B )
{
 // ScalarProduct( g1{B} , g2 )

 T1 t = 0;
 for( Index h ; ( h = *(B++) ) < Inf<Index>() ; )
  t += g1[ h ] * (*(g2++));

 return( t );
 }

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

template<class T1, class T2>
inline T1 ScalarProductBB( const T1 *g1 , const T2 *g2 , cIndex_Set B )
{
 // ScalarProduct( g1{B} , g2{B} )

 T1 t = 0;
 for( Index h ; ( h = *(B++) ) < Inf<Index>() ; )
  t += g1[ h ] * g2[ h ];

 return( t );
 }

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*-- Vector operations: sum/difference of vectors, multiplication by a    --*/
/*-- scalar, assignment, setting, scaling ...                             --*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

template<class T>
inline void VectAssign( T *const g , const T x , cIndex n )
{
 // g[ i ] = x for each i = 0 .. n - 1

 for( T *tg = g + n ; tg > g ; )
  *(--tg) = x;
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

template<class T>
inline void VectAssign( T *const g , const T x , cIndex_Set B )
{
 // g{B} = x, all other entries of g unchanged

 for( Index h ; ( h = *(B++) ) < Inf<Index>() ; )
  g[ h ] = x;
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

template<class T1, class T2>
inline void VectAssign( T1 *g1 , const T2 *g2 , Index n )
{
 // g1 := g2

 for( ; n-- ; )
  *(g1++) = *(g2++);
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

inline Index_Set VectAssign( Index_Set g1 , cIndex_Set g2 )
{
 // special version of the above for g2 an Inf<Index>()-terminated vector of indices
 // do not write the terminating Inf<Index>(), but returns the pointer to the
 // position where it should be written

 for( Index h ; ( h = *(g2++) ) < Inf<Index>() ; )
  *(g1++) = h;

 return( g1 );
 }  

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

template<class T1, class T2>
inline void VectAssign( T1 *g1 , const T2 *g2 , cIndex_Set B )
{
 // g1 = g2{B}

 for( Index h ; ( h = *(B++) ) < Inf<Index>() ; )
  *(g1++) = g2[ h ];
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

template<class T1, class T2>
inline void VectAssign( T1 *g1 , cIndex_Set B1 , const T2 *g2 , cIndex_Set B2 )
{
 // g1{B} = g2{B}, B = intersection of B1 and B2, all other unchanged

 Index h = *B1;
 Index k = *B2;
 if( ( h < Inf<Index>() ) && ( k < Inf<Index>() ) )
  for(;;)
   if( k < h )
   {
    if( ( k = *(++B2) ) == Inf<Index>() )
     break;
    g2++;
    }
   else
    if( h < k )
    {
     if( ( h = *(++B1) ) == Inf<Index>() )
      break;
     g1++;
     }
    else
    {
     (*(g1++)) = (*(g2++));
     if( ( k = *(++B2) ) == Inf<Index>() )
      break;
     if( ( h = *(++B1) ) == Inf<Index>() )
      break;
     }

 }  // end( VectAssign( g1 , B1 , g2 , B2 )

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

template<class T1, class T2>
inline void VectAssign( T1 *g1 , cIndex_Set B1 , const T2 *g2 ,
			cIndex_Set B2 ,	const T1 gg )
{
 // g1{B} = g2{B}, B = intersection of B1 and B2, g1{B1 / B} = gg 

 Index k = *B2;
 for( Index h ; ( h = *(B1++) ) < Inf<Index>() ; ) {
  while( k < h ) {
   k = *(++B2);
   g2++;
   }

  if( k == Inf<Index>() )
   break;

  if( h == k )
   *(g1++) = *(g2++);
  else
   *(g1++) = gg;
  }
 }  // end( VectAssign( g1 , B1 , g2 , B2 , gg )

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

template<class T1, class T2, class T3>
inline void VectAssign( T1 *g1 ,  const T2 *g2 , const T3 x , Index n )
{
 // g1 := x * g2

 for( ; n-- ; )
  *(g1++) = x * (*(g2++));
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

template<class T1, class T2, class T3>
inline void VectAssign( T1 *g1 , const T2 *g2 , const T3 x , cIndex_Set B )
{
 // g1 = x * g2{B}

 for( Index h ; ( h = *(B++) ) < Inf<Index>() ; )
  *(g1++) = x * g2[ h ];
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

template<class T1, class T2, class T3>
inline void VectAssign( T1 *g1 , cIndex_Set B1 , const T3 x , const T2 *g2 ,
			cIndex_Set B2 )
{
 // g1{B} = x * g2{B}, B = intersection of B1 and B2, all other unchanged

 Index h = *B1;
 Index k = *B2;
 if( ( h < Inf<Index>() ) && ( k < Inf<Index>() ) )
  for(;;)
   if( k < h ) {
    if( ( k = *(++B2) ) == Inf<Index>() )
     break;
    g2++;
    }
   else
    if( h < k ) {
     if( ( h = *(++B1) ) == Inf<Index>() )
      break;
     g1++;
     }
    else {
     (*(g1++)) = x * (*(g2++));
     if( ( k = *(++B2) ) == Inf<Index>() )
      break;
     if( ( h = *(++B1) ) == Inf<Index>() )
      break;
     }

 }  // end( VectAssign( g1 , B1 , x , g2 , B2 )

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

template<class T1, class T2, class T3>
inline void VectAssign( T1 *g1 , cIndex_Set B1 , const T3 x , const T2 *g2 ,
			cIndex_Set B2 , const T1 gg )
{
 // g1{B} = x * g2{B}, B = intersection of B1 and B2, g1{B1 / B} = gg 

 Index k = *B2;
 for( Index h ; ( h = *(B1++) ) < Inf<Index>() ; ) {
  while( k < h ) {
   k = *(++B2);
   g2++;
   }

  if( k == Inf<Index>() )
   break;

  if( h == k )
   *(g1++) = x * (*(g2++));
  else
   *(g1++) = gg;
  }
 }  // end( VectAssign( g1 , B1 , x , g2 , B2 , gg )

/*--------------------------------------------------------------------------*/

template<class T>
inline void VectAssignB( T *const g , const T x , cIndex_Set B )
{
 // g{B} = x, all other entries of g unchanged

 for( Index h ; ( h = *(B++) ) < Inf<Index>() ; )
  g[ h ] = x;
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

template<class T>
inline void VectAssignB( T *g , const T x , cIndex_Set B , cIndex n ,
			 const T gg = 0 )
{
 // g{B} = x, all other entries (0 .. n - 1) of g1 = gg

 Index h = *(B++);
 for( Index i = 0 ; i < n ; )
  if( h == i++ )
  {
   *(g++) = x;
   h = *(B++);
   }
  else
   *(g++) = gg;
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

template<class T1, class T2>
inline void VectAssignB( T1 *g1 , const T2 *g2 , cIndex_Set B )
{
 // g1{B} = g2, all other entries of g1 unchanged

 for( Index h ; ( h = *(B++) ) < Inf<Index>() ; )
  g1[ h ] = *(g2++);
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

template<class T1, class T2>
inline void VectAssignB( T1 *g1 , const T2 *g2 , cIndex_Set B , cIndex n ,
			 const T1 gg = 0 )
{
 // g1{B} = g2, all other entries (0 .. n - 1) of g1 = gg

 Index h = *(B++);
 for( Index i = 0 ; i < n ; )
  if( h == i++ ) {
   *(g1++) = *(g2++);
   h = *(B++);
   }
  else
   *(g1++) = gg;
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

template<class T1, class T2>
inline void VectAssignB( T1 *g1 , const T2 *g2 ,
			 const T1 x , cIndex_Set B )
{
 // g1{B} = g2, all other entries of g1 = unchanged

 for( Index h ; ( h = *(B++) ) < Inf<Index>() ; )
  g1[ h ] = x * (*(g2++));
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

template<class T1, class T2>
inline void VectAssignB( T1 *g1 , const T2 *g2 , const T1 x , cIndex_Set B ,
			 cIndex n , const T1 gg = 0 )
{
 // g1{B} = x * g2, all other entries (0 .. n - 1) of g1 = gg

 Index h = *(B++);
 for( Index i = 0 ; i < n ; )
  if( h == i++ )
  {
   *(g1++) = x * (*(g2++));
   h = *(B++);
   }
  else
   *(g1++) = gg;
 }

/*--------------------------------------------------------------------------*/

template<class T1, class T2>
inline void VectAssignBB( T1 *g1 , const T2 *g2 , cIndex_Set B )
{
 // g1{B} = g2{B}, all other entries unchanged

 for( Index h ; ( h = *(B++) ) < Inf<Index>() ; )
  g1[ h ] = g2[ h ];
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

template<class T1, class T2>
inline void VectAssignBB( T1 *g1 , const T2 *g2 , const T1 x , cIndex_Set B )
{
 // g1{B} = x * g2{B}, all other entries unchanged

 for( Index h ; ( h = *(B++) ) < Inf<Index>() ; )
  g1[ h ] = x * g2[ h ];
 }

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

template<class T1, class T2>
inline void VectMAssign( T1 *g1 , const T2 *g2 , Index n )
{
 // g1 := - g2

 for( ; n-- ; )
  *(g1++) = - *(g2++);
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

template<class T1, class T2>
inline void VectMAssign( T1 *g1 , const T2 *g2 , cIndex_Set B )
{
 // g1 = - g2{B}

 for( Index h ; ( h = *(B++) ) < Inf<Index>() ; )
  *(g1++) = - g2[ h ];
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

template<class T1, class T2>
inline void VectMAssign( T1 *g1 , cIndex_Set B1 , T2 *g2 , cIndex_Set B2 )
{
 // g1{B} = - g2{B}, B = intersection of B1 and B2, all other unchanged

 Index h = *B1;
 Index k = *B2;
 if( ( h < Inf<Index>() ) && ( k < Inf<Index>() ) )
  for(;;)
   if( k < h ) {
    if( ( k = *(++B2) ) == Inf<Index>() )
     break;
    g2++;
    }
   else
    if( h < k ) {
     if( ( h = *(++B1) ) == Inf<Index>() )
      break;
     g1++;
     }
    else {
     (*(g1++)) = - (*(g2++));
     if( ( k = *(++B2) ) == Inf<Index>() )
      break;
     if( ( h = *(++B1) ) == Inf<Index>() )
      break;
     }

 }  // end( VectMAssign( g1 , B1 , g2 , B2 )

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

template<class T1, class T2>
inline void VectMAssign( T1 *g1 , cIndex_Set B1 , const T2 *g2 ,
			 cIndex_Set B2 , const T1 gg )
{
 // g1{B} = - g2{B}, B = intersection of B1 and B2, g1{B1 / B} = gg 

 Index k = *B2;
 for( Index h ; ( h = *(B1++) ) < Inf<Index>() ; ) {
  while( k < h ) {
   k = *(++B2);
   g2++;
   }

  if( k == Inf<Index>() )
   break;

  if( h == k )
   *(g1++) = - *(g2++);
  else
   *(g1++) = gg;
  }
 }  // end( VectMAssign( g1 , B1 , g2 , B2 , gg )

/*--------------------------------------------------------------------------*/

template<class T1, class T2>
inline void VectMAssignB( T1 *g1 , const T2 *g2 , cIndex_Set B )
{
 // g1{B} = - g2, all other entries of g1 = unchanged

 for( Index h ; ( h = *(B++) ) < Inf<Index>() ; )
  *(g1++) = - g2[ h ];
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

template<class T1, class T2>
inline void VectMAssignB( T1 *g1 , const T2 *g2 , cIndex_Set B , cIndex n ,
			  const T1 gg = 0 )
{
 // g1{B} = - g2, all other entries (0 .. n - 1) of g1 = gg

 Index h = *(B++);
 for( Index i = 0 ; i < n ; )
  if( h == i++ ) {
   *(g1++) = - *(g2++);
   h = *(B++);
   }
  else
   *(g1++) = gg;
 }

/*--------------------------------------------------------------------------*/

template<class T1, class T2>
inline void VectMAssignBB( T1 *g1 , const T2 *g2 , cIndex_Set B )
{
 // g1{B} = - g2{B}, all other entries unchanged

 for( Index h ; ( h = *(B++) ) < Inf<Index>() ; )
  g1[ h ] = - g2[ h ];
 }

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

template<class T>
inline void VectSum( T *const g , const T x , cIndex n )
{
 // g[ i ] += x for each i = 0 .. n - 1

 for( T *tg = g + n ; tg > g ; )
  *(--tg) += x;
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

inline void VectSum( Index_Set g , cIndex x )
{
 // special version of the above for g an Inf<Index>()-terminated vector of indices

 while( *g < Inf<Index>() )
  *(g++) += x;
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

template<class T>
inline void VectSum( T *const g , const T x , cIndex_Set B )
{
 // g{B} += x, all other entries of g unchanged

 for( Index h ; ( h = *(B++) ) < Inf<Index>() ; )
  g[ h ] += x;
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

template<class T1, class T2>
inline void VectSum( T1 *g1 , const T2 *g2 , Index n )
{
 // g1 += g2

 for( ; n-- ; )
  *(g1++) += *(g2++);
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

template<class T1, class T2>
inline void VectSum( T1 *g1 , const T2 *g2 , cIndex_Set B )
{
 // g1 += g2{B}

 for( Index h ; ( h = *(B++) ) < Inf<Index>() ; )
  *(g1++) += g2[ h ];
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

template<class T1, class T2>
inline void VectSum( T1 *g1 , cIndex_Set B1 , const T2 *g2 , cIndex_Set B2 )
{
 // g1{B} += g2{B}, B = intersection of B1 and B2, all other unchanged

 Index h = *B1;
 Index k = *B2;
 if( ( h < Inf<Index>() ) && ( k < Inf<Index>() ) )
  for(;;)
   if( k < h ) {
    if( ( k = *(++B2) ) == Inf<Index>() )
     break;
    g2++;
    }
   else
    if( h < k ) {
     if( ( h = *(++B1) ) == Inf<Index>() )
      break;
     g1++;
     }
    else {
     (*(g1++)) += (*(g2++));
     if( ( k = *(++B2) ) == Inf<Index>() )
      break;
     if( ( h = *(++B1) ) == Inf<Index>() )
      break;
     }

 }  // end( VectSum( g1 , B1 , g2 , B2 )

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

template<class T1, class T2, class T3>
inline void VectSum( T1 *g1 , const T2 *g2 , const T3 x , Index n )
{
 // g1 += x * g2

 for( ; n-- ; )
  *(g1++) += x * (*(g2++));
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

template<class T1, class T2, class T3>
inline void VectSum( T1 *g1 , const T2 *g2 , const T3 x , cIndex_Set B )
{
 // g1 += x * g2{B}

 for( Index h ; ( h = *(B++) ) < Inf<Index>() ; )
  *(g1++) += g2[ h ];
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

template<class T1, class T2, class T3>
inline void VectSum( T1 *g1 , cIndex_Set B1 , const T3 x , const T2 *g2 ,
		     cIndex_Set B2 )
{
 // g1{B} += x * g2{B}, B = intersection of B1 and B2, all other unchanged

 Index h = *B1;
 Index k = *B2;
 if( ( h < Inf<Index>() ) && ( k < Inf<Index>() ) )
  for(;;)
   if( k < h ) {
    if( ( k = *(++B2) ) == Inf<Index>() )
     break;
    g2++;
    }
   else
    if( h < k ) {
     if( ( h = *(++B1) ) == Inf<Index>() )
      break;
     g1++;
     }
    else {
     (*(g1++)) += x * (*(g2++));
     if( ( k = *(++B2) ) == Inf<Index>() )
      break;
     if( ( h = *(++B1) ) == Inf<Index>() )
      break;
     }

 }  // end( VectSum( g1 , B1 , x , g2 , B2 )

/*--------------------------------------------------------------------------*/

template<class T>
inline void VectSumB( T *g , const T x , cIndex_Set B , cIndex n , const T gg )
{
 // g{B} += x, all other entries (0 .. n - 1) of g1 += gg

 Index h = *(B++);
 for( Index i = 0 ; i < n ; )
  if( h == i++ ) {
   *(g++) += x;
   h = *(B++);
   }
  else
   *(g++) += gg;
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

template<class T1, class T2>
inline void VectSumB( T1 *g1 , const T2 *g2 , cIndex_Set B )
{
 // g1{B} += g2, all other entries of g1 unchanged

 for( Index h ; ( h = *(B++) ) < Inf<Index>() ; )
  g1[ h ] += *(g2++);
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

template<class T1, class T2>
inline void VectSumB( T1 *g1 , const T2 *g2 , cIndex_Set B , cIndex n ,
		      const T1 gg )
{
 // g1{B} += g2, all other entries (0 .. n - 1) of g1 += gg

 Index h = *(B++);
 for( Index i = 0 ; i < n ; )
  if( h == i++ ) {
   *(g1++) += *(g2++);
   h = *(B++);
   }
  else
   *(g1++) += gg;
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

template<class T1, class T2, class T3>
inline void VectSumB( T1 *g1 , const T2 *g2 , const T3 x , cIndex_Set B )
{
 // g1{B} += x * g2, all other entries of g1 unchanged

 for( Index h ; ( h = *(B++) ) < Inf<Index>() ; )
  g1[ h ] += x * (*(g2++));
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

template<class T1, class T2, class T3>
inline void VectSumB( T1 *g1 , const T2 *g2 , const T3 x , cIndex_Set B ,
		      cIndex n , const T1 gg )
{
 // g1{B} += x * g2, all other entries (0 .. n - 1) of g1 += gg

 Index h = *(B++);
 for( Index i = 0 ; i < n ; )
  if( h == i++ ) {
   *(g1++) += x * (*(g2++));
   h = *(B++);
   }
  else
   *(g1++) += gg;
 }

/*--------------------------------------------------------------------------*/

template<class T1, class T2>
inline void VectSumBB( T1 *g1 , const T2 *g2 , cIndex_Set B )
{
 // g1{B} += g2{B}, all other entries unchanged

 for( Index h ; ( h = *(B++) ) < Inf<Index>() ; )
  g1[ h ] += g2[ h ];
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

template<class T1, class T2, class T3>
inline void VectSumBB( T1 *g1 , const T2 *g2 , const T3 x , cIndex_Set B )
{
 // g1{B} += x * g2{B}, all other entries unchanged

 for( Index h ; ( h = *(B++) ) < Inf<Index>() ; )
  g1[ h ] += x * g2[ h ];
 }

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

template<class T>
inline void VectSubtract( T *const g , const T x , cIndex n )
{
 // g[ i ] -= x for each i = 0 .. n - 1; useful for *unsigned* data types
 // where VectSum( g , -x , n ) would not work

 for( T *tg = g + n ; tg > g ; )
  *(--tg) -= x;
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

inline void VectSubtract( Index_Set g , cIndex x )
{
 // special version of the above for g an Inf<Index>()-terminated vector of indices
 // (note that indices are generally unsigned)

 while( *g < Inf<Index>() )
  *(g++) -= x;
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

template<class T>
inline void VectSubtract( T *const g , const T x , cIndex_Set B )
{
 // g{B} -= x, all other entries of g unchanged

 for( Index h ; ( h = *(B++) ) < Inf<Index>() ; )
  g[ h ] -= x;
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

template<class T1, class T2>
inline void VectSubtract( T1 *g1 , const T2 *g2 , Index n )
{
 // g1 -= g2 (element-wise)

 for( ; n-- ; )
  *(g1++) -= *(g2++);
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

template<class T1, class T2>
inline void VectSubtract( T1 *g1 , const T2 *g2 , cIndex_Set B )
{
 // g1 -= g2{B} (element-wise)

 for( Index h ; ( h = *(B++) ) < Inf<Index>() ; )
  *(g1++) -= g2[ h ];
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

template<class T1, class T2>
inline void VectSubtract( T1 *g1 , cIndex_Set B1 , const T2 *g2 ,
			  cIndex_Set B2 )
{
 // g1{B} -= g2{B}, B = intersection of B1 and B2, all other unchanged

 Index h = *B1;
 Index k = *B2;
 if( ( h < Inf<Index>() ) && ( k < Inf<Index>() ) )
  for(;;)
   if( k < h ) {
    if( ( k = *(++B2) ) == Inf<Index>() )
     break;
    g2++;
    }
   else
    if( h < k ) {
     if( ( h = *(++B1) ) == Inf<Index>() )
      break;
     g1++;
     }
    else {
     (*(g1++)) -= (*(g2++));
     if( ( k = *(++B2) ) == Inf<Index>() )
      break;
     if( ( h = *(++B1) ) == Inf<Index>() )
      break;
     }

 }  // end( VectSubtract( g1 , B1 , g2 , B2 )

/*--------------------------------------------------------------------------*/

template<class T>
inline void VectSubtractB( T *g , const T x , cIndex_Set B , cIndex n ,
			   const T gg )
{
 // g{B} -= x, all other entries (0 .. n - 1) of g1 -= gg

 Index h = *(B++);
 for( Index i = 0 ; i < n ; )
  if( h == i++ ) {
   *(g++) -= x;
   h = *(B++);
   }
  else
   *(g++) -= gg;
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

template<class T1, class T2>
inline void VectSubtractB( T1 *g1 , const T2 *g2 , cIndex_Set B )
{
 // g1{B} -= g2, all other entries of g1 unchanged

 for( Index h ; ( h = *(B++) ) < Inf<Index>() ; )
  g1[ h ] -= *(g2++);
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

template<class T1, class T2>
inline void VectSubtractB( T1 *g1 , const T2 *g2 , cIndex_Set B , cIndex n ,
			   const T1 gg )
{
 // g1{B} -= g2, all other entries (0 .. n - 1) of g1 -= gg

 Index h = *(B++);
 for( Index i = 0 ; i < n ; )
  if( h == i++ ) {
   *(g1++) -= *(g2++);
   h = *(B++);
   }
  else
   *(g1++) -= gg;
 }

/*--------------------------------------------------------------------------*/

template<class T1, class T2>
inline void VectSubtractBB( T1 *g1 , const T2 *g2 , cIndex_Set B )
{
 // g1{B} -= g2{B}, all other entries unchanged

 for( Index h ; ( h = *(B++) ) < Inf<Index>() ; )
  g1[ h ] -= g2[ h ];
 }

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

template<class T>
inline void VectScale( T *g , const T x , Index n )
{
 // g *= x, x a scalar

 for( ; n-- ; )
  *(g++) *= x;
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

template<class T>
inline void VectScale( T *g , const T x , cIndex_Set B )
{
 // g{B} *= x, x a scalar

 for( Index h ; ( h = *(B++) ) < Inf<Index>() ; )
  g[ h ] *= x;
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

template<class T1, class T2>
inline void VectScale( T1 *g1 , const T2 *g2 , Index n )
{
 // g1 *= g2

 for( ; n-- ; )
  *(g1++) *= *(g2++);
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

template<class T1, class T2>
inline void VectScale( T1 *g1 , const T2 *g2 , cIndex_Set B )
{
 // g1 *= g2{B}

 for( Index h ; ( h = *(B++) ) < Inf<Index>() ; )
  *(g1++) *= g2[ h ];
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

template<class T1, class T2>
inline void VectScale( T1 *g1 , cIndex_Set B1 , const T2 *g2 , cIndex_Set B2 )
{
 // g1{B} *= g2{B}, B = intersection of B1 and B2, all other unchanged

 Index h = *B1;
 Index k = *B2;
 if( ( h < Inf<Index>() ) && ( k < Inf<Index>() ) )
  for(;;)
   if( k < h ) {
    if( ( k = *(++B2) ) == Inf<Index>() )
     break;
    g2++;
    }
   else
    if( h < k ) {
     if( ( h = *(++B1) ) == Inf<Index>() )
      break;
     g1++;
     }
    else {
     (*(g1++)) *= (*(g2++));
     if( ( k = *(++B2) ) == Inf<Index>() )
      break;
     if( ( h = *(++B1) ) == Inf<Index>() )
      break;
     }

 }  // end( VectScale( g1 , B1 , g2 , B2 )

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

template<class T1, class T2, class T3>
inline void VectScale( T1 *g1 , const T2 *g2 , const T3 x , Index n )
{
 // g1 *= x * g2

 for( ; n-- ; )
  *(g1++) *= x * (*(g2++));
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

template<class T1, class T2, class T3>
inline void VectScale( T1 *g1 , const T2 *g2 , const T3 x , cIndex_Set B )
{
 // g1 *= x * g2{B}

 for( Index h ; ( h = *(B++) ) < Inf<Index>() ; )
  *(g1++) *= g2[ h ];
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

template<class T1, class T2, class T3>
inline void VectScale( T1 *g1 , cIndex_Set B1 , const T3 x , const T2 *g2 ,
		       cIndex_Set B2 )
{
 // g1{B} *= x * g2{B}, B = intersection of B1 and B2, all other unchanged

 Index h = *B1;
 Index k = *B2;
 if( ( h < Inf<Index>() ) && ( k < Inf<Index>() ) )
  for(;;)
   if( k < h ) {
    if( ( k = *(++B2) ) == Inf<Index>() )
     break;
    g2++;
    }
   else
    if( h < k ) {
     if( ( h = *(++B1) ) == Inf<Index>() )
      break;
     g1++;
     }
    else {
     (*(g1++)) *= x * (*(g2++));
     if( ( k = *(++B2) ) == Inf<Index>() )
      break;
     if( ( h = *(++B1) ) == Inf<Index>() )
      break;
     }

 }  // end( VectScale( g1 , B1 , x , g2 , B2 )

/*--------------------------------------------------------------------------*/

template<class T1, class T2>
inline void VectScaleB( T1 *g1 , const T2 *g2 , cIndex_Set B )
{
 // g1{B} *= g2, all other entries of g1 unchanged

 for( Index h ; ( h = *(B++) ) < Inf<Index>() ; )
  g1[ h ] *= *(g2++);
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

template<class T1, class T2, class T3>
inline void VectScaleB( T1 *g1 , const T2 *g2 , const T3 x , cIndex_Set B )
{
 // g1{B} *= x * g2, all other entries of g1 unchanged

 for( Index h ; ( h = *(B++) ) < Inf<Index>() ; )
  g1[ h ] *= x * (*(g2++));
 }

/*--------------------------------------------------------------------------*/

template<class T1, class T2>
inline void VectScaleBB( T1 *g1 , const T2 *g2 , cIndex_Set B )
{
 // g1{B} *= g2{B}, all other entries unchanged

 for( Index h ; ( h = *(B++) ) < Inf<Index>() ; )
  g1[ h ] *= g2[ h ];
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

template<class T1, class T2, class T3>
inline void VectScaleBB( T1 *g1 , const T2 *g2 , const T3 x , cIndex_Set B )
{
 // g1{B} *= x * g2{B}, all other entries unchanged

 for( Index h ; ( h = *(B++) ) < Inf<Index>() ; )
  g1[ h ] *= x * g2[ h ];
 }

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

template<class T>
inline void VectIScale( T *g , const T x , Index n )
{
 // g := g / x, x a scalar; useful e.g. for *integer* types where
 // VectScale( g , 1 / x , n ) would not work

 for( ; n-- ; )
  *(g++) /= x;
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

template<class T>
inline void VectIScale( T *g , const T x , cIndex_Set B )
{
 // g{B} := g{B} / x, x a scalar

 for( Index h ; ( h = *(B++) ) < Inf<Index>() ; )
  g[ h ] /= x;
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

template<class T1, class T2>
inline void VectIScale( T1 *g1 , const T2 *g2 , Index n )
{
 // g1 /= g2

 for( ; n-- ; )
  *(g1++) /= *(g2++);
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

template<class T1, class T2>
inline void VectIScale( T1 *g1 , const T2 *g2 , cIndex_Set B )
{
 // g1 /= g2{B}

 for( Index h ; ( h = *(B++) ) < Inf<Index>() ; )
  *(g1++) /= g2[ h ];
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

template<class T1, class T2>
inline void VectIScale( T1 *g1 , cIndex_Set B1 , const T2 *g2 , cIndex_Set B2 )
{
 // g1{B} /= g2{B}, B = intersection of B1 and B2, all other unchanged

 Index h = *B1;
 Index k = *B2;
 if( ( h < Inf<Index>() ) && ( k < Inf<Index>() ) )
  for(;;)
   if( k < h ) {
    if( ( k = *(++B2) ) == Inf<Index>() )
     break;
    g2++;
    }
   else
    if( h < k ) {
     if( ( h = *(++B1) ) == Inf<Index>() )
      break;
     g1++;
     }
    else {
     (*(g1++)) /= (*(g2++));
     if( ( k = *(++B2) ) == Inf<Index>() )
      break;
     if( ( h = *(++B1) ) == Inf<Index>() )
      break;
     }

 }  // end( VectIScale( g1 , B1 , g2 , B2 )

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

template<class T1, class T2, class T3>
inline void VectIScale( T1 *g1 , const T2 *g2 , const T3 x , Index n )
{
 // g1 /= x * g2

 for( ; n-- ; )
  *(g1++) /= ( x * (*(g2++)) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

template<class T1, class T2, class T3>
inline void VectIScale( T1 *g1 , const T2 *g2 , const T3 x , cIndex_Set B )
{
 // g1 /= x * g2{B}

 for( Index h ; ( h = *(B++) ) < Inf<Index>() ; )
  *(g1++) /= ( x * g2[ h ] );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

template<class T1, class T2, class T3>
inline void VectIScale( T1 *g1 , cIndex_Set B1 , const T3 x , const T2 *g2 ,
			cIndex_Set B2 )
{
 // g1{B} /= x * g2{B}, B = intersection of B1 and B2, all other unchanged

 Index h = *B1;
 Index k = *B2;
 if( ( h < Inf<Index>() ) && ( k < Inf<Index>() ) )
  for(;;)
   if( k < h ) {
    if( ( k = *(++B2) ) == Inf<Index>() )
     break;
    g2++;
    }
   else
    if( h < k ) {
     if( ( h = *(++B1) ) == Inf<Index>() )
      break;
     g1++;
     }
    else {
     (*(g1++)) /= ( x * (*(g2++)) );
     if( ( k = *(++B2) ) == Inf<Index>() )
      break;
     if( ( h = *(++B1) ) == Inf<Index>() )
      break;
     }

 }  // end( VectIScale( g1 , B1 , x , g2 , B2 )

/*--------------------------------------------------------------------------*/

template<class T1, class T2>
inline void VectIScaleB( T1 *g1 , const T2 *g2 , cIndex_Set B )
{
 // g1{B} /= g2, all other entries of g1 unchanged

 for( Index h ; ( h = *(B++) ) < Inf<Index>() ; )
  g1[ h ] /= *(g2++);
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

template<class T1, class T2, class T3>
inline void VectIScaleB( T1 *g1 , const T2 *g2 , const T3 x , cIndex_Set B )
{
 // g1{B} /= x * g2, all other entries of g1 unchanged

 for( Index h ; ( h = *(B++) ) < Inf<Index>() ; )
  g1[ h ] /= x * (*(g2++));
 }

/*--------------------------------------------------------------------------*/

template<class T1, class T2>
inline void VectIScaleBB( T1 *g1 , const T2 *g2 , cIndex_Set B )
{
 // g1{B} /= g2{B}, all other entries unchanged

 for( Index h ; ( h = *(B++) ) < Inf<Index>() ; )
  g1[ h ] *= g2[ h ];
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

template<class T1, class T2, class T3>
inline void VectIScaleBB( T1 *g1 , const T2 *g2 , const T3 x , cIndex_Set B )
{
 // g1{B} /= x * g2{B}, all other entries unchanged

 for( Index h ; ( h = *(B++) ) < Inf<Index>() ; )
  g1[ h ] *= x * g2[ h ];
 }

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

template<class T1, class T2>
inline void VectAdd( T1 *g , const T2 *g1 , const T1 x , Index n )
{
 // g := g1 + x

 for( ; n-- ; )
  *(g++) = *(g1++) + x;
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

template<class T1, class T2, class T3>
inline void VectAdd( T1 *g , const T2 *g1 , const T3 *g2 , Index n )
{
 // g := g1 + g2

 for( ; n-- ; )
  *(g++) = *(g1++) + *(g2++);
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

template<class T1, class T2, class T3>
inline void VectAdd( T1 *g , const T2 *g1 , const T3 *g2 , cIndex_Set B ) 
{ 
 // g := g1 + g2{B}

 for( Index h ; ( h = *(B++) ) < Inf<Index>() ; )
  *(g++) = *(g1++) + g2[ h ];
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

template<class T1, class T2, class T3>
inline void VectAdd( T1 *g , const T2 *g1 , const T3 *g2 , const T3 x ,
                     Index n )
{
 // g := g1 + x * g2

 for( ; n-- ; )
  *(g++) = *(g1++) + x * (*(g2++));
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

template<class T1, class T2>
inline void VectAdd( T1 *g , const T2 *g1 , const T2 x1 , const T1 x , Index n )
{
 // g := x1 * g1 + x

 for( ; n-- ; )
  *(g++) = x1 * (*(g1++)) + x;
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

template<class T1, class T2, class T3>
inline void VectAdd( T1 *g , const T2 *g1 , const T2 x1 , const T3 *g2 ,
		     const T3 x2 , Index n )
{
 // g := x1 * g1 + x2 * g2

 for( ; n-- ; )
  *(g++) = x1 * (*(g1++)) + x2 * (*(g2++));
 }

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

template<class T1, class T2, class T3>
inline void VectDiff( T1 *g , const T2 *g1 , const T3 *g2 , Index n )
{
 // g := g1 - g2

 for( ; n-- ; )
  *(g++) = *(g1++) - *(g2++);
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

template<class T1, class T2, class T3>
inline void VectDiff( T1 *g , const T2 *g1 , const T3 *g2 , cIndex_Set B ) 
{ 
 // g := g1 - g2{B}

 for( Index h ; ( h = *(B++) ) < Inf<Index>() ; )
  *(g++) = *(g1++) - g2[ h ];
 }

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

template<class T1, class T2>
inline void VectMult( T1 *g , const T2 *g1 , const T1 x , Index n )
{
 // g := g1 * x

 for( ; n-- ; )
  *(g++) = *(g1++) * x;
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

template<class T1, class T2, class T3>
inline void VectMult( T1 *g , const T2 *g1 , const T3 *g2 , Index n )
{
 // g := g1 * g2

 for( ; n-- ; )
  *(g++) = *(g1++) * (*(g2++));
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

template<class T1, class T2, class T3>
inline void VectMult( T1 *g , const T2 *g1 , const T3 *g2 , cIndex_Set B ) 
{ 
 // g := g1 * g2{B}

 for( Index h ; ( h = *(B++) ) < Inf<Index>() ; )
  *(g++) = *(g1++) * g2[ h ];
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

template<class T1, class T2, class T3, class T4>
inline void VectMult( T1 *g , const T2 x , const T3 *g1 , const T4 *g2 ,
		      Index n )
{
 // g := x * g1 * g2

 for( ; n-- ; )
  *(g++) = x * (*(g1++)) * (*(g2++));
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

template<class T1, class T2, class T3, class T4>
inline void VectMult( T1 *g , const T2 x , const T3 *g1 , const T4 *g2 ,
		      cIndex_Set B ) 
{ 
 // g := x * g1 * g2{B}

 for( Index h ; ( h = *(B++) ) < Inf<Index>() ; )
  *(g++) = x * (*(g1++)) * g2[ h ];
 }

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

template<class T1, class T2, class T3>
inline void VectDivide( T1 *g , const T2 *g1 , const T3 *g2 , Index n )
{
 // g := g1 / g2

 for( ; n-- ; )
  *(g++) = *(g1++) / (*(g2++));
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

template<class T1, class T2, class T3>
inline void VectDivide( T1 *g , const T2 *g1 , const T3 *g2 , cIndex_Set B ) 
{ 
 // g := g1 / g2{B}

 for( Index h ; ( h = *(B++) ) < Inf<Index>() ; )
  *(g++) = *(g1++) / g2[ h ];
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

template<class T1, class T2, class T3, class T4>
inline void VectDivide( T1 *g , const T2 x , const T3 *g1 , const T4 *g2 ,
			Index n )
{
 // g := x * ( g1 / g2 )

 for( ; n-- ; )
  *(g++) = x * ( (*(g1++)) / (*(g2++)) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

template<class T1, class T2, class T3, class T4>
inline void VectDivide( T1 *g , const T2 x , const T3 *g1 , const T4 *g2 ,
			cIndex_Set B ) 
{ 
 // g := x * ( g1 / g2{B} )

 for( Index h ; ( h = *(B++) ) < Inf<Index>() ; )
  *(g++) = x * ( (*(g1++)) / g2[ h ] );
 }

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

template<class T>
inline void VectXcg( T *g1 , T *g2 , Index n )
{
 // swap of g1 and g2 (element-wise)

 for( ; n-- ; g1++ , g2++ )
  Swap( *g1 , *g2 );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

template<class T>
inline void VectXcg( T *g1 , T *g2 , cIndex_Set B )
{
 // swap of g1 and g2{B}

 for( Index h ; ( h = *(B++) ) < Inf<Index>() ; g1++ )
  Swap( *g1 , g2[ h ] );
 }

/*--------------------------------------------------------------------------*/

template<class T>
inline void VectXcgB( T *g1 , T *g2 , cIndex_Set B )
{
 // swap of g1{B} and g2

 for( Index h ; ( h = *(B++) ) < Inf<Index>() ; g2++ )
  Swap( g1[ h ] , *g2 );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

template<class T>
inline void VectXcgBB( T *g1 , T *g2 , cIndex_Set B )
{
 // swap of g1{B} and g2{B}, all other entries are unchanged

 for( Index h ; ( h = *(B++) ) < Inf<Index>() ; )
  Swap( g1[ h ] , g2[ h ] );
 }

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*-- Sparse/dense vector transformations: turn "sparse" vectors into      --*/
/*-- "dense" ones and vice-versa.                                         --*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

template<class T>
inline Index_Set Sparsify( T* g , Index_Set B , Index n , Index Bs = 0 )
{
 // turns g from a "dense" n-vector to a "sparse" one, eliminating all items
 // that are exactly == 0; writes the set of nonzero items in B, with names
 // from Bs onwards, ordered in increasing sense; returns a pointer to the
 // first element in B after the last index vritten: this can be used for
 // computing the number of nonzeroes in the "sparsified" vector and/or for
 // Inf<Index>()-terminating the set

 for( ; n ; n-- , g++ )
  if( *g )
   *(B++) = Bs++;
  else
   break;

 if( n ) {
  T* tg = g++;
  for( Bs++ ; --n ; g++ , Bs++ )
   if( *g ) {
    *(tg++) = *g;
    *(B++) = Bs;
    }
  }

 return( B );

 }  // end( Sparsify )

/*--------------------------------------------------------------------------*/

template<class T>
inline Index_Set SparsifyT( T* g , Index_Set B , Index n , const T eps ,
			    Index Bs = 0 )
{
 // as Sparsify(), but elements are considered nonzero only if they are >=
 // eps (the idea is that all elements are >= 0)

 for( ; n ; n-- , g++ )
  if( *g >= eps )
   *(B++) = Bs++;
  else
   break;

 if( n ) {
  T* tg = g++;
  for( Bs++ ; --n ; g++ , Bs++ )
   if( *g >= eps ) {
    *(tg++) = *g;
    *(B++) = Bs;
    }
  }

 return( B );

 }  // end( SparsifyT )

/*--------------------------------------------------------------------------*/

template<class T>
inline Index_Set SparsifyAT( T* g , Index_Set B , Index n , const T eps ,
			     Index Bs = 0 )
{
 // as SparsifyT(), but elements are considered nonzero only if their ABS()
 // is >= eps

 for( ; n ; n-- , g++ )
  if( ABS( *g ) >= eps )
   *(B++) = Bs++;
  else
   break;

 if( n ) {
  T* tg = g++;
  for( Bs++ ; --n ; g++ , Bs++ )
   if( ABS( *g ) >= eps ) {
    *(tg++) = *g;
    *(B++) = Bs;
    }
  }

 return( B );

 }  // end( SparsifyAT )

/*--------------------------------------------------------------------------*/

template<class T>
inline void Densify( T* g , cIndex_Set B , Index m , Index n , cIndex k = 0 )
{
 // turns g from a "sparse" m-vector, whose set of nonzero elements is B, to
 // a "dense" n-vector padded with zeroes where necessary; B has to be ordered
 // in increasing sense, but does not need to be Inf<Index>()-terminated (m gives the
 // same information); note that the function will write in g[ n - 1 ], hence
 // the vector has to have been properly allocated
 // if k > 0, the function only works for the subvector of g between k and
 // n - 1, that is, g[ 0 ] .. g[ k - 1 ] are left intact, while the
 // subvector is densified; it is required that B[] only contains indices
 // >= k (and < n, of course)

 // since the vector is densified in place, processing has to be done from
 // the end backwards to avoid overriding the data
 // note that there are two different cases:
 // - if B[ 0 ] > k, all the elements g[ k ] .. g[ B[ 0 ] - 1 ] are zero, and
 //   the beginning of B[] will be hit before hitting the g + k
 // - if instead B[ 0 ] == k, then for some h >= 0 we have that
 //   g[ k ] .. g[ k + h ] are "already in the right place" (the beginning
 //   of the "sparse" g[] is in fact dense already), and therefore one can
 //   stop copying at g + h, before hitting the beginning of B[]
 // we distinguish between the two cases in order to optimize the code

 if( m ) {        // there is at least a nonzero element
  T *g2 = g + n;
  T *g1 = g + k + m - 1;

  if( *B > k ) {  // first case - - - - - - - - - - - - - - - - - - - - - - -

   cIndex_Set B1 = B + m;
   for( Index h = *(--B1) ;; )
    if( h == --g2 - g ) {
     *g2 = *(g1--);
     if( B1 > B )
      h = *(--B1);
     else
      break;
     }
    else
     *g2 = 0;

   for( ; --g2 > g1 ; )
    *g2 = 0;
   }
  else {          // second case- - - - - - - - - - - - - - - - - - - - - - -

   B += m;
   for( Index h = *(--B) ; --g2 > g1 ; )
    if( h == g2 - g ) {
     *g2 = *(g1--);
     h = *(--B);
     }
    else
     *g2 = 0;

   }              // end second case- - - - - - - - - - - - - - - - - - - - -
  }
 else             // m == 0: just zero-out the vector
  VectAssign( g + k , T( 0 ) , n - k );

 } // end( Densify )

/*--------------------------------------------------------------------------*/

template<class T>
inline void Compact( T* g , cIndex_Set B , Index n )
{
 // takes a "dense" n-vector g and "compacts" it deleting the elements whose
 // indices are in B (all elements of B[] must be in the range 0 .. n, B[]
 // must be ordered in increasing sense and Inf<Index>()-terminated)
 // the remaining entries in g[] are shifted left of the minimum possible
 // amount in order to fill the holes left by the deleted ones

 Index i = *(B++);  // current position where to write
 Index j = i + 1;   // element to copy

 for( Index h = *(B++) ; h < Inf<Index>() ; j++ ) {
  while( j < h )
   g[ i++ ] = g[ j++ ];

  h = *(B++);
  }

 VectAssign( g + i , g + j , n - j );

 }  // end( Compact )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*--                                                                      --*/
/*-- Functions for array manipulation (merging, shifting ...). The types  --*/
/*-- involved in these operations must support ordering relations '=='    --*/
/*-- and '>', as well as the assigment '='.                               --*/
/*--                                                                      --*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

template<class T>
inline void Merge( T *g , const T *g2 , const T *g3 , const T Stp )
{
 // merges the two ordered and "Stp-terminated" vectors g2 and g3 (that is,
 // there must be an element in both g2 and g3 that is >= Stp, and that is
 // is considered to be the terminator) into g1, that will therefore be
 // ordered and Stp-terminated

 T i = *g2;
 T j = *g3;

 for(;;)
  if( i < j ) {
   if( ( *(g++) = i ) >= Stp )
    break;

   i = *(++g2);
   }
  else {
   if( ( *(g++) = j ) >= Stp )
    break;

   if( i == j )
    i = *(++g2);

   j = *(++g3);
   }
 }  // end( Merge )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

template<class T>
inline void ShiftVect( T *g , Index n )
{
 // g[ 0 ] = g[ 1 ], g[ 1 ] = g[ 2 ] ... g[ n - 1 ] = g[ n ]

 for( T *t = g ; n-- ; t = g )
  *t = *(++g);
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

template<class T>
inline void ShiftVect( T *g , Index n , cIndex k )
{
 // g[ i ] = g[ i + k ], i = 0 ... n - 1

 for( T *t = g + k ; n-- ; )
  *(g++) = *(t++);
 }

/*--------------------------------------------------------------------------*/

template<class T>
inline void RotateVect( T *g , Index n )
{
 // as ShiftVect(), but g[ 0 ] is moved to g[ n ]

 const T tmp = *g;

 for( T *t = g ; n-- ; t = g )
  *t = *(++g);

 *g = tmp;
 }

/*--------------------------------------------------------------------------*/

template<class T>
inline void ShiftRVect( T *g , Index n )
{
 // g[ n ] = g[ n - 1 ], g[ n - 1 ] = g[ n - 2 ], ..., g[ 1 ] = g[ 0 ]

 for( T *t = ( g += n ) ; n-- ; t = g )
  *t = *(--g);
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

template<class T>
inline void ShiftRVect( T *g , Index n , cIndex k )
{
 // g[ i ] = g[ i - k ], i = n - 1 ... 0

 g += n;
 for( T *t = g + k ; n-- ; )
  *(--t) = *(--g);
 }

/*--------------------------------------------------------------------------*/

template<class T>
inline void RotateRVect( T *g , Index n )
{
 // as ShiftRVect(), but g[ n ] is moved to g[ 0 ]

 T *t = ( g += n );
 const T tmp = *t;

 for( ; n-- ; t = g )
  *t = *(--g);

 *g = tmp;
 }

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*--                                                                      --*/
/*-- Searching/ordering functions. The types involved in these operations --*/
/*-- must support ordering relations '==' and '>', as well as the         --*/
/*-- assigment '='.                                                       --*/
/*--                                                                      --*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

template<class T>
inline Index Match( const T *g , const T x , cIndex n )
{
 // if v contains elements identical to x, then the smallest index among all
 // such elements is reported (0 .. n - 1), else n is reported

 Index i = 0;
 for( ; ( i < n ) && ( *(g++) != x ) ; )
  i++;

 return( i );
 }

/*--------------------------------------------------------------------------*/

template<class T>
inline bool EqualVect( const T *g1 ,  const T *g2 , Index n )
{
 // returns true <=> the two n-vectors g1 and g2 are element-wise identical

 for( ; n-- ; )
  if( *(g1++) != *(g2++) )
   return( false );

 return( true );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

template<class T>
inline bool EqualVect( const T *g1 ,  const T *g2 , cIndex n , cIndex_Set B )
{
 // returns true <=> the two n-vectors g1 and g2 are element-wise identical,
 // where g2 is given in sparse form, i.e., g2[ i ] is the B[ i ]-th element

 Index i = 0;
 for( Index h ; ( h = *(B++) ) < Inf<Index>() ; i++ ) {
  for( ; i < h ; i++ )
   if( *(g1++) )
    return( false );

  if( *(g1++) != *(g2++) )
   return( false );
  }

 for( ; i < n ; i++ )
  if( *(g1++) )
   return( false );

 return( true );
 }

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

template<class T>
inline Index BinSearch( const T *Set , Index Stop , const T What )
{
 // perform a binary search on the ordered set of Stop elements (of type T and
 // without replications) contained in the array Set: Set must be "infinity-
 // terminated", i.e. Set[ Stop ] > What. Searches for the element What, that
 // may or may not be in Set: if What is found, its position is reported,
 // otherwise the position of the smallest element > What in Set is reported.
 // If there are no elements > What in Set, then Stop is reported

 Index i = Stop / 2;
 for( Index Strt = 0 ; Strt != Stop ; ) {
  if( Set[ i ] == What )
   break;
  else
   if( Set[ i ] > What )
    Stop = i;
   else
    Strt = i + 1;

  i = ( Strt + Stop ) / 2;
  }

 return( i );

 }  // end( BinSearch )

/*--------------------------------------------------------------------------*/

template<class T>
inline Index BinSearch1( const T *Set , Index Stop , const T What )
{
 // like BinSearch() above, but What *must be* in Set (=> Set is nonempty)

 Index Strt = 0;
 Index i = Stop / 2;
 for( Index h ; ( h = Set[ i ] ) != What ; ) {
  if( h > What )
   Stop = i - 1;
  else
   Strt = i + 1;

  i = ( Strt + Stop ) / 2;
  }

 return( i );

 }  // end( BinSearch1 )

/*--------------------------------------------------------------------------*/

template<class T>
inline Index BinSearch2( const T *Set , Index Stop , const T What )
{
 // like BinSearch() above, but What must *not* be in Set (=> the position
 // of the smallest element > What in Set will be returned)

 Index i = Stop / 2;
 for( Index Strt = 0 ; Strt != Stop ; ) {
  if( Set[ i ] > What )
   Stop = i;
  else
   Strt = i + 1;

  i = ( Strt + Stop ) / 2;
  }

 return( i );

 }  // end( BinSearch2 )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

template<class T>
inline void HeapIns( T *H , const T x , Index n )
{
 // H is a binary heap of elements of type T, ordered in increasing sense
 // (i.e. the root of the heap is the smallest element) and containing n
 // elements: the new element x is inserted in H

 for( H-- , n++ ;;) {
  Index p = n / 2;
  if( ( ! p ) || ( H[ p ] <= x ) )
   break;

  H[ n ] = H[ p ];
  n = p;
  }

 H[ n ] = x;

 }  // end( HeapIns )

/*--------------------------------------------------------------------------*/

template<class T>
inline Index HeapDel( T *H , cIndex n )
{
 // H is a binary heap of elements of type T, as above: returns the smallest
 // element (the root) deleting it from H. n is now intended to be the
 // position of the last element in H, i.e. | H | - 1, i.e. | H | *after*
 // the deletion

 const T h = *H;
 *H = H[ n ];

 for( Index i = 0 ; i < n  ; ) {
  Index j = 2 * i;

  if( ++j >= n )
   break;

  Index k = j;

  if( ++k < n )
   if( H[ k ] < H[ j ] )
    k = j;

  if( H[ i ] <= H[ j ] )
   break;

  Swap( H[ i ] , H[ j ] );
  i = j;
  }

 return( h );
 }

/*--------------------------------------------------------------------------*/

#if( OPT_USE_NAMESPACES )
 };  // end( namespace OPTtypes_di_unipi_it )
#endif

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif  /* OPTvect.h included */

/*--------------------------------------------------------------------------*/
/*---------------------- End File OPTvect.h --------------------------------*/
/*--------------------------------------------------------------------------*/

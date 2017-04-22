/*--------------------------------------------------------------------------*/
/*----------------------------- File SubGrad.C -----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 *
 * Implementation of the <tt>SubGrad</tt> class, which implements the
 * NDOSolver interface for NonDifferentiable Optimization Solvers, as
 * described in NDOSlver.h, using a (deflected, projected, incremental)
 * subgradient-type algorithm.
 *
 * \version 1.00
 *
 * \date 22 - 10 - 2015
 *
 * \author Antonio Frangioni \n
 *         Operations Research Group \n
 *         Department of Informatics \n
 *         University of Pisa \n\n
 *
 * \author Enrico Gorgone \n
 *         Graphs and Mathematical Optimization (GOM) Group \n
 *         Department of Informatics \n
 *         Free University of Brussels \n\n
 *
 * Copyright &copy 2001 - 2015 by Antonio Frangioni, Enrico Gorgone.
 */

/*--------------------------------------------------------------------------*/
/*--------------------------- IMPLEMENTATION -------------------------------*/
/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "SubGrad.h"

#include "math.h"
#include "OPTvect.h"

/*--------------------------------------------------------------------------*/
/*-------------------------------- MACROS ----------------------------------*/
/*--------------------------------------------------------------------------*/

#if( LOG_SG )
 #define VLOG( l , x ) if( NDOLLvl > l ) *NDOLog << x
 #define VLOG2( l , c , x ) if( ( NDOLLvl > l ) && c ) *NDOLog << x

#else
 #define VLOG( l , x )
 #define VLOG2( l , c , x )
#endif

/*--------------------------------------------------------------------------*/
/*-------------------------------- USING -----------------------------------*/
/*--------------------------------------------------------------------------*/

#if( OPT_USE_NAMESPACES )
 using namespace NDO_di_unipi_it;
#endif

/*--------------------------------------------------------------------------*/
/*----------------------------- CONSTANTS ----------------------------------*/
/*--------------------------------------------------------------------------*/

static cHpNum Multplier = 1.0;

static const unsigned char PrjSG =  1;   // project of Gi[wFi] = \f$g_i\f$
static const unsigned char PrjDIRM1 = 2; // project of DirM1 = \f$d_{i-1}\f$
static const unsigned char PrjDIR = 4;   // project of DirM1 = \f$d_i\f$

static const unsigned char SafeRULE = 1;  // safe rule
static const unsigned char StepRST = 2;   // stepsize-restricted

static const unsigned char RstAlg =  1;  // don't reset algorithmic parameters
static const unsigned char RstCrr =  2;  // don't reset current point
static const unsigned char RstSbg =  4;  // don't reset subgradients
static const unsigned char RstCnt =  8;  // don't reset constraints
static const unsigned char RstFiV = 16;  // don't reset FiVals

static const bool first = true;          // indicate the position of the
static const bool second = false;        // direction

/*--------------------------------------------------------------------------*/
/*----------------------------- FUNCTIONS ----------------------------------*/
/*--------------------------------------------------------------------------*/

int myrandom( int i ) {
 return rand() % i;
 }

/*--------------------------------------------------------------------------*/
/*---------------------- IMPLEMENTATION OF SubGrad -------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- PUBLIC METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/
/*---------------------------- CONSTRUCTOR ---------------------------------*/
/*--------------------------------------------------------------------------*/

SubGrad::SubGrad( istream *iStrm )
         :
         NDOSolver( iStrm )
{
 // initialize algorithmic parameters - - - - - - - - - - - - - - - - - - - -

 DfltdSfInpt( iStrm , SGPar1 , Index( 0 ) );
 DfltdSfInpt( iStrm , SGPar2 , HpNum( 0 ) );
 DfltdSfInpt( iStrm , SGPar3 , Index( 1 ) );
 DfltdSfInpt( iStrm , SGPar4 , bool( true ) );
 DfltdSfInpt( iStrm , SGPar5 , Index( 0 ) );
 srand( SGPar5 );

 if( SGPar2 < 0 )
  SGPar2 = 0;

 // some initializations - - - - - - - - - - - - - - - - - - - - - - - - - - -

 KpBstL = false;
 LHasChgd = true;  // LHasChgd will be again set to true by FormLambda1
    
 Q2KNP = 0;
 }  // end( SubGrad::SubGrad() )

/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

void SubGrad::SetStepsize( Stepsize *STP )
{
 stepsize = STP;
 if( stepsize && Oracle )
  InitStepsize();

 }  // end( SubGrad::SetStepsize() )

/*--------------------------------------------------------------------------*/

void SubGrad::SetDeflection( Deflection *Vol )
{
 deflection = Vol;
 if( Oracle )
  InitDeflection();

 }  // end( SubGrad::SetDeflection() )

/*--------------------------------------------------------------------------*/

void SubGrad::SetQKNP( CQKnPClass *KNP )
{
 Q2KNP = KNP;

 }  // end( SubGrad::SetQKNP() )

/*--------------------------------------------------------------------------*/

void SubGrad::SetFiOracle( FiOracle *Fi )
{
 if( Oracle )  // changing from a previous oracle- - - - - - - - - - - - - - -
  MemDealloc();   // deallocate memory

 if( Fi ) {    // setting a new oracle - - - - - - - - - - - - - - - - - - - -

  // "throw" the method of the base class- - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  NDOSolver::SetFiOracle( Fi );

  // read information about Fi - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  MaxNumVar = Oracle->GetMaxNumVar();
  MaxNConst = Oracle->MaxNConst();

  if( NrFi > 1 )
   for( Index k = 0 ; k++ < NrFi ; )
    if( Oracle->GetBNC( k ) )
     throw( NDOException(
		   "SubGrad::SetFiOracle: easy component not allowed yet" ) );

  // allocate memory - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  LambdaBar = new LMNum[ MaxNumVar ];                // stability center
  VectAssign( LambdaBar , LMNum( 0 ) , MaxNumVar );  // default = 0

  Lambda = new LMNum[ MaxNumVar ];                   // tentative point
  VectAssign( Lambda , LMNum( 0 ) , MaxNumVar );     // default = 0

  if( KpBstL )
   LambdaBest = new LMNum[ MaxNumVar ];              // best point so far
  else
   LambdaBest = 0;

  LowerBound = -Inf<HpNum>();               // lower bound on the function
  TrueLB = false;

  FiBest = FiLambda = Inf<HpNum>();         // Fi( Lambda ) is not known

  // box constraints  - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  ub = new LMNum[ MaxNumVar ];           // upper bounds on the variables
  lb = new LMNum[ MaxNumVar ];           // lower bounds on the variables

  BoxConst = false;
  for( Index i = 0 ; i < NumVar ; i++ ) {
   if( Oracle->GetUC( i ) )
    lb[ i ] = - Inf<LMNum>();
   else {
    lb[ i ] = 0;
    BoxConst = true;
    }

   ub[ i ] = Oracle->GetUB( i );
   if( ub[ i ] < Inf<LMNum>() )
    BoxConst = true;
   }

  // allocate memory for subgradient Gi[wFi] and direction dir- - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  NrmGi = 0;                     // norm of Gi[wFi]
  Gi = new SgNum[ MaxNumVar ];   // Gi[wFi]
  dir = new SgNum[ MaxNumVar ];  // the direction will be allocated after
                                 // a call to SetDeflection()

  // definition of disjunctive simplex constraints, namely each variable is
  // associated to only one knapsack constraint - - - - - - - - - - - - - - -

  if( MaxNConst ) {
   CnstBeg = new Index[ MaxNConst ];
   CnstVol = new HpNum[ MaxNConst ];
   CnstNxt = new SIndex[ MaxNumVar ];
   VectAssign( CnstBeg , Index( Inf<Index>() ) , MaxNConst );
   VectAssign( CnstVol , HpNum( Inf<HpNum>() ) , MaxNConst );
   VectAssign( CnstNxt , SIndex(-1) , MaxNumVar );
   MaxName = MaxNConst + 3;
   }
  else
   MaxName = 3;

  // 0 = the subgradient Gi[wFi], 1 = the direction Dir,
  // 2 = the direction DirM1

  ReSetAlg( RstCrr | RstSbg | RstCnt );

  // tell the oracle about the NDOSolver and its settings - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  Oracle->SetNDOSolver( this );
  Oracle->SetMaxName( MaxName );
  Oracle->SetPrecision( EInit );

  // asking the projection of both DirM1 and and Gi implies the automatic
  // projection of dir- - - - - - - - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( ( SGPar1 & PrjDIRM1 ) && ( SGPar1 & PrjSG ) )
   SGPar1 &= ~PrjDIR;

  // warning: the following things can only be done *after* that
  // Oracle->SetMaxName() has been invoked, because they use methods of the
  // oracle which depends on knowledge of the MaxName to work properly
  // insert the constant subgradient of the 0-th component - - - - - - - - - -

  MultBse = new Index[2];
  MultBse[ 1 ] = Inf<Index>();

  // checking for 0-th component - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  cIndex_Set SG0Bse;
  cIndex SG0BDm = Oracle->GetGi( Gi , SG0Bse , MaxName );

  ZeroComp = SG0BDm > 0;
  SubGrad::SetPar( SubGrad::kSGPar2 , SGPar2 );
  // initialization of data structures for incremental done there inside

  // initialize HatLambda  - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( SGPar4 ) {
   LambdaHat = new LMNum[ MaxNumVar ];
   VectAssign( LambdaHat , LMNum( 0 ) , MaxNumVar );
   }
  else
   LambdaHat = 0;

  FiHat = 0;     // Fi( HatLambda ) is not known

  // initialize both stepsize and deflection solvers - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( stepsize )
   InitStepsize( );
  InitDeflection( );

  // other initializations- - - - - - - - - - - - - - - - - - - - - - - - - -

  fs = FiOracle::kFiNorm;
  EmptySet = false;       // the feasible set is not empty

  } // end if( Fi )
 else {
  if( deflection )
   deflection->Format();

  if( stepsize )
   stepsize->Format();
  }
 }  // end( SubGrad::SetFiOracle() )

/*--------------------------------------------------------------------------*/

void SubGrad::SetLambda( cLMRow tLambda )
{
 // some exceptions- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( ! Oracle )
  throw( NDOException( "SubGrad::SetLambda: Oracle == 0" ) );

 // assign tLambda to Lambda - - - - - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( tLambda )
  VectAssign( LambdaBar , tLambda , NumVar );
 else
  VectAssign( LambdaBar , LMNum(0) , NumVar );

 FiBest = FiBar = Inf<HpNum>();  // so far Fi( Lambda ) is not known

 } // end ( SubGrad::SetLambda() )

/*--------------------------------------------------------------------------*/

void SubGrad::KeepBestLambda( const bool KBL )
{
 if( KpBstL != KBL ) {
  if( KpBstL ) {
   delete[] LambdaBest;
   LambdaBest = 0;
   }
  else
   LambdaBest = new LMNum[ MaxNumVar ];

  KpBstL = KBL;
  }
 }  // end( SubGrad::KeepBestLambda() )

/*--------------------------------------------------------------------------*/

void SubGrad::SetPar( const int wp , const int value )
{
 switch( wp ) {
  case( kSGPar1 ):  // - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   SGPar1 = value;
   if( ( SGPar1 & PrjDIRM1 ) && ( SGPar1 & PrjSG ) )
    SGPar1 &= ~PrjDIR;
   if( ! deflection )
    SGPar1 &= ~( PrjDIRM1 | PrjDIR );
   break;

  case( kSGPar3 ):  // - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   SGPar3 = value;
   if( ! deflection )
    SGPar3 |= StepRST;
   break;

  case( kSGPar5 ):  // - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   SGPar5 = value;
   break;

  default:          //- - - - - - - - - - - - - - - - - - - - - - - - - - - -
   NDOSolver::SetPar( wp , value );
  }
 }  // end( SubGrad::SetPar( Index ) )

/*--------------------------------------------------------------------------*/

void SubGrad::SetPar( const int wp , cHpNum value )
{
 switch( wp ) {
  case( kEInit ):  // - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   EInit = value;
   if( Oracle )
    Oracle->SetPrecision( EInit );
   break;

  case( kSGPar2 ):  //- - - - - - - - - - - - - - - - - - - - - - - - - - - -

   SGPar2 = value;

   InnIter = false;
   if( ZeroComp ) {  // 0-th component is present, so it has to be considered
    NItIncr = ceil( double( NrFi + 1 ) * SGPar2 );
    for( Index i = 0 ; i <= NrFi ; i++ )
     Seq.push_back( i );
    }
   else {
    NItIncr = ceil( double( NrFi ) * SGPar2 );
    for( Index i = 0 ; i++ < NrFi ; )
     Seq.push_back( i );
    }

   if( NItIncr ) {
    SGPar4 = false;
    if( deflection )
     throw( NDOException( "SubGrad::SetPar: deflection is not allowed" ) );
    }
   break;

  default:          //- - - - - - - - - - - - - - - - - - - - - - - - - - - -
   NDOSolver::SetPar( wp , value );
  }
 }  // end( SubGrad::SetPar( HpNum ) )

/*--------------------------------------------------------------------------*/

void SubGrad::SetPar( const int wp , const bool value ) {
 switch( wp ) {

  case( kSGPar4 ):   //- - - - - - - - - - - - - - - - - - - - - - - - - - - -

   if( NItIncr )     // forbidden for incremental version
    SGPar4 = false;
   else
    SGPar4 = value;
   break;

  default:           //- - - - - - - - - - - - - - - - - - - - - - - - - - - -
   throw( NDOException( "SubGrad::SetPar: unknown parameter" ) );
  }
 }  // end( SubGrad::SetPar( bool ) )

/*--------------------------------------------------------------------------*/

void SubGrad::SetNDOLog( ostream *outs , const char lvl )
{
 NDOSolver::SetNDOLog( outs , lvl );

 #if( LOG_SG )
  if( NDOLLvl > 1 ) {
   *NDOLog << endl << "SubGrad Solver: Vars = " << NumVar << " (" << MaxNumVar
	   << ") ~  NItIncr = " << NItIncr <<  " ~ EpsLin = " << EpsLin
	   << " ~ FiEps = " << EInit << " ~ t* = " << tStar;

   *NDOLog << endl << " SGParameters = ( " << SGPar1 << " ~ " << SGPar2
		   << " ~ " << SGPar3 << " ~ " << ( SGPar4? "true" : "false")
		   <<  " ~ "  << SGPar5 << " )" << endl;
   }
 #endif

 }  // end( SubGrad::SetNDOLog() )

/*--------------------------------------------------------------------------*/
/*-------------------- METHODS FOR SOLVING THE PROBLEM ---------------------*/
/*--------------------------------------------------------------------------*/

NDOSolver::NDOStatus SubGrad::Solve( void )
{
 // initializations- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 Result = kOK;
 SCalls++;
 if( NDOt )
  NDOt->Start();

 Index comp; // a component

 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // main cycle starts here- - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 do {

  // function log: print the information relative to the stability center
  // and the best point - - - - - - - - - - - - - - - - - - - - - - - - - - -

  Log1();

  // check for running time - - - - - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( MaxTime && NDOt )
   if( NDOt->Read() > MaxTime ) {
    Result = kStpTime;
   break;
   }

  // To compute the next point Lambda1, select the direction Dir and the
  // stepsize Nu using, respectively, Deflection and Stepsize objects
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( FiBar < Inf<HpNum>() ) {
   // if the function value is unknown  skip the direction and stepsize
   // computation and go straight to FiAndGi()

   // direction log: print the information relative to DIRM1  - - - - - - - -
   // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   Log2();

   // FormD call: a new direction is provided, but Dir is not updated
   // yet because AddVariables() and RemoveVariables() need to work
   // with the prior vector Dir - - - - - - - - - - - - - - - - - - - - - - -
   // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   FormD( );

   // check the status of the FiOracle and take the action, however
   // because we do not succeed in handling the variables generation
   // combined with incremental approach, GetFiStatus is invoked only during
   // outer iterations  - - - - - - - - - - - - - - - - - - - - - - - - - - -
   // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   if( !InnIter )
    fs = Oracle->GetFiStatus();

   if( fs == FiOracle::kFiStop ) {
    VLOG( 1 , " ~ FiOracle:STOP" << endl );
    Result = kStopped;
    break;
    }

   if( fs == FiOracle::kFiError ) {
    VLOG( 1 , " ~ Error in the FiOracle" << endl );
    Result = kError;
    break;
    }

   // if something is changed go back to FormD(): the deflection coefficient
   // Alfa, and in turn Dir, must be again computed, and, in the
   // deflection-restricted approach, even the stepsize Nu- - - - - - - - - -
   // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   if( fs == FiOracle::kFiChgd ) {
    VLOG( 1 , " ~ Fi changed: loop" << endl );
    continue;
    }

   // save the search direction Dir: note that the previous DIR is not saved
   // as DirM1, but it is merely discharged; only and one vector is kept
   // for storing the direction- - - - - - - - - - - - - - - - - - - - - - - -
   //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   SaveDir();

   // if the Stepsize object needs the scalar product between DirM1 and Gi,
   // then the direction vector is kept- - - - - - - - - - - - - - - - - - - -

   if( stepsize->NeedsdkM1Gk() && ( ! deflection ) )
    VectAssign( dir , Gi , NumVar );

   // direction log: print the information relative to DIR   - - - - - - - - -
   //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   Log2();

   // a stop condition regarding consecutive "short" stepsizes - - - - - - - -
   //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   if( step <= 1e-8 * tStar ) {
    CSmallStep++;
    if( CSmallStep > ( ( NItIncr > 0 )? 100 * NrFi : 100 ) ) {
     Result = kStopped;
     break;
     }
    }
   else
    CSmallStep = 0;

   }  // end( search direction computation ) - - - - - - - - - - - - - - - - -
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  // check for optimality  - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   if( ( ! InnIter ) && IsOptimal() )
   break;

  // calculate the tentative point Lambda1 (not projected yet)   - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  FormLambda1();

  VLOG( 1 , endl << "           " << " |Lambda1| = "
	         << Norm( Lambda , NumVar )  );

  // a real iteration (iterations where Fi() is not evaluated count as well) -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  ParIter++;

  // is an inner iteration?- - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  InnIter = NItIncr? GiEvaltns % ( NItIncr + 1 ): 0;
  LHasProj = false;

  // pick up a component - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( InnIter ) {  // if it is the first component, re-shuffle the vector
   comp = ( ( GiEvaltns % ( NItIncr + 1 ) ) - 1 ) %
            ( ZeroComp? NrFi + 1 : NrFi );
   if( comp == 0 )
    random_shuffle( Seq.begin() , Seq.end() , myrandom );
   comp = Seq.at( comp );
   }
  else
   comp = Inf<Index>();

  // project Lambda1 onto the feasible set - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  for( ;; ) {
   if( ProjFsb( Lambda , EmptySet ) ) {  // go to the feasible point - - - - -
    LHasProj = LHasChgd = true;
    VLOG( 1 , " ~ |Lambda1| = " << Norm( Lambda , NumVar )  );
    if( EmptySet )
     break;
    }

   // never evaluate the objective function at the same point - - - - - - - -
   // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   if( ! LHasChgd )
    throw( NDOException( "SubGrad::Solve: this should not happen" ) );

   // compute Fi( Lambda1 ) and a subgradient / constraint Gi[wFi]- - - - - -
   // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   FiAndGi( comp );

   // if the Lambda1 is infeasible, go back to the projection - - - - - - - -
   // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   if( FiLambda < Inf<HpNum>() )
    break;

   } // end function evaluation - - - - - - - - - - - - - - - - - - - - - - -
     // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  // project the subgradient Gi[wFi] over the tangent cone at tentative point
  // Lambda1- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( SGPar1 & PrjSG )
   ProjTng( Gi , Lambda );
  NrmGi = Norm( Gi , NumVar );

  // compute the scalar product between Dir and Gi[wFi] at Lambda1- - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( stepsize->NeedsdkM1Gk() ) {
   dM1Gk = ScalarProduct( dir , Gi , NumVar );
   dM1GkDone = true;
   }
  else
   dM1GkDone = false;

  // print the logs regarding Fi and Gi[wFi] at Lambda1 - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  Log3();
  Log4();

  // check for feasibility- - - - - - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( EmptySet ) {
   Result = kUnfsbl;
   break;
   }

  if( FiLambda == -Inf<HpNum>() ) {
   Result = kUnbndd;
   break;
   }

  // check whether or not the Lower Bounds have changed - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  UpdtLowerBound();
  VLOG2( 1 , LowerBound > - Inf<HpNum>() , " ~ UB = " << -LowerBound );

  // check for Lower Bound stopping condition - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( TrueLB )
   if( FiBest - EpsLin * ABS( FiBest ) <= LowerBound )
    break;

  if( FiBest <= LowerBound  * ( 1.0 - ( LowerBound > 0.0 ? EpsLin : - EpsLin )
				) ) {
   Result = kUnbndd;
   break;
   }

  // if FiBar (the function at LambdaBar) is unknown, or no deflection is
  // used, then move to Lambda1 - - - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( ( FiBar == Inf<HpNum>() ) || ( ! deflection ) ) {
   GotoLambda1();

   if( deflection ) {
    HpNum Mlt= 0;
    Index NmSt = MaxNConst;
    Oracle->Aggregate( &Mlt , &NmSt , 1 , MaxNConst + 2 );
    }

   DoSS = true;

   VLOG( 1 , endl );
   continue;
   }

  // Null Step (NS) vs Serious Step (SS)- - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( ( DoSS = deflection->DoSS() ) ) {  // SS- - - - - - - - - - - - - - - -

   VLOG2( 1 , deflection->Delta() < Inf<HpNum>() , endl << " SS[" << CSSCntr
	  );

   GotoLambda1();

   CSSCntr++;
   CNSCntr = 0;
   }
  else {  // NS- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   VLOG2( 1 , deflection->Delta() > - Inf<HpNum>() , endl << " NS["
	      << CNSCntr );

   UpdateSigma();

   CNSCntr++;
   CSSCntr = 0;

   }  // end SS/NS - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  VLOG( 1 , "]: DFi (" << FiBar - FiLambda << ") < m1 * Dv ("
	     << deflection->Delta() << ")" );
  VLOG( 1 , endl );

  } while( ( ! MaxIter ) || ( ParIter < MaxIter ) );

 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // main cycle ends here- - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( MaxIter && ( ParIter >= MaxIter ) && ( ! Result ) )
  Result = kStpIter;

 if( NDOt )
  NDOt->Stop();

 return( Result );

 }  // end( SubGrad::Solve() )

/*--------------------------------------------------------------------------*/

void SubGrad::ReSetAlg( unsigned char RstLvl )
{
 if( ! ( RstLvl & RstAlg ) ) {  // reset algorithmic parameters - - - - - - -
  ParIter = 0;             // reset iterations counerst, comprised
  CSmallStep = 0;          // consecutive NS/SS count
  CSSCntr = CNSCntr = 0;
  SetPar( NDOSolver::kEInit , HpNum( 0 ) );  // reset Fi precision
  }

 if( ! ( RstLvl & RstCrr ) )  // reset the current point to all-0- - - - - - -
  SetLambda();

 if( ! ( RstLvl & ( RstSbg | RstCnt ) ) ) {  // reset everything - - - - - - -
  if( Oracle )
   Oracle->Deleted();   // tell the oracle (if any) about it

  VectAssign( CnstBeg , Index( Inf<Index>() ) , MaxNConst );
  VectAssign( CnstVol , HpNum( Inf<HpNum>() ) , MaxNConst );
  VectAssign( CnstNxt , SIndex(-1) , MaxNumVar );
  }
 else
  if( ! ( RstLvl & RstSbg ) ) {  // reset the subgradient and the directions
	                         // (but not the knapsack constraints)- - - -
   for( Index i = MaxName ; i-- > MaxNConst ; )
    Oracle->Deleted( i );
   }
  else
   if( ! ( RstLvl & RstCnt ) ) {  // reset the knapsack constraints - - - - -
    for( Index i = MaxNConst ; i-- ; )
     if( CnstBeg[ i ] < Inf<Index>() )
      Oracle->Deleted( i );
    VectAssign( CnstBeg , Index( Inf<Index>() ) , MaxNConst );
    VectAssign( CnstVol , HpNum( Inf<HpNum>() ) , MaxNConst );
    VectAssign( CnstNxt , SIndex(-1) , MaxNumVar );
    }

 if( ! ( RstLvl & RstFiV ) )  // reset Fi( LambdaBar )- - - - - - - - - - - -
  FiBar = Inf<HpNum>();

 }  // end( SubGrad::ReSetAlg() )

/*--------------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/

cLMRow SubGrad::ReadBestSol( cIndex_Set &I , Index &D )
{
 I = 0;
 D = NumVar;
 return( KpBstL ? LambdaBest : LambdaBar );

 }  // end( SubGrad::ReadBestSol() )

/*--------------------------------------------------------------------------*/

HpNum SubGrad::ReadFiVal( cIndex wFi )
{
 if( wFi!= Inf<Index>() )
  throw( NDOException( "SubGrad::ReadFiVal: wFi != INF" ) );

 return( FiBar );

 }  // end( SubGrad::ReadFiVal() )

/*--------------------------------------------------------------------------*/

cLMRow SubGrad::ReadSol( cIndex_Set &I , Index &D )
{
 I = 0;
 D = NumVar;
 return( LambdaBar );

 }  // end( SubGrad::ReadSol() )

/*--------------------------------------------------------------------------*/

HpNum SubGrad::ReadBestFiVal( cIndex wFi )
{
 if( wFi!= Inf<Index>() )
  throw( NDOException( "SubGrad::ReadBestFiVal: wFi != INF" ) );

 return( FiBest );

 }  // end( SubGrad::ReadBestFiVal() )

/*--------------------------------------------------------------------------*/

HpNum SubGrad::ReadHatFiVal( void )
{
 if( SGPar4 )
  return( FiHat );
 else
  return( Inf<HpNum>() );

 }  // end( SubGrad::ReadHatFiVal() )

/*--------------------------------------------------------------------------*/

bool SubGrad::IsOptimal( HpNum eps )
{
 HpNum FiL = FiBest;

 if( FiL == Inf<HpNum>() )
  return( false );
 else {
  if( FiL < 0 ) FiL = - FiL;
  if( FiL < 1 ) FiL = 1;

  if( eps <= 0 )
   eps = EpsLin;

  if( SGPar4 )
   return( tStar * sqrt( NrmDir ) + min( Epsilon , HatEpsilon )
	   <= eps * FiL );
  else
   return( tStar * sqrt( NrmDir ) + Epsilon <= eps * FiL );
  }
 }  // end( SubGrad::IsOptimal() )

/*--------------------------------------------------------------------------*/

cHpRow SubGrad::ReadMult( cIndex_Set &I , Index &D , cIndex wFi )
{
 if( wFi!= Inf<Index>() )
  throw( NDOException( "SubGrad::ReadBestFiVal: wFi != INF" ) );

 if( FiBar < Inf<HpNum>() ) {
  if( deflection )
   *MultBse = MaxNConst + ( DirPos == first ? 2 : 1 );
  else
   *MultBse = MaxNConst;

  I = MultBse;
  D = 1;
  return( & Multplier );
  }
 else {
  I = 0;
  D = 0;
  return( 0 );
  }
 }  // end( SubGrad::ReadMult() )

/*--------------------------------------------------------------------------*/

HpNum SubGrad::ReadLBMult( cIndex wFi )
{
 if( wFi!= Inf<Index>() )
  throw( NDOException( "SubGrad::ReadBestFiVal: wFi != INF" ) );

 return( 0 );

 }  // end( SubGrad::ReadLBMult() )

/*--------------------------------------------------------------------------*/
/*------------- METHODS FOR ADDING / REMOVING / CHANGING DATA --------------*/
/*--------------------------------------------------------------------------*/

void SubGrad::AddVariables( Index NNwVrs , cLMRow IVs )
{
 SgRow Cnstr = 0;  // the i-th constraint
 Index_Set SGBse1 = 0;
 Index_Set SGBse2;
 Index h;

 // some exception - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( NumVar >= MaxNumVar ) // if no space for any new variable, return
  return;

 if( !NNwVrs )             // no variables to be added, just return
  return;

 if( NumVar + NNwVrs > MaxNumVar )  // not enough space for all the new vars
  NNwVrs = MaxNumVar - NumVar;      // put in only the first ones

 if( IVs )  { // check if any of the variables is strictly nonzero- - - - - -
  VectAssign( LambdaBar + NumVar , IVs , NNwVrs );
  SetLambda();
  }
 else {
  VectAssign( LambdaBar + NumVar , LMNum( 0 ) , NNwVrs );
  if( KpBstL )
   VectAssign( LambdaBest + NumVar , LMNum( 0 ) , NNwVrs );
  }

 cIndex_Set SGBse;    // dictionary vector
 Index SGBDim;        // dimension of direction or the subgradient

 // recover the restriction of the constraints to the new subspace - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( MaxNConst ) {  // the case of knapsack constraints is treated separately
	            // - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  Cnstr = new SgNum[ NNwVrs ];
  for( Index i = 0 ; i < MaxNConst ; i++ )
   if( CnstBeg[ i ] < Inf<Index>() ) {
    // get (a part of) the knapsack constraint  - - - - - - - - - - - - - - -
    //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    SGBDim = Oracle->GetGi( Cnstr , SGBse , i , NumVar , NumVar + NNwVrs );

    if( SGBDim ) {
     SGBse1 = new Index[ NNwVrs + 1 ];
     if( !SGBse ) {
      SGBse2 = Sparsify( Cnstr , SGBse1 , NNwVrs );
      *SGBse2 = Inf<Index>();
      SGBDim = SGBse2 - SGBse1;
      }
     else
      *VectAssign( SGBse1 , SGBse ) = Inf<Index>();

     for( Index j = 0 ; j < SGBDim ; j++ )
      if( ABS( Cnstr[ j ] - 1 ) >= Oracle->GetBndEps() )
       throw( NDOException(
	       "SubGrad::AddVariables: this is not a knapsack constraint.") );

     // write the constraint - - - - - - - - - - - - - - - - - - - - - - - - -
     //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

     for( h = CnstBeg[ i ] ; CnstNxt[ h ] < MaxNumVar ;  )
      h = CnstNxt[ h ];

     SGBse2 = SGBse1;
     for( Index h1 ; ( h1 = *(SGBse2++) ) < Inf<Index>() ; ) {
      CnstNxt[ h ] = h1;
      h = h1;
      }
     CnstNxt[ h ] = MaxNumVar;

     delete[] SGBse1;
     }  // end( SGBDim )
    }  // end scan of knapsack constraints - - - - - - - - - - - - - - - - - -
  }  // end knapsack part- - - - - - - - - - - - - - - - - - - - - - - - - - -
     //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 delete[] Cnstr;
 delete[] SGBse1;

 // get the restriction of DirM1 to new subspace - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 Index strt;
 if( deflection ) {  // if no deflection is done, Dir does not exist - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // if the projection is demanded, DirM1 needs all the unprojected variables;
  // however, the new variables can be projected independently from the other
  // ones if the feasible region consists of the box constraints

  if( ( SGPar1 & PrjDIRM1 ) && MaxNConst )
   strt = 0;
  else
   strt = NumVar;

  // DirPos == first => DirM1 is at the position MaxNConst + 1
  // DirPos == second => DirM1 is at the position MaxNConst + 2
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  SGBDim = Oracle->GetGi( dir + strt , SGBse , MaxNConst + ( first ? 1 : 2 ) ,
			  strt , NumVar + NNwVrs );

  if( SGBse )   // if it is given in sparse format, densify it
   Densify( dir + strt , SGBse , SGBDim , NumVar + NNwVrs , strt );

  if( SGBDim ) {
   if( SGPar1 & PrjDIRM1 )
    ProjTng( dir, LambdaBar , strt , NumVar + NNwVrs );
   NrmDir = Norm( dir , NumVar + NNwVrs );
   dM1GkDone = false;
   }
  }  // end( direction update )- - - - - - - - - - - - - - - - - - - - - - - -

 // get the restriction of Gi[wFi] to new subspace - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // if the subgradient projection is demanded, Gi[wFi] needs all the
 // unprojected variables; however, the new variables can be projected
 // independently from the other ones if the feasible region only consists of
 // the box constraints

 if( ( SGPar1 & PrjSG ) && MaxNConst )
  strt = 0;
 else
  strt = NumVar;

 SGBDim = Oracle->GetGi( Gi + strt , SGBse , MaxNConst , strt ,
			 NumVar + NNwVrs );
 if( SGBDim ) {
  if( SGBse )   // if it is given in sparse format, densify it
   Densify( Gi + strt , SGBse , SGBDim , NumVar + NNwVrs , strt );
  if( SGPar1 & PrjSG )
   ProjTng( Gi , LambdaBar , strt , NumVar + NNwVrs );
  NrmGi = Norm( Gi , NumVar + NNwVrs );
  dM1GkDone = false;
  }

 // get the new upper and lower bounds - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 for( register Index i = NumVar ; i < NumVar + NNwVrs ; i++ ) {
  lb[ i ] = Oracle->GetUC( i ) ? - Inf<LMNum>() : 0;
  ub[ i ] = Oracle->GetUB( i );
  }

 // update NumVar - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 NumVar += NNwVrs;

 UpdtLowerBound();

 }  // end( SubGrad::AddVariables() )

/*--------------------------------------------------------------------------*/

void SubGrad::RemoveVariables( cIndex_Set whch , Index hwmny )
{
 Index cnstr, rr;
 bool FoundVar;

 // if the set of variables to be removed is not appropriately set, throw
 // exception  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( ! ( whch && whch[ hwmny ] == Inf<Index>() ) &&
     ! ( ( ! whch ) && ( ! hwmny ) ) )
  throw( NDOException( "SubGrad::RemoveVariables: " ) );

 // check if any of the variables is strictly nonzero- - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 bool RmvNnz = false;
 if( whch ) {
  cIndex_Set tw = whch;
  for( Index h ; ( h = *(tw++) ) < Inf<Index>() ; )
   if( ABS( LambdaBar[ h ] ) > Oracle->GetBndEps() ) {
    RmvNnz = true;
    break;
    }
  }
 else
  for( Index h = 0 ; h < NumVar ; h++ )
   if( ABS( LambdaBar[ h ] ) > Oracle->GetBndEps() ) {
    RmvNnz = true;
    break;
    }

 if( RmvNnz ) {  // if so, first make a change of current point- - - - - - - -
  LMRow NewL = new LMNum[ NumVar ];
  if( whch ) {
   VectAssign( NewL , LambdaBar , NumVar );
   VectAssign( NewL , LMNum( 0 ) , whch );
   }
  else
   VectAssign( NewL , LMNum( 0 ) , NumVar );
  SetLambda();
  delete[] NewL;
  }

 // remove the knapsack constraints- - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( MaxNConst ) {
  if( whch ) {
   cIndex_Set tw = whch;
   for( Index h ; ( h = *(tw++) ) < Inf<Index>() ; ) {
    // find the constraint - - - - - - - - - - - - - - - - - - - - - - - - - -
    for( cnstr = 0 ; cnstr < MaxNConst ; cnstr++ ) {
     if( CnstBeg[ cnstr ] < Inf<Index>() )
      for( rr = CnstBeg[ cnstr ] ; CnstNxt[ rr ] < MaxNumVar ;  ) {
       if( rr == h ) {
    	FoundVar = true;
        break;
        }
       rr = CnstNxt[ rr ];
       }

     if( FoundVar )
      break;
     }

    // remove variable - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    if( CnstBeg[ cnstr ] == h ) {
     if( CnstNxt[ h ] >= MaxNumVar ) {
      CnstBeg[ cnstr ] = Inf<Index>();
      Oracle->Deleted( cnstr );
      }
     else
      CnstBeg[ cnstr ] = CnstNxt[ h ];
     }
    else {
     rr = CnstBeg[ h ];
     while( CnstNxt[ rr ] != h )
      rr = CnstNxt[ rr ];
     CnstNxt[ rr ] = CnstNxt[ h ];
     }
    CnstNxt[ h ] = -1;
    }
   }  // end( if( whch ) )
  else {
   for( Index i = MaxNConst ; i-- ; )
    if( CnstBeg[ i ] < Inf<Index>() )
     Oracle->Deleted( i );
   VectAssign( CnstBeg , Index( Inf<Index>() ) , MaxNConst );
   VectAssign( CnstVol , HpNum( Inf<HpNum>() ) , MaxNConst );
   VectAssign( CnstNxt , SIndex(-1) , MaxNumVar );

  }  // end( total remotion ) - - - - - - - - - - - - - - - - - - - - - - - -
  }  // Knapsack - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 // remove components in subgradient Gi[wFi], direction Dir, LambdaBar,
 // LambdaBest, l and u- - - - - - - - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( whch ) {
  Compact( Gi , whch , NumVar );   // subgradient
  NrmGi = Norm( Gi , NumVar );

  if( deflection ) {
   Compact( dir , whch , NumVar);  // Dir
   NrmDir = Norm( dir , NumVar );
   }

  dM1GkDone = false;

  Compact( ub , whch , NumVar );   // lower and upper bound
  Compact( lb , whch , NumVar );

  Compact( LambdaBar , whch , NumVar );    // LambdaBar
  if( KpBstL )
   Compact( LambdaBest , whch , NumVar );  // LambdaBest

  NumVar -= hwmny;
  }
 else
  NumVar = 0;

 UpdtLowerBound();

 }  // end( SubGrad::RemoveVariables() )

/*--------------------------------------------------------------------------*/

void SubGrad::ChgFiV( cIndex wFi )
{
 throw( NDOException( "SubGrad::ChgFiV: not implemented yet" ) );

 }  // end( SubGrad::ChgFiV() )

/*--------------------------------------------------------------------------*/

void SubGrad::ChgSbG( cIndex strt , Index stp , cIndex wFi )
{
 throw( NDOException( "SubGrad::ChgSbG: not implemented yet" ) );

 }  // end( SubGrad::ChgSbG() )

/*--------------------------------------------------------------------------*/
/*------------------------------ DESTRUCTOR --------------------------------*/
/*--------------------------------------------------------------------------*/

SubGrad::~SubGrad()
{
 // output times and statistics - - - - - - - - - - - - - - - - - - - - - - -

 VLOG( 1 , endl << "Total Fi() evaluations = " << FiEvaltns );
 VLOG2( 1 , NDOt , endl << "Tot. time (s): " << NDOt->Read() );
 VLOG2( 1 , Oracle , endl << "Fi() time (s): " << Oracle->FiTime() );
 VLOG( 1 , endl );

 // memory deallocation  - - - - - - - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( Oracle ) {
  MemDealloc();

  if( deflection )
   deflection->Format( );

  if( stepsize )
   stepsize->Format();
  }
 }  // end( ~SubGrad )

/*--------------------------------------------------------------------------*/
/*-------------------------- PROTECTED METHODS -----------------------------*/
/*--------------------------------------------------------------------------*/

void SubGrad::FiAndGi( cIndex wFi )
{
 // pass Lambda to the Oracle - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( LHasChgd ) {  // if Lambda1 has changed, set the new one to the Oracle
  Oracle->SetLambda( Lambda );
  LHasChgd = false;
  }

 // if outer iteration, evaluate Fi in Lambda1  - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( wFi == Inf<Index>() ) {  // call Fi()- - - - - - - - - - - - - - - - - -
  FiLambda = Oracle->Fi( );
  FiEvaltns++;

  if( FiLambda == -Inf<HpNum>() )    // Fi() unbounded
   return;

  // update FiBest, if necessary - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( FiLambda < FiBest ) {
   FiBest = FiLambda;
   if( KpBstL )
    VectAssign( LambdaBest , Lambda , NumVar );
   }
  }  // end( function evaluation - - - - - - - - - - - - - - - - - - - - - - -

 // if Lambda is unfeasible, a new subgradient is computed; otherwise a new
 // constraint is added - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( FiLambda == Inf<HpNum>() )
  EvalConst( );  // the item is a constraint
 else
  EvalGi( wFi ); // the item is a subgradient

 }  // end( SubGrad::FiAndGi() )

/*--------------------------------------------------------------------------*/

void SubGrad::FormD( void ) {
 // project the subgradient and compute its norm - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // take in account some comments and details about FormD():
 // (i) the method computes the sugradient in Lambda, if not already computed;
 // (ii) the method computes the deflection coefficient Alpha, and a fortiori
 //      the direction Dir;
 // (iii) the method computes the stepsize Nu using the deflection-restricted
 //       scheme;
 // (iv) if the variables-set has been changed after a call to GetFiStatus,
 //      even the projection of both Gi and Dir is done in AddVariables or
 //      RemoveVariables.
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( ( fs == FiOracle::kFiNorm ) && ( ! DoSS ) && ( SGPar1 & PrjSG ) ) {
  // the projection of Gi[wFi] is already performed but relatively to the
  // tangent cone at the tentative point of the previous iteration: if a NS
  // occurred, the projection has to be recomputed

  cIndex_Set SGBse;
  Index SGBDim = Oracle->GetGi( Gi , SGBse );
  if( SGBse )
   Densify( Gi , SGBse , SGBDim , NumVar );

  ProjTng( Gi , LambdaBar );
  NrmGi = Norm( Gi , NumVar );

  dM1GkDone = false;
  }

 // subgradient log: print the information relative to Gi[wFi] - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 Log4();

 // project DirM1 and compute its norm - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if ( ( SGPar1 & PrjDIRM1 ) && ( fs == FiOracle::kFiNorm ) )  {
  if( SGPar1 & PrjDIR ) {  // DirM1 is projected over the tangent cone at the
	                   // previous LambdaBar - - - - - - - - - - - - - - -
   if ( DoSS )  {  // the projection is no longer valid, re-project DirM1
    ProjTng( dir , LambdaBar );
    NrmDir = Norm( dir , NumVar );
    dM1GkDone = false;
    }
   }
  else {  // project DirM1 - - - - - - - - - - - - - - - - - - - - - - - - - -
   ProjTng( dir , LambdaBar );
   NrmDir = Norm( dir , NumVar );

   dM1GkDone = false;
   }
  }  // end DirM1 projection - - - - - - - - - - - - - - - - - - - - - - - - -
     //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 // (re)compute the scalar product between DirM1 and Gi[wFi] at Lambda - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( deflection ) {
  if( !dM1GkDone )
   dM1Gk = ScalarProduct( Gi , dir , NumVar );
  }
 else {  // if no deflection is done, Dir <==> Gi[wFi]
  dGk = NrmDir = NrmGi;
  if( !dM1GkDone )
   dM1Gk = dGk;
  }

 dM1GkDone = true;

 // stepsise and deflection computation  - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( SGPar3 & StepRST ) {  // stepsize-restricted approach - - - - - - - - - -

  if( deflection ) {  // deflection coefficient computation
   deflection->NewDEF(  );
   alpha = deflection->GetDFLCoeff();
   }

  VLOG( 1 , " ~ alpha = " << alpha );
  }
 else {  // deflection-restricted approach - - - - - - - - - - - - - - - - - -
         //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  dGk = dM1Gk;  // this scheme doesn't distinguish between dGk and dM1Gk

  if( ! InnIter )  // stepsize computation - - - - - - - - - - - - - - - - - -
   stepsize->NewStep( );

  step = stepsize->GetStepsize( InnIter );
  VLOG( 1 , " ~ step = " << step );

  deflection->NewDEF( );  // no deflection is not allowed
  alpha = deflection->GetDFLCoeff();
  VLOG( 1 , " ~ alpha = " << alpha );

  // safe rule - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( ( SGPar3 & SafeRULE )  && !InnIter ) {
   HpNum target =  ( step * NrmDir ) /
                   ( FiBar - stepsize->GetLev() + step * NrmDir );
   alpha = min( 1.0 , max( alpha , target ) );
   VLOG( 1 , " ~ fixed alpha = " << alpha );
   }
  }  // end (deflection-restricted approach) - - - - - - - - - - - - - - - - -
     //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 // tell the oracle to keep track the direction- - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // DirPos == first => DirM1 is at the position MaxNConst + 1
 // and Dir is registered in the free entry MaxNConst + 2,
 // and vice versa for the case DirPos == second

 if( deflection )
  if(  alpha == 1 ) {  // pass the subgradient to the FiOracle object
   HpNum Mlt = 1;
   Index NmSt = MaxNConst;
   Oracle->Aggregate( &Mlt , &NmSt , 1 ,
		      MaxNConst + ( DirPos == first ? 2 : 1 ) );
   }
  else {
   HpNum Mlt[ 2 ] = { alpha , 1.0 - alpha };
   if( DirPos == first ) {
    Index NmSt[ 2 ] = { MaxNConst , MaxNConst + 1 };
    Oracle->Aggregate( Mlt , NmSt , 2 , MaxNConst + 2 );
    }
   else {
    Index NmSt[ 2 ] = { MaxNConst , MaxNConst + 2 };
    Oracle->Aggregate( Mlt , NmSt , 2 , MaxNConst + 1 );
    }
   }

 // print Dir norm - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 VLOG( 1 , " ~ trial |D| = " << NrmDir );

 } // end( SubGrad::FormD() )

/*--------------------------------------------------------------------------*/

void SubGrad::SaveDir( void )
{
 // compute the linearization error of Dir at LambdaBar- - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( alpha == 1 )
  Epsilon = Sigma;
 else
  Epsilon = alpha * Sigma + ( HpNum( 1 ) - alpha ) * Epsilon;

 // change LambdaHat - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( SGPar4 )
  ChgLambdaHat();

 // save the search direction Dir  - - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( deflection ) {
  // Dir is the convex combination of DirM1 and Gi[wFi]- - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( alpha == 1 )
   VectAssign( dir , Gi , NumVar );
  else
   VectAdd( dir , Gi , alpha , dir , HpNum( 1 ) - alpha , NumVar );

  // swap the position for DirM1 and Dir - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( DirPos == first )
   DirPos = second;
  else
   DirPos = first;

  // project Dir on tangent cone at LambdaBar  - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( SGPar1 & PrjDIR )
   ProjTng( dir , LambdaBar );

  // compute the (squared) norm of the direction and its scalar product with
  // the subgradient   - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  NrmDir = Norm( dir , NumVar );
  dGk = ScalarProduct( Gi , dir , NumVar );

  }  // end( deflection )  - - - - - - - - - - - - - - - - - - - - - - - - - -
     //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 // stepsize computation - - - - - - - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( SGPar3 & StepRST ) {
  if( ! InnIter )
   stepsize->NewStep( );

  step = stepsize->GetStepsize( InnIter );
  VLOG( 1 , " ~ step = " << step );

  // safe rule - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( SGPar3 & SafeRULE ) {
   stepsize->SetMaxBeta( alpha );
   step = stepsize->GetStepsize( InnIter );
   VLOG( 1 , " ~ (fixed) step = " << step );
   }

  }  // end stepsize-restricted scheme - - - - - - - - - - - - - - - - - - - -
 }  // end( SubGrad::SaveDir() )

/*--------------------------------------------------------------------------*/

void SubGrad::GotoLambda1( void )
{
 // compute the linearization error of Dir at the new LambdaBar  - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( FiBar < Inf<HpNum>() && deflection ) {
  HpNum SclPrd;
  if( ! LHasProj )
   SclPrd = step * NrmDir;
  else {
   SclPrd = 0;
   for( Index i = 0 ; i < NumVar ; i++ )
    SclPrd += dir[ i ] * ( LambdaBar[ i ] - Lambda[ i ] );
   }

  Epsilon += FiLambda - FiBar + SclPrd;
  }

 // move Lambda to LambdaBar= - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 VectAssign( LambdaBar , Lambda , NumVar );
 if( ! InnIter )
  FiBar = FiLambda;

 }  // end( SubGrad::GotoLambda1() )

/*--------------------------------------------------------------------------*/

void SubGrad::FormLambda1( void )
{
 if( deflection )
  for( Index i = 0 ; i < NumVar ; i++ )
   Lambda[ i ] = LambdaBar[ i ] - step * dir[ i ];
 else
  for( Index i = 0 ; i < NumVar ; i++ )
   Lambda[ i ] = LambdaBar[ i ] - step * Gi[ i ];

 LHasChgd = true;

 }  // end( SubGrad::FormLambda1() )

/*--------------------------------------------------------------------------*/
/*--------------------------- PRIVATE METHODS ------------------------------*/
/*--------------------------------------------------------------------------*/

void SubGrad::EvalConst( void )
{
 // some exception - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( ! Oracle->NewGi( 0 ) )
  throw( NDOException(
		   "SubGrad::EvalCons: the constraint cannot be computed" ) );

 if( ! MaxNConst )
  throw( NDOException("SubGrad::EvalCons: there are not constraints" ) );

 // fetch the item (constraint) from the Oracle- - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 cIndex_Set SGBse;
 Index SGBDim = Oracle->GetGi( Gi , SGBse );

 if( SGBDim ) {  // if it is not a zero constraint - - - - - - - - - - - - - -
  Index wh = 0;
  while( CnstBeg[ wh ] < Inf<Index>() )
   wh++;

  if( wh == MaxNConst )
   throw( NDOException(
		   "SubGrad::EvalCons: max number of constraints reached" ) );

  // sparsify the item - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  Index_Set SGBse1 = new Index[ NumVar + 1 ];
  Index_Set SGBse2;
  if( !SGBse ) {
   SGBse2 = Sparsify( Gi , SGBse1 , NumVar );
   *SGBse2 = Inf<Index>();
   SGBDim = SGBse2 - SGBse1;
   }
  else
   *VectAssign( SGBse1 , SGBse ) = Inf<Index>();

  // only knapsack constraints are allowed - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  for( Index i = 0 ; i < SGBDim ; i++ )
   if( ABS( Gi[ i ] - 1 ) >= Oracle->GetBndEps() )
    throw( NDOException(
		   "SubGrad::EvalCons: this is not a knapsack constraint" ) );

  // write the constraint in CnstNxt - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  SGBse2 = SGBse1;
  Index h = *(SGBse2++);
  CnstBeg[ wh ] = h;
  for( Index h1 ; ( h1 = *(SGBse2++) ) < Inf<Index>() ; ) {
   CnstNxt[ h ] = h1;
   h = h1;
   }

  CnstNxt[ h ] = MaxNumVar;
  CnstVol[ wh ] = Oracle->GetVal();
  Oracle->SetGiName( wh );

  delete[] SGBse1;

  }  // end( if( SGBDim ) )- - - - - - - - - - - - - - - - - - - - - - - - - -
 }  // end( SubGrad::EvalConst() )

/*--------------------------------------------------------------------------*/

void SubGrad::EvalGi( Index wFi )
{
 // some exception - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( wFi > NrFi )
  wFi = Inf<Index>();

 if( wFi && ! Oracle->NewGi( wFi ) )
  throw( NDOException(
		    "SubGrad::EvalGi: the subgradient cannot be computed" ) );

 // fetch the item (subgradient) from the Oracle - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 cIndex_Set SGBse;
 Index SGBDim = Oracle->GetGi( Gi , SGBse , wFi? Inf<Index>() : MaxName );

 // densify the item  - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( SGBse )
  Densify( Gi , SGBse , SGBDim , NumVar );

 // compute the linearization error of Gi[wFi] at Lambda1 - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 Sigma = Oracle->GetVal();
 if( SGPar4 )
  SigmaHat = Sigma;

 // write the item  - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( wFi && !InnIter )
  Oracle->SetGiName( MaxNConst );
 GiEvaltns++;

 }  // end( SubGrad::EvalGi() )

/*--------------------------------------------------------------------------*/
/*------------------------ Projection methods ------------------------------*/
/*--------------------------------------------------------------------------*/

bool SubGrad::ProjFsb( LMRow Lpt , bool & EmptySet )
{
 // project Lpt in the feasible set, and if some changes occurs in Lpt
 // return true - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( ( ! BoxConst ) && ( ! MaxNConst ) )
  return( false );
    
 if( MaxNConst && !Q2KNP )
  throw( NDOException( "SubGrad::ProjFsb: the pointer to SetQKNP object is not assigned" ) );

 bool LptHasChgd = false;  // it is true, if Lpt changes
 EmptySet = false;         // it is true, if the feasible set is empty

 LMRow PjLpt = 0; // the projected point

 // the case with knapsack constraints is treated apart - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( MaxNConst ) {  // there are knapsack constraints - - - - - - - - - - - -
                    //- - - - - - - - - - - - - - - - - - - - - - - - - - - -
  Index Nvar;       // number of variables in the i-th knapsack constraint
  Row pC, pD, pA, pB;  // problem description
  Number pV;

  PjLpt = new LMNum[ NumVar ];  // reset the projected point
  VectAssign( PjLpt , HpNum( Inf<HpNum>() ) , NumVar );

  for( Index i = 0 ; i < MaxNConst ; i++ )
   if( CnstBeg[ i ] < Inf<Index>() ) {
    // count the number of variables of the i-th knapsack constraint- - - - -

    Nvar = 0;
    for( Index j = CnstBeg[ i ] ; j < MaxNumVar ; j = CnstNxt[ j ] )
     Nvar++;

    // allocate memory for the i-th knapsack problem- - - - - - - - - - - - -

    pC = new Number[ Nvar ];  // linear costs
    pD = new Number[ Nvar ];  // quadratic costs,
    pA = new Number[ Nvar ];  // lower bounds
    pB = new Number[ Nvar ];  // upper bounds
    pV = CnstVol[ i ];        // volume

    // construct the data of the i-th quadratic knapsack problem- - - - - - -

    Nvar = 0;
    for( Index j = CnstBeg[ i ] ; j < MaxNumVar ; j = CnstNxt[ j ]  ) {
     pB[ Nvar ] = ub[ j ];
     pA[ Nvar ] = lb [ j ];
     pC[ Nvar ] = -2.0 * Lpt[ j ];
     pD[ Nvar++ ] = 1;
     }

    Q2KNP->LoadSet( Nvar , pC , pD , pA , pB , pV , false );
    delete[] pC;
    delete[] pD;
    delete[] pA;
    delete[] pB;

    // solve the i-th quadratic knapsack problem- - - - - - - - - - - - - - -
    CQKnPClass::CQKStatus status = Q2KNP->SolveKNP();

    cRow Lpt1 = 0;  // auxiliary array
    switch( status ) {
    case( CQKnPClass::kOK ):       // optimal solution found
     Q2KNP->KNPGetFO();
    case( CQKnPClass::kStopped ):  // optimization process stopped
     // get the projected point

      Lpt1 = Q2KNP->KNPGetX( );
      for( Index j = CnstBeg[ i ] ; j < MaxNumVar ; j = CnstNxt[ j ]  )
       PjLpt[ j ] = *(Lpt1++);

      break;

     case( CQKnPClass::kUnfeasible ):
      EmptySet = true;
      break;

     default:
      throw( NDOException( "SubGrad::ProjFsb: solution not found" ) );
     }
    }  // end( solution of quadratic knapsack problems )

  if( ! EmptySet ) {
   // project the components that are not in any knapsack problem - - - - - -

   for( Index j = 0 ; j < NumVar ; j++ )
    if( PjLpt[ j ] == Inf<HpNum>() ) {
     if( Lpt[ j ] > ub[ j ] )
      PjLpt[ j ] = ub[ j ];
     else
      if( Lpt[ j ] < lb[ j ] )
       PjLpt[ j ] = lb[ j ];
      else
       PjLpt[ j ] = Lpt[ j ];
     }

   // compare Lpt and PjLpt - - - - - - - - - - - - - - - - - - - - - - - - -

   if( ! EqualVect( Lpt , PjLpt , NumVar ) ) {
    VectAssign( Lpt , PjLpt , NumVar );
    LptHasChgd = true;
    }
   }
  }
 else {  // only box constraints- - - - - - - - - - - - - - - - - - - - - - -
         // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  for( Index i = 0 ; i < NumVar ; i++ )
   if( ub[ i ] - lb[ i ] < - Oracle->GetBndEps() ) {
    EmptySet = true;
    break;
    }

  // project Lpt in the box constraints - - - - - - - - - - - - - - - - - - -

  if( ! EmptySet )
   for( Index i = 0 ; i < NumVar ; i++ )
    if( Lpt[ i ] > ub[ i ] ) {
     Lpt[ i ] = ub[ i ];
     LptHasChgd = true;
     }
    else
     if( Lpt[ i ] < lb[ i ] ) {
      Lpt[ i ] = lb[ i ];
      LptHasChgd = true;
      }

  }  // end( else ) - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 delete[] PjLpt;
 return( LptHasChgd );

 }  // end( SubGrad::ProjFsb() )

/*--------------------------------------------------------------------------*/

void SubGrad::ProjTng( SgRow Gpt , cLMRow Lpt , cIndex strt , Index stp )
{
 // project Gpt in the tangent cone at Lpt  - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( ( ! BoxConst ) && ( ! MaxNConst ) )
  return;
    
 if( MaxNConst && !Q2KNP )
  throw( NDOException( "SubGrad::ProjTng: the pointer to SetQKNP object is not assigned" ) );

 if( stp > NumVar )
  stp = NumVar;

 LMRow PjGpt = 0; // the projected point
 HpNum lhs;

 if( MaxNConst ) {  // there are knapsack constraints - - - - - - - - - - - -
                    //- - - - - - - - - - - - - - - - - - - - - - - - - - - -
  Index Nvar;       // number of variables in the i-th knapsack constraint
  Row pC, pD, pA, pB;  // problem description
  Number pV;

  PjGpt = new LMNum[ NumVar ];  // reset the projected point
  VectAssign( PjGpt , HpNum( Inf<HpNum>() ) , NumVar );

  for( Index i = 0 ; i < MaxNConst ; i++ )
   if( CnstBeg[ i ] < Inf<Index>() ) {
    // count the number of variables of the i-th knapsack constraint- - - - -

    lhs = 0;
    Nvar = 0;
    for( Index j = CnstBeg[ i ] ; j == MaxNumVar ; j = CnstNxt[ j ] ) {
     Nvar++;
     lhs += Lpt[ j ];
     }

    HpNum Eps = max( HpNum(1) , ABS( CnstVol[ i ] ) ) * Oracle->GetBndEps();
    if( ABS( lhs - CnstVol[ i ] ) <= Eps ) {
     // allocate memory for the i-th knapsack problem - - - - - - - - - - - -
     pC = new Number[ Nvar ];  // linear costs
     pD = new Number[ Nvar ];  // quadratic costs,
     pA = new Number[ Nvar ];  // lower bounds
     pB = new Number[ Nvar ];  // upper bounds
     pV = 0;                   // volume

     // construct the data of the i-th quadratic knapsack problem - - - - - -

     Nvar = 0;
     for( Index j = CnstBeg[ i ] ; j == MaxNumVar ; j = CnstNxt[ j ] ) {
      if( Lpt[ j ] == ub [ j ] )
       pB[ Nvar ] = 0;
      else
       pB[ Nvar ] = Inf<HpNum>();
      if( Lpt[ j ] == lb [ j ] )
       pA[ Nvar ] = 0;
      else
       pA[ Nvar ] = -Inf<HpNum>();
  	  pC[ Nvar ] = -2.0 * Gpt[ j ];
  	  pD[ Nvar++ ] = 1;
  	  }

     // solve the i-th quadratic knapsack problem - - - - - - - - - - - - - -
     Q2KNP->LoadSet( Nvar , pC , pD , pA , pB , pV , false );
     delete[] pC;
     delete[] pD;
     delete[] pA;
     delete[] pB;

     if( Q2KNP->SolveKNP() == CQKnPClass::kOK ) {
      cRow Gpt1 = Q2KNP->KNPGetX( );
      for( Index j = CnstBeg[ i ] ; j == MaxNumVar ; j = CnstNxt[ j ] )
       PjGpt[ j ] = *(Gpt1++);
      }
     else
      throw( NDOException( "SubGrad::ProjTng: solution not found" ) );
     }
    }  // end i-th knapsack constraint

  //  project the components that are not in any knapsack problem - - - - - -

  for( Index j = 0 ; j < NumVar ; j++ )
   if( PjGpt[ j ] == Inf<HpNum>() ) {
    if( ( Lpt[ j ] != ub[ j ] || Gpt[ j ] <= 0 ) &&
	( Lpt[ j ] != lb[ j ] || Gpt[ j ] >= 0 ) )
     PjGpt[ j ] = Gpt[ j ];
    else
     PjGpt[ j ] = 0;
    }
  }
 else {  // only box constraints- - - - - - - - - - - - - - - - - - - - - - -
         // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  for( Index j = strt ; j < stp ; j++ )
   if( ( Lpt[ j ] >= ub[ j ] && Gpt[ j ] <= 0 ) ||
       ( Lpt[ j ] <= lb[ j ] && Gpt[ j ] >= 0 ) )
    Gpt[ j ] = 0;

  }  // end( else ) - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 delete[] PjGpt;

 }  // end( SubGrad::ProjTng() )

/*--------------------------------------------------------------------------*/

void SubGrad::InitStepsize( void )
{
 step = 0;
 stepsize->Format( );

 }  // end( SubGrad::InitStepsize() )

/*--------------------------------------------------------------------------*/

void SubGrad::InitDeflection( void )
{
 // set to 0 all the information connected with the direction- - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 NrmDir = dGk = dM1Gk = 0;
 Epsilon = HatEpsilon = 0;
 alpha = 1;
 dM1GkDone = false;

 // initialize the direction dir to zero - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 VectAssign( dir , HpNum( 0 ) , NumVar );

 // if the pointer to Deflection object is not 0, invoke Format() in
 // Deflection,  otherwise some conditions have to be taken into account:
 // (i) the stepsize-restricted scheme must be used and (ii) the projection
 // of the direction is useless
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( deflection ) {
  DirPos = second;
  deflection->Format( );
  }
 else {
  SGPar3 |= StepRST;
  SGPar1 &= ~( PrjDIRM1 | PrjDIR );
  }
 }  // end( SubGrad::InitDeflection() )

/*--------------------------------------------------------------------------*/

void SubGrad::ChgLambdaHat( void )
{
 // compute the LambdaHat- - - - - - - - - - - - - - - -  - - - - - - - - - - 
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( alpha == 1 ) {
  HatEpsilon = SigmaHat;
  VectAssign( LambdaHat , Lambda , NumVar );
  FiHat = FiLambda;
  }
 else {
  // compute the linearization error of DirM1 in LambdaHat- - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  HpNum SclPrd = 0;
  for( Index i = 0 ; i < NumVar ; i++ )
   SclPrd += dir[ i ] * ( LambdaHat[ i ] - Lambda[ i ] );

  HatEpsilon += alpha * ( FiLambda - FiHat + SclPrd );

  // compute the linearization error of Gi[wFi] in LambdaHat  - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  SclPrd = 0;
  for( Index i = 0 ; i < NumVar ; i++ )
   SclPrd += Gi[ i ] * ( Lambda[ i ] - LambdaHat[ i ] );

  SigmaHat += (1.0 - alpha ) * ( FiHat - FiLambda + SclPrd );

  // compute the LambdaHat and FiHat - - - - - - - - - -  - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  VectAdd( LambdaHat , Lambda , alpha , LambdaHat , 1.0 - alpha , NumVar );
  FiHat = alpha * FiLambda + (1.0 - alpha) * FiHat;

  // compute the linearization error of Dir in LambdaHat  - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  HatEpsilon = alpha * SigmaHat + ( HpNum(1) - alpha ) * HatEpsilon;
  }  // end( alpha < 1 )- - - - - - - - - - - - - - - - - - - - - - - - - - -

 VLOG( 1 , " ~ HFi = " << -FiHat );

 }  // end( SubGrad::ChgHatlambda() )

/*--------------------------------------------------------------------------*/

void SubGrad::UpdateSigma( void )
{
 // some exception - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( ! deflection )
  throw( NDOException( "SubGrad::UpdateSigma: this is not allowed" ) );

 // compute the linearization error of Gi[wFi] at Lambdabar - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 Sigma += FiBar - FiLambda;
 if( ! LHasProj )
  Sigma -= step * ScalarProduct( Gi , dir , NumVar );
 else {
  HpNum SclPrd = 0;
  for( Index i = 0 ; i < NumVar ; i++ )
   SclPrd +=  Gi[ i ] * ( Lambda[ i ] - LambdaBar[ i ] );
  Sigma += SclPrd;
  }

 VLOG( 1 , " ~ Sigma1 = " << Sigma );

 }  // end( SubGrad::UpdateSigma() )

/*--------------------------------------------------------------------------*/

void SubGrad::UpdtLowerBound( void )
{
 cHpNum LwrBnd = Oracle->GetLowerBound();
 if( LwrBnd > LowerBound ) {
  LowerBound = LwrBnd;
  TrueLB = ( LwrBnd > Oracle->GetMinusInfinity() );
  }
 }  // end( SubGrad::UpdtLowerBound() )

/*--------------------------------------------------------------------------*/

void SubGrad::MemDealloc( void )
{
 // the point - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 delete[] LambdaBar;
 delete[] Lambda;
 delete[] LambdaBest;
 delete[] LambdaHat;

 // box and simplex constraints - - - - - - - - - - - - - - - - - - - - - - -

 delete[] lb;
 delete[] ub;

 if( MaxNConst ) {
  delete[] CnstBeg;
  delete[] CnstNxt;
  delete[] CnstVol;
  }

 // the subgradients and direction- - - - - - - - - - - - - - - - - - - - - -

 delete[] Gi;
 delete[] dir;

 Deflection *StpDef = dynamic_cast<Deflection*>( stepsize );
 if( StpDef )
  delete stepsize;
 else {
  delete stepsize;
  if( deflection )
   delete deflection;
  }

 if( Q2KNP )
  delete Q2KNP;

 delete[] MultBse;

 }  // end( SubGrad::MemDealloc() )

/*--------------------------------------------------------------------------*/
/*-------------------------- End CLASS SubGrad -----------------------------*/
/*--------------------------------------------------------------------------*/
/*-------------------------- CLASS Stepsize --------------------------------*/
/*--------------------------------------------------------------------------*/
/** These methods provide derived classes of Stepsize (who are not friends of
    SubGrad, while Stepsize is) access to protected data of SubGrad. */

/*--------------------------------------------------------------------------*/

FiOracle* Stepsize::GetOracle( void )
{
 return( Solver->Oracle );
 }

/*--------------------------------------------------------------------------*/

HpNum Stepsize::GetCoeffDefl( void )
{
 return( Solver->alpha );
 }

/*--------------------------------------------------------------------------*/

cHpNum Stepsize::GetGiNorm( void )
{
 return( Solver->NrmGi );
 }

/*--------------------------------------------------------------------------*/

cHpNum Stepsize::GetDNorm( void )
{
 return( Solver->NrmDir );
 }

/*--------------------------------------------------------------------------*/

cHpNum Stepsize::GetdGk( void )
{
 return( Solver->dGk );
 }

/*--------------------------------------------------------------------------*/

HpNum Stepsize::GetdkM1Gk( void )
{
 return( Solver->dM1Gk );
 }

/*--------------------------------------------------------------------------*/

Index Stepsize::GetNItIcr( void )
{
 return( Solver->NItIncr );
 }

/*--------------------------------------------------------------------------*/

cHpNum Stepsize::ReadFkVal( void )
{
 return( Solver->FiLambda );
 }

/*--------------------------------------------------------------------------*/

cHpNum Stepsize::ReadFiBar( void )
{
 return( Solver->FiBar );
 }

/*--------------------------------------------------------------------------*/
/*-------------------------- CLASS Deflection ------------------------------*/
/*--------------------------------------------------------------------------*/
/** These methods provide derived classes of Deflection (who are not friends
    of SubGrad, while Deflection is) access to protected data of SubGrad. */

FiOracle* Deflection::GetOracle( void )
{
 return( Solver->Oracle );
 }

/*--------------------------------------------------------------------------*/

HpNum Deflection::GetStepsize( void )
{
 return( Solver->step );
 }

/*--------------------------------------------------------------------------*/

cHpNum Deflection::GetGiNorm( void )
{
 return( Solver->NrmGi );
 }

/*--------------------------------------------------------------------------*/

cHpNum Deflection::GetDNorm( void )
{
 return( Solver->NrmDir );
 }

/*--------------------------------------------------------------------------*/

cHpNum Deflection::GetdGk( void )
{
 return( Solver->dGk );
 }

/*--------------------------------------------------------------------------*/

cHpNum Deflection::GetSigma( void )
{
 return( Solver->Sigma );
 }

/*--------------------------------------------------------------------------*/

cHpNum Deflection::GetEpsilon( void )
{
 return( Solver->Epsilon );
 }

/*--------------------------------------------------------------------------*/

cHpNum Deflection::ReadFVal( void )
{
 return( Solver->FiLambda );
 }

/*--------------------------------------------------------------------------*/
/*------------------------- End File SubGrad.C -----------------------------*/
/*--------------------------------------------------------------------------*/

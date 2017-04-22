
\mainpage <H1>The CNDSM Project Documentation</H1>

\section intro Introduction

This file is a short user manual for the \c CNDSM project, a subgradient
method (SM) that solves two Lagrangian dual problems of the Fixed-Charge
Multicommodity Capacitated Network Design (FC-MCND) problem). Two kinds of
relaxations are provided: <em>Flow relaxation</em> (FR) and <em>Knapsack
relaxation</em> (KR). The method has the potentialities and possibilities
to deal with a lot of variants of SM such as the deflected, incremental and
projected one,  because it is implemented according the <em>object
oriented</em> paradigm. 

\subsection disclaimer Standard Disclaimer

This code is provided "as is", without any explicit or implicit warranty
that it will properly behave or it will suit you needs. Although codes
reaching the distribution phase have usually been extensively tested, we
cannot guarantee that they are absolutely bug-free (who can?). Any use of
the codes is at you own risk: in no case we could be considered liable for
any damage or loss you would eventually suffer, either directly or indirectly,
for having used this code. More details about the non-warranty attached to
this code are available in the license description file.

The code also comes with a "good will only" support: feel free to contact us
for any comments/critics/bug report/request help you may have, we will be
happy to try to answer and help you. But we cannot spend much time solving
your problems, just the time to read a couple of e-mails and send you fast
suggestions if the problem is easily solvable. Apart from that, we can't offer
you any support.

\subsection release Release

Current version is: 1.00

Current date is: April 22, 2017

\section howto Package description

This release comes out with the following files:

- \b doc/LGPL.txt: Description of the LGPL license.

- \b doc/refman.pdf: Pdf version of the manual.

- \b doc/html/*: Html version of the  manual.

- \b NDO/NDOSlver/NDOSlver.h: Contains the declaration of class \c NDOSolver.
  It is an <em>abstract class</em> with <em>pure virtual methods</em>, so that
  you cannot declare objects of type \c NDOSolver, but only of derived classes
  of its. \c NDOSolver offers a general interface for NonDifferentiable
  Optimization Solvers codes. The function to be minimized must be derived
  from the <em> FiOracle interface </em> .

- \b /MMCF/extlib/NDO/SubGrad/SubGrad.h:  Contains the declaration of class \c SubGrad, 
 implementing a <em> subgradient method </em> (SM) and offering multiple variants thereof. The 
 method is in fact  based on abstract formulae for the <em> stepsize rule</em>, (SR) and the 
 <em> deflection rule</em> (DR), namely the method relies on the objects of the friend classes 
 \a Stepsize  and \a Deflection, which have to return, respectively,
 the stepsize and the deflection coefficient. This class derives from \c NDOSolver, hence most of its 
 interface is defined and discussed in NDOSolver.h. Furthermore,
 a few implementation-dependent details (and compile-time switches) which are worth knowing 
 are described in this file.
   
- \b /MMCF/extlib/NDO/SubGrad/SubGrad.C: Contains the implementation of the \c SubGrad
 class. You should not need to read it.
  
- \b /MMCF/extlib/NDO/SubGrad/Stepsize.h:  Contains the declaration of class \c Stepsize. It is an
 <em>abstract class</em> with <em>pure virtual methods</em>, so that you
 cannot declare objects of type \c Stepsize, but only of derived classes of
 its. The class defines the interface for SR and interacts with SubGrad class. 
 
- \b /MMCF/extlib/NDO/SubGrad/Stepsize/ColorTV.h: Contains the declaration of class \c ColorTV, 
 implementing a target value SR, the one introduced in the original  <em>Volume algorithm </em>.
 Furthermore, a few implementation-dependent details (and compile-time switches) 
 which are worth knowing are described in this file.

- \b  /MMCF/extlib/NDO/SubGrad/Stepsize/FumeroTV.h: Contains the declaration of class \c FumeroTV, 
 implementing a target value SR. The rule is characterized  by an exponential formula for the level.
 Furthermore, a few implementation-dependent details (and compile-time switches) 
 which are worth knowing are described in this file.
 
- \b  /MMCF/extlib/NDO/SubGrad/Stepsize/Polyak.h: Contains the declaration of class \c Polyak, 
 implementing  the prior target value SR. Both the level and the factor beta are constants. Furthermore, 
 a few implementation-dependent details (and compile-time switches) which are 
 worth knowing are described in this file.

- \b /MMCF/extlib/NDO/SubGrad/Deflection.h:  Contains the declaration of class \c Deflection. It is an
 <em>abstract class</em> with <em>pure virtual methods</em>, so that you
 cannot declare objects of type \c Deflection, but only of derived classes of
 its. The class defines the interface for DR and interacts with SubGrad class. 
 
-  \b  /MMCF/extlib/NDO/SubGrad/Deflection/PrimalDual.h: Contains the declaration of class \c PrimalDual, 
 implementing two variants of the <em>Primal-Dual subgradient method</em> (PDSM): the simple and 
 weighted averages. Furthermore, a few implementation-dependent details 
 (and compile-time switches) which are worth knowing are described in this file.
  
-  \b /MMCF/extlib/NDO/SubGrad/Deflection/Volume.h: Contains the declaration of class \c Volume, 
 implementing a revisited <em>Volume algorithm</em>. The deflection coefficient is found solving a simple 
 quadratic problem.  Furthermore, a few implementation-dependent details 
 (and compile-time switches) which are worth knowing are described in this file.
   
- \b MMCF/extlib/NDO/NDOSlver/FiOracle.h:  Contains the declaration of class \c FiOracle. It is an
 <em>abstract class</em> with <em>pure virtual methods</em>, so that you
 cannot declare objects of type \c FiOracle, but only of derived classes of
 its. FiOracle sets the interface for the "black boxes" (oracles) computing 
 [approximate] function values and [epsilon-]subgradients for (constrained) NonDifferentiable
 Optimization algorithms.
  
- \b /MMCF/FiLagrRelax/Flow/FlwFiOrcl.h:  Contains the declaration of class \c FlwFiOrcl, implementing the 
 <em> Flow relaxation </em> (FR) of the FC-MCND problem. This class derives 
 from \c FiOracle, hence most of its interface is defined and discussed in FiOracle.h. 
 Furthermore, a few implementation-dependent details (and compile-time switches) 
 which are worth knowing are described in this file.

- \b /MMCF/FiLagrRelax/Flow/FlwFiOrcl.C: Contains the implementation of the \c FlwFiOrcl class. 
 You should not need to read it.

- \b /MMCF/FiLagrRelax/Knapsack/KnpsFiOrcl.h:  Contains the declaration of class \c KnpsFiOrcl, implementing the 
 <em> Knapsack relaxation </em> (FR) of the FC-MCND problem. This class derives 
 from \c FiOracle, hence most of its interface is defined and discussed in FiOracle.h. 
 Furthermore, a few implementation-dependent details (and compile-time switches) 
 which are worth knowing are described in this file.

- \b /MMCF/FiLagrRelax/Knapsack/KnpsFiOrcl.C: Contains the implementation of the \c KnpsFiOrcl class. 
You should not need to read it.
     
- \b Main/Main.C: Contains an example of use of the provided SubGrad solver
 to solve the Lagrangian dual of the FC-MCND problem based on both the FR 
 and the KR.  The user need to set a few macro in order to manage the variant of 
 SM to be used. <b> The project uses  other classes like the solvers for the Min Cost Flow (MCF)
 problem. You should not need to read it. </b>

\section install Software Installation

 -#  The configure script is responsible to build the software on your specific system. You can adjust some configuration options by editing \a build/config before the execution thereof. This script needs the full paths to the CPLEX libraries and the include directory, respectively,  setting  --with-cplex-libdir= \< <em>path to CPLEX include dir </em> \>  and --with-cplex-includedir= \< <em>path to CPLEX library</em> \> ;
 
 -# To build the CNDSM project, type \a make in the \a build directory;
 
 -# To clean the project,  type <em> make clean </em> and <em> make distclean </em>  in the \a build directory.



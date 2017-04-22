##############################################################################
################################# makefile ###################################
##############################################################################
#                                                                            #
#   makefile of libNDO                                                       #
#                                                                            #
#   all modules are given in input                                           #
#                                                                            #
#   Input:  $(CC)        = compiler command		                     #
#           $(SW)        = compiler options                                  #
#           $(OPTH)      = the OPT****.h files to include	             #
#           $(OPTINC)    = the -I$(OPT****.h directory)		             #
#           $(libOSIINC) = the -I$(include directories) for libOSI           #
#           $(libNDODIR) = the basic directory for the library	             #
#                                                                            #
#   Output: $(NDOH)      = the include files for for NDOSolver/FiOracle      #
#           $(FIOH)      = the include files for FiOracle (only)             #
#           $(NDOINC)    = the -I$(NDOSolver/FiOracle directory)             #
#           $(libNDOOBJ) = the library (that must be built)	             #
#           $(libNDOLIB) = the library + external libreries + -L<libdirs>    #
#           $(libNDOH)   = all the .h files to include		             #
#           $(libNDOINC) = all the -I$(include directories)                  #
#                                                                            #
#                                VERSION 2.01	                             #
#                                                                            #
#                	        28 - 02 - 2014	                             #
#                                                                            #
#			        Antonio Frangioni		             #
#   			     Operations Research Group		             #
#                           Dipartimento di Informatica		             #
#   	                        Universita' di Pisa	                     #
#                                                                            #
##############################################################################

# module name
libNDOOBJ = $(libNDODIR)lib/libNDO.a

# define & include the necessary modules- - - - - - - - - - - - - - - - - - -
# if a module is not used in the current configuration, just comment out the
# corresponding include line
# each module outputs some macros to be used here:
# *OBJ is the final object/library
# *LIB external libreries + -L<libdirs>
# *H   is the list of all include files
# *INC is the -I<include directories>

# NDOSlver
NDODIR = $(libNDODIR)NDOSlver/
include $(NDODIR)makefile

# SubGradient and its all stuff (volume and stpsize)- - - - - - - - - - - - -

SGRDIR = $(libNDODIR)SubGrad/
include $(SGRDIR)makefile

VOLSDIR = $(SGRDIR)Deflection/
include $(VOLSDIR)makefile 

STPSDIR = $(SGRDIR)Stepsize/
include $(STPSDIR)makefile 

# main module (archiving phase) - - - - - - - - - - - - - - - - - - - - - - -

# object files (this is local)
lcllibNDOOBJ =	$(SGROBJ)

# libreries
libNDOLIB =	$(SGRLIB)

# includes
libNDOH   =	$(NDOH)\
		$(SGRH)\
		$(VOLSH)\
		$(STPSH)

# include dirs
libNDOINC =	$(NDOINC)\
		$(SGRINC)\
		$(VOLSINC)\
		$(STPSINC)

$(libNDOOBJ): $(lcllibNDOOBJ)
	ar -ruv $(libNDOOBJ) $(lcllibNDOOBJ)

# clean - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

clean::
	rm -f $(libNDOOBJ)

############################ End of makefile #################################

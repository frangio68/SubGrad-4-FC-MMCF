##############################################################################
################################# makefile ###################################
##############################################################################
#                                                                            #
#   makefile of libCQK - do not build libCQK.a, just $(libCQKOBJ)            #
#   all modules are given in input                                           #
#                                                                            #
#   Input:  $(CC)        = compiler command                                  #
#           $(SW)        = compiler options                                  #
#           $(libCQKDIR) = the basic directory for the library               #
#                                                                            #
#   Output: $(libCQKOBJ) = the library                                       #
#           $(libCQKLIB) = external libreries + -L<libdirs>                  #
#           $(libCQKH)   = all the .h files of the CQK library               #
#           $(libCQKINC) = all the -I$(include dirs) of the CQK library	     #
#                                                                            #
#                                VERSION 1.00	                             #
#                                17 - 05 - 2012                              #
#                                                                            #
#                              Antonio Frangioni                             #
#                            Operations Research Group                       #
#                           Dipartimento di Informatica	                     #
#                            Universita' di Pisa                             #
#                                                                            #
##############################################################################

# define & include the necessary modules- - - - - - - - - - - - - - - - - - -
# if a module is not used in the current configuration, just comment out the
# corresponding include line
# each module outputs some macros to be used here:
# *OBJ is the final object/library
# *LIB external libreries + -L<libdirs>
# *H   is the list of all include files
# *INC is the -I<include directories>

# CQKnPClass
CQKDIR = $(libCQKDIR)CQKnPClass/
include $(CQKDIR)makefile

# DualCQKnP
DQKDIR = $(libCQKDIR)DualCQKnP/
include $(DQKDIR)makefile

# ExDualCQKnP
EQKDIR = $(libCQKDIR)ExDualCQKnP/
include $(EQKDIR)makefile

# object files- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
libCQKOBJ =	$(CPKOBJ)\
		$(DQKOBJ)\
		$(EQKOBJ)

# libreries - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
libCQKLIB =	$(CPKLIB)\
		$(DQKLIB)\
		$(EQKLIB)

# includes- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
libCQKH =	$(CQKH)\
		$(DQKH)\
		$(EQKH)

# include dirs- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
libCQKINC =	$(CQKINC)\
		$(DQKINC)\
		$(EQKINC)

# clean - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

clean::
	rm -f $(libCQKOBJ)

############################ End of makefile #################################

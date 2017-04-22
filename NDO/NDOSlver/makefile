##############################################################################
################################ makefile ####################################
##############################################################################
#									     #
#   makefile of NDOSlver						     #
#									     #
#   Input:  $(NDODIR) = the directory where the NDOSolver/FiOracle source is #
#									     #
#   Output: $(NDOH)   = the .h files to include for NDOSolver/FiOracle	     #
#           $(FIOH)   = the .h files to include for FiOracle (only)	     #
#           $(NDOINC) = the -I$(NDOSolver/FiOracle directory)		     #
#									     #
#                                VERSION 3.00				     #
#                	        09 - 05 - 2012				     #
#									     #
#			      Antonio Frangioni				     #
#   			   Operations Research Group			     #
#			  Dipartimento di Informatica			     #
#   			     Universita' di Pisa			     #
#									     #
##############################################################################

# clean - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

clean::
	rm -f $(NDODIR)*~

# macroes to be exported- - - - - - - - - - - - - - - - - - - - - - - - - - -

FIOH   = $(NDODIR)FiOracle.h
NDOH   = $(NDODIR)NDOSlver.h $(FIOH)
NDOINC = -I$(NDODIR)

# dependencies: every .o from its .C + every recursively included .h- - - - -

########################## End of makefile ###################################

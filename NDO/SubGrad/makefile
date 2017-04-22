##############################################################################
################################ makefile ####################################
##############################################################################
#                                                                            #
#   makefile of SubGrad                                                      #
#                                                                            #
#   Input:  $(CC)     = compiler command		                     #
#           $(SW)     = compiler options                                     #
#           $(OPTH)   = the OPT****.h files to include		             #
#           $(OPTINC) = the -I$(OPT****.h directory)		             #
#   	    $(NDOH)   = the .h files to include for SubGrad		     #
#           $(NDOINC) = the -I$(NDOSolver/FiOracle directories)	             #
#           $(SGRDIR) = the directory where the SubGrad source is            #
#                                                                            #
#   Output: $(SGROBJ) = the final object(s) / library		             #
#           $(SGRLIB) = external libreries + -L<libdirs>	             #
#           $(SGRH)   = the .h files to include for SubGrad                  #
#           $(SGRINC) = the -I$(Bundle/MPSolver directory)	             #
#                                                                            #
#                                VERSION 2.00                                #
#                	        28 - 02 - 2014 		                     #
#                                                                            #
#                              Antonio Frangioni                             #
#                           Operations Research Group                        #
#                          Dipartimento di Informatica                       #
#                              Universita' di Pisa                           #
#                                                                            #
##############################################################################

# clean - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

clean::
	rm -f $(SGRDIR)*.o $(SGRDIR)*~

# define & include the necessary modules- - - - - - - - - - - - - - - - - - -
# if a module is not used in the current configuration, just comment out the
# corresponding include line
# each module outputs some macros to be used here:
# *OBJ is the final object/library
# *LIB external libreries + -L<libdirs>
# *H   is the list of all include files
# *INC is the -I<include directories>

# the whole libCQK, with all required external libreries
libCQKDIR = $(SGRDIR)../../CQKnP/
include $(libCQKDIR)lib/makefile

# macroes to be exported- - - - - - - - - - - - - - - - - - - - - - - - - - -

SGROBJ = $(SGRDIR)SubGrad.o $(libCQKOBJ)
SGRLIB = $(libCQKLIB)
SGRINC = -I$(SGRDIR) $(libCQKINC)
SGRH   = $(SGRDIR)SubGrad.h $(SGRDIR)Stepsize.h $(SGRDIR)Deflection.h \
	 $(libCQKH)

# dependencies: every .o from its .C + every recursively included .h- - - - -

$(SGRDIR)SubGrad.o: $(SGRDIR)SubGrad.C $(SGRH) $(NDOH) $(OPTH)
	$(CC) -c $*.C -o $@ $(NDOINC) $(CQKINC) $(OPTINC) $(SW)

########################## End of makefile ###################################

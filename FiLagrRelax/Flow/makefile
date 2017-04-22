##############################################################################
################################# makefile ###################################
##############################################################################
#                                                                            #
#   makefile of FlwFiOrcl                                                    #
#   all modules are given in input                                           #
#                                                                            #
#   Input:  $(CC)        = compiler command                                  #
#           $(SW)        = compiler options                                  #
#           $(OPTH)      = the OPT****.h files to include		     #
#           $(OPTINC)    = the -I$(OPT****.h directory)			     #
#           $(FIOH)      = the include files for FiOracle                    #
#           $(NDOINC)    = the -I$(NDOSolver/FiOracle directories)           #
#           $(libMCFH)   = all the .h files of the MCF library		     #
#           $(libMCFINC) = all the -I$(include dirs) of the MCF library	     #
#           $(MMFH)      = the MMCFClas.h file(s) to include                 #
#           $(MMFINC)    = the -I$(MMCFClas directory)			     #
#           $(MMGH)      = the Graph.h file(s) to include                    #
#           $(MMGINC)    = the -I$(Graph directory)                          #
#           $(MFBH)      = the include files for MMCFFlwB                    #
#           $(MFBINC)    = the -I$(MMCFFlwB directories)                     #
#           $(FLWDIR)    = the directory where the source is                 #
#                                                                            #
#   Output: $(FLWOBJ)    = the final object file(s)                          #
#           $(FLWLIB)    = external libreries + -L<libdirs>                  #
#           $(FLWBH)     = the .h files to include                           #
#           $(FLWINC)    = the -I$(FlowFiOracle directory)                   #
#                                                                            #
#                                VERSION 3.00                                #
#                               10 - 05 - 2012                               #
#                                                                            #
#                              Antonio Frangioni                             #
#                           Operations Research Group                        #
#                          Dipartimento di Informatica                       #
#                              Universita' di Pisa                           #
#                                                                            #
##############################################################################

# clean - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

clean::
	rm -f $(FLWDIR)*.o $(FLWDIR)*~

# macroes to be exported- - - - - - - - - - - - - - - - - - - - - - - - - - -

FLWOBJ = $(FLWDIR)FlwFiOrcl.o
FLWLIB = 
FLWH   = $(FLWDIR)FlwFiOrcl.h
FLWINC = -I$(FLWDIR)

# dependencies: every .o from its .C + every recursively included .h- - - - -

$(FLWOBJ): $(FLWDIR)FlwFiOrcl.C $(FLWH) $(MFBH) $(MMGH) $(FIOH) $(libMCFH) \
	$(OPTH)
	$(CC) -c $*.C -o $@ $(MMFINC) $(MMGINC) $(MFBINC) $(NDOINC) \
	$(libMCFINC) $(OPTINC) $(SW)

############################ End of makefile #################################

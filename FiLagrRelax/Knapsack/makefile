##############################################################################
################################# makefile ###################################
##############################################################################
#                                                                            #
#   makefile of KnpsFiOrcl                                                   #
#   all modules are given in input                                           #
#                                                                            #
#   Input:  $(CC)     = compiler command                                     #
#           $(SW)     = compiler options                                     #
#           $(OPTH)   = the OPT****.h files to include			     #
#           $(OPTINC) = the -I$(OPT****.h directory)			     #
#           $(MMGH)   = the Graph.h file(s) to include                       #
#           $(MMGINC) = the -I$(Graph directory)                             #
#           $(FIOH)   = the include files for FiOracle                       #
#           $(NDOINC) = the -I$(NDOSolver/FiOracle directories)              #
#           $(KFLDIR) = the directory where the source is                    #
#                                                                            #
#   Output: $(KFLOBJ) = the final object file(s)                             #
#           $(KFLLIB) = external libreries + -L<libdirs>                     #
#           $(KFLBH)  = the .h files to include                              #
#           $(KFLINC) = the -I$(FlowFiOracle directory)                      #
#                                                                            #
#                                VERSION 3.00                                #
#                               16 - 05 - 2012                               #
#                                                                            #
#                              Antonio Frangioni                             #
#                           Operations Research Group                        #
#                          Dipartimento di Informatica                       #
#                              Universita' di Pisa                           #
#                                                                            #
##############################################################################

# clean - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

clean::
	rm -f $(KFLDIR)*.o $(KFLDIR)*~

# macroes to be exported- - - - - - - - - - - - - - - - - - - - - - - - - - -

KFLOBJ = $(KFLDIR)KnpsFiOrcl.o
KFLLIB = 
KFLH   = $(KFLDIR)KnpsFiOrcl.h
KFLINC = -I$(KFLDIR) 

# dependencies: every .o from its .C + every recursively included .h- - - - -

$(KFLOBJ): $(KFLDIR)KnpsFiOrcl.C $(KFLH) $(MMFH) $(FIOH) $(MMGH) $(OPTH)
	$(CC) -c $*.C -o $@ $(MMGINC) $(NDOINC) $(MMFINC) $(OPTINC) $(SW)

############################ End of makefile #################################

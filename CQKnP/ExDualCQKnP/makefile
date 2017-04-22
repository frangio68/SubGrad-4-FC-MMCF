##############################################################################
#                                                                            #
#   makefile of ExDualCQKnP                                                  #
#   all modules are given in input                                           #
#                                                                            #
#   Input:  $(CC)     = compiler command                                     #
#           $(SW)     = compiler options                                     #
#           $(CQKH)   = the include file(s) for CQKnPClass.h                 #
#           $(CQKINC) = the -I$(CQKnPClass directory)                        #
#           $(DQKH)   = the .h files to include for DualCQKnP.h              #
#           $(DQKINC) = the -I$(DualCQKnP directory)                         #
#           $(EQKDIR) = the directory where the source is                    #
#                                                                            #
#   Output: $(EQKLIB) = the final object file(s)                             #
#           $(EQKH)   = the .h files to include                              #
#           $(EQKINC) = the -I$(DualCQKnP directory)                         #
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

# clean - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

clean::
	rm -f $(EQKDIR)*.o $(EQKDIR)*~

# macroes to be exported- - - - - - - - - - - - - - - - - - - - - - - - - - -

EQKOBJ = $(EQKDIR)ExDualCQKnP.o
EQKLIB = 
EQKH   = $(EQKDIR)ExDualCQKnP.h
EQKINC = -I$(EQKDIR)

# dependencies: every .o from its .C + every recursively included .h- - - - -

$(EQKOBJ): $(EQKDIR)ExDualCQKnP.C $(EQKH) $(DQKH) $(CQKH)
	$(CC) -c $*.C -o $@ $(CQKINC) $(DQKINC) $(SW)

############################ End of makefile #################################
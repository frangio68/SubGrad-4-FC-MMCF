
##############################################################################
#                                                                            #
#   makefile of DualCQKnP                                                    #
#   all modules are given in input                                           #
#                                                                            #
#   Input:  $(CC)     = compiler command                                     #
#           $(SW)     = compiler options                                     #
#           $(CQKH)   = the include file(s) for CQKnPClass.h                 #
#           $(CQKINC) = the -I$(CQKnPClass directory)                        #
#           $(DQKDIR) = the directory where the source is                    #
#                                                                            #
#   Output: $(DQKLIB) = the final object file(s)                             #
#           $(DQKH)   = the .h files to include                              #
#           $(DQKINC) = the -I$(DualCQKnP directory)                         #
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
	rm -f $(DQKDIR)*.o $(DQKDIR)*~

# macroes to be exported- - - - - - - - - - - - - - - - - - - - - - - - - - -

DQKOBJ = $(DQKDIR)DualCQKnP.o
DQKLIB = 
DQKH   = $(DQKDIR)DualCQKnP.h
DQKINC = -I$(DQKDIR)

# dependencies: every .o from its .C + every recursively included .h- - - - -

$(DQKOBJ): $(DQKDIR)DualCQKnP.C $(DQKH) $(CQKH)
	$(CC) -c $*.C -o $@ $(CQKINC) $(SW)

############################ End of makefile #################################
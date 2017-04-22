##############################################################################
################################# makefile ###################################
##############################################################################
#                                                                            #
#   makefile of Graph                                                        #
#                                                                            #
#   Input:  $(CC)     = compiler command				     #
#           $(SW)     = compiler options                                     #
#           $(MMGDIR) = the directory where the source is                    #
#									     #
#   Output: $(MMGOBJ) = the final object(s) / library			     #
#           $(MMGLIB) = external libreries + -L<libdirs>		     #
#           $(MMGH)   = the .h files to include	(apart from the OPT* ones)   #
#           $(MMGINC) = the -I$(include directories)			     #
#                                                                            #
#                                VERSION 3.00				     #
#                	        07 - 05 - 2012				     #
#									     #
#			      Antonio Frangioni				     #
#   			   Operations Research Group			     #
#			  Dipartimento di Informatica			     #
#   			     Universita' di Pisa			     #
#									     #
##############################################################################

# clean - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

clean::
	rm -f $(MMGDIR)*.o $(MMGDIR)*~

# macroes to be exported- - - - - - - - - - - - - - - - - - - - - - - - - - -

MMGOBJ = $(MMGDIR)Graph.o
MMGLIB =
MMGH   = $(MMGDIR)Graph.h
MMGINC = -I$(MMGDIR)

# dependencies: every .o from its .C + every recursively included .h- - - - -

$(MMGDIR)Graph.o: $(MMGDIR)Graph.C $(MMGH)
	$(CC) -c $*.C -o $@ $(SW)  

############################ End of makefile #################################

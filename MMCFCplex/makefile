##############################################################################
#									     #
#   makefile of MMCFCple       						     #
#   all modules are given in input                                           #
#									     #
#   Input:  $(CC)        = compiler command				     #
#           $(SW)        = compiler options                                  #
#           $(OPTH)      = the OPT****.h files to include		     #
#           $(OPTINC)    = the -I$(OPT****.h directory)			     #
#           $(MMCFH)     = the include files for MMCFClass		     #
#           $(MMFINC)    = the -I$(MMCFClass.h directory)		     #
#           $(MMGH)      = the Graph.h file(s) to include                    #
#           $(MMGINC)    = the -I$(Graph directory)			     #
#           $(libCPXINC) = the -I$(include directories)	for libCPX	     #
#           $(MMXDIR)    = the directory where the source is                 #
#									     #
#   Output: $(MMXLIB)    = the final object file(s)			     #
#           $(MMXH)      = the .h files to include			     #
#           $(MMXINC)    = the -I$(MMCFCple directory)			     #
#									     #
#                                VERSION 4.00				     #
#                	        12 - 12 - 2009				     #
#									     #
#			      Antonio Frangioni				     #
#   			   Operations Research Group			     #
#			  Dipartimento di Informatica			     #
#   			     Universita' di Pisa			     #
#									     #
##############################################################################

# clean - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

clean::
	rm -f $(MMXDIR)*.o $(MMXDIR)*~

# macroes to be exported- - - - - - - - - - - - - - - - - - - - - - - - - - -

MMXOBJ = $(MMXDIR)MMCFCple.o
MMXLIB = 
MMXH   = $(MMXDIR)MMCFCple.h
MMXINC = -I$(MMXDIR)

# dependencies: every .o from its .C + every recursively included .h- - - - -

$(MMXDIR)MMCFCple.o: $(MMXDIR)MMCFCple.C $(MMXH) $(MMGH) $(MMCFH) $(OPTH)
	$(CC) -c $*.C -o $@ $(SW) $(MMGINC) $(MMFINC) $(libCPXINC) $(OPTINC)

############################ End of makefile #################################

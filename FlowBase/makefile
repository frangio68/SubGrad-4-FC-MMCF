##############################################################################
################################# makefile ###################################
##############################################################################
#                                                                            #
#   makefile of MMCFFlwB                                                     #
#   all modules are given in input                                           #
#                                                                            #
#   Input:  $(CC)        = compiler command                                  #
#           $(SW)        = compiler options                                  #
#           $(OPTH)      = the OPT****.h files to include		     #
#           $(OPTINC)    = the -I$(OPT****.h directory)			     #
#           $(MMFH)      = the MMCFClas.h file(s) to include                 #
#           $(MMFINC)    = the -I$(MMCFClas directory)			     #
#           $(MMGH)      = the Graph.h file(s) to include                    #
#           $(MMGINC)    = the -I$(Graph directory)                          #
#           $(MFBDIR)    = the directory where the source is                 #
#                                                                            #
#   Output: $(MFBOBJ) = the final object file(s)                             #
#           $(MFBLIB) = external libreries + -L<libdirs>                     #
#           $(MFBH)   = the .h files to include                              #
#           $(MFBINC) = the -I$(include directories)                         #
#                                                                            #
#                                VERSION 3.00                                #
#                               09 - 05 - 2012                               #
#                                                                            #
#                              Antonio Frangioni                             #
#                           Operations Research Group                        #
#                          Dipartimento di Informatica                       #
#                              Universita' di Pisa                           #
#                                                                            #
##############################################################################

# clean - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

clean::
	rm -f $(MFBDIR)*.o $(MFBDIR)*~

# define & include the necessary modules- - - - - - - - - - - - - - - - - - -
# if a module is not used in the current configuration, just comment out the
# corresponding include line
# each module outputs some macros to be used here:
# *LIB external libreries + -L<libdirs>
# *OBJ the final object(s) / library
# *H   is the list of all include files
# *INC is the -I<include directories>

libMCFDIR = $(MFBDIR)../MCF/
include $(libMCFDIR)lib/makefile-c

# macroes to be exported- - - - - - - - - - - - - - - - - - - - - - - - - - -

MFBOBJ = $(MFBDIR)MMCFFlwB.o $(libMCFOBJ)
MFBLIB = $(libMCFLIB)
MFBH   = $(MFBDIR)MMCFFlwB.h $(libMCFH)
MFBINC = -I$(MFBDIR) $(libMCFINC)

# dependencies: every .o from its .C + every recursively included .h- - - - -

$(MFBDIR)MMCFFlwB.o: $(MFBDIR)MMCFFlwB.C $(MFBH) $(MMGH) $(MMFH) \
	$(libMCFH) $(OPTH)
	$(CC) -c $*.C -o $@ $(MMGINC) $(MFBINC) $(MMFINC) $(OPTINC) \
	$(libMCFINC) $(SW)

############################ End of makefile #################################

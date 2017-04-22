##############################################################################
#									     #
#   main makefile of CpxSolver						     #
#									     #
##############################################################################

# module
NAME = CpxSolver

# basic directory
DIR = ./

# debug - production switches
SW = -O3

# libreries
LIB = -lm

# compiler
CC = g++

# default target- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

all default: $(NAME)

# clean - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

clean::
	rm -f $(DIR)*.o $(DIR)*~ $(NAME)

# distclean - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

distclean: 
	-rm -f makefile

# define & include the necessary modules- - - - - - - - - - - - - - - - - - -
# if a module is not used in the current configuration, just comment out the
# corresponding include line
# each module outputs some macros to be used here:
# *LIB external libreries + -L<libdirs>
# *OBJ the final object(s) / library
# *H   is the list of all include files
# *INC is the -I<include directories>

# OPT****.h
OPTH	= $(DIR)OPTUtils.h $(DIR)OPTvect.h $(DIR)OPTtypes.h
OPTINC	= -I$(DIR)

# Graph
MMGDIR = $(DIR)../Graph/
include $(MMGDIR)makefile

# MMCFClass
MMFDIR =  $(DIR)../MMCFClas/
include $(MMFDIR)makefile

# Cplex
MMXDIR = ../MMCFCplex/
include $(MMXDIR)makefile-c

# main module (linking phase) - - - - - - - - - - - - - - - - - - - - - - - -

# final object files
MOBJ =	$(MMGOBJ)\
	$(MMXOBJ)

# final libreries + -L<libdirs>
MLIB =	$(MMGLIB)\
	$(MMXLIB)

# creating the module

$(NAME): $(MOBJ) $(DIR)Main.o
	$(CC) -o $(NAME) $(DIR)Main.o $(MOBJ) $(MLIB) $(LIB) $(SW)  

# dependencies: every .o from its .C + every recursively included .h- - - - -

# final includes
MH =    $(MMFH)\
        $(MMGH)\
	$(MMXH)\
	$(OPTH)

# final include paths
MINC =  $(MMFINC)\
        $(MMGINC)\
	$(MMXINC)\
	$(OPTINC)

$(DIR)Main.o: $(DIR)Main.C $(MH)
	$(CC) -c $*.C -o $@ $(MINC) $(SW)

############################ End of makefile #################################

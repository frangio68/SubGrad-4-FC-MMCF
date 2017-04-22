##############################################################################
#									     #
#   main makefile of CNDSM						     #
#									     #
##############################################################################

# module
NAME = SgSolver

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

# the whole libNDO, with all required external libreries
libNDODIR = $(DIR)../NDO/
include $(libNDODIR)lib/makefile

# Flow
FLWDIR = $(DIR)../FiLagrRelax/Flow/
include $(FLWDIR)makefile-c

# knapsack
KFLDIR = $(DIR)../FiLagrRelax/Knapsack/
include $(KFLDIR)makefile

# main module (linking phase) - - - - - - - - - - - - - - - - - - - - - - - -

# object files (this is local)
MOBJ =	$(MMGOBJ)\
	$(libNDOOBJ)\
	$(FLWOBJ)\
	$(KFLOBJ)

# libreries
MLIB =  $(MMGLIB)\
	$(libNDOLIB)\
	$(FLWLIB)\
	$(KFLLIB)

# includes
MH =    $(MMFH)\
        $(MMGH)\
	$(libNDOH)\
	$(FLWH)\
	$(KFLH)

# include dirs
MINC =  $(MMFINC)\
        $(MMGINC)\
	$(libNDOINC)\
	$(FLWINC)\
	$(KFLINC)

# main module (linking phase) - - - - - - - - - - - - - - - - - - - - - - - -

$(NAME): $(MOBJ) $(DIR)Main.o
	$(CC) -o $(NAME) $(DIR)Main.o $(MOBJ) $(MLIB) $(LIB) $(SW) 

# dependencies: every .o from its .C + every recursively included .h- - - - -

$(DIR)Main.o: $(DIR)Main.C $(MH) 
	$(CC) -c $*.C -o $@ $(MINC) $(SW)	

############################ End of makefile #################################

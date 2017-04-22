##############################################################################
################################# makefile ###################################
##############################################################################
#                                                                            #
#   makefile of Stepsize algorithms                                          #
#                                                                            #
#   Input:  $(CC)     = compiler command                                     #
#           $(SW)     = compiler options                                     #
#           $(STPSDIR) = the directory where the step algorithms source are  #
#           $(SGRINC)   = the -I$(Stepsize.h directory)	                     #
#                                                                            #
#   Output: $(STPSOBJ) = the final object(s) / library                       #
#           $(STPSLIB) = external libreries + -L<libdirs>                    #
#           $(STPSH)   = the .h files to include                             #
#           $(STPSINC) = the -I$(Step directory)                             #
#                                                                            #
#                                VERSION 1.00                                #
#                               28 - 02 - 2014                               #
#                                                                            #
#                              Antonio Frangioni                             #
#                               Enrico Gorgone                               #
#                           Operations Research Group                        #
#                          Dipartimento di Informatica                       #
#                             Universita' di Pisa                            #
#                                                                            #
##############################################################################

# define & include the necessary modules- - - - - - - - - - - - - - - - - - -

# Fumero version
FUMH = $(STPSDIR)FumeroTV.h

# Polyak version
PLYH = $(STPSDIR)Polyak.h

# Diminishing version
DIMH = $(STPSDIR)ColorTV.h

# macroes to be exported- - - - - - - - - - - - - - - - - - - - - - - - - - -

# final include paths
STPSINC = -I$(STPSDIR)\
         $(SGRINC)

# final includes
STPSH = $(FUMH)\
       $(PLYH)\
       $(DIMH)

########################### End of makefile ##################################

















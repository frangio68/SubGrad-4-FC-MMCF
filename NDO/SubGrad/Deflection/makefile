##############################################################################
################################# makefile ###################################
##############################################################################
#                                                                            #
#   makefile of Deflection algorithms                                        #
#                                                                            #
#   Input:  $(CC)     = compiler command                                     #
#           $(SW)     = compiler options                                     #
#           $(VOLSDIR) = the directory where VolumKND algorithms source are  #
#           $(SGRINC)   = the -I$(Deflection.h directory)                    #
#                                                                            #
#   Output: $(VOLSOBJ) = the final object(s) / library                       #
#           $(VOLSLIB) = external libreries + -L<libdirs>                    #
#           $(VOLSH)   = the .h files to include                             #
#           $(VOLSINC) = the -I$(Deflection algorithms directory)            #
#                                                                            #
#                                VERSION 1.00                                #
#                               28 - 12 - 2014                               #
#                                                                            #
#                              Antonio Frangioni                             #
#                               Enrico Gorgone                               #
#                           Operations Research Group                        #
#                          Dipartimento di Informatica                       #
#                             Universita' di Pisa                            #
#                                                                            #
##############################################################################

# define & include the necessary modules- - - - - - - - - - - - - - - - - - -

# StandardVolume algorithm
STVH = $(VOLSDIR)Volume.h

# Primal-Dual simple and weighted averages
PDSWH = $(VOLSDIR)PrimalDual.h

# macroes to be exported- - - - - - - - - - - - - - - - - - - - - - - - - - -

# final include paths
VOLSINC = -I$(VOLSDIR)\
         $(SGRINC)\
         $(OPTINC)

# final includes
VOLSH = $(STVH)\
       $(PDSWH)\
       $(OPTH)

########################### End of makefile ##################################

















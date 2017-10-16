SHELL = /bin/bash

# 'make'        build executable file
# 'make clean'  removes all .o and executable files

NO_COLOR=\x1b[0m
OK_COLOR=\x1b[32;01m
ERROR_COLOR=\x1b[31;01m
WARN_COLOR=\x1b[33;01m

OK_STRING=$(OK_COLOR)[OK]$(NO_COLOR)
ERROR_STRING=$(ERROR_COLOR)[ERRORS]$(NO_COLOR)
WARN_STRING=$(WARN_COLOR)[WARNINGS]$(NO_COLOR)

AWK_CMD = awk '{ printf "%-30s %-10s\n",$$1, $$2; }'
PRINT_ERROR = printf "$@ $(ERROR_STRING)\n" | $(AWK_CMD) && printf "$(CMD)\n$$LOG\n" && false
PRINT_WARNING = printf "$@ $(WARN_STRING)\n" | $(AWK_CMD) && printf "$(CMD)\n$$LOG\n"
PRINT_OK = printf "$@ $(OK_STRING)\n" | $(AWK_CMD)
BUILD_CMD = LOG=$$($(CMD) 2>&1) ; if [ $$? -eq 1 ]; then $(PRINT_ERROR); elif [ "$$LOG" != "" ] ; then $(PRINT_WARNING); else $(PRINT_OK); fi;

override CXXFLAGS += --std=c++11 -Wno-multichar
ifeq ($(shell echo "int main(){}" | $(CXX) --stdlib=libc++ -x c - -o /dev/null >& /dev/null; echo $$?), 0)
	override CXXFLAGS += --stdlib=libc++
endif

# Apple clang 5.0 requires -fcolor-diagnostics
# GCC 4.9 requires -fdiagnostics-color
ifeq ("$(shell echo "int main(){}" | $(CXX) -fdiagnostics-color -x c - -o /dev/null 2>&1)", "")
	override CXXFLAGS += -fdiagnostics-color
else ifeq ("$(shell echo "int main(){}" | $(CXX) -fcolor-diagnostics -x c - -o /dev/null 2>&1)", "")
	override CXXFLAGS += -fcolor-diagnostics
endif


# define the C compiler to use
CC = g++

MV = mv

CRFLAGS = `root-config --cflags`

GC = g++ -g -std=c++11 -w -fPIC

# define any compile-time flags
CFLAGS = -w -g -std=c++11

# define any directories containing header files other than /usr/include
INCLUDES = -I/usr/include/eigen3/ -I ~/TRACS_Concurrency_master/include/  -I/usr/local/root/include/ -I/usr/include/qt4/QtCore/ -I/usr/include/qt4/ -I/usr/include/qt4/QtGui/

# define library paths in addition to /usr/lib
LFLAGS = -L/usr/local/root/lib -L ~/TRACS_Concurrency_master/lib -L/usr/lib/x86_64-linux-gnu/ -L/usr/lib 

# define any libraries to link into executable:
LIBS = -ldolfin -lCore -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMultiProc -pthread -lm -ldl -rdynamic -lTreePlayer -lTreeViewer -lHistPainter -lQtCore -lTable -lFFTW -lFITSIO -lGX11TTF -lPyROOT -lMinuit2 -lMathMore -lCling -lRooFit -lRooFitCore -lMatrix -lTMeasHeader -lTMeas -lTWaveform -lboost_system

PRINT = @echo -e "\e[1;34mBuilding $<\e[0m"


# define the C source files
SDIR = src/
SRCS = $(SDIR)DoTRACSFit.cpp $(SDIR)TRACSFit.cpp $(SDIR)CarrierCollection.cpp $(SDIR)Carrier.cpp $(SDIR)CarrierMobility.cpp $(SDIR)CarrierTransport.cpp  $(SDIR)Global.cpp $(SDIR)SMSDetector.cpp $(SDIR)SMSDSubDomains.cpp $(SDIR)Threading.cpp $(SDIR)TRACSInterface.cpp $(SDIR)H1DConvolution.C $(SDIR)Utilities.cpp $(SDIR)TMeas.cpp $(SDIR)TWaveform.cpp $(DIR)TMeasHeader.cpp

ODIR = obj/
OBJ_ = DoTRACSFit.o TRACSFit.o CarrierCollection.o Carrier.o CarrierMobility.o CarrierTransport.o Global.o SMSDetector.o SMSDSubDomains.o Threading.o TRACSInterface.o H1DConvolution.o Utilities.o TMeas.o TWaveform.o TMeasHeader.o TMeasDict.o TMeasHeaderDict.o TWaveDict.o
OBJB_ = DoTracsOnly.o TRACSFit.o CarrierCollection.o Carrier.o CarrierMobility.o CarrierTransport.o Global.o SMSDetector.o SMSDSubDomains.o Threading.o TRACSInterface.o H1DConvolution.o Utilities.o TMeas.o TWaveform.o TMeasHeader.o TMeasDict.o TMeasHeaderDict.o TWaveDict.o
OBJC_ = MfgTRACSFit.o TRACSFit.o CarrierCollection.o Carrier.o CarrierMobility.o CarrierTransport.o Global.o SMSDetector.o SMSDSubDomains.o Threading.o TRACSInterface.o H1DConvolution.o Utilities.o TMeas.o TWaveform.o TMeasHeader.o TMeasDict.o TMeasHeaderDict.o TWaveDict.o
OBJEDGE_ = Edge_tree.o TMeas.o TWaveform.o TMeasHeader.o TMeasDict.o TMeasHeaderDict.o TWaveDict.o

OBJ := $(patsubst %,$(ODIR)%,$(OBJ_))
OBJB := $(patsubst %,$(ODIR)%,$(OBJB_))
OBJC := $(patsubst %,$(ODIR)%,$(OBJC_))
OBJEDGE := $(patsubst %,$(ODIR)%,$(OBJEDGE_))

all: DoTracsOnly DoTRACSFit MfgTRACSFit Edge_tree 

MAIN = DoTRACSFit

MAINB = DoTracsOnly

MAINC = MfgTRACSFit

EDGE = Edge_tree

DoTRACSFit: $(OBJ)
	@echo +++Compilation OK!
	@echo +++Linking in progress...
	@$(CC) $(CFLAGS) $(INCLUDES) -o myApp/$(MAIN) $(OBJ) $(LFLAGS) $(LIBS)
	@$(BUILD_CMD)
	@echo +++Building SUCCESSFUL!
	
DoTracsOnly: $(OBJB)
	@echo +++Compilation OK!
	@echo +++Linking in progress...
	@$(CC) $(CFLAGS) $(INCLUDES) -o myApp/$(MAINB) $(OBJB) $(LFLAGS) $(LIBS)
	@$(BUILD_CMD)
	@echo +++Building SUCCESSFUL!

MfgTRACSFit: $(OBJC)
	@echo +++Compilation OK!
	@echo +++Linking in progress...
	@$(CC) $(CFLAGS) $(INCLUDES) -o myApp/$(MAINC) $(OBJC) $(LFLAGS) $(LIBS)
	@$(BUILD_CMD)
	@echo +++Building SUCCESSFUL!
	
Edge_tree: $(OBJEDGE)
	@echo +++Compilation OK!
	@echo +++Linking in progress...
	@$(CC) $(CFLAGS) $(INCLUDES) -o myApp/$(EDGE) $(OBJEDGE) $(LFLAGS) $(LIBS)
	@$(BUILD_CMD)
	@echo +++Building SUCCESSFUL!
	

$(ODIR)DoTRACSFit.o: src/DoTRACSFit.cpp
	@$(PRINT)
	@$(CC) $(CFLAGS) $(INCLUDES) -c $(SDIR)DoTRACSFit.cpp -o $@
	@$(BUILD_CMD)
	
$(ODIR)DoTracsOnly.o: src/DoTracsOnly.cpp
	@$(PRINT)
	@$(CC) $(CFLAGS) $(INCLUDES) -c $(SDIR)DoTracsOnly.cpp -o $@
	@$(BUILD_CMD)

$(ODIR)MfgTRACSFit.o: src/MfgTRACSFit.cpp
	@$(PRINT)
	@$(CC) $(CFLAGS) $(INCLUDES) -c $(SDIR)MfgTRACSFit.cpp -o $@
	@$(BUILD_CMD)
	
$(ODIR)Edge_tree.o: src/Edge_tree.cpp
	@$(PRINT)
	@$(CC) $(CFLAGS) $(INCLUDES) -c $(SDIR)Edge_tree.cpp -o $@
	@$(BUILD_CMD)
	
$(ODIR)TRACSFit.o: $(SDIR)TRACSFit.cpp
	@$(PRINT)
	@$(CC) $(CFLAGS) $(INCLUDES) -c $(SDIR)TRACSFit.cpp -o $@ 
	@$(BUILD_CMD)

$(ODIR)CarrierCollection.o: $(SDIR)CarrierCollection.cpp
	@$(PRINT)
	@$(CC) $(CFLAGS) $(INCLUDES) -c $(SDIR)CarrierCollection.cpp -o $@
	@$(BUILD_CMD)

$(ODIR)Carrier.o: $(SDIR)Carrier.cpp
	@$(PRINT)
	@$(CC) $(CFLAGS) $(INCLUDES) -c $(SDIR)Carrier.cpp -o $@
	@$(BUILD_CMD)

$(ODIR)CarrierMobility.o: $(SDIR)CarrierMobility.cpp
	@$(PRINT)
	@$(CC) $(CFLAGS) $(INCLUDES) -c $(SDIR)CarrierMobility.cpp -o $@
	@$(BUILD_CMD)

$(ODIR)CarrierTransport.o: $(SDIR)CarrierTransport.cpp
	@$(PRINT)
	@$(CC) $(CFLAGS) $(INCLUDES) -c $(SDIR)CarrierTransport.cpp -o $@
	@$(BUILD_CMD)

$(ODIR)Global.o: $(SDIR)Global.cpp
	@$(PRINT)
	@$(CC) $(CFLAGS) $(INCLUDES) -c $(SDIR)Global.cpp -o $@
	@$(BUILD_CMD)

$(ODIR)SMSDetector.o: $(SDIR)SMSDetector.cpp
	@$(PRINT)
	@$(CC) $(CFLAGS) $(INCLUDES) -c $(SDIR)SMSDetector.cpp -o $@
	@$(BUILD_CMD)

$(ODIR)SMSDSubDomains.o: $(SDIR)SMSDSubDomains.cpp
	@$(PRINT)
	@$(CC) $(CFLAGS) $(INCLUDES) -c $(SDIR)SMSDSubDomains.cpp -o $@
	@$(BUILD_CMD)

$(ODIR)Threading.o: $(SDIR)Threading.cpp
	@$(PRINT)
	@$(CC) $(CFLAGS) $(INCLUDES) -c $(SDIR)Threading.cpp -o $@
	@$(BUILD_CMD)

$(ODIR)TRACSInterface.o: $(SDIR)TRACSInterface.cpp
	@$(PRINT)
	@$(CC) $(CFLAGS) $(INCLUDES) -c $(SDIR)TRACSInterface.cpp -o $@
	@$(BUILD_CMD)

$(ODIR)H1DConvolution.o: $(SDIR)H1DConvolution.C
	@$(PRINT)
	@$(CC) $(CFLAGS) $(INCLUDES) -c $(SDIR)H1DConvolution.C -o $@
	@$(BUILD_CMD)

$(ODIR)Utilities.o: $(SDIR)Utilities.cpp
	@$(PRINT)
	@$(CC) $(CFLAGS) $(INCLUDES) -c $(SDIR)Utilities.cpp -o $@	
	@$(BUILD_CMD)
	
$(ODIR)TMeas.o: $(SDIR)TMeas.cpp
	$(PRINT)
	$(CC) $(CFLAGS) $(INCLUDES) -c $(SDIR)TMeas.cpp -o $@
	rootcling -f $(SDIR)TMeasDict.C -c -p -I~/TRACS_Concurrency_master/include/ ~/TRACS_Concurrency_master/include/TMeas.h ~/TRACS_Concurrency_master/include/LinkDef.h
	$(BUILD_CMD)

$(ODIR)TWaveform.o: $(SDIR)TWaveform.cpp
	$(PRINT)
	$(CC) $(CFLAGS) $(INCLUDES) -c $(SDIR)TWaveform.cpp -o $@
	rootcling -f $(SDIR)TWaveDict.C -c -p -I~/TRACS_Concurrency_master/include/ ~/TRACS_Concurrency_master/include/TWaveform.h
	$(BUILD_CMD)
	
$(ODIR)TMeasHeader.o: $(SDIR)TMeasHeader.cpp
	$(PRINT)
	$(CC) $(CFLAGS) $(INCLUDES) -c $(SDIR)TMeasHeader.cpp -o $@
	rootcling -f $(SDIR)TMeasHeaderDict.C -c -p -I~/TRACS_Concurrency_master/include/ ~/TRACS_Concurrency_master/include/TMeasHeader.h
	$(BUILD_CMD)

$(ODIR)TMeasDict.o: $(SDIR)TMeasDict.C
	$(GC) $(CRFLAGS) -c $(SDIR)TMeasDict.C -o  $@


$(ODIR)TMeasHeaderDict.o: $(SDIR)TMeasHeaderDict.C
#	@echo nothing.............
	$(GC) $(CRFLAGS) -c $(SDIR)TMeasHeaderDict.C -o $@


$(ODIR)TWaveDict.o: $(SDIR)TWaveDict.C
#	@echo nothing............
	$(GC) $(CRFLAGS) -c $(SDIR)TWaveDict.C -o $@



# define the C object files 
#
# This uses Suffix Replacement within a macro:
#   $(name:string1=string2)
#         For each word in 'name' replace 'string1' with 'string2'
# Below we are replacing the suffix .c of all words in the macro SRCS
# with the .o suffix
#
#OBJS = $(SRCS:.c=.o)

# define the executable file 
#MAIN = DoTRACSFit

#
# The following part of the makefile is generic; it can be used to 
# build any executable just by changing the definitions above and by
# deleting dependencies appended to the file from 'make depend'
#

#.PHONY: depend clean

#all:    $(MAIN)
##	@echo Compilation SUCCESS!!. Executable file has been created.

#$(MAIN): $(OBJS) 
#	$(CC) $(CFLAGS) $(INCLUDES) -o $(MAIN) $(OBJS) $(LFLAGS) $(LIBS)

# this is a suffix replacement rule for building .o's from .c's
# it uses automatic variables $<: the name of the prerequisite of
# the rule(a .c file) and $@: the name of the target of the rule (a .o file) 
# (see the gnu make manual section about automatic variables)
#.c.o:
#	$(CC) $(CFLAGS) $(INCLUDES) -c "$<"  -o "$@"

clean:
	@$(RM) $(ODIR)*.o *~ myApp/$(MAIN) *~ myApp/$(MAINB) *~ myApp/$(MAINC) *~ myApp/$(EDGE)
	@$(RM) $(SDIR)TMeasDict.C $(SDIR)TMeasHeaderDict.C $(SDIR)TWaveDict.C
	@$(RM) $(SDIR)*.pcm
	@$(BUILD_CMD)

depend: $(SRCS)
	makedepend $(INCLUDES) $^


# DO NOT DELETE THIS LINE -- make depend needs it

SHELL       = /bin/sh
.SUFFIXES:
.DELETE_ON_ERROR:

CC          = gcc
CXX         = g++
FC          = gfortran

# Environment independent flags
CXXINCPATH :=
LIBPATH    :=
CXXFLAGS    = -O2 -std=c++0x $(CXXINCPATH) -I. -Wall -W -Wformat=2 -Wcast-align -Wwrite-strings -Wconversion -Wfloat-equal -Wpointer-arith
FFLAGS      = -O2 -Wall -Wextra -Wcharacter-truncation -Wunderflow -Wimplicit-interface
F95FLAGS    = -O2 -std=f95 --free-form -Wall -Wextra -Wcharacter-truncation -Wunderflow -Wimplicit-interface

FLYMAKE_CXXFLAGS = -std=c++0x -Wall -W -Wformat=2 -Wcast-align -Wwrite-strings -Wconversion -Wfloat-equal -Wpointer-arith -fsyntax-only -fmessage-length=0 -fno-diagnostics-show-option
FLYMAKE_FFLAGS   = -fsyntax-only -fmessage-length=0 -fno-diagnostics-show-option
FLYMAKE_F95FLAGS = -fsyntax-only -fmessage-length=0 -fno-diagnostics-show-option

# GCC dependent flags
ifeq ($(shell expr `gcc -dumpversion | sed -e 's/\.\([0-9]\)\b/0\1/g' -e 's/\.//g'` \>= 40900), 1)
	CXXFLAGS += -fdiagnostics-color
	FFLAGS   += -fdiagnostics-color
	F95FLAGS += -fdiagnostics-color
	FLYMAKE_CXXFLAGS += -fno-diagnostics-color -fno-diagnostics-show-caret
	FLYMAKE_FFLAGS   += -fno-diagnostics-color -fno-diagnostics-show-caret
	FLYMAKE_F95FLAGS += -fno-diagnostics-color -fno-diagnostics-show-caret
endif

# Environment dependent flags
ifneq (,$(shell find /usr/local/gcc -maxdepth 0 2>/dev/null))
	LIBPATH    := $(LIBPATH)    -L/usr/local/gcc/lib
	CXXINCPATH := $(CXXINCPATH) -I/usr/local/gcc/include -I/usr/local/gcc/include/boost-numeric-bindings
endif
ifneq (,$(shell find /Local -maxdepth 0 2>/dev/null))
	LIBPATH     := $(LIBPATH)    -L/Local/boost/lib
	CXXINCPATH  := $(CXXINCPATH) -I/Local/boost/include
endif

### Main part
targets = gen_tree_slha

.PHONY: all
all: $(targets)
gen_tree_slha: LDFLAGS += -lblas -llapack -lgfortran

-include .depend
.depend: $(wildcard *.c) $(wildcard *.cpp)
	@mkdep $^ $(CXXINCPATH)

$(targets) : % : .depend %.o
	@$(if $(DOERASE),rm -f $@)
	$(CXX) $(filter-out .depend,$^) -o $@ $(LIBPATH) $(CXXFLAGS) $(LDFLAGS)

%.o : %.cpp
	@$(if $(DOERASE),rm -f $*)
	@cp $< $<.`date "+%H%M%S"`.bak~
	$(CXX) $(CXXFLAGS) -c $< -o $@ ${MAKEARG}

%.o : %.f
	$(FC) $(FFLAGS) -c $< -o $@

%.o : %.f95
	$(FC) $(FFLAGS) -c $< -o $@


.PHONY: clean rm
clean:
	@rm -vf *~ *.bak tmp* .depend $(targets) *.exe *.so *.o
rm:
	@rm -vf *~ *.bak tmp*


$(targets:%=%!) : DOERASE = do
$(targets:%=%!) : %! : %
	@./$* $(ARG)

$(targets:%=%_) : %_ : %_clean %
$(targets:%=%_clean) : %_clean :
	-@rm -f $* $*.o


### Flymake ###
.PHONY: check-syntax-cpp check-syntax-f
check-syntax-cpp:
	$(CXX) $(FLYMAKE_CXXFLAGS) -S ${CHK_SOURCES} $(CXXINCPATH)
check-syntax-f:
	$(FC) $(FLYMAKE_FFLAGS) -S ${CHK_SOURCES}
check-syntax-f95:
	$(FC) $(FLYMAKE_F95FLAGS) -S ${CHK_SOURCES}


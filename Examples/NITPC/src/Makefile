Target = nitpc

ObjSuf = .o
SrcSuf = .cxx

OBJDIR = $(GARFIELD_HOME)/Object
SRCDIR = $(GARFIELD_HOME)/Source
INCDIR = $(GARFIELD_HOME)/Include
HEEDDIR = $(GARFIELD_HOME)/Heed
LIBDIR = $(GARFIELD_HOME)/Library

CFLAGS = -Wall -Wextra -Wno-long-long \
	-O3 -fno-common -c \
	-I$(INCDIR) -I$(HEEDDIR)
CFLAGS += $(shell root-config --cflags)
CFLAGS += -g

LDFLAGS = $(shell root-config --glibs)
LDFLAGS += -lGeom -lgfortran -lm
LDFLAGS += -L$(LIBDIR) -lGarfield

TargetSrc = $(Target:%=%$(SrcSuf))
TargetObj = $(Target:%=%$(ObjSuf))

all	: $(Target)

$(Target) : $(TargetObj)
	$(CXX) -o $@ $@.o $(LDFLAGS)
	mv $@ ../

.cxx.o  :
	$(CXX) $(CFLAGS) -c -o $@ $<

.SUFFIXES: .o .cxx

clean   :
	rm -f *.o
	rm -f *~
	rm -f ../$(Target)



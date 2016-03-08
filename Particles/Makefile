SHELL = /bin/bash

.DELETE_ON_ERROR:

.PHONY: all clean


ROOTCONFIG  := root-config
ROOTCFLAGS  := $(shell $(ROOTCONFIG) --cflags)
ROOTLDFLAGS := $(shell $(ROOTCONFIG) --ldflags)
ROOTLIBS    := $(shell $(ROOTCONFIG) --libs)
ROOTINCDIR  := $(shell $(ROOTCONFIG) --incdir)

CXX       := g++
CXXFLAGS  += -std=c++11 -O2 -Wall -fPIC $(ROOTCFLAGS)
LD        = g++
LDFLAGS   = -O2 $(ROOTLDFLAGS)

INCLUDES  := -I/$(ROOTINCDIR)
LIBS      := $(ROOTLIBS)

##############################################################################
all: particlesToRoot

particlesToRoot:particlesToRoot.o
	$(LD) $(LDFLAGS) $^ $(LIBS) -o $@

particlesToRoot.o:particlesToRoot.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

clean:
	@rm -rf particlesToRoot particlesToRoot.o

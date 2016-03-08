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
all: makePt2BPlot

makePt2BPlot:makePt2BPlot.o
	$(LD) $(LDFLAGS) $^ $(LIBS) -o $@

makePt2BPlot.o:makePt2BPlot.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

clean:
	@rm -rf makePt2BPlot makePt2BPlot.o

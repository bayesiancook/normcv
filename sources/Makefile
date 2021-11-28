CC=g++
CPPFLAGS= -Wall -O3 -std=c++11
LDFLAGS= -O3
SRCS=  Random.cpp Chrono.cpp

OBJS=$(patsubst %.cpp,%.o,$(SRCS))
ALL_SRCS=$(wildcard *.cpp)
ALL_OBJS=$(patsubst %.cpp,%.o,$(ALL_SRCS))

PROGSDIR=../data
ALL= simu truefit sstruefit cponormcv ssnormcv isnormcv 
PROGS=$(addprefix $(PROGSDIR)/, $(ALL))

.PHONY: all clean
all: $(PROGS)

# Rules for generate the dependencies automatically

%.d: %.cpp
	@echo "Generating dependencies for $<..."; \
	 set -e; rm -f $@; $(CC) -MM $(CPPFLAGS) $< > $@.$$$$; \
	 sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; rm -f $@.$$$$


# Rules to create .o files from .cpp files
%.o: %.cpp %.d
	$(CC) -c $(CPPFLAGS) $*.cpp

# Include the dependencies unless the request was to clean eveything up
ifneq ($(MAKECMDGOALS),clean)
-include $(ALL_OBJS:.o=.d)
endif

# Targets

$(PROGSDIR)/simu: SimulateData.o $(OBJS)
	$(CC) SimulateData.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

$(PROGSDIR)/sstruefit: SampleSizeTrueFit.o $(OBJS)
	$(CC) SampleSizeTrueFit.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

$(PROGSDIR)/truefit: TrueFit.o $(OBJS)
	$(CC) TrueFit.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

$(PROGSDIR)/cponormcv: CPONormalCV.o $(OBJS)
	$(CC) CPONormalCV.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

$(PROGSDIR)/ssnormcv: SSNormalCV.o $(OBJS)
	$(CC) SSNormalCV.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

$(PROGSDIR)/isnormcv: ISNormalCV.o $(OBJS)
	$(CC) ISNormalCV.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

clean:
	-rm -f *.o *.d *.d.*
	-rm -f $(PROGS)


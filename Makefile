CC	= g++
LD	= g++

#CCFLAGS	= -g -O0 `gsl-config --cflags` `root-config --cflags` -D_FINISH_DEBUG_
CCFLAGS	= -g -O0 `gsl-config --cflags` `root-config --cflags`
LDFLAGS	= -g -O0 -lgsl `gsl-config --libs` `root-config --libs`

all: CfBackgroundFlow

%.o: %.cc
	$(CC) $^ -o $@ $(CCFLAGS) -c

CorrFctnDirectYlm.o: CorrFctnDirectYlm.cxx CorrFctnDirectYlm.h sf.h
	$(CC) $< -o $@ $(CCFLAGS) -c

FemtoFlowDatabase.o: FemtoFlowDatabase.cxx FemtoFlowDatabase.h
	$(CC) $< -o $@ $(CCFLAGS) -c

cfylmFLOW.o: cfylmFLOW.cxx
	$(CC) $< -o $@ $(CCFLAGS) -c

CfBackgroundFlow: cfylmFLOW.o CorrFctnDirectYlm.o FemtoFlowDatabase.o ylm.o
	$(LD) -o $@ $^ $(LDFLAGS)

clean:
		rm -f CfBackgroundFlow *.o

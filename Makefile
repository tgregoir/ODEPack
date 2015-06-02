FSPLIT := $(shell pwd)/utils/fsplit

OPKSA1_FILES := rumach rumsum scfode sewset sintdy sprepj ssolsy ssrcom sstode svnorm
OPKSA1_FLAGS := $(foreach dir,$(OPKSA1_FILES),-e$(dir)) \

OPKSA2_FILES := isamax iumach ixsav saxpy sdot sgbfa sgbsl sgefa sgesl sscal xerrwv xsetf xsetun
OPKSA2_FLAGS := $(foreach dir,$(OPKSA2_FILES),-e$(dir)) \

.PHONY: all extract clean

all: clean $(FSPLIT) extract

$(FSPLIT): utils/fsplit.c
	gcc -Dlint utils/fsplit.c -o $(FSPLIT)

extract: original/opksa1.f original/opksa2.f original/opksmain.f
	mkdir -p extracted
	mkdir -p extracted/sub_opksa1
	mkdir -p extracted/sub_opksa2
	mkdir -p extracted/sub_opksmain
	cd extracted/sub_opksa1 && $(FSPLIT) $(OPKSA1_FLAGS) ../../original/opksa1.f
	cd extracted/sub_opksa2 && $(FSPLIT) $(OPKSA2_FLAGS) ../../original/opksa2.f
	cd extracted/sub_opksmain && $(FSPLIT) -eslsode ../../original/opksmain.f
	#cd extracted/sub_opksa1 && gfortran -c *.f
	#cd extracted/sub_opksa2 && gfortran -c *.f
	#cd extracted/sub_opksmain && gfortran -c *.f

clean:
	rm -rf extracted
	rm -f $(FSPLIT)

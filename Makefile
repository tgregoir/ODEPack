FSPLIT := $(shell pwd)/utils/fsplit

OPKSA1_FILES := rumach rumsum scfode sewset sintdy sprepj ssolsy ssrcom sstode svnorm
OPKSA1_FLAGS := $(foreach dir,$(OPKSA1_FILES),-e$(dir)) \

OPKSA2_FILES := isamax iumach ixsav saxpy sdot sgbfa sgbsl sgefa sgesl sscal xerrwv xsetf xsetun
OPKSA2_FLAGS := $(foreach dir,$(OPKSA2_FILES),-e$(dir)) \

.PHONY: all extract clean

all: clean $(FSPLIT) extract

$(FSPLIT): utils/fsplit.c
	@echo "Building fsplit..."
	@gcc -Dlint utils/fsplit.c -o $(FSPLIT)

extract: original/opksa1.f original/opksa2.f original/opksmain.f
	@echo "Creating directories..."
	@mkdir -p extracted
	@mkdir -p extracted/sub_opksa1
	@mkdir -p extracted/sub_opksa2
	@mkdir -p extracted/sub_opksmain
	@echo "Extracting Fortran subroutines..."
	@cd extracted/sub_opksa1 && $(FSPLIT) $(OPKSA1_FLAGS) ../../original/opksa1.f
	@cd extracted/sub_opksa2 && $(FSPLIT) $(OPKSA2_FLAGS) ../../original/opksa2.f
	@cd extracted/sub_opksmain && $(FSPLIT) -eslsode ../../original/opksmain.f
	@echo "Converting into C code..."
	@cd extracted/sub_opksa1 && f2c *.f
	@cd extracted/sub_opksa2 && f2c *.f
	@cd extracted/sub_opksmain && f2c *.f
	@echo "Moving C files to a proper directory..."
	@mkdir -p cfiles
	@mv extracted/sub_opksa1/*.c cfiles
	@mv extracted/sub_opksa2/*.c cfiles
	@mv extracted/sub_opksmain/*.c cfiles
	#@echo "Patching code for C++..."
	#@sed -i '1i #ifdef __cplusplus\nextern "C" {\n#endif\n' cfiles/*.c
	#@sed -i '$$a#ifdef __cplusplus\n}\n#endif\n' cfiles/*.c
	@echo "Compiling libf2c..."
	@make -C libf2c
	@echo "int MAIN__(void) {}\nint xargc;\nchar **xargv;" > cfiles/fix.c
	@cd cfiles && gcc -c *.c

clean:
	@echo "Cleaning up..."
	@rm -rf extracted
	@rm -rf cfiles
	@rm -f $(FSPLIT)
	@make -C libf2c clean
	@make -C driver clean

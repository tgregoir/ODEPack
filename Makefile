FSPLIT := $(shell pwd)/utils/fsplit

.PHONY: all extract clean

all: clean $(FSPLIT) extract

$(FSPLIT): utils/fsplit.c
	gcc -Dlint utils/fsplit.c -o $(FSPLIT)

extract: original/opksa1.f original/opksa2.f original/opksmain.f
	mkdir -p extracted
	mkdir -p extracted/sub_opksa1
	mkdir -p extracted/sub_opksa2
	mkdir -p extracted/sub_opksmain
	cd extracted/sub_opksa1 && $(FSPLIT) ../../original/opksa1.f
	cd extracted/sub_opksa2 && $(FSPLIT) ../../original/opksa2.f
	cd extracted/sub_opksmain && $(FSPLIT) ../../original/opksmain.f

clean:
	rm -rf extracted
	rm -f $(FSPLIT)

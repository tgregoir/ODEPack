LIBF2C_FILES := s_copy fmt err sig_die open util close \
                endfile sfe wsfe wrtfmt wref ctype fmtlib \
                s_stop pow_dd pow_ri r_sign
OFILES := $(foreach f,$(LIBF2C_FILES), ../libf2c/$(f).o)

all: demo

demo: test.c lsode.c
	@gcc -c test.c
	@gcc -c lsode.c
	@gcc test.o lsode.o ../cfiles/*.o $(OFILES) -o demo -lm

clean:
	@rm -f test.o demo

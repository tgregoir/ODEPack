LIBF2C_FILES := s_copy fmt err sig_die open util close \
                endfile sfe wsfe wrtfmt wref ctype fmtlib \
                s_stop pow_dd pow_ri r_sign
OFILES := $(foreach f,$(LIBF2C_FILES), ../libf2c/$(f).o)

all: demo

demo: test.c odesolve.c
	@gcc -c test.c
	@gcc -c odesolve.c
	@gcc test.o odesolve.o ../cfiles/*.o $(OFILES) -o demo -lm

clean:
	@rm -f test.o odesolve.o demo

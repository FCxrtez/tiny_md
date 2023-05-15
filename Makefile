CC      = clang $(USER_DEFINES)
CFLAGS  = -Ofast -march=native -Rpass=loop-vectorize -Rpass-missed=loop-vectorize -Rpass-analysis=loop-vectorize -ffp-contract=fast
WFLAGS	= -std=c11 -Wall -Wextra -g
LDFLAGS	= -lm

# TARGETS	= tiny_md viz
TARGETS	= tiny_md
SOURCES	= $(shell echo *.c)
OBJECTS = core.o wtime.o

all: $(TARGETS)

viz: viz.o $(OBJECTS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS) -lGL -lGLU -lglut

tiny_md: tiny_md.o $(OBJECTS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

%.o: %.c
	$(CC) $(WFLAGS) $(CPPFLAGS) $(CFLAGS) -c $<

clean:
	rm -f $(TARGETS) *.o *.xyz *.log .depend

.depend: $(SOURCES)
	$(CC) -MM $^ > $@

test:
	python3 test/test_main.py

rebuild-and-test:
	make clean
	make
	make test

benchmark:
	python3 benchmark/main.py $(name) --n-values 256 500 864 1372 2048 --num-runs 3

plot:
	python3 benchmark/main.py $(name) --mode plot

rebuild-and-benchmark:
	make clean
	make
	make benchmark

-include .depend

.PHONY: clean all test benchmark rebuild-and-test rebuild-and-benchmark plot

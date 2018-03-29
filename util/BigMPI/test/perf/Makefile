CC     = mpicc
CFLAGS = -O2 -Wall -std=c99

TESTS = reductions typecontig typepiggy

all: $(TESTS)

%: %.c
	$(CC) $(CFLAGS) $< -o $@

clean:
	-rm -f $(TESTS)


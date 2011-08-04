#MAKEFILE for ccsim ammonia-water ccombined cycle simulator
#By: Mark Taurence
#Created 1-1-10


target: twophasedyn.o nh3propslib.o
	cc -o $@ $<

twophasedyn.o: twophasedyn.cc
	cc -c -o $@ $<

nh3propslib.o: nh3propslib.cc
	cc -c -o $@ $<

.PHONY: clean

clean:
	rm -f *.o



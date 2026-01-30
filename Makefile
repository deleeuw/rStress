
SSRC = smacofIsotone.c smacofMPInverseV.c smacofSort.c \
	smacofSSRStressEngine.c smacofSSRStressMajorize.c smacofSSRStressMonotone.c

%.o: %.c smacofSSRStress.h
	clang -c $@

shlib: smacofSSRStress.h $(SSRC)
	R CMD SHLIB -o smacofSSRStress.so $(SSRC)

clean:
	rm -f *.o

pristine:
	rm -f *.o *.so


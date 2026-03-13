
SSRC = smacofIsotone.c smacofMPInverseV.c smacofSort.c \
	smacofSSFStressEngine.c smacofSSFStressMajorize.c smacofSSFStressMonotone.c \
	smacofSSPrint.c smacofSSFStressFlist.c

%.o: %.c smacofSSFStress.h
	clang -c $@

shlib: smacofSSFStress.h $(SSRC)
	R CMD SHLIB -o smacofSSFStress.so $(SSRC)

clean:
	rm -f *.o

pristine:
	rm -f *.o *.so


LIBS= -Wl,-rpath,${TACC_MKL_LIB} -L${TACC_MKL_LIB} -lmkl_intel_lp64 -lmkl_sequential -lmkl_lapack -lmkl_core

fopts = -I${TACC_MKL_INC} -O3 -fp-model precise

OBJECTS_GIVENS = test_givens.o \
                 sub.abort2.o \
                 set_h_y_i.o \
                 givens.o


test_givens: $(OBJECTS_GIVENS)
	ifort -o test_givens ${fopts} ${LIBS} $(OBJECTS_GIVENS)

tidy :
	rm -f $(OBJECTS_GIVENS)

clean : tidy
	rm -f test_givens

.f.o:
	ifort -c ${fopts} $<

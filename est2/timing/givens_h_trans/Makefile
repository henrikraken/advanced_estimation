LIBS= -Wl,-rpath,${TACC_MKL_LIB} -L${TACC_MKL_LIB} -lmkl_intel_lp64 -lmkl_sequential -lmkl_lapack -lmkl_core

fopts = -I${TACC_MKL_INC} -O3 -fp-model precise

OBJECTS_GIVENS_h_trans = test_givens_h_trans.o \
                 sub.abort2.o \
                 set_h_y_i.o \
                 givens_h_trans.o

test_givens_h_trans: $(OBJECTS_GIVENS_h_trans)
	ifort -o test_givens_h_trans ${fopts} ${LIBS} $(OBJECTS_GIVENS_h_trans)

tidy :
	rm -f $(OBJECTS_GIVENS_h_trans)

clean : tidy
	rm -f test_givens_h_trans

.f.o:
	ifort -c ${fopts} $<

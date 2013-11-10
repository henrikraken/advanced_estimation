LIBS= -Wl,-rpath,${TACC_MKL_LIB} -L${TACC_MKL_LIB} -lmkl_intel_lp64 -lmkl_sequential -lmkl_lapack -lmkl_core

fopts = -I${TACC_MKL_INC} -O3 -fp-model precise

OBJECTS_test_house_blas2 = test_house_blas2.o \
                sub.abort2.o \
                house_blas2.o \
                set_h_y_i.o

test_house_blas2: $(OBJECTS_test_house_blas2)
	ifort -o test_house_blas2 ${fopts} ${LIBS} $(OBJECTS_test_house_blas2)

test: test.o
	ifort -o test ${fopts} ${BLAS_LIBS} ${LAPACK_LIBS} test.o

tidy :
	rm -f $(OBJECTS_test_house_blas2)

clean : tidy
	rm -f test_house_blas2

.f.o:
	ifort -c ${fopts} $<

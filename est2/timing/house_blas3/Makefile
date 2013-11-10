LIBS= -Wl,-rpath,${TACC_MKL_LIB} -L${TACC_MKL_LIB} -lmkl_intel_lp64 -lmkl_sequential -lmkl_lapack -lmkl_core

fopts = -I${TACC_MKL_INC} -O3 -fp-model precise

OBJECTS_test_house_blas3 = test_house_blas3.o \
                sub.abort2.o \
                house_blas3.o \
                set_h_y_i.o

test_house_blas3: $(OBJECTS_test_house_blas3)
	mpif90 -o test_house_blas3 ${fopts} ${LIBS} $(OBJECTS_test_house_blas3)

test: test.o
	mpif90 -o test ${fopts} ${BLAS_LIBS} ${LAPACK_LIBS} test.o

tidy :
	rm -f $(OBJECTS_test_house_blas3)

clean : tidy
	rm -f test_house_blas3

.f.o:
	mpif90 -c ${fopts} $<

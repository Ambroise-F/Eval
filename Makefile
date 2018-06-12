CC=/opt/rh/devtoolset-7/root/usr/bin/gcc


CFLAGS=-O3 -fopenmp -mfma -fno-trapping-math

eval : evalnext2.c
	$(CC) -S $(CFLAGS) -mavx2 $< -o $@.s; grep $@.s -i -e mul -e fma -e vfn -e vfm
	$(CC) $(CFLAGS) -mavx2 $< -o $@

eval2 : evalnext2.c
	$(CC) -S $(CFLAGS) -march=knl $< -o $@.s; grep $@.s -i -e mul -e fma -e vfn -e vfm
	$(CC) $(CFLAGS) -march=knl $< -o $@
 

clean :
	@rm -f eval eval2 a.out *.s *~

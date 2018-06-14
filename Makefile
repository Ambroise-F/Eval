CC=/opt/rh/devtoolset-7/root/usr/bin/gcc


CFLAGS=-O3# -fopenmp -mfma -fno-trapping-math -mavx2


all : eval og verif

eval : evalnext2.c
	$(CC) -S $(CFLAGS) $< -o $@.s; grep $@.s -i -e fma -e vfn -e vfm
	$(CC) $(CFLAGS) $< -o $@

og : evalnext2_og.c
	$(CC) -S $(CFLAGS) $< -o $@.s; grep $@.s -i -e fma -e vfn -e vfm
	$(CC) $(CFLAGS) $< -o $@

eval2 : evalnext2.c
	$(CC) -S $(CFLAGS) -march=knl $< -o $@.s; grep $@.s -i -e mul -e fma -e vfn -e vfm
	$(CC) $(CFLAGS) -march=knl $< -o $@
 


#--- VERIF 

evalverif : evalnext2.c
	$(CC) -S $(CFLAGS) -DVERIF=1 $< -o $@.s; grep $@.s -i -e fma -e vfn -e vfm
	$(CC) $(CFLAGS) -DVERIF=1 $< -o $@

ogverif : evalnext2_og.c
	$(CC) -S $(CFLAGS) -DVERIF=1 $< -o $@.s; grep $@.s -i -e fma -e vfn -e vfm
	$(CC) $(CFLAGS) -DVERIF=1 $< -o $@

verif : evalnext2.c evalnext2_og.c evalverif ogverif 
	bash verif.sh





clean :
	@rm -f eval eval2 a.out *.s *~ og.txt eval.txt ogverif evalverif og

CC=/opt/rh/devtoolset-7/root/usr/bin/gcc


CFLAGS=-O3 -fopenmp -march=knl -mfma -fno-trapping-math#-fopt-info-vec=vec.out#-fopenmp -fno-trapping-math -mavx2 -march=knl# -mavx2


all : eval og verif

eval : evalnext2.c
	$(CC) -S -fverbose-asm $(CFLAGS) $< -o $@.s -lm; grep $@.s -i -e fma -e vfn -e vfm
	$(CC) $(CFLAGS) $< -o $@ -lm

og : evalnext2_og.c
	$(CC) -S $(CFLAGS) $< -o $@.s; grep $@.s -i -e fma -e vfn -e vfm
	$(CC) $(CFLAGS) $< -o $@

eval2 : evalnext2.c
	$(CC) -S $(CFLAGS) -march=knl $< -o $@.s -lm; grep $@.s -i -e mul -e fma -e vfn -e vfm
	$(CC) $(CFLAGS) -march=knl $< -o $@ -lm
 


#--- VERIF 

evalverif : evalnext2.c
	$(CC) -S $(CFLAGS) -DVERIF=1 $< -o $@.s -lm; grep $@.s -i -e fma -e vfn -e vfm
	$(CC) $(CFLAGS) -DVERIF=1 $< -o $@ -lm

ogverif : evalnext2_og.c
	$(CC) -S $(CFLAGS) -DVERIF=1 $< -o $@.s -lm; grep $@.s -i -e fma -e vfn -e vfm
	$(CC) $(CFLAGS) -DVERIF=1 $< -o $@ -lm

verif : evalnext2.c evalnext2_og.c evalverif ogverif 
	bash verif.sh


clean :
	@rm -f eval eval2 a.out *.s *~ og.txt eval.txt ogverif evalverif og

#opt
#macro
#align

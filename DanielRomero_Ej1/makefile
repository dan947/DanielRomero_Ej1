placas: placas.pdf  
	evince placas.pdf 
placas.pdf: grafica.py data.dat	
	python grafica.py
data.dat: placas.c
	mpicc -o placas placas.c 
	./a.out > data.dat
.PHONY: clean
clean:
rm -f a.out a.dat placas.pdf

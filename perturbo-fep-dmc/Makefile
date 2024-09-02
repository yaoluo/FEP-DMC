
default: all

all: qe2pert perturbo

qe2pert: pw ph
	cd qe2pert-src && make
	(mkdir -p bin)
	(cd bin; ln -sf ../qe2pert-src/qe2pert.x .)

perturbo: pw ph
	cd pert-src && make
	(mkdir -p bin)
	(cd bin; ln -sf ../pert-src/perturbo.x .)

pw: 
	cd .. && make pw

ph: 
	cd .. && make ph

clean:
	cd pert-src && make clean
	cd qe2pert-src && make clean

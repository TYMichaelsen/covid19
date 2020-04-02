all: 
	bash ./prep_build.sh
	singularity build --fakeroot singularity-covid19.sif Covid19.Singularity

clean:
	rm -f singularity-covid19.sif
	rm -rf build
	rm -f singularity-r_3.6.3.sif

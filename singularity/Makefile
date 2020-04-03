all: 
	bash ./prep_build.sh
	singularity build --fakeroot -F singularity-covid19.sif covid19.singularity

clean:
	rm -f singularity-covid19.sif
	rm -rf build
	rm -f singularity-r_3.6.3.sif

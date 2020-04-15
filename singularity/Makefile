all:
	singularity build --fakeroot -F covid19_1.5.sif covid19_base.singularity
	singularity build --fakeroot -F covid19_latest.sif covid19.singularity

clean:
	rm -f covid19_1.5.sif
	rm -f covid19_latest.sif

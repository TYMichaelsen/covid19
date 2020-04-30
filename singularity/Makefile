all: covid19_1.5.sif
	./wrapbuild.sh covid19_latest.sif covid19.singularity

covid19_1.5.sif:
	./wrapbuild.sh covid19_1.5.sif covid19_base.singularity

clean:
	rm -f covid19_1.5.sif
	rm -f covid19_latest.sif

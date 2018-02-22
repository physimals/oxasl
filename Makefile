PROJNAME = oxford_asl

# Pass Git revision details
GIT_VERSION:=$(shell git describe --dirty)
GIT_DATE:=$(shell git log -1 --format=%ad --date=local)

# Always rebuild scripts
.PHONY: FORCE

meta.yaml: FORCE
	sed -e "s/@GIT_VERSION@/${GIT_VERSION}/" -e "s/@GIT_DATE@/${GIT_DATE}/" <meta.yaml.in >meta.yaml

clean:
	rm -f meta.yaml

FORCE:


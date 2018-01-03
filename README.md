# nf-repeatmasking
Nextflow workflow for automatic repeat detection, classification and masking.

# Introduction
This pipeline is a copy
of [MAKER](http://www.yandell-lab.org/software/maker.html)
recommendation for
'[advanced repeat library construction](http://weatherby.genetics.utah.edu/MAKER/wiki/index.php/Repeat_Library_Construction-Advanced)'.

# Prerequisites
The only requirement is [Nextflow](https://www.nextflow.io), which can
be installed on any POSIX system (Linux, Solaris, OS X) with Java 7 or 8:

```sh
curl -s https://get.nextflow.io | bash
```

The pipeline itself requires a great many different pieces of
software. There are two ways to install all of the software packages
and scripts required for this pipeline - natively, or by using docker
(recommended).

## Docker
This is certainly the easier (and more reproducible) method. I've
already built a docker image that includes almost all of the software
necessary to run the pipeline. I'm unable to include the RepBase
repeat database due to licencing restrictions, but you can register
yourself for free (academic use only) and bundle it into a new docker
image yourself:

```sh
# We don't need the tar.gz once the image is built, so we'll use a temporary directory.
cd `mktemp -d`

# Download the RepBase repeat library (replace RB_USERNAME and RB_PASSWORD with your username and password)
wget --user $RB_USERNAME \
	 --password $RB_PASSWORD \
	 -O repeatmaskerlibraries.tar.gz \
     http://www.girinst.org/server/RepBase/protected/repeatmaskerlibraries/RepBaseRepeatMaskerEdition-20170127.tar.gz

# Make an (almost) empty Dockerfile.
#  All of the important instructions are in the repeatmasker-onbuild image. You can either grab the pre-build copy from Dockerhub, 
#  or you can build it from the Dockerfiles/nf-repeatmasking-onbuild directory of this repository.
echo "FROM robsyme/nf-repeatmasking-onbuild" > Dockerfile

# Build a new docker images called 'repeats'.
#  When building, this image looks for a file called 'repeatmaskerlibraries.tar.gz' which it pulls into the image.
docker build -t robsyme/nf-repeatmasking .
```

## Natively
You'll need the following pieces of software

- Bioperl
- Hmmer
- MITE Hunter
- Genometools
- RepeatMasker
- Blast+ (v2.4.0)
- RepeatModeler
- R
- ggplot2
- dplyr
- tidyr
- magrittr

## Citing

[![DOI](https://zenodo.org/badge/1134750.svg)](https://zenodo.org/badge/latestdoi/1134750)


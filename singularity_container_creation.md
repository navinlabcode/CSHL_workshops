SIF Files for workshops located here:

```bash
/volumes/USR2/Ryan/singularity/cshl_course/scdna.sif
/volumes/USR2/Ryan/singularity/cshl_course/scrna.sif
/volumes/USR2/Ryan/singularity/cshl_course/spatial.sif
```

TODO:
* Add data sets into sif files?
* See if I can get the Rstudio integration to work
* All SIFs need to be tested for all requested packages still.

## single cell RNA
Requested packages:
* seurat 5
* harmony
* copyKAT
* inferCNV
* singleR
* fastMNN


scrna.def
```bash
Bootstrap: docker
From: ubuntu:latest

%environment
    # set up all essential environment variables
    export LC_ALL=C
    export PATH=/opt/miniconda3/bin:$PATH
    export PYTHONPATH=/opt/miniconda3/lib/python3.9/:$PYTHONPATH

%post
 
	# update and install essential dependencies
	apt-get -y update
	apt-get update && apt-get install -y automake \
	build-essential \
	bzip2 \
	wget \
	git \
	default-jre \
	unzip \
	zlib1g-dev \
	parallel

	# download, install, and update miniconda3
	wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
	bash Miniconda3-latest-Linux-x86_64.sh -b -f -p /opt/miniconda3/
	rm Miniconda3-latest-Linux-x86_64.sh
	chmod --recursive a+rw /opt/miniconda3

	# install dependencies via conda
	export PATH="/opt/miniconda3/bin:$PATH"

	#build full conda environment for sif
	conda install -y -c conda-forge mamba 
	conda config --add channels bioconda
	conda config --add channels conda-forge
	
	# denotes installed in sandbox without error
	#install additional tools
	mamba install -y -f bioconda::bwa #
	mamba install -y -f bioconda::samtools #
	mamba install -y -f bioconda::bedtools #
	mamba install -y -f bioconda::fastqc #
	mamba install -y -f bioconda::multiqc #
	mamba install -y -f anaconda::graphviz #
	mamba install -y -f conda-forge::parallel #
	conda install -y -f conda-forge::ncurses #

	#install base R packages
	conda install -y -f r-base=4.2 #
	mamba install -y -f conda-forge::r-devtools #
	mamba install -y -f conda-forge::r-biocmanager=1.30.19 #
	mamba install -y -f conda-forge::r-rlang #
	mamba install -y -f conda-forge::r-ggplot2 #

	#R utility libraries
	R --slave -e 'install.packages("remotes", repos="http://cran.us.r-project.org")' #
	R --slave -e 'install.packages("circlize", repos="http://cran.us.r-project.org")' #
	R --slave -e 'install.packages("optparse", repos="http://cran.us.r-project.org")' #
	R --slave -e 'install.packages("patchwork", repos="http://cran.us.r-project.org")' #
	R --slave -e 'install.packages("plyr", repos="http://cran.us.r-project.org")' #
	R --slave -e 'install.packages("stringr", repos="http://cran.us.r-project.org")' #
	R --slave -e 'install.packages("tidyverse", repos="http://cran.us.r-project.org")' #
	R --slave -e 'install.packages("RColorBrewer", repos="http://cran.us.r-project.org")' #

	#Bioconductor packages through conda
	mamba install -y -f bioconda::bioconductor-biocparallel #
	mamba install -y -f bioconda::bioconductor-bsgenome.hsapiens.ucsc.hg38 #
	mamba install -y -f bioconda::bioconductor-ensdb.hsapiens.v86 #
	mamba install -y -f bioconda::bioconductor-org.hs.eg.db #
	mamba install -y -f bioconda::bioconductor-txdb.hsapiens.ucsc.hg38.knowngene #
	mamba install -y -f bioconda::bioconductor-decoupler #
	mamba install -y -f bioconda::bioconductor-scran #
	mamba install -y -f bioconda::bioconductor-infercnv #
	mamba install -y -f bioconda::bioconductor-complexheatmap #
	mamba install -y -f bioconda::bioconductor-biovizbase #
	mamba install -y -f bioconda::bioconductor-singler

	#Funner stuff!
	R --slave -e 'install.packages("rliger", repos="http://cran.us.r-project.org")' #
	R --slave -e 'devtools::install_github("navinlabcode/copykat")'
	R --slave -e 'install.packages("Seurat", repos="http://cran.us.r-project.org")' #
	R --slave -e 'remotes::install_github("satijalab/seurat-wrappers")' #
	R --slave -e 'install.packages("harmony")'

%labels
    Author Ryan Mulqueen
    Version v0.2
    MyLabel scRNA


```
Build on GEO with fakeroot.

```bash
singularity build --fakeroot scrna.sif scrna.def
```

## Spatial transcriptomics SIF DONE
Requested packages:
* CellTrek
* Seurat 

spatial.def
```bash
Bootstrap: docker
From: ubuntu:latest

%environment
	# set up all essential environment variables
	export LC_ALL=C
	export PATH=/opt/miniconda3/bin:$PATH
	export PYTHONPATH=/opt/miniconda3/lib/python3.9/:$PYTHONPATH
	export LC_ALL=C.UTF-8
%post
	# update and install essential dependencies
	apt-get -y update
	apt-get update && apt-get install -y automake \
	build-essential \
	bzip2 \
	wget \
	git \
	default-jre \
	unzip \
	zlib1g-dev \
	parallel \
	libglpk40 \
	gfortran

	# download, install, and update miniconda3
	wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
	bash Miniconda3-latest-Linux-x86_64.sh -b -f -p /opt/miniconda3/
	rm Miniconda3-latest-Linux-x86_64.sh

	# install dependencies via conda
	export PATH="/opt/miniconda3/bin:$PATH"
	conda install -y -c conda-forge mamba 

	#install R packages
	mamba install -y -f conda-forge::r-base=4.1
	mamba install -y -f conda-forge::r-devtools
	mamba install -y -f conda-forge::r-essentials
	#mamba install -y -f bioconda::bioconductor-bsgenome.hsapiens.ucsc.hg38 #
	#mamba install -y -f bioconda::bioconductor-ensdb.hsapiens.v86 #

	#R utility libraries
	R --slave -e 'install.packages("remotes", repos="http://cran.us.r-project.org")' #
	R --slave -e 'install.packages("circlize", repos="http://cran.us.r-project.org")' #
	R --slave -e 'install.packages("optparse", repos="http://cran.us.r-project.org")' #
	R --slave -e 'install.packages("patchwork", repos="http://cran.us.r-project.org")' #
	R --slave -e 'install.packages("plyr", repos="http://cran.us.r-project.org")' #
	R --slave -e 'install.packages("stringr", repos="http://cran.us.r-project.org")' #
	R --slave -e 'install.packages("tidyverse", repos="http://cran.us.r-project.org")' #
	R --slave -e 'install.packages("RColorBrewer", repos="http://cran.us.r-project.org")' #

	R --slave -e 'install.packages("BiocManager", repos="http://cran.us.r-project.org")'
	R --slave -e 'BiocManager::install("biovizBase")'
	R --slave -e 'install.packages("Signac", repos="http://cran.us.r-project.org")'
	R --slave -e 'install.packages("Seurat", repos="http://cran.us.r-project.org")'
	R --slave -e 'install.packages("stringr", repos="http://cran.us.r-project.org")'
	R --slave -e 'install.packages("devtools", repos="http://cran.us.r-project.org")'
	R --slave -e 'devtools::install_github("navinlabcode/CellTrek")'
	R --slave -e 'BiocManager::install("EnsDb.Hsapiens.v86")'

%labels
	Author Ryan Mulqueen
	Version v0.1
	MyLabel SpatialTranscriptomics

```

```bash
singularity build --fakeroot spatial.sif spatial.def
singularity build --fakeroot spatial.def

singularity build --fakeroot --sandbox spatial/ docker://ubuntu:latest
singularity shell --fakeroot --writable spatial/

```


# Single cell DNA SIF DONE
Requested packages:
* copyKIT
* R version 4.2
* R packages: tidyverse

### Define the image to be created

scdna.def
```bash
Bootstrap: docker
From: ubuntu:latest

%environment
	# set up all essential environment variables
	export LC_ALL=C
	export PATH=/opt/miniconda3/bin:$PATH
	export PYTHONPATH=/opt/miniconda3/lib/python3.9/:$PYTHONPATH
	export LC_ALL=C.UTF-8
%post
	# update and install essential dependencies
	apt-get -y update
	apt-get update && apt-get install -y automake \
	build-essential \
	bzip2 \
	wget \
	git \
	default-jre \
	unzip \
	zlib1g-dev \
	parallel \
	libglpk40 \
	gfortran

	# download, install, and update miniconda3
	wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
	bash Miniconda3-latest-Linux-x86_64.sh -b -f -p /opt/miniconda3/
	rm Miniconda3-latest-Linux-x86_64.sh

	# install dependencies via conda
	export PATH="/opt/miniconda3/bin:$PATH"
	conda install -y -c conda-forge mamba 
	mamba install -y -f bioconda::samtools #
	mamba install -y -f bioconda::bedtools #
	mamba install -y -f conda-forge::parallel #

	#install R packages
	mamba install -y -f conda-forge::r-base #=4.2
	mamba install -y -f conda-forge::r-devtools
	mamba install -y -f conda-forge::r-tidyverse
	#mamba install -y -f bioconda::bioconductor-ensdb.hsapiens.v86


	#R utility libraries
	R --slave -e 'install.packages("remotes", repos="http://cran.us.r-project.org")' #
	R --slave -e 'install.packages("circlize", repos="http://cran.us.r-project.org")' #
	R --slave -e 'install.packages("optparse", repos="http://cran.us.r-project.org")' #
	R --slave -e 'install.packages("patchwork", repos="http://cran.us.r-project.org")' #
	R --slave -e 'install.packages("plyr", repos="http://cran.us.r-project.org")' #
	R --slave -e 'install.packages("stringr", repos="http://cran.us.r-project.org")' #
	R --slave -e 'install.packages("tidyverse", repos="http://cran.us.r-project.org")' #
	R --slave -e 'install.packages("RColorBrewer", repos="http://cran.us.r-project.org")' #

	R --slave -e 'devtools::install_github("navinlabcode/copykit")'
	conda install -y -f --no-deps conda-forge::r-igraph
	conda install -y -f --no-deps bioconda::bioconductor-bluster
	conda install -y -f --no-deps bioconda::bioconductor-copynumber
	conda install -y -f --no-deps bioconda::bioconductor-ggtree
	wget https://github.com/navinlabcode/copykit/releases/download/v.0.1.2/copykit_0.1.2.tar.gz
	R --slave -e 'install.packages("copykit_0.1.2.tar.gz", repos = NULL)' # the install_github is broken so pulling from archive

	R --slave -e 'install.packages("optparse", repos="http://cran.us.r-project.org")' #


%labels
	Author Ryan Mulqueen
	Version v0.1
	MyLabel Copykit 

```

```bash
singularity build --fakeroot scdna.sif scdna.def
```


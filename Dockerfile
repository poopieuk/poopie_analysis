FROM bioconductor/bioconductor_docker:RELEASE_3_18

# Install system dependencies for Python & pip
RUN apt-get update && apt-get install -y \
    python3 \
    python3-pip \
    python3-setuptools \
    python3-dev \
    && apt-get clean

# Install R packages
RUN R -e "install.packages(c('jsonlite','remotes','ggplot2','tidyverse','R.utils','optparse', 'dplyr', 'tidyr'), repos='https://cloud.r-project.org/')"
RUN R -e "BiocManager::install(c('dada2', 'phyloseq', 'ShortRead', 'DECIPHER'), ask=FALSE)"
# Optional extra packages you may need
RUN R -e "install.packages(c('data.table','reshape2','plyr','stringr','gridExtra'), repos='https://cloud.r-project.org')"

# Install Python packages for reporting
RUN pip3 install \
    pandas \
    matplotlib \
    seaborn \
    fpdf \
    fuzzywuzzy \
    reportlab \
    python-Levenshtein

WORKDIR /pipeline
COPY . /pipeline

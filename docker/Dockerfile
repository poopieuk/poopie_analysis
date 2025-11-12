FROM bioconductor/bioconductor_docker:RELEASE_3_18

# Install CRAN packages
RUN R -e "install.packages(c( \
    'tidyverse', \
    'data.table', \
    'optparse', \
    'jsonlite', \
    'ggtext', \
    'RColorBrewer', \
    'patchwork', \
    'ggpubr', \
    'glue', \
    'remotes' \
), repos='https://cloud.r-project.org')" 

# Install Bioconductor packages (including dada2, phyloseq, DECIPHER)
RUN R -e "BiocManager::install(c( \
    'dada2', \
    'phyloseq', \
    'DECIPHER', \
    'Biostrings', \
    'ShortRead', \
    'GenomicRanges' \
))"

# Optional: Install any GitHub packages you need
# RUN R -e "remotes::install_github('your/repo')"

# Clean up
RUN apt-get clean && \
    rm -rf /var/lib/apt/lists/* /tmp/downloaded_packages

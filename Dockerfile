FROM bioconductor/bioconductor_docker:RELEASE_3_18

# Install needed R packages
RUN R -e "install.packages(c('optparse'), repos='https://cloud.r-project.org/')"

RUN R -e "BiocManager::install(c('dada2', 'phyloseq', 'ShortRead', 'DECIPHER'), ask=FALSE)"

# Install Python packages for your report script
RUN pip install pandas matplotlib seaborn fpdf

# Set working directory
WORKDIR /pipeline
COPY . /pipeline

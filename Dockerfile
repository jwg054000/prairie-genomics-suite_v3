# Prairie Genomics Suite R Shiny - Docker Deployment
FROM rocker/shiny-verse:4.3.0

# Install system dependencies
RUN apt-get update && apt-get install -y \
    libxml2-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    && rm -rf /var/lib/apt/lists/*

# Install R packages
RUN install2.r --error \
    shiny \
    shinydashboard \
    shinyWidgets \
    DT \
    plotly \
    ggplot2 \
    dplyr \
    readr \
    readxl \
    RColorBrewer \
    pheatmap \
    ggrepel \
    viridis

# Install Bioconductor packages
RUN R -e "if (!require('BiocManager')) install.packages('BiocManager'); BiocManager::install('DESeq2')"

# Copy app files
COPY . /srv/shiny-server/prairie-genomics-suite/

# Set permissions
RUN chown -R shiny:shiny /srv/shiny-server/prairie-genomics-suite/

# Expose port
EXPOSE 3838

# Run app
CMD ["/usr/bin/shiny-server"]
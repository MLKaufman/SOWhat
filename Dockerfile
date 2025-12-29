# Dockerfile for SOWhat (Seurat Object Annotation Shiny App)

# Base image
FROM rocker/r-ver:4.4.2

# Set maintainer
LABEL maintainer="MLKaufman"

# Install system dependencies
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    zlib1g-dev \
    libfontconfig1-dev \
    libfreetype6-dev \
    libcairo2-dev \
    pkg-config \
    git \
    build-essential \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libglpk-dev \
    pandoc \
    && apt-get clean && rm -rf /var/lib/apt/lists/*
# Set the working directory
WORKDIR /app

# Copy renv.lock and renv infrastructure
COPY renv.lock renv.lock

# Install renv and restore the environment
RUN R -e "install.packages('renv', repos='https://cloud.r-project.org')" && \
    R -e "renv::init()" && \
    R -e "renv::install('BiocManager', prompt = FALSE)" && \
    R -e "renv::install('gitcreds', prompt = FALSE)" && \
    R -e "options(repos = BiocManager::repositories())" && \
    R -e "renv::restore(prompt = FALSE)"

# Copy the rest of the application files
COPY . .

# Expose the Shiny app port
EXPOSE 3838

# Run the Shiny app
CMD ["R", "-e", "shiny::runApp('/app', port=3838, host='0.0.0.0')"]

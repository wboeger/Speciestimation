FROM rocker/r-ver:4.4.2

# --- System dependencies -----------------------------------------------
# build-essential/g++: compile rstan/StanHeaders C++ code
# libcurl/libssl/libxml2: httr/curl-family R package deps
# pandoc: required by rmarkdown::render()
# libfontconfig/libfreetype/libpng/libtiff/libjpeg/libharfbuzz/libfribidi: ragg/textshaping (ggplot2 graphics device)
# cmake: required by RcppParallel (a genuine Imports/LinkingTo dependency of rstan)
# libuv1-dev: required by `fs` (a genuine Imports dependency of shiny/bslib/rmarkdown/DT)
# wget/perl/xzdec: required by tinytex::install_tinytex() (PDF manual rendering)
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    cmake \
    pandoc \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libfontconfig1-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    zlib1g-dev \
    libuv1-dev \
    wget \
    perl \
    xzdec \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

# --- R package installation (cached as its own layer) -------------------
COPY build/install_packages.R build/install_packages.R
RUN Rscript build/install_packages.R

# --- TinyTeX (minimal LaTeX distribution) for PDF manual rendering -------
RUN Rscript -e "tinytex::install_tinytex()"


# --- App source -----------------------------------------------------------
COPY . .

# --- Pre-compile Stan models into this image (platform-specific; cannot be
# baked on a developer machine and copied in) -----------------------------
RUN Rscript build/compile_models.R

# --- Render the static user manual to www/ (Shiny auto-serves www/ as
# static assets, so the download link needs no server-side handler) -------
RUN Rscript build/render_manual.R

EXPOSE 8080
CMD ["Rscript", "run.R"]

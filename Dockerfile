# RNAflow — reproducible container
#
# Pins the exact platform that RNAflow is developed against (R 4.5 /
# Bioconductor 3.22) so the heavy Bioconductor dependency stack resolves
# identically on any machine. This is the primary reproducibility guarantee
# for the project — it fixes the R/Bioc release and system libraries, which a
# renv lockfile alone cannot (see dev/make_renv_lock.R for the optional
# CRAN/Bioc version pin on top of this).
#
# Build:  docker build -t rnaflow .
# Run:    docker run --rm -p 8080:8080 rnaflow
# Open:   http://localhost:8080
FROM bioconductor/bioconductor_docker:RELEASE_3_22

LABEL org.opencontainers.image.title="RNAflow" \
      org.opencontainers.image.description="Interactive Shiny platform for end-to-end bulk RNA-seq analysis" \
      org.opencontainers.image.source="https://github.com/KmBioChemo/RNAflow" \
      org.opencontainers.image.licenses="MIT"

WORKDIR /rnaflow

# 1) Dependencies first, so Docker's layer cache only reinstalls the (slow)
#    Bioconductor stack when DESCRIPTION changes — not on every source edit.
COPY DESCRIPTION ./DESCRIPTION
RUN Rscript -e 'if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes"); \
                remotes::install_deps(dependencies = TRUE, upgrade = "never")'

# 2) Install RNAflow itself.
COPY . .
RUN R CMD INSTALL --no-multiarch --with-keep.source .

# Shiny must listen on all interfaces to be reachable from outside the container.
EXPOSE 8080
CMD ["Rscript", "-e", "options(shiny.host = '0.0.0.0'); RNAflow::run_app(port = 8080L, launch_browser = FALSE)"]

# Install Bioconductor packages needed for BSmooth analysis

module --force purge
module load StdEnv/2023
module load r/4.4.0
module load r-bundle-bioconductor/3.20

export R_LIBS_USER="$HOME/R/%p-library/%v"
mkdir -p "$R_LIBS_USER"

R --vanilla -q <<'EOF'
if (!requireNamespace("BiocManager", quietly=TRUE)) {
  install.packages("BiocManager", repos="https://cloud.r-project.org")
}
BiocManager::install("bsseq", update=FALSE, ask=FALSE)
q(save="no")
EOF


##then verify installation:
module --force purge
module load StdEnv/2023
module load r/4.4.0
module load r-bundle-bioconductor/3.20

R --vanilla -q <<'EOF'
cat(.libPaths(), sep="\n")
cat("bsseq installed? ", requireNamespace("bsseq", quietly=TRUE), "\n", sep="")
q(save="no")
EOF

#note the submission script 2_test.sh is updated to specify the Rlibrary path for the package.
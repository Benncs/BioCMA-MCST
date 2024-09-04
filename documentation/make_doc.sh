#!/bin/bash

# Set up directory paths
BUILD_DIR=./build
CLI_DIR=./cli
CLI_PDF_NAME=cli
MODEL_DIR=./models
EXAMPLE_DIR=./examples
MODEL_PDF_NAME=model
EXAMPLE_PDF_NAME=examples

DOC_DIR=docs
EXAMPLES_SCRIPT_DIR=examples

# Function to build PDFs
build_pdf() {
    # Change to the specified directory or exit if it fails
    cd "$2" || { echo "Error: Could not change to directory $2"; exit 1; }

    # Compile LaTeX document and redirect output to /dev/null to suppress messages
    pdflatex main.tex || { echo "Error: pdflatex failed"; exit 1; }
    bibtex main  || { echo "Error: bibtex failed"; exit 1; }
    pdflatex main.tex  || { echo "Error: pdflatex failed"; exit 1; }

    # Remove auxiliary files and redirect output to /dev/null to suppress messages
    rm -f *.aux *.bbl *.blg *.log *.lot *.out *.xml *.toc *.lof *-blx.bib  2>&1 || { echo "Error: Failed to remove auxiliary files"; exit 1; }

    # Go back to the original directory
    cd ..

    # Copy the generated PDF to the specified destination
    cp "$2/main.pdf" "$3/$1.pdf" || { echo "Error: Failed to copy PDF file"; exit 1; }
}

# Store the current working directory
current_pwd=$(pwd)

# Copy PDF files to the specified destination or print an error message
# cp ./${EXAMPLES_SCRIPT_DIR}/*.pdf ./${DOC_DIR}/${EXAMPLE_DIR}/assets || { echo "Error: Failed to copy PDF files"; exit 1; }

# Change to the documentation directory
# cd ${DOC_DIR}

# Create the build directory if it doesn't exist
mkdir -p "${BUILD_DIR}" || { echo "Error: Failed to create build directory"; exit 1; }

# Build the model and example PDFs
build_pdf "${CLI_PDF_NAME}" "${CLI_DIR}" "${BUILD_DIR}" || { echo "Error: Failed to build model PDF"; exit 1; }

# build_pdf "${MODEL_PDF_NAME}" "${MODEL_DIR}" "${BUILD_DIR}" || { echo "Error: Failed to build model PDF"; exit 1; }
# build_pdf "${EXAMPLE_PDF_NAME}" "${EXAMPLE_DIR}" "${BUILD_DIR}" || { echo "Error: Failed to build example PDF"; exit 1; }

# Return to the original working directory
cd ${current_pwd}

# If the script reaches this point, it has executed successfully
echo "Script executed successfully."

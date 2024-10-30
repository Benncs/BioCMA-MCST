import os
import subprocess
import shutil

# Directory paths
BUILD_DIR = './build'
MODEL_DIR = './models'
EXAMPLE_DIR = './examples'
MODEL_PDF_NAME = 'model'
EXAMPLE_PDF_NAME = 'examples'

DOC_DIR = 'docs'
EXAMPLES_SCRIPT_DIR = 'examples'

def run_command(command, cwd=None):
    try:
        result = subprocess.run(command, cwd=cwd, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        return result.stdout.decode('utf-8')
    except subprocess.CalledProcessError as e:
        print(f"Error: {e.stderr.decode('utf-8')}")
        exit(1)

def build_pdf(pdf_name, source_dir, destination_dir):
    # Change to the specified directory
    os.chdir(source_dir)
    
    # Compile LaTeX document
    run_command("pdflatex main.tex")
    run_command("bibtex main")
    run_command("pdflatex main.tex")

    # Remove auxiliary files
    aux_files = ['*.aux', '*.bbl', '*.blg', '*.log', '*.lot', '*.out', '*.xml', '*.toc', '*.lof', '*-blx.bib']
    for aux_file in aux_files:
        run_command(f'rm -f {aux_file}')
    
    # Return to the original directory
    os.chdir('..')

    # Copy the generated PDF to the specified destination
    src_pdf = os.path.join(source_dir, "main.pdf")
    dest_pdf = os.path.join(destination_dir, f"{pdf_name}.pdf")
    try:
        shutil.copy(src_pdf, dest_pdf)
    except Exception as e:
        print(f"Error: Failed to copy PDF file: {e}")
        exit(1)

if __name__== '__main__':

    #TODO Continue script 

    # Store the current working directory
    current_pwd = os.getcwd()

    # Create the build directory if it doesn't exist
    os.makedirs(BUILD_DIR, exist_ok=True)

    # Build the model and example PDFs
    build_pdf(MODEL_PDF_NAME, MODEL_DIR, BUILD_DIR)
    # build_pdf(EXAMPLE_PDF_NAME, EXAMPLE_DIR, BUILD_DIR)

    # Return to the original working directory
    os.chdir(current_pwd)

    # If the script reaches this point, it has executed successfully
    print("Script executed successfully.")

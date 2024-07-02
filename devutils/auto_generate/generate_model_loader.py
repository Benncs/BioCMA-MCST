import os
import sys

def list_model_files(directory):
    """
    Lists model source files and headers in the given directory.
    """
    try:
        src_files = os.listdir(os.path.join(directory, "src"))
        
        header_files = os.listdir(os.path.join(directory, "public/models"))

        model_files = [file[:-4] for file in src_files if file.startswith("model_") and file.endswith(".cpp")]
        
        model_headers = [header for header in header_files if header.endswith(".hpp")]

        return model_files, model_headers

    except FileNotFoundError as e:
        print(f"Error: {e}")
        return [], []
    except PermissionError as e:
        print(f"Error: {e}")
        return [], []

def generate_includes(model_headers):
    """
    Generates include directives for the model headers.
    """
    return "\n".join([f"#include <models/{header}>" for header in model_headers])

def generate_function_body(model_files):
    """
    Generates the body of the load_model_ function.
    """
    function_body = ""
    for model_name in model_files:
        function_body += f"""    if (name == "{model_name}")
    {{
        return get_{model_name}();
    }}
    else """
    # Add the default else statement to handle unknown models
    function_body += """{
        throw std::runtime_error("Model not found: " + name);
    }"""

    return function_body

def generate_cpp_file(template_path, output_path, includes, body):
    """
    Generates the C++ file by replacing the placeholders in the template.
    """
    try:
        with open(template_path, 'r') as template_file:
            template_content = template_file.read()

        # Replace the placeholders
        content = template_content.replace("@INCLUDES@", includes)
        content = content.replace("@BODY@", body)

        # Write the modified content to the output file
        with open(output_path, 'w') as output_file:
            output_file.write(content)

    except FileNotFoundError as e:
        print(f"Error: {e}")
    except IOError as e:
        print(f"Error: {e}")

def generate_header(template_path,output_path):
    content = ""
    with open(template_path,'r') as file:
        content = file.read()

    with open(output_path,'w') as file:
         file.write(content)

if __name__ == "__main__":
   
    # Read command-line arguments
    args = sys.argv
    if len(args)!=6:
        raise Exception("Bad argument")

    models_path = args[1]
    source_template_path = args[2] 
    source_output_path = args[3] 

    header_template_path = args[4]
    header_output_path = args[5]

    files, headers = list_model_files(models_path)

    includes = generate_includes(headers)
    body = generate_function_body(files)

    # Generate the C++ file
    generate_cpp_file(source_template_path, source_output_path, includes, body)

    generate_header(header_template_path,header_output_path)
import os
import sys
from typing import Tuple,List

def to_camel_case(snake_str:str)->str:
    components = snake_str.split("_")
    return "".join(x.capitalize() for x in components[1:])


def list_model_files(directory:str)->Tuple[List[str],List[str]]:
    """
    Lists model source files and headers in the given directory.
    """
    try:
        src_files = os.listdir(os.path.join(directory, "src"))

        header_files = os.listdir(os.path.join(directory, "public/models"))

        model_files = [
            file[:-4]
            for file in src_files
            if file.startswith("model_") and file.endswith(".cpp")
        ]

        model_headers = [header for header in header_files if header.endswith(".hpp")]

        return model_files, model_headers

    except FileNotFoundError as e:
        print(f"Error: {e}")
        return [], []
    except PermissionError as e:
        print(f"Error: {e}")
        return [], []


def generate_includes(model_headers:List[str])->str:
    """
    Generates include directives for the model headers.
    """
    return "\n".join([f"#include <models/{header}>" for header in model_headers])


def generate_loader_body(model_files:List[str])->str:
    """
    Generates the body of the load_model_ function.
    """

    function_body = ""
    for i, model_name in enumerate(model_files):
        function_body += f"""
    case {i}:
    {{
        return MC::init<Models::{to_camel_case(model_name)}>(
            info, number_particle, liq_volume, liquid_neighbors, x0);
    }}
    """
    # Add the default else statement to handle unknown models
    # function_body += """{
    #     throw std::runtime_error("Model not found: " + name);
    # }"""

    return function_body


def generate_list_body(model_files:List[str])->str:
    function_body = ""
    for model_name in model_files:
        function_body += f'list.emplace_back("{model_name}");\r\n'
    return function_body

def generate_selection_body(model_files: List[str]) -> str:
    body = ""

    map_elements = ", ".join([f'{{ "{name}", {index} }}' for index, name in enumerate(model_files)])

    body += "static std::unordered_map<std::string_view, int> map{ {\"default\", -1},"+f"{map_elements}"+"};"

    return body

def generate_cpp_file(template_path, output_path, includes, body, map_selection):
    """
    Generates the C++ file by replacing the placeholders in the template.
    """
    try:
        with open(template_path, "r") as template_file:
            template_content = template_file.read()

        # Replace the placeholders
        content = template_content.replace("@INCLUDES@", includes)
        content = content.replace("@SWITCH_BODY@", body)
        # content = content.replace("@AM_BODY@", list_body)
        content = content.replace("@MODEL_INDEX_MAP@",map_selection)

        # Write the modified content to the output file
        with open(output_path, "w") as output_file:
            output_file.write(content)

    except FileNotFoundError as e:
        print(f"Error: {e}")
    except IOError as e:
        print(f"Error: {e}")


def generate_header(template_path, output_path):
    content = ""
    with open(template_path, "r") as file:
        content = file.read()

    with open(output_path, "w") as file:
        file.write(content)


def generate_variant(template_path: str, output_path: str, includes, model_files):
    try:
        with open(template_path, "r") as template_file:
            template_content = template_file.read()
            content = template_content.replace("@INCLUDES@", includes)
            body = "MC::ParticlesContainer<DefaultModel>,"
            for model in model_files:
                body += f"MC::ParticlesContainer<Models::{to_camel_case(model)}>,"

            body = body[:-1]

            content = content.replace("@VARIANT_TYPE@", body)

            

            with open(output_path, "w") as output_file:
                output_file.write(content)
    except FileNotFoundError as e:
        print(f"Error: {e}")
    except IOError as e:
        print(f"Error: {e}")


if __name__ == "__main__":
    # Read command-line arguments
    args = sys.argv
    if len(args) != 8:
        raise Exception("Bad argument")

    models_path = args[1]
    source_template_path = args[2]
    source_output_path = args[3]

    header_template_path = args[4]
    header_output_path = args[5]

    variant_template_path = args[6]
    variant_output_path = args[7]

    files, headers = list_model_files(models_path)

    includes = generate_includes(headers)
    loader_body = generate_loader_body(files)
    # list_body = generate_list_body(files)
    map_selection = generate_selection_body(files)

    # Generate the C++ file
    generate_cpp_file(
        source_template_path, source_output_path, includes, loader_body, map_selection
    )

    generate_header(header_template_path, header_output_path)

    generate_variant(variant_template_path, variant_output_path, includes, files)

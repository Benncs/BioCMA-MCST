import os
import sys
from typing import Tuple, List


def to_camel_case(snake_str: str) -> str:
    components = snake_str.split("_")
    return "".join(x.capitalize() for x in components)


def generate_includes(model_list: List[str]) -> str:
    """
    Generates include directives for the model headers.
    """
    return "\n".join([f"#include <models/{header}.hpp>" for header in model_list])


def generate_loader_body(model_files: List[str]) -> str:
    """
    Generates the body of the load_model_ function.
    """

    function_body = ""
    for i, model_name in enumerate(model_files):
        function_body += f"""
    case {i}:
    {{
        return MC::init<Models::{to_camel_case(model_name)}>(
             number_particle, liq_volume, liquid_neighbors,total_mass);
    }}
    """
    # Add the default else statement to handle unknown models
    # function_body += """{
    #     throw std::runtime_error("Model not found: " + name);
    # }"""

    return function_body


def generate_list_body(model_files: List[str]) -> str:
    function_body = ""
    for model_name in model_files:
        function_body += f'list.emplace_back("{model_name}");\r\n'
    return function_body


def generate_selection_body(model_files: List[str]) -> str:
    body = ""

    map_elements = ", ".join(
        [f'{{ "{name}", {index} }}' for index, name in enumerate(model_files)]
    )

    body += (
        'static std::unordered_map<std::string_view, int> map{ {"default", -1},'
        + f"{map_elements}"
        + "};"
    )

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
        content = content.replace("@MODEL_INDEX_MAP@", map_selection)

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


def generate_variant(
    template_path: str, output_path: str, includes, model_files, add_py_variant
):
    try:
        with open(template_path, "r") as template_file:
            template_content = template_file.read()
            content = template_content.replace("@INCLUDES@", includes)
            body = "MC::ParticlesContainer<DefaultModel>,"
            for model in model_files:
                body += f"MC::ParticlesContainer<Models::{to_camel_case(model)}>,"

            body = body[:-1]
            if add_py_variant:
                body += ",MC::ParticlesContainer<PythonWrap::PimpModel>"
            
     
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
    if len(args) < 11:
        print(args)
        raise Exception("Bad argument")



    models_folder_root = args[1]

    template_model_loader_path = args[2]
    template_model_loader_header = args[3]
    template_variant_model = args[4]
    source_output_path = args[5]
    header_output_path = args[6]
    variant_output_path = args[7]
    add_py_variant = args[8] != ""
    add_udf_variant = args[9] != ""

    model_list_name = args[10:]


    # files, headers = list_model_files(models_folder_root,add_udf_variant)

    includes = generate_includes(model_list_name)
    loader_body = generate_loader_body(model_list_name)
    # # list_body = generate_list_body(files)
    map_selection = generate_selection_body(model_list_name)

    # # Generate the C++ file
    generate_cpp_file(
        template_model_loader_path, source_output_path, includes, loader_body, map_selection
    )

    generate_header(template_model_loader_header, header_output_path)

    generate_variant(
        template_variant_model, variant_output_path, includes, model_list_name, add_py_variant
    )

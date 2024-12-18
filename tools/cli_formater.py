#!/usr/bin/python3

from lxml import etree
import sys

import os

current_file_path = os.path.abspath(__file__)
current_directory = os.path.dirname(current_file_path)


def create_xml_parser():
    XSD_PATH = current_directory + "/../devutils/datamodel/input_scheme.xsd"
    # Load XSD schema
    with open(XSD_PATH, "r") as xsd_file:
        xsd_content = xsd_file.read()

    # Create XML schema object
    schema = etree.XMLSchema(etree.XML(xsd_content))

    # Create parser with schema validation
    parser = etree.XMLParser(schema=schema)
    return parser


def read_xml_values(xml_path, parser, target_name):
    # Read XML content
    with open(xml_path, "rb") as xml_file:
        xml_content = xml_file.read()
    cli_args = ""

    try:
        # Parse the XML input
        root = etree.fromstring(xml_content, parser)

        # Find all control elements under cases
        controls = root.findall("control")

        for control in controls:
            name = control.get("name")
            if name == target_name:
                cma_case_path = control.findtext("cma_case_path")
                recursive = control.find("cma_case_path").get("recursive")
                final_time = float(control.findtext("final_time"))
                numper_particle = int(control.findtext("numper_particle"))
                delta_time_element = control.find("delta_time")
                delta_time = (
                    float(delta_time_element.text)
                    if delta_time_element is not None
                    else None
                )
                results_file_name_element = control.find("results_file_name")
                results_file_name = (
                    results_file_name_element.text
                    if results_file_name_element is not None
                    else None
                )
                number_exported_result = int(control.findtext("number_exported_result"))
                model_name = control.findtext("model_name")
                initialiser_path = control.findtext("initialiser_path")
                cli_args = ""

                if recursive:
                    cli_args += "-r 1 "
                cli_args += f"-f {cma_case_path} "

                cli_args += f"-d {final_time} "
                cli_args += f"-np {numper_particle} "
                if delta_time > 0:
                    cli_args += f"-dt {delta_time} "

                if results_file_name != "" and results_file_name is not None:
                    cli_args += f"-er {results_file_name} "

                if number_exported_result > 0:
                    cli_args += f"-nex {number_exported_result} "

                if initialiser_path != "" and initialiser_path is not None:
                    cli_args += f"-fi {initialiser_path} "

                cli_args += f"-mn {model_name} "
                break
        if cli_args == "":
            raise Exception("Case not found")
        return cli_args, results_file_name

    except etree.XMLSyntaxError as e:
        print(f"XML Syntax Error: {str(e)}")
    except etree.DocumentInvalid as e:
        print(f"Invalid XML: {str(e)}")


def format_cli(args, file="cases.xml"):
    if (len(args)) == 2:
        name = args[1]
        xml_file_path = current_directory + "/" + file

        # Create XML parser with schema validation
        parser = create_xml_parser()

        # Read XML and extract values
        return read_xml_values(xml_file_path, parser, name)
    else:
        return "h 1"


if __name__ == "__main__":
    args = sys.argv

    cli = format_cli(args)[0]
    print(cli)

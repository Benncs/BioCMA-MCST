#!/usr/bin/python3

from lxml import etree


def create_xml_parser():
    XSD_PATH = "devutils/datamodel/input_scheme.xsd"
    # Load XSD schema
    with open(XSD_PATH, 'r') as xsd_file:
        xsd_content = xsd_file.read()

    # Create XML schema object
    schema = etree.XMLSchema(etree.XML(xsd_content))

    # Create parser with schema validation
    parser = etree.XMLParser(schema=schema)
    return parser

def read_xml_values(xml_path, parser):
    # Read XML content
    with open(xml_path, 'rb') as xml_file:
        xml_content = xml_file.read()

    try:
        # Parse the XML input
        root = etree.fromstring(xml_content, parser)

        # Read values from XML
        cma_case_path = root.findtext("cma_case_path")
        recursive = root.find("cma_case_path").get("recursive")
        final_time = float(root.findtext("final_time"))
        numper_particle = int(root.findtext("numper_particle"))
        delta_time = float(root.findtext("delta_time"))
        results_file_name = root.findtext("results_file_name")
        number_exported_result = int(root.findtext("number_exported_result"))
        model_name = root.findtext("model_name")

        cli_args = ""
        if recursive:
            cli_args+="-r 1 "
        cli_args+= f"-f { cma_case_path} "

        cli_args+= f"-d {final_time} "
        cli_args+= f"-np {numper_particle} "
        if(delta_time>0):
            cli_args+= f"-dt {delta_time} "
        if(results_file_name!=""):
            cli_args+= f"-er {results_file_name} "
        if(number_exported_result>0):
            cli_args += f"-nex {number_exported_result} "
        cli_args+=f"-mn {model_name}"

        print(cli_args)

    except etree.XMLSyntaxError as e:
        print(f"XML Syntax Error: {str(e)}")
    except etree.DocumentInvalid as e:
        print(f"Invalid XML: {str(e)}")


xml_file_path = "/mnt/c/Users/casale/Documents/code/cpp/biomc/devutils/datamodel/input_scheme.xml"  

# Create XML parser with schema validation
parser = create_xml_parser()

# Read XML and extract values
read_xml_values(xml_file_path, parser)

import re
import os 
current_file_path = os.path.abspath(__file__)
current_directory = os.path.dirname(current_file_path)


header_template_path = f"{current_directory}/template_model.hxx"
source_template_path = f"{current_directory}/template_model.cxx"
model_name = "model_test"




def validate_model_name(name: str) -> bool:
    
    pattern = r'^model_\w+$'
    
    if re.match(pattern, name):
        return True
    else:
        return False


def generate_model_header(name):
    content = ""
    with open(header_template_path,'r') as file:
        content = file.read()

    content = content.replace('@__model__name__@',name)
    
    return content

def generate_model_source(name):
    content = ""
    with open(source_template_path,'r') as file:
        content = file.read()

    content = content.replace('@__model__name__@',name)
    return content



if __name__=="__main__":

    if validate_model_name(model_name):
        header = generate_model_header(model_name)
        source = generate_model_source(model_name)

        with open(f"{model_name}.cpp",'w') as file:
            file.write(source)

        with open(f"{model_name}.hpp",'w') as file:
            file.write(header)
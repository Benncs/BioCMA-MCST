import os

def check_correct_name(name:str)->bool:
    print("TODO")
    return True

def get_root_sim():
    if os.environ.get("BIOMC_SIMFOLDER"):
        print("ee")
        return os.environ["BIOMC_SIMFOLDER"]
    else:
#        homedir=os.path.expanduser('~')
#        return f"{homedir}/biomc_simulation"
        return "/tmp/biomc_simulation"

if __name__ == "__main__":
    simulation_name=input("Simulation name: ")
    simulation_folder_root = get_root_sim()
    if not check_correct_name(simulation_name):
        print("Incorrect name, exiting..")
        exit(-1)

    out_folder = f"{simulation_folder_root}/{simulation_name}"
    os.makedirs(out_folder, exist_ok=False)

    os.chdir(out_folder)
    
    

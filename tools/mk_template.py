import os
import shutil
import tarfile
from datetime import datetime
from string import Template
import subprocess
import platform

ARCHIVE_PATH = "./tool/template.tar"


def get_user_name():
    def get_git_username():
        try:
            return (
                subprocess.check_output(["git", "config", "user.name"])
                .decode("utf-8")
                .strip()
            )
        except (subprocess.CalledProcessError, FileNotFoundError):
            return None

    def get_system_username():
        if platform.system() == "Windows":
            return os.getlogin()
        else:
            import pwd

            return pwd.getpwuid(os.getuid()).pw_name

    return get_git_username() or get_system_username()


def check_correct_name(name: str) -> bool:
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


def mk_template(data):
    with tarfile.open("template.tar", "r:*") as tar:
        tar.extractall(path="./")
    os.remove("template.tar")

    with open("run_simulation.py", "r", encoding="utf-8") as f:
        content = f.read()

    template = Template(content)
    result = template.substitute(data)
    print(result)

    with open("run_simulation.py", "w", encoding="utf-8") as f:
        f.write(result)


if __name__ == "__main__":
    simulation_folder = input("Simulation folder: ")
    simulation_folder_root = get_root_sim()
    if not check_correct_name(simulation_folder):
        print("Incorrect name, exiting..")
        exit(-1)

    simulation_name = input("Simulation name (default same as folder): ")
    if not simulation_name:
        simulation_name = simulation_folder
    user_description = input("Simulation description: ")
    out_folder = f"{simulation_folder_root}/{simulation_folder}"
    os.makedirs(out_folder, exist_ok=False)
    shutil.copy("./tools/template.tar", out_folder)
    os.chdir(out_folder)
    data = {
        "SIMULATION_NAME": simulation_name,
        "DESCRIPTION": user_description,
        "DATE": datetime.today().strftime("%Y-%m-%d %H:%M:%S"),
        "AUTHOR": get_user_name(),
        "SIMFILENAME": simulation_folder,
    }

    mk_template(data)

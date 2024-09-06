import multiprocessing
from cli_formater import format_cli
from runner import get_executable, exec
from biomc.__main__ import main
import subprocess
import os
import sys
from contextlib import contextmanager

@contextmanager
def silence_stdout():
    old_target = sys.stdout
    try:
        with open(os.devnull, "w") as new_target:
            sys.stdout = new_target
            yield new_target
    finally:
        sys.stdout = old_target

file_header = """<?xml version="1.0" encoding="UTF-8"?>
<cases>
"""


duration = [8*3600]

rd = ""

for i, d in enumerate(duration):
    raw_case = f"""
    <control name="0d_{i}">
        <cma_case_path>./cma_data/0d/</cma_case_path>
        <final_time>{d}</final_time>
        <numper_particle>1000</numper_particle>
        <delta_time>0.1</delta_time>
        <results_file_name>batch_fast_{i}</results_file_name>
        <number_exported_result>30</number_exported_result>
        <model_name>model_monod</model_name> 
        <initialiser_path>./cma_data/0d_init.h5</initialiser_path>
    </control>"""
    rd += raw_case


full_file = file_header + rd + " </cases>"

with open("./tools/batch.xml", "w") as file:
    file.write(full_file)

def run_command(i):

    cli,res_path = format_cli(["_", f"0d_{i}"], "batch.xml")
    print(f"Executing: {cli}")
    with silence_stdout():
        command = get_executable("release", False) + " " + cli + " -nt 6"
        exec(command, '6',stdout=subprocess.PIPE)
        main([res_path])


def _main():
    with multiprocessing.Pool(processes=2) as pool:
        pool.map(run_command, range(len(duration)))

if __name__=="__main__":
        _main()
        print("All commands executed.")
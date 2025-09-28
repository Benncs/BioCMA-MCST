import subprocess
import re
import time
import os
import signal


def format_rhs(match):
    """
    Function to format RHS values in green.

    Args:
        match (re.Match): Match object containing the regex match.

    Returns:
        str: Formatted string with ANSI escape codes for green text.
    """
    return f"{match.group(1)}\033[92m{match.group(2)}\033[0m"


def wrap_timer(f, do_measure: bool):
    try:

        old_action = signal.signal(signal.SIGINT, signal.SIG_IGN)
        if do_measure:
            start_time = time.perf_counter()
            process = f()
            return_code = process.wait()
            end_time = time.perf_counter()
            elapsed_time = end_time - start_time
            print(f"Command executed in \033[92m{elapsed_time:.6f}\033[0m seconds")
            return return_code
        else:
            try:
                process = f()
            finally:
                # Restore original SIGINT handling in the parent
                signal.signal(signal.SIGINT, old_action)
            return process.wait()
    except KeyboardInterrupt:
        pass


def exec(command, n_thread, do_measure: bool = True, do_kokkos_measure=False, **kwargs):
    env_var = os.environ.copy()
    env_var["OMP_PLACES"] = "threads"
    env_var["OMP_PROC_BIND"] = "spread"
    env_var["OMP_NUM_THREADS"] = n_thread
    if do_kokkos_measure:
        # env_var["KOKKOS_TOOLS_LIBS"] = "/usr/local/lib64/libkp_memory_events.so"
        # env_var["KOKKOS_TOOLS_LIBS"]="/usr/local/lib64/libkp_kernel_timer.so"
        env_var["KOKKOS_TOOLS_LIBS"]="/usr/local/lib/libkp_kernel_timer.so" #libkp_memory_usage

    result = command.replace("-", "\n-")
    pattern = re.compile(r"(-\w+\s)(\S+)")
    formatted_command = pattern.sub(format_rhs, result)
    print("\r\n")
    print(formatted_command)
    print("\n")
    return wrap_timer(
        lambda: subprocess.Popen(command, shell=True, env=env_var, **kwargs), do_measure
    )

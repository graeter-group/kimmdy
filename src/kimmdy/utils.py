import subprocess as sp

def run_shell_cmd(s):
    sp.run(s, shell=True)

def get_shell_stdout(s):
    process = sp.run(s, shell=True, capture_output=True, encoding='utf-8')
    return process.stdout

def check_gmx_version():
    return get_shell_stdout('gmx --quiet --version')

import subprocess as sp
import os

def get_git_hash(path):
    if os.path.isfile(path):
        path = os.path.dirname(path)
    try:
        stdout, stderr = sp.Popen('git rev-parse HEAD', cwd=path, 
                                  shell=True, stdout=sp.PIPE).communicate()
    except OSError:
        return ''
    else:
        return stdout.strip()
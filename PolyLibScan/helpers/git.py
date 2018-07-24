import subprocess as sp
import os
from functools import partial

def git_interface(path, cli):
    if os.path.isfile(path):
        path = os.path.dirname(path)
    try:
        stdout, stderr = sp.Popen(cli, cwd=path,
                                  shell=True, stdout=sp.PIPE).communicate()
    except OSError:
        return ''
    else:
        return stdout.strip()

get_git_diff = partial(git_interface, cli='git diff')
get_git_hash = partial(git_interface, cli='git rev-parse HEAD')

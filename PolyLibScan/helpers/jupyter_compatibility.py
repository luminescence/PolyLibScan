from tqdm import tqdm, tqdm_notebook


def isnotebook():
    try:
        shell = get_ipython().__class__.__name__
        # Jupyter notebook or qtconsole
        if shell == 'ZMQInteractiveShell':
            return True
        # Terminal running IPython
        elif shell == 'TerminalInteractiveShell':
            return False
        # Other type (?)
        else:
            return False
    except NameError:
        # Probably standard Python interpreter
        return False


def agnostic_tqdm(*args, **kwargs):
    """show progress bar based on tqdm. Choose which tqdm to use based on Kernel (e.g. python, ipython or
    jupyter-notebook and pass all arguments"""

    if isnotebook():
        sel_func = tqdm_notebook    # only works within notebook
    else:
        sel_func = tqdm     # only works outside notebook

    return sel_func(*args, **kwargs)   # all args and kwargs are preserved

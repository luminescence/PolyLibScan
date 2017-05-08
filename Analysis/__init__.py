import misc.git as _git

__git_hash__ = _git.get_git_hash(__file__)
__version__ = "0.1.0"

from job import Project, Job, Run
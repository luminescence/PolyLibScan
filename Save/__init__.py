import misc.git as _git
import sys

this = sys.modules[__name__]
this.__git_hash__ = _git.get_git_hash(__file__)

from JobSave import JobSave
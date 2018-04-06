import unittest as ut

from PolyLibScan.helpers.git import get_git_hash, git_interface


class Test_git(ut.TestCase):
    """Test git interface"""

    def __init__(self, *args, **kwargs):
        super(Test_git, self).__init__(*args, **kwargs)

    def test_git_hash(self):
        git_hash = get_git_hash('.')
        self.assertEqual(len(git_hash), 40)

    def test_interface(self):
        """check general function of interface, e.g. if branch master is present"""
        branches = git_interface('.', 'git branch').split('\n')
        self.assertTrue('  master' in branches)

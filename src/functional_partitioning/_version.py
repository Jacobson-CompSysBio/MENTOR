# # Error during installation (pip install -e .):
# # AttributeError: module 'functional_partitioning' has no attribute '__version__'. Did you mean: '_version'?
# def get_version():
#     import os
#     import subprocess
#     cwd = os.path.getcwd()
#     os.chdir(os.path.dirname(os.path.abspath(__file__)))
#     version = subprocess.check_output(
#         'git describe --tags --dirty --always'.split()
#     ).decode().strip()
#     os.chdir(cwd)
#     return version

def get_version():
    return 'v0.5.1'

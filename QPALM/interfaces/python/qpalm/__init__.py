"""Proximal Augmented Lagrangian method for Quadratic Programs"""

__version__ = "1.2.6"

import os
import typing

def _is_truthy(s: typing.Optional[str]):
    if s is None: return False
    return not s.lower() in ('', 'false', 'no', 'off', '0')

if not typing.TYPE_CHECKING and _is_truthy(os.getenv('QPALM_PYTHON_DEBUG')):
    from qpalm._qpalm_d import *
    from qpalm._qpalm_d import __version__ as c_version
else:
    from qpalm._qpalm import *
    from qpalm._qpalm import __version__ as c_version
assert __version__ == c_version

del _is_truthy, typing, os
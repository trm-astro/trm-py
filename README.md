# Python Package trm-py

This is the python wrapper around trm-subs and some other trm c++ software.

This is the update to many of the repositories found [here](https://github.com/trmrsh?tab=repositories&q=&type=&language=python&sort=). The code will create a drop in replacement via the form `from trm_py import trm_<name>` or `from trm_py.trm_<name> import <name>`, `import trm_py` will also be supported with tool access via `trm_py.<name>.<tool>`.

This code will be installable via `pip install trm-py` from PyPi (and potentially via Conda). This code will not require additional C++ libaray installations as it will come formost as a `bdist` for most common platforms and a `sdist` to support build and installation on uncommon platforms.

Until further notice this is a WIP and should not be considered usable feor scientific use.

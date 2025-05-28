# This file is thin wrappers around the lcurve executables
# to allow them to be installed via the trm-py package

import sys
from .._cpp import _cpp_lcurve


# def lprofile():
#     call_tool("lprofile")


# def levmarq():
#     call_tool("levmarq")


def lroche():
    _cpp_lcurve.lroche(sys.argv[1:])


# def picture():
#     call_tool("picture")


# def simplex():
#     call_tool("simplex")


# def visualise():
#     call_tool("visualise")


# This file is thin wrappers around the lcurve executables
# to allow them to be installed via the trm-py package

import sys
import trm_py
from pathlib import Path

# Get the location of the trm_py package
trm_py_path = Path(trm_py.__path__[0], "_cpp")


def call_tool(tool_name, tool_path=trm_py_path):
    """
    Call a tool of a given name with the system arguments
    tool_name: name of the tool
    tool_path: path to the tool (default is the trm_py module path with trm_py/_cpp)
    """
    import os
    import subprocess

    # Get the path to the tool
    tool_path = tool_path/tool_name

    # Check if the tool exists
    if not os.path.exists(tool_path):
        print(f"Tool {tool_name} does not exist in {tool_path}")
        return

    # Call the tool with the system arguments minus the python script name
    args = sys.argv[1:]
    # Check if the tool is executable
    if not os.access(tool_path, os.X_OK):
        print(f"Tool {tool_name} is not executable")
        return
    # make the tool path absolute using pathlib
    tool_path = tool_path.resolve()

    # convert to a string and add the arguments
    tool_path = str(tool_path)
    command = [tool_path, *args]

    subprocess.check_call(command)


def lprofile():
    call_tool("lprofile")


def levmarq():
    call_tool("levmarq")


def lroche():
    call_tool("lroche")


def lroches():
    call_tool("lroches")


def picture():
    call_tool("picture")


def simplex():
    call_tool("simplex")


def visualise():
    call_tool("visualise")


# System to build the trm c++ package required by trm_py
from pybind11.setup_helpers import Pybind11Extension, build_ext


def build(setup_kwargs):
    ext_modules = [
        Pybind11Extension("trm_py", ["src/main.cpp"]),
    ]
    setup_kwargs.update(
        {
            "ext_modules": ext_modules,
            "cmdclass": {"build_ext": build_ext},
            "zip_safe": False,
        }
    )

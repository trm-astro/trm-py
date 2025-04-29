import os
import subprocess
import site
import platform
import sys


def run_command(command, cwd=None):
    """Run a shell command and raise an error if it fails."""
    print(f"Running command: {' '.join(command)}")
    subprocess.check_call(command, cwd=cwd)


def run_conan(source_dir):
    # Use Conan to install the dependencies
    print(f"Installing dependencies in {source_dir}")
    run_command([
        'conan',
        'install',
        source_dir,
        '-build=missing',
    ])
    print("Dependencies installed")

    # Use Conan to create a build profile
    print("Forcing conan profile")
    run_command([
        'conan',
        'profile',
        'detect',
        '--force',
    ])
    print("Conan profile created")


def run_cmake(system, arch, source_dir, build_dir):
    # Get the Python site-packages directory
    python_site_dir = site.getsitepackages()[0]

    # Configure the project using CMake
    print("Configuring and building the project")

    cmake_toolchain_file = os.path.join(build_dir,
                                        'Release',
                                        'generators',
                                        'conan_toolchain.cmake')

    cmake_args = [
        'cmake',
        '--preset', 'conan-release',
        f'-DCMAKE_TOOLCHAIN_FILE={cmake_toolchain_file}',
        '-DPLPLOT_BUILD_TYPE=4',
        f'-DPYTHON_SITE_PACKAGES={python_site_dir}',
        f'-S {source_dir}',
        f'-B {build_dir}',
    ]

    if system == 'windows':
        cmake_args.append('-G')
        cmake_args.append('NMake Makefiles')  # Use NMake for Windows builds

    run_command(cmake_args)
    print("Project configured")

    print("Building the project")
    # Build the project using CMake
    run_command([
        'cmake',
        '--build', build_dir,
        '--config', 'Release',
    ])
    print("Project built")

    # Install the built shared object file to the py_subs directory
    print("Installing the built shared object file")
    run_command([
        'cmake',
        '--install', build_dir,
        '--prefix', source_dir,
    ])
    print("Installation complete")


def main():
    # Detect the platform and architecture
    system = platform.system().lower()
    arch = platform.machine().lower()
    print(f"Building on {system} ({arch})")

    # Define the source and build directories
    source_dir = os.path.dirname(os.path.abspath(__file__))
    build_dir = os.path.join(source_dir, 'build')
    # Create the build directory if it doesn't exist
    os.makedirs(build_dir, exist_ok=True)

    # Run Conan to install dependencies
    run_conan(source_dir)

    # Run CMake to configure and build the project
    run_cmake(system, arch, source_dir, build_dir)


if __name__ == "__main__":
    try:
        main()
    except subprocess.CalledProcessError as e:
        print(f"An error occurred: {e}")
        sys.exit(1)

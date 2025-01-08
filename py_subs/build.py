import os
import subprocess
import site


def main():
    # Define the source and build directories
    source_dir = os.path.dirname(os.path.abspath(__file__))
    build_dir = os.path.join(source_dir, 'build')
    # Create the build directory if it doesn't exist
    os.makedirs(build_dir, exist_ok=True)

    # Use Conan to install the dependencies
    print(f"Installing dependencies in {source_dir}")
    subprocess.check_call([
        'conan',
        'install',
        source_dir,
        '-build=missing',
    ])
    print("Dependencies installed")

    # Use Conan to create a build profile
    print("Forcing conan profile")
    subprocess.check_call([
        'conan',
        'profile',
        'detect',
        '--force',
    ])
    print("Conan profile created")

    python_site_dir = site.getsitepackages()[0]

    # Configure the project using CMake
    print("Configuring and building the project")
    subprocess.check_call([
        'cmake',
        '--preset', 'conan-release',
        '-DCMAKE_TOOLCHAIN_FILE=' + os.path.join(build_dir, 'Release/generators/conan_toolchain.cmake'),
        '-DPLPLOT_BUILD_TYPE=4',
        '-DPYTHON_SITE_PACKAGES=' + python_site_dir,
        '-S', source_dir,
        '-B', build_dir,
    ])
    print("Project configured")

    # Build the project using CMake
    subprocess.check_call([
        'cmake',
        '--build',
        build_dir,
    ])
    print("Project built")

    # Install the built shared object file to the py_subs directory
    print("Installing the built shared object file")
    subprocess.check_call([
        'cmake',
        '--install', build_dir,
        '--prefix', f'{source_dir}/',
    ])
    print("Installation complete")


if __name__ == "__main__":
    main()

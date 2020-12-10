import os
import sys
import subprocess
import argparse
import shutil


def build(sourcedir, builddir, target=[], debug=False, toolchain=None, cflags=[]):
    cfg = 'Debug' if debug else 'Release'
    cmake_args = ['-DCMAKE_BUILD_TYPE=' + cfg]
    if not toolchain is None:
        cmake_args += ['-DCMAKE_TOOLCHAIN_FILE=' + toolchain]
    build_args = ['--config', cfg]
    build_args += ['--', '-j']

    for cflag in cflags:
        cmake_args += ["-D"+cflag[0]+"="+cflag[1]]
    print(cmake_args)

    os.makedirs(builddir, exist_ok=True)

    subprocess.check_call(['cmake', sourcedir] +
                          cmake_args, cwd=builddir)
    subprocess.check_call(['cmake', '--build', '.'] +
                          build_args + [target], cwd=builddir)  # basically make



# get command line arguments
ap = argparse.ArgumentParser()
ap.add_argument("targets", nargs='?',
                help="...")
ap.add_argument("-d", "--debug", action='store_true',
                help="...")
ap.add_argument("-b", "--build", default="1",
                help="...")
ap.add_argument("-p", "--with_parallel", default="1",
                help="...")
ap.add_argument("-c", "--clean", action='store_true',
                help="...")
ap.add_argument("-r", "--run", default="1",
                help="...")
# ap.add_argument("-t", "--target", default="micro",
#                 help="...")
ap.add_argument('extra_args', nargs='*')
# TODO add a no compilation flag
args = vars(ap.parse_args())
args['build'] = args['build'] != "0"
args['run'] = args['run'] != "0"
args['with_parallel'] = args['with_parallel'] != "0"

if args['targets'] is None or len(args['targets']) == 0:
    target = "mesh2yarns"
else:
    target = args['targets'] # not tested on multiple targets

# BUILD
sourcedir = os.path.join(os.getcwd(),"src")  # where the main CMakeLists.txt is
builddir = os.path.join(
    os.getcwd(), "build", "build-Debug" if args['debug'] else "build-Release")

if not args['with_parallel']:
    builddir += "_NP"

if args['clean']:
    shutil.rmtree(builddir)

cflags=[]
if not args['with_parallel']:
    cflags.append(["NO_PARALLEL", "TRUE"])

toolchain = "../vcpkg/scripts/buildsystems/vcpkg.cmake" # relative to build dir
if args['build']:
    build(sourcedir=sourcedir, builddir=builddir, target=target,
        debug=args['debug'], toolchain=toolchain, cflags=cflags)

# RUN
if args['run']:
    workdir = os.getcwd() # run from parent
    executable = [os.path.join(builddir, target)] + args['extra_args']

    print("Executing:", executable)
    try:
        subprocess.run(executable, cwd=workdir)
    except subprocess.CalledProcessError as exc:
        print("PY: CalledProcessError")
        print(exc.returncode, exc.output, exc.stderr)
    except KeyboardInterrupt:
        print("PY: Aborting execution")

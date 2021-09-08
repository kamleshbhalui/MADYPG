# ensure the script runs from the MADYPG directory, even if called from within the script subdirectory
# by moving to the parent directory of this script
cd "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")/.."
python exec.py mesh2yarns 9
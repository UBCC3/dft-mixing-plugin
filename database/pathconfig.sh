
PSIPATH="/usr/local/psi4/lib"
PLUGINPATH=$(dirname "${PWD}")

export PYTHONPATH="${PLUGINPATH}:${PSIPATH}:${PWD}/scripts:${PYTHONPATH}"

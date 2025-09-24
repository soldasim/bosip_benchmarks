
This repo contains benchmarks of the BOSIP.jl method.

# Installation

**OUTDATED: This complex installation is no longer needed. Python is no longer used anywhere. :)** Just install julia and use the /src environment.

The installation and the startup are quite complex due to using PyCall.jl and PythonCall.jl simultaneously. PythonCall/JuliaCall are neede for some SBIBM tasks. PyCall.jl is much nicer to use for me.

PC:
- `pyenv install 3.10.16`
- `pyenv local 3.10.16`

Cluster:
- `ml Python/3.10`

Then:
- `python -m venv venv`
- `. venv/bin/activate`
- `pip install --upgrade pip`
- `pip install -r requirements.txt`
- manually edit `/venv/lib/python3.10/site-packages/diffeqtorch/diffeqtorch.py` - see the TODOs
- python: `import juliacall` -> creates `/venv/julia_env` used by PythonCall/JuliaCall
- (use the created `/venv/julia_env` environment as the main Julia environment for the project)
- `julia`
- `] activate venv/julia_env`
- `ENV["PYTHON"] = Sys.which("python")` # make sure PyCall.jl uses the correct python to during build
- `] add PyCall, PythonCall`
- `] build`
- `] add ...`

# Startup (on every run)

**OUTDATED: This complex startup is no longer needed. Python is no longer used anywhere. :)** Just run julia and use the /src environment.

Evaluate the script `sh julia.sh`, which will set up everything (including activating all environments) and start a Julia REPL.

Alternatively, do the set up manually, in which case do the following (in this order):
- `export PYTHON_JULIACALL_HANDLE_SIGNALS=yes` - something to avoid segfaults with JuliaCall
- `. venv/bin/activate` - activate python env
- `julia`
- `] activate venv/julia_env` - activate julia env
- `ENV["PYCALL_JL_RUNTIME_PYTHON"] = Sys.which("python")` - set python for PyCall.jl
- `using PyCall` - only then import
- (check `PyCall.pyimport("sys").executable`)
- `ENV["JULIA_CONDAPKG_BACKEND"] = "Null"` - tell PythonCall.jl to not use custom conda env
- `ENV["JULIA_PYTHONCALL_EXE"] = "@PyCall"` - tell PythonCall.jl to use the same python as PyCall.jl
- `using PythonCall`
- (check `PythonCall.pyimport("sys").executable`)
- `PyCall.pyimport(\"juliacall\")` - import juliacall (into python)

# MCP service to run [`PLMC`](https://github.com/debbiemarkslab/plmc)

## Installation

### Quick Setup (Recommended)

Run the automated setup script:

```bash
cd plmc_mcp
bash quick_setup.sh
```

The script will create the conda environment, clone the PLMC repository, build the binaries, install all dependencies, and display the Claude Code configuration. See `quick_setup.sh --help` for options like `--skip-env` or `--skip-repo`.

### Manual Installation (Alternative)

```bash
# Create environment
mamba env create -p ./env python=3.10 pip
mamba activate ./env

# EVcouplings
pip install https://github.com/debbiemarkslab/EVcouplings/archive/develop.zip

# Install plmc
git clone https://github.com/debbiemarkslab/plmc.git
cd plmc

# or cd repo/plmc
make all-openmp
mamba install -c conda-forge -c bioconda hhsuite

pip install fastmcp pydantic
```

## Local usage
```shell
PLMC_DIR=./repo
DATASET=notebooks/example

# reformat a3m to a2m (a3m is the mmseqs2 output format)
reformat.pl a3m a2m  $DATASET/DHFR.a3m $DATASET/DHFR.a2m
python notebooks/rm_a2m_query_gaps.py $DATASET/DHFR.a2m $DATASET/alignment.a2m

$PLMC_DIR/plmc/bin/plmc \
    -o $DATASET/plmc/uniref100.model_params \
    -c $DATASET/plmc/uniref100.EC \
    -f $PROTEIN \
    -le 16.2 -lh 0.01 -m 200 -t 0.2 \
    -g $DATASET/alignment.a2m

```

## MCP usage
### Install mcp
```shell
# Install `plmc` mcp
fastmcp install claude-code tool-mcps/plmc_mcp/src/server.py --python tool-mcps/plmc_mcp/env/bin/python
```
## Call MCP

```markdown
I have created a a3m file for subtilisin BPN' in file @examples/case2.1_subtilisin/subtilisin.a3m. Can you help build a ev model using plmc mcp and create it to @examples/case2.1_subtilisin/plmc directory. The wild-type sequence is @examples/case2.1_subtilisin/wt.fasta.

Please convert the relative path to absolution path before calling the MCP servers.
```

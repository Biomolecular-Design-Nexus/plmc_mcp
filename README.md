# MCP service to run [`PLMC`](https://github.com/debbiemarkslab/plmc)

## Create environment
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
# reformat a3m to a2m (a3m is the mmseqs2 output format)
reformat.pl a3m a2m  example/DHFR.a3m example/DHFR.a2m

```

## MCP usage
```markdown
- a2m path: /home/xux/Desktop/ProteinMCP/ProteinMCP/mcp-servers/plmc_mcp/notebooks/example/DHFR.a2m
- focus_seq_id: query
- output_dir:/home/xux/Desktop/ProteinMCP/ProteinMCP/mcp-servers/plmc_mcp/notebooks/example/plmc
```

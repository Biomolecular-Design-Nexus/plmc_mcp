"""
Model Context Protocol (MCP) for plmc_mcp

Pseudolikelihood Maximization for Coevolution Analysis (PLMC) tools for EV+Onehot fitness prediction modeling.

This MCP Server provides tools to generate PLMC model parameters and evolutionary couplings
needed as inputs for EV+Onehot fitness prediction models.

Tools:
1. plmc_generate_model: Generate PLMC model parameters and evolutionary couplings from A2M alignment
   - Outputs: model_params (binary), EC file (text)
   - Uses preset parameters optimized for EV+Onehot modeling (lambda_e=16.2, lambda_h=0.01, m=200, theta=0.2)

2. plmc_convert_a3m_to_a2m: Convert A3M alignment to A2M format and clean query gaps
   - Required preprocessing step as PLMC requires A2M format
   - Removes positions where query sequence has gaps to prevent downstream modeling issues

Reference: https://github.com/debbiemarkslab/plmc
"""

from fastmcp import FastMCP
import sys
from pathlib import Path

# Add current directory to path for imports
current_dir = Path(__file__).parent
sys.path.insert(0, str(current_dir))

# Import statements (alphabetical order)
from tools.readme import readme_mcp

# Server definition and mounting
mcp = FastMCP(name="plmc_mcp")
mcp.mount(readme_mcp)

if __name__ == "__main__":
    mcp.run()
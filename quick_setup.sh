#!/bin/bash
#===============================================================================
# PLMC MCP Quick Setup Script
#===============================================================================
# This script sets up the complete environment for PLMC MCP server.
# Evolutionary coupling analysis using PLMC for protein fitness prediction.
#
# After cloning the repository, run this script to set everything up:
#   cd plmc_mcp
#   bash quick_setup.sh
#
# Once setup is complete, register in Claude Code with the config shown at the end.
#
# Options:
#   --skip-env        Skip conda environment creation
#   --skip-plmc       Skip PLMC compilation
#   --skip-hhsuite    Skip HHsuite installation
#   --help            Show this help message
#===============================================================================

set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ENV_DIR="${SCRIPT_DIR}/env"
PYTHON_VERSION="3.10"
REPO_DIR="${SCRIPT_DIR}/repo"
PLMC_REPO="https://github.com/debbiemarkslab/plmc.git"

# Print banner
echo -e "${BLUE}"
echo "=============================================="
echo "       PLMC MCP Quick Setup Script           "
echo "=============================================="
echo -e "${NC}"

# Helper functions
info() { echo -e "${BLUE}[INFO]${NC} $1"; }
success() { echo -e "${GREEN}[SUCCESS]${NC} $1"; }
warn() { echo -e "${YELLOW}[WARNING]${NC} $1"; }
error() { echo -e "${RED}[ERROR]${NC} $1"; }

# Check for conda/mamba
check_conda() {
    if command -v mamba &> /dev/null; then
        CONDA_CMD="mamba"
        info "Using mamba (faster package resolution)"
    elif command -v conda &> /dev/null; then
        CONDA_CMD="conda"
        info "Using conda"
    else
        error "Neither conda nor mamba found. Please install Miniconda or Mambaforge first."
        exit 1
    fi
}

# Parse arguments
SKIP_ENV=false
SKIP_PLMC=false
SKIP_HHSUITE=false

while [[ $# -gt 0 ]]; do
    case $1 in
        --skip-env) SKIP_ENV=true; shift ;;
        --skip-plmc) SKIP_PLMC=true; shift ;;
        --skip-hhsuite) SKIP_HHSUITE=true; shift ;;
        -h|--help)
            echo "Usage: ./quick_setup.sh [OPTIONS]"
            echo ""
            echo "Options:"
            echo "  --skip-env        Skip conda environment creation"
            echo "  --skip-plmc       Skip PLMC compilation"
            echo "  --skip-hhsuite    Skip HHsuite installation"
            echo "  -h, --help        Show this help message"
            exit 0
            ;;
        *) warn "Unknown option: $1"; shift ;;
    esac
done

# Check prerequisites
info "Checking prerequisites..."
check_conda

if ! command -v git &> /dev/null; then
    error "git is not installed. Please install git first."
    exit 1
fi
success "Prerequisites check passed"

# Step 1: Create conda environment
echo ""
echo -e "${BLUE}Step 1: Setting up conda environment${NC}"

if [ "$SKIP_ENV" = true ]; then
    info "Skipping environment creation (--skip-env)"
elif [ -d "$ENV_DIR" ] && [ -f "$ENV_DIR/bin/python" ]; then
    info "Environment already exists at: $ENV_DIR"
else
    info "Creating conda environment with Python ${PYTHON_VERSION}..."
    $CONDA_CMD create -p "$ENV_DIR" python=${PYTHON_VERSION} -y

    info "Installing Python dependencies..."
    "${ENV_DIR}/bin/pip" install pydantic fastmcp
fi

# Step 2: Install EVcouplings
echo ""
echo -e "${BLUE}Step 2: Installing EVcouplings${NC}"

if [ "$SKIP_ENV" = true ]; then
    info "Skipping EVcouplings installation (--skip-env)"
else
    info "Installing EVcouplings..."
    "${ENV_DIR}/bin/pip" install https://github.com/debbiemarkslab/EVcouplings/archive/develop.zip
fi

# Step 3: Clone and compile PLMC
echo ""
echo -e "${BLUE}Step 3: Setting up PLMC${NC}"

if [ "$SKIP_PLMC" = true ]; then
    info "Skipping PLMC compilation (--skip-plmc)"
elif [ -d "$REPO_DIR/plmc" ] && [ -f "$REPO_DIR/plmc/bin/plmc" ]; then
    info "PLMC already exists and compiled"
else
    if [ ! -d "$REPO_DIR/plmc" ]; then
        info "Cloning PLMC repository..."
        mkdir -p "$REPO_DIR"
        git clone "$PLMC_REPO" "$REPO_DIR/plmc"
    fi

    info "Compiling PLMC with OpenMP support..."
    cd "$REPO_DIR/plmc"
    make all-openmp || warn "PLMC compilation failed"
    cd "$SCRIPT_DIR"
fi

# Step 4: Install HHsuite
echo ""
echo -e "${BLUE}Step 4: Installing HHsuite${NC}"

if [ "$SKIP_HHSUITE" = true ]; then
    info "Skipping HHsuite installation (--skip-hhsuite)"
else
    info "Installing HHsuite from conda-forge..."
    $CONDA_CMD install -c conda-forge -c bioconda -p "$ENV_DIR" hhsuite -y || warn "Cannot install hhsuite without conda/mamba"
fi

# Step 5: Verify installation
echo ""
echo -e "${BLUE}Step 5: Verifying installation${NC}"

"${ENV_DIR}/bin/python" -c "import fastmcp; import pydantic; print('Core packages OK')" && success "Core packages verified" || error "Package verification failed"

if [ -f "$REPO_DIR/plmc/bin/plmc" ]; then
    success "PLMC binary found"
else
    warn "PLMC binary not found"
fi

# Print summary
echo ""
echo -e "${GREEN}=============================================="
echo "           Setup Complete!"
echo "==============================================${NC}"
echo ""
echo "Environment: $ENV_DIR"
echo "PLMC:        $REPO_DIR/plmc"
echo ""
echo -e "${YELLOW}Claude Code Configuration:${NC}"
echo ""
cat << EOF
{
  "mcpServers": {
    "plmc": {
      "command": "${ENV_DIR}/bin/python",
      "args": ["${SCRIPT_DIR}/src/plmc_mcp.py"]
    }
  }
}
EOF
echo ""
echo "To add to Claude Code:"
echo "  claude mcp add plmc -- ${ENV_DIR}/bin/python ${SCRIPT_DIR}/src/plmc_mcp.py"
echo ""

#!/bin/bash

# Foldseek Pocket Miner - Web Application Launcher
# =================================================

set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

echo -e "${GREEN}Starting Foldseek Pocket Miner Web Application${NC}"
echo "=============================================="

# Get the directory of this script
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$SCRIPT_DIR"

# Check if foldseek is available
if ! command -v ~/bin/foldseek &> /dev/null && ! command -v foldseek &> /dev/null; then
    echo -e "${YELLOW}Warning: foldseek not found in PATH. Make sure it's installed.${NC}"
fi

# Check if PyMOL is available (optional)
if ! command -v pymol &> /dev/null; then
    echo -e "${YELLOW}Note: PyMOL not found. Some features may be limited.${NC}"
fi

# Create virtual environment if it doesn't exist
if [ ! -d "venv" ]; then
    echo -e "${YELLOW}Creating virtual environment...${NC}"
    python3 -m venv venv
fi

# Activate virtual environment
source venv/bin/activate

# Install dependencies
echo "Installing dependencies..."
pip install -q -r web/requirements.txt
pip install -q -e .

# Create necessary directories
mkdir -p web/uploads web/results

echo ""
echo -e "${GREEN}Starting web server...${NC}"
echo -e "Open your browser at: ${GREEN}http://localhost:5000${NC}"
echo ""
echo "Press Ctrl+C to stop the server"
echo ""

# Run the Flask app
cd web
python app.py

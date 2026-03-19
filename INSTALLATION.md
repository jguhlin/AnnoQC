# Installation Guide

This guide covers installing AnnoQC using various methods. Choose the one that best fits your workflow.

## Option 1: Pixi (Recommended)

Pixi provides a reproducible environment with all dependencies managed for you.

### Prerequisites

None! Pixi will manage Rust and all dependencies.

### Installation

```bash
# Install pixi
curl -fsSL https://pixi.sh/install.sh | bash

# Clone the repository
git clone https://github.com/jguhlin/AnnoQC.git
cd AnnoQC

# Run annoqc
pixi run annoqc --help
```

### First Run

```bash
pixi run annoqc prepare
pixi run annoqc analyze --config config.example.toml
```

## Option 2: Precompiled Binary

Download a precompiled binary for your platform from [GitHub Releases](https://github.com/jguhlin/AnnoQC/releases).

### Linux

```bash
wget https://github.com/jguhlin/AnnoQC/releases/latest/download/annoqc-linux-x86_64.tar.gz
tar -xzf annoqc-linux-x86_64.tar.gz
chmod +x annoqc
./annoqc --version
```

### macOS

```bash
wget https://github.com/jguhlin/AnnoQC/releases/latest/download/annoqc-macos-x86_64.tar.gz
tar -xzf annoqc-macos-x86_64.tar.gz
chmod +x annoqc
./annoqc --version
```

### macOS Apple Silicon (ARM64)

```bash
wget https://github.com/jguhlin/AnnoQC/releases/latest/download/annoqc-macos-aarch64.tar.gz
tar -xzf annoqc-macos-aarch64.tar.gz
chmod +x annoqc
./annoqc --version
```

### Windows

```powershell
# Download from https://github.com/jguhlin/AnnoQC/releases/latest
# Extract the zip file
# Open PowerShell in the extracted directory
.\annoqc.exe --version
```

### Optional: Add to PATH

To run `annoqc` from anywhere, add it to your PATH:

**Linux/macOS:**
```bash
sudo mv annoqc /usr/local/bin/
```

**Windows:**
Add the directory containing `annoqc.exe` to your system PATH.

## Option 3: Cargo

Install from source using Cargo (Rust package manager).

### Prerequisites

- Rust toolchain (1.70+)
- C compiler and linker
- OpenSSL development libraries (on some platforms)

### Install Rust

```bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
source $HOME/.cargo/env
```

### Install AnnoQC

```bash
cargo install AnnoQC --locked
```

The `--locked` flag ensures you use the exact dependency versions from `Cargo.lock`.

### Verify Installation

```bash
annoqc --version
```

## Option 4: Docker

AnnoQC is available as a Docker image with all dependencies pre-installed.

### Pull the Image

```bash
docker pull ghcr.io/jguhlin/annoqc:latest
```

### Run with Docker

```bash
docker run --rm -v /path/to/data:/data ghcr.io/jguhlin/annoqc:latest \
  annoqc analyze --config /data/config.toml
```

### Docker Compose Example

```yaml
version: '3.8'
services:
  annoqc:
    image: ghcr.io/jguhlin/annoqc:latest
    volumes:
      - ./data:/data
    working_dir: /data
```

## External Dependencies

AnnoQC requires the following external tools for full functionality:

### DIAMOND (Required)

DIAMOND is used for homology searches.

**Install via conda:**
```bash
conda install -c bioconda diamond
```

**Install via apt (Ubuntu/Debian):**
```bash
sudo apt-get install diamond-aligner
```

**Install from source:**
```bash
wget http://github.com/bbuchfink/diamond/releases/download/v2.0.15/diamond-linux64.tar.gz
tar -xzf diamond-linux64.tar.gz
sudo mv diamond/diamond /usr/local/bin/
```

### HMMER (Optional)

HMMER is used for domain architecture analysis.

**Install via conda:**
```bash
conda install -c bioconda hmmer
```

**Install via apt (Ubuntu/Debian):**
```bash
sudo apt-get install hmmer
```

### MAFFT (Optional)

MAFFT is used for multiple sequence alignment.

**Install via conda:**
```bash
conda install -c bioconda mafft
```

**Install via apt (Ubuntu/Debian):**
```bash
sudo apt-get install mafft
```

## Platform-Specific Notes

### Ubuntu/Debian

```bash
# System dependencies
sudo apt-get update
sudo apt-get install -y \
    build-essential \
    pkg-config \
    libssl-dev \
    curl
```

### macOS

```bash
# Install Xcode command line tools
xcode-select --install

# Install Homebrew (if not already installed)
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"

# Install dependencies
brew install openssl
```

### Windows

Install:
1. [Visual Studio Build Tools](https://visualstudio.microsoft.com/visual-cpp-build-tools/)
2. [Rustup](https://rustup.rs/)
3. [Git for Windows](https://git-scm.com/download/win)

## Troubleshooting

### Error: "command not found: annoqc"

The `annoqc` binary is not in your PATH. Either:
- Add the binary location to your PATH, or
- Use the full path to the binary: `/path/to/annoqc --help`

### Error: "DIAMOND not found"

Install DIAMOND using one of the methods above, or ensure it's in your PATH.

### Error: "OpenSSL headers not found"

Install OpenSSL development headers:

**Ubuntu/Debian:**
```bash
sudo apt-get install libssl-dev pkg-config
```

**macOS:**
```bash
brew install openssl
export OPENSSL_DIR=$(brew --prefix openssl)
```

### Error: "Permission denied" when running binary

Make the binary executable:
```bash
chmod +x annoqc
```

### Build fails on macOS Apple Silicon

If you're building from source on Apple Silicon:

```bash
# Install Rust with arch64 support
rustup default stable-aarch64-apple-darwin

# Build
cargo build --release
```

## Verification

After installation, verify everything works:

```bash
# Check version
annoqc --version

# Check help
annoqc --help

# Verify external tools (optional)
diamond version
hmmscan --version  # if using HMMER
mafft --version    # if using MAFFT
```

## Next Steps

Once installed, see the [Quickstart Guide](https://jguhlin.github.io/AnnoQC/quickstart.html) for your first run.

For configuration options, see [analyze.md](https://jguhlin.github.io/AnnoQC/analyze.html).

## Getting Help

If you encounter issues not covered here:

1. Check [GitHub Issues](https://github.com/jguhlin/AnnoQC/issues)
2. Search existing discussions
3. Open a new issue with:
   - Your OS and version
   - Installation method used
   - Error message (full output)
   - Steps to reproduce

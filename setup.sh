#!/bin/bash
# setup.sh

echo "Creating necessary directories..."
mkdir -p references
mkdir -p tmp
mkdir -p inputs

echo "Downloading reference data from Zenodo..."
wget -O references.tar.gz "https://zenodo.org/records/24052169/files/references.tar.gz?download=1"

echo "Extracting reference data..."
# The -C flag forces all the loose files inside the tarball into the references/ folder
tar -xzf references.tar.gz -C references/

echo "Cleaning up compressed archive..."
rm references.tar.gz

echo "Setup complete! You can now run: docker compose up"
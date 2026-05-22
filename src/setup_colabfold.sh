#!/bin/bash

# Move to the root of the project to keep things organized
cd /mnt/c/Users/maxim/OrthCom

echo "Downloading LocalColabFold installer..."
wget https://raw.githubusercontent.com/YoshitakaMo/localcolabfold/main/install_colabbatch_linux.sh

echo "Running the installer (this will take a few minutes)..."
bash install_colabbatch_linux.sh

echo "Cleaning up installer file..."
rm install_colabbatch_linux.sh

echo "ColabFold installation complete! The engine is in the 'localcolabfold' directory."

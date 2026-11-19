# Use Miniforge so conda/mamba is globally available for Snakemake
FROM condaforge/miniforge3:latest

# Set the working directory inside the container
WORKDIR /prs-app

# Copy the environment file
COPY app/environment.yml /prs-app/app/environment.yml

# Create the environment for Streamlit/Snakemake
RUN mamba env create -f /prs-app/app/environment.yml

# Activate the environment globally for subsequent commands
ENV PATH=/opt/conda/envs/prs-app/bin:$PATH

EXPOSE 8501

# Command to run the app (relative to WORKDIR)
CMD ["streamlit", "run", "app/main.py", "--server.address=0.0.0.0"]
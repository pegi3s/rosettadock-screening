FROM pegi3s/docker:29.0.1

# Install Python 3 and necessary packages
RUN apt-get update && apt-get install -y \
    python3 \
    python3-pip \
    && rm -rf /var/lib/apt/lists/*
RUN apt update && apt install -y python3-pandas python3-openpyxl

# Set working directory for the pipeline and copied scripts
WORKDIR /opt

# Copy Python scripts to the runtime image
COPY pipeline.py /opt/
COPY script_mat.py /opt/
COPY script_pre.py /opt/

# Ensure config and data can be mounted at /data
# The config file should be mounted as a volume at runtime

# Run the pipeline from the copied image scripts
ENTRYPOINT ["python3", "/opt/pipeline.py"]

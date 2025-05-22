# Base image with Python 3.10 (as used by Hugging Face by default)
FROM python:3.10-slim-buster

# Install system dependencies needed for RDKit and other libraries
# RDKit needs libstdc++6, sometimes other common libs
RUN apt-get update && apt-get install -y \
    build-essential \
    libstdc++6 \
    git \
    git-lfs \
    ffmpeg \
    libsm6 \
    libxext6 \
    cmake \
    rsync \
    libgl1-mesa-glx \
    # Clean up APT cache to reduce image size
    && rm -rf /var/lib/apt/lists/* \
    && git lfs install

# Set working directory inside the container
WORKDIR /app

# Copy your requirements.txt to the working directory
# We're manually defining requirements here for absolute control
COPY requirements.txt .

# Install Python dependencies from requirements.txt
# This time, we will manually define them in requirements.txt (next step)
# This will ensure pip gets the versions we specify.
RUN pip install --no-cache-dir -r requirements.txt

# Copy your application code into the container
COPY app.py .

# Expose the port Gradio runs on
EXPOSE 7860

# Command to run the application when the container starts
CMD ["python", "app.py"]

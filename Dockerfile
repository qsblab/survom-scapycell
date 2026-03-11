FROM python:3.11-slim

# Install system-level dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    libgl1 \
    libhdf5-dev \
    curl \
    git \
    && rm -rf /var/lib/apt/lists/*

# Set working directory inside the container
WORKDIR /app

# Upgrade pip and install build tools early
RUN pip install --upgrade pip && pip install numpy Cython

# Copy Python requirements and install them
COPY requirements.txt .
RUN pip install -r requirements.txt

# Copy all app source code into the container
COPY . .

# Ensure necessary internal directories exist
RUN mkdir -p /app/output /app/data /app/uploaded_data /root/.celloracle \
    && chmod -R 777 /app/output /app/data /app/uploaded_data

# Expose Dash app port
EXPOSE 8050

# Set the default command to run your app
CMD ["python", "app.py"]

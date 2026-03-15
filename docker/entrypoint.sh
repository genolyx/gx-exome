#!/bin/bash
set -e

echo "==================================="
echo "Carrier Screening Pipeline Container"
echo "==================================="

# Create necessary directories
mkdir -p /data/fastq /data/analysis /data/output /data/log
mkdir -p /var/log/supervisor

# Check if environment variables are set
echo "Checking environment..."
echo "BASE_DIR: ${BASE_DIR:-/data}"
echo "PORTAL_API_URL: ${PORTAL_API_URL:-Not set}"

# Initialize Nextflow if needed
if [ ! -f /root/.nextflow/history ]; then
    echo "Initializing Nextflow..."
    nextflow info
fi

# Start services
echo "Starting services..."
echo "- Dashboard: http://localhost:5000"
echo "- Daemon API: http://localhost:8080"
echo "==================================="

exec "$@"

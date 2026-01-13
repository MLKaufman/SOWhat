#!/usr/bin/env bash

set -euo pipefail

# Configuration
IMAGE_NAME="mlkaufman/sowhat:latest"
CONTAINER_NAME="sowhat_app"
SHINY_PORT=3838

echo "=== SoWhat Docker Runner ==="

# Check if Docker is installed
if ! command -v docker &> /dev/null; then
    echo "Error: Docker is not installed or not in PATH"
    exit 1
fi

# Check if Docker daemon is running
if ! docker info &> /dev/null; then
    echo "Error: Docker daemon is not running"
    exit 1
fi

# Pull the latest image
echo "Pulling latest image: ${IMAGE_NAME}..."
docker pull "${IMAGE_NAME}"

# Stop and remove existing container if it exists
if docker ps -a --format '{{.Names}}' | grep -q "^${CONTAINER_NAME}$"; then
    echo "Removing existing container: ${CONTAINER_NAME}..."
    docker rm -f "${CONTAINER_NAME}"
fi

# Run the container
echo "Starting Shiny app container: ${CONTAINER_NAME}..."
docker run \
    --name "${CONTAINER_NAME}" \
    -p ${SHINY_PORT}:3838 \
    -e SHINY_PORT=3838 \
    "${IMAGE_NAME}"

echo "To view logs: docker logs -f ${CONTAINER_NAME}"
echo "To stop: docker stop ${CONTAINER_NAME}"
echo "To remove: docker rm ${CONTAINER_NAME}"
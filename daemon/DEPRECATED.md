# DEPRECATED: Legacy Daemon

This directory contains the legacy daemon implementation that was originally built for the Dark Gene Pipeline. It has been superseded by the **genolyx/service-daemon** project, which provides a unified, multi-service daemon framework.

## Migration Guide

The `service-daemon` project replaces both `daemon.py` (file watcher) and `api_server.py` (Flask API) with a single FastAPI-based service that supports multiple analysis services through a plugin architecture.

| Legacy (this directory) | New (service-daemon) | Description |
|---|---|---|
| `daemon.py` | `app/queue_manager.py` | Job queue management and monitoring |
| `api_server.py` | `app/main.py` | REST API endpoints |
| `PortalAPI` class | `app/platform_client.py` | Platform/Portal communication |
| Flask + watchdog | FastAPI + asyncio | Web framework |
| Single service | Multi-service plugins | Architecture |

## How to Use service-daemon

1. Clone the service-daemon repository:
   ```bash
   git clone https://github.com/genolyx/service-daemon.git
   ```

2. Configure for Carrier Screening:
   ```bash
   cp .env.local .env
   # Edit .env with actual paths
   ```

3. Run:
   ```bash
   ENABLED_SERVICES=carrier_screening python3 -m uvicorn app.main:app --host 0.0.0.0 --port 8000
   ```

4. Or use Docker Compose (from `docker/` directory):
   ```bash
   docker compose up -d
   ```

## Files in This Directory

These files are kept for reference only and should not be used in production:

- `daemon.py` - Legacy file watcher daemon
- `api_server.py` - Legacy Flask API server
- `Dockerfile` - Legacy Docker build (use `docker/Dockerfile` instead)
- `docker-compose.yml` - Legacy compose (use `docker/docker-compose.yml` instead)

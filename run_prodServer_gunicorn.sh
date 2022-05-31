#!/usr/bin/bash

sudo env PATH=$PATH gunicorn src.main:app --workers 4  --worker-class uvicorn.workers.UvicornWorker -b 0.0.0.0:1010 --access-logfile '-' -t 120

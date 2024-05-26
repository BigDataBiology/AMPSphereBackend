from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from src.router import amp_router, family_router, default_router
import os


os.environ['worker_class'] = 'uvicorn.workers.UvicornH11Worker'


with open('description.md', 'r') as f:
    description = f.read()


app = FastAPI(
    title="AMPSphere API",
    description=description,
    version="0.1.0",
    contact={
        "name": "Luis Pedro Coelho",
        "url": "https://luispedro.org",
        "email": "luis@luispedro.org",
    },
    license_info={
        "name": "MIT",
        "url": "https://github.com/BigDataBiology/AMPSphereBackend/blob/main/COPYING.MIT",
    },
    docs_url='/',
    redoc_url='/redoc',
    openapi_url='/openapi.json'
)

if 'ADD_CORS_HEADERS' in os.environ:
    app.add_middleware(
         CORSMiddleware,
         allow_origins=["*"],
         allow_methods=["*"],
         allow_headers=["*"],
         allow_credentials=True,
     )


for router in [amp_router, family_router, default_router]:
    app.include_router(router=router)

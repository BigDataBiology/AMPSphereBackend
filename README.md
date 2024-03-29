# AMPSphereBackend

This repository contains source code for the backend of [AMPSphere website](https://ampsphere.big-data-biology.org)

## Components and description

```txt
├── config/*                # Define default configuration of the backend here.
├── conftest.py             # This file only exists for unittest purpose. 
├── database/*              # AMPSphere database and sequence search databases. 
├── description.md          # Markdown content to be displayed on the Swagger API page.
├── __init__.py   
├── README.md  
├── requirements.txt
├── setup.cfg
├── src
│   ├── crud.py             # Define database-related functions here: including Create, Read, Update, and Delete
│   ├── database.py         # Generate a database session and a base class (used to define models of the database) here
│   ├── __init__.py         
│   ├── main.py             # Main FastAPI application.
│   ├── models.py           # Define models (classes) of the database
│   ├── __pycache__
│   ├── router.py           # Mount functions to URLs here
│   ├── schemas.py          # Define basic classes to return by the backend
│   └── utils.py            # Define database-unrelated functions here.
├── taxdump.tar.gz
└── tests
    ├── coverage_html       # Coverage report (HTML).
    ├── __init__.py         
    ├── inputs.tsv          # Define test cases and expected http response code here.
    ├── __pycache__
    └── testing.py          # Main testing script.
```

## Deploy instructions

### Prepare necessary data
Copy the latest AMPSphere data (should be generated based on Celio's script) to data/original_data. The folder will be look like.
```text
data/original_data
├── AMPSphere_generation_v.2022-03
│   ├── analysis
│   ├── data
│   ├── main.py
│   ├── README.md
│   ├── seqlogo
│   └── utils
├── metadata_analysis
│   ├── data
│   ├── main.py
│   ├── outputs
│   ├── README.md
│   └── utils
└── zenodo_repo
    ├── AMPSphere_v.2022-03.faa.gz
    ├── AMPSphere_v.2022-03.fna.xz
    ├── AMPSphere_v.2022-03.general_geneinfo.tsv.gz
    ├── AMPSphere_v.2022-03.quality_assessment.tsv.gz
    ├── README.md
    └── SPHERE_v.2022-03.levels_assessment.tsv.gz
```

### Prepare the environment

```bash
conda install -c default -c conda-forge -c bioconda \
    fastapi gunicorn sqlalchemy pandas pytest pydantic biopython requests
    mmseqs2 hmmer uvicorn

# For testing
conda install -c default -c conda-forge httpx
```

### Generate the backend data and test on dev server

**Important**: Before executing the data generation step, make sure you have the right VERSION_CODE argument in the script (edit it to modify).

```bash
./generate_data.sh  
python -m pytest tests
```

This will cost hours to half a day depends on the available computing resources.


### Transfer all the data to the production server

```bash
tar zcvf data.tgz --exclude "data/original_data" data
scp data.tgz ampsphere-api.big-data-biology.org:/AMPSphere/website/AMPSphereBackend/
```

### Uncompress the data and do unittest before production

```bash
tar zxvf data.tgz
python -m pytest tests/
```

If everything is Okay, please run the backend in production mode.


Go to its API page and test it manually if necessary: [https://ampsphere-api.big-data-biology.org/](https://ampsphere-api.big-data-biology.org/)

## Contact

- [Luis Pedro Coelho](https://luispedro.org)

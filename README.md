# AMPSphereBackend

This repository contains source code for the backend of [AMPSphere website](https://ampsphere.big-data-biology.org)


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

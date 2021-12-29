from sqlalchemy.orm import Session
from src.database import SessionLocal
from typing import List, Dict
from fastapi import Depends
from fastapi.responses import FileResponse
from src import schemas
from src import utils
from src import crud
from fastapi import APIRouter


def get_db():
    db = SessionLocal()
    try:
        yield db
    finally:
        db.close()


# Change here.
default_route_summary = ' '
version = '/v1'

amp_router = APIRouter(
    prefix=version + "/amps",
    tags=['amp']
)


# TODO define consistent schema for AMP object.
@amp_router.get(path="",
                response_model=schemas.PagedAMPs,
                summary=default_route_summary)
def amps(db: Session = Depends(get_db),
         family: str = None,
         habitat: str = None,
         sample: str = None,
         microbial_source: str = None,
         pep_length_interval: str = None,
         mw_interval: str = None,
         pI_interval: str = None,
         charge_interval: str = None,
         page_size: int = 20,
         page: int = 0):
    """
    Available filters:
    - :param family: AMP family accession (e.g., AMP10.000_000).
    - :param habitat: Human-readable habitat name (e.g., human gut).
    - :param sample: Metagenome sample accession / Progenomes2 genome accession.
    - :param microbial_source: Taxonomy name (e.g., Escherichia coli).
    - :param pep_length_interval: Peptide length interval (format: `min_len,max_len`, e.g., `40,50`).
    - :param mw_interval: Molecular weight interval (format: `min_mw,max_mw`, e.g., `813,1000`).
    - :param pI_interval: Isoelectric point interval (format: `min_pI,max_pI`, e.g., `4,12`).
    - :param charge_interval: Charge at pH 7 interval (format: `min_charge,max_charge`, e.g., `-57,44`).
    """
    return crud.get_amps(db, page=page, page_size=page_size,
                         family=family, habitat=habitat, microbial_source=microbial_source, sample=sample,
                         pep_length_interval=pep_length_interval, mw_interval=mw_interval,
                         pI_interval=pI_interval, charge_interval=charge_interval)


@amp_router.get(path="/{accession}",
                response_model=schemas.AMP,
                summary=default_route_summary)
def amp(accession: str,
        db: Session = Depends(get_db)):
    return crud.get_amp(accession, db)


@amp_router.get(path="/{accession}/helicalwheel",
                # response_class=FileResponse,
                summary=default_route_summary)
def amp_helicalwheel(accession: str):
    path = crud.get_amp_helicalwheel(accession)
    print(path)
    return FileResponse(path)


@amp_router.get(path="/{accession}/features",
                response_model=schemas.AMPFeatures,
                summary=default_route_summary)
def amp_features(accession: str = 'AMP10.000_000',
                 db: Session = Depends(get_db)):
    # TODO get sequence here.
    return crud.get_amp_features(accession, db)


@amp_router.get(path="/{accession}/distributions",
                response_model=schemas.Distributions,
                summary=default_route_summary)
def distributions(accession: str = 'AMP10.000_000',
                  db: Session = Depends(get_db)):
    return crud.get_distributions(accession=accession, db=db)


@amp_router.get(path="/{accession}/metadata",
                response_model=schemas.PagedMetadata,
                summary=default_route_summary)
def metadata(accession: str = 'AMP10.000_000',
             db: Session = Depends(get_db),
             page: int = 0,
             page_size: int = 20):
    return crud.get_amp_metadata(accession=accession, db=db, page=page, page_size=page_size)


family_router = APIRouter(
    prefix=version + '/families',
    tags=['family']
)


@family_router.get(path="",
                   response_model=schemas.PagedFamilies,
                   summary=default_route_summary)
def families(db: Session = Depends(get_db),
             habitat: str = None,
             sample: str = None,
             microbail_source: str = None,
             page_size: int = 5,
             page: int = 0):
    families = crud.get_families(
        db, page=page, page_size=page_size,
        habitat=habitat, microbail_source=microbail_source, sample=sample
    )
    return families


@family_router.get(path="/{accession}",
                   response_model=schemas.Family,
                   summary=default_route_summary)
def family(accession: str = 'SPHERE-III.001_396', db: Session = Depends(get_db)):
    families = crud.get_family(accession, db)
    return families


@family_router.get(path="/{accession}/features",
                   response_model=Dict[str, schemas.AMPFeatures],
                   summary=default_route_summary)
def fam_features(accession: str = 'SPHERE-III.001_396', db: Session = Depends(get_db)):
    # TODO get sequence here.
    return crud.get_fam_features(accession, db=db)


@family_router.get(path="/{accession}/distributions",
                   response_model=schemas.Distributions,
                   summary=default_route_summary)
def fam_distributions(accession: str = 'SPHERE-III.001_396', db: Session = Depends(get_db)):
    return crud.get_distributions(accession=accession, db=db)


@family_router.get("/{accession}/downloads",
                   response_model=schemas.FamilyDownloads,
                   summary=default_route_summary)
def fam_downloads(accession: str = 'SPHERE-III.001_396',
                  db: Session = Depends(get_db)):
    return crud.get_fam_downloads(accession=accession, db=db)


@family_router.get("/{accession}/downloads/{file}",
                   # response_model=schemas.FamilyDownloads,
                   summary=default_route_summary)
def fam_download_file(accession: str, file: str):
    return FileResponse(utils.fam_download_file(accession=accession, file=file))


default_router = APIRouter(
    prefix=version + '',
    tags=['default']
)


@default_router.get(path="/statistics",
                    response_model=schemas.Statistics,
                    summary=default_route_summary)
def get_statistics(db: Session = Depends(get_db)):
    return crud.get_statistics(db)


@default_router.get(path='/available_filters',
                    #response_model=schemas.Filters,
                    summary=default_route_summary
                    )
def get_filters(db: Session = Depends(get_db)):
    return crud.get_filters(db)


@default_router.get(path="/downloads",
                    response_model=List[str],
                    summary=default_route_summary)
def get_downloads():
    downloads = utils.get_downloads()
    return downloads


@default_router.get(path="/downloads/{file}",
                    # response_class=FileResponse,
                    summary=default_route_summary)
async def download_file(file: str):
    return FileResponse(utils.download(file))


@default_router.get(path="/search/text",
                    response_model=schemas.PagedAMPs,
                    summary=default_route_summary)
def text_search(db: Session = Depends(get_db),
                query: str = 'AMP10.000_000',
                page: int = 0,
                page_size: int = 20):
    return crud.search_by_text(db, text=query, page=page, page_size=page_size)


@default_router.get(path="/search/mmseqs",
                    response_model=List[schemas.mmSeqsSearchResult],
                    summary=default_route_summary)
def mmseqs_search(query: str = 'KKVKSIFKKALAMMGENEVKAWGIGIK', db: Session = Depends(get_db)):
    return crud.mmseqs_search(query, db)


@default_router.get(path="/search/hmmer",
                    response_model=List[schemas.HMMERSearchResult],
                    summary=default_route_summary)
def hmmscan_search(query: str = 'KKVKSIFKKALAMMGENEVKAWGIGIK'):
    return crud.hmmscan_search(query)
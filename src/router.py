from sqlalchemy.orm import Session
from src.database import SessionLocal
from typing import List, Dict
from fastapi import Depends
from fastapi.responses import FileResponse, JSONResponse, StreamingResponse
from src import schemas
from src import utils
from src import crud
from fastapi import APIRouter, Request


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
                response_class=JSONResponse,
                summary=default_route_summary)
def amps(db: Session = Depends(get_db),
         exp_evidence: str = None,
         antifam: str = None,
         RNAcode: str = None,
         coordinates: str = None,
         family: str = None,
         habitat: str = None,
         sample_genome: str = None,
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
                         exp_evidence=exp_evidence, antifam=antifam, RNAcode=RNAcode, coordinates=coordinates,
                         sample_genome=sample_genome,
                         family=family, habitat=habitat, microbial_source=microbial_source,
                         pep_length_interval=pep_length_interval, mw_interval=mw_interval,
                         pI_interval=pI_interval, charge_interval=charge_interval)


@amp_router.get(path="/{accession}",
                response_model=schemas.AMP,
                response_class=JSONResponse,
                summary=default_route_summary)
def amp(accession: str,
        db: Session = Depends(get_db)):
    return crud.get_amp(accession, db)


@amp_router.get(path="/{accession}/features",
                response_model=schemas.AMPFeatures,
                response_class=JSONResponse,
                summary=default_route_summary)
def amp_features(accession: str = 'AMP10.000_000',
                 db: Session = Depends(get_db)):
    # TODO get sequence here.
    return crud.get_amp_features(accession, db)


@amp_router.get(path="/{accession}/distributions",
                response_model=schemas.Distributions,
                response_class=JSONResponse,
                summary=default_route_summary)
def distributions(accession: str = 'AMP10.000_000',
                  db: Session = Depends(get_db)):
    return crud.get_distributions(accession=accession, db=db)


@amp_router.get(path="/{accession}/metadata",
                response_model=schemas.PagedMetadata,
                response_class=JSONResponse,
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
                   response_class=JSONResponse,
                   summary=default_route_summary)
def families(request: Request,
             db: Session = Depends(get_db),
             habitat: str = None,
             sample: str = None,
             microbail_source: str = None,
             page_size: int = 5,
             page: int = 0):
    families = crud.get_families(
        db=db, request=request, page=page, page_size=page_size,
        habitat=habitat, microbail_source=microbail_source, sample=sample)
    return families


@family_router.get(path="/{accession}",
                   response_model=schemas.Family,
                   response_class=JSONResponse,
                   summary=default_route_summary)
def family(accession: str, request: Request, db: Session = Depends(get_db)):
    families = crud.get_family(accession, db=db, request=request)
    return families


@family_router.get(path="/{accession}/features",
                   response_model=Dict[str, schemas.AMPFeatures],
                   response_class=JSONResponse,
                   summary=default_route_summary)
def fam_features(accession: str = 'SPHERE-III.001_396', db: Session = Depends(get_db)):
    # TODO get sequence here.
    return crud.get_fam_features(accession, db=db)


@family_router.get(path="/{accession}/distributions",
                   response_model=schemas.Distributions,
                   response_class=JSONResponse,
                   summary=default_route_summary)
def fam_distributions(accession: str = 'SPHERE-III.001_396', db: Session = Depends(get_db)):
    return crud.get_distributions(accession=accession, db=db)


@family_router.get("/{accession}/downloads",
                   response_model=schemas.FamilyDownloads,
                   response_class=JSONResponse,
                   summary=default_route_summary)
def fam_downloads(accession: str,
                  request: Request,
                  db: Session = Depends(get_db)):
    return crud.get_fam_downloads(accession=accession, db=db, request=request)


@family_router.get("/{accession}/downloads/{file}",
                   # response_model=schemas.FamilyDownloads,
                   response_class=FileResponse,
                   summary=default_route_summary)
def fam_download_file(accession: str, file: str):
    return FileResponse(utils.fam_download_file(accession=accession, file=file))


default_router = APIRouter(
    prefix=version + '',
    tags=['default']
)


@default_router.get(path="/statistics",
                    response_model=schemas.Statistics,
                    response_class=JSONResponse,
                    summary=default_route_summary)
def get_statistics(db: Session = Depends(get_db)):
    return crud.get_statistics(db)


@default_router.get(path='/current_available_options',
                    # response_model=schemas.Filters,
                    response_class=JSONResponse,
                    summary=default_route_summary
                    )
def get_filtered_options(db: Session = Depends(get_db),
                         quality: str = None,
                         family: str = None,
                         habitat: str = None,
                         sample: str = None,
                         microbial_source: str = None,
                         pep_length_interval: str = None,
                         mw_interval: str = None,
                         pI_interval: str = None,
                         charge_interval: str = None):
    return crud.get_filtered_options(db, quality=quality, family=family, habitat=habitat,
                                     microbial_source=microbial_source, sample=sample,
                                     pep_length_interval=pep_length_interval, mw_interval=mw_interval,
                                     pI_interval=pI_interval, charge_interval=charge_interval)


@default_router.get(path='/all_available_options',
                    # response_model=schemas.Filters,
                    response_class=JSONResponse,
                    summary=default_route_summary
                    )
def get_all_options(db: Session = Depends(get_db)):
    return crud.get_all_options(db)


@default_router.get(path='/in_db/{entity_type}/{accession}',
                    summary=default_route_summary)
def in_db(db: Session = Depends(get_db), entity_type: str = 'family', accession: str = 'SPHERE-III.000_428'):
    return crud.entity_in_db(db=db, entity_type=entity_type, accession=accession)


@default_router.get(path="/downloads",
                    response_model=List[str],
                    response_class=JSONResponse,
                    summary=default_route_summary)
def get_downloads():
    downloads = utils.get_downloads()
    return downloads


@default_router.get(path="/downloads/{file}",
                    response_class=StreamingResponse,
                    summary=default_route_summary)
def download_file(file: str):
    def iterfile():
        with open(utils.download(file), mode="r" if not file.endswith('sqlite') else 'rb') as f:
            yield from f
    media_types = dict(
        sqlite="application/vnd.sqlite3",
        tsv="text/tab-separated-values",
        mmseqsdb="application/octet-stream",
        hmm="application/octet-stream",

    )
    return StreamingResponse(iterfile(), media_type=media_types[file.split('.')[-1]])



@default_router.get(path="/search/mmseqs",
                    response_model=List[schemas.mmSeqsSearchResult],
                    response_class=JSONResponse,
                    summary=default_route_summary)
def mmseqs_search(query: str = 'KKVKSIFKKALAMMGENEVKAWGIGIK', db: Session = Depends(get_db)):
    return crud.mmseqs_search(query, db)


@default_router.get(path="/search/hmmer",
                    response_model=List[schemas.HMMERSearchResult],
                    response_class=JSONResponse,
                    summary=default_route_summary)
def hmmscan_search(query: str = 'KKVKSIFKKALAMMGENEVKAWGIGIK'):
    return crud.hmmscan_search(query)

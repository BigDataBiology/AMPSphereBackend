from typing import List, Dict
from fastapi import Depends
from fastapi.responses import FileResponse, JSONResponse
from src import schemas, crud, database, utils
from fastapi import APIRouter


# Change here.
default_route_summary = ' '
version = '/v1'

amp_router = APIRouter(
    prefix=version + "/amps",
    tags=['amp']
)


@amp_router.get(path="",
                response_model=schemas.PagedAMPs,
                response_class=JSONResponse,
                summary=default_route_summary)
def amps(
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
    return crud.get_amps(page=page, page_size=page_size,
                         exp_evidence=exp_evidence, antifam=antifam, RNAcode=RNAcode, coordinates=coordinates,
                         sample_genome=sample_genome,
                         family=family, habitat=habitat, microbial_source=microbial_source,
                         pep_length_interval=pep_length_interval, mw_interval=mw_interval,
                         pI_interval=pI_interval, charge_interval=charge_interval)


@amp_router.get(path="/{accession}",
                response_model=schemas.AMP,
                response_class=JSONResponse,
                summary=default_route_summary)
def amp(accession: str):
    return crud.get_amp(accession)


@amp_router.get(path="/{accession}/features",
                response_model=schemas.AMPFeatures,
                response_class=JSONResponse,
                summary=default_route_summary)
def amp_features(accession: str = 'AMP10.000_000'):
    # TODO get sequence here.
    return crud.get_amp_features(accession)


@amp_router.get(path="/{accession}/distributions",
                response_model=schemas.Distributions,
                response_class=JSONResponse,
                summary=default_route_summary)
def distributions(accession: str = 'AMP10.000_000'):
    return crud.get_distributions(accession=accession)


@amp_router.get(path="/{accession}/metadata",
                response_model=schemas.PagedMetadata,
                response_class=JSONResponse,
                summary=default_route_summary)
def metadata(accession: str = 'AMP10.000_000',
             page: int = 0,
             page_size: int = 20):
    return crud.get_amp_metadata(accession=accession, page=page, page_size=page_size)


@amp_router.get(path="/{accession}/coprediction",
                response_class=JSONResponse,
                summary=default_route_summary)
def coprediction(accession: str = 'AMP10.000_000'):
    r = database.coprediction.get(accession, None)
    if r is None:
        raise HTTPException(status_code=400, detail='invalid accession received.')
    return [{'predictor': k, 'value': v} for k, v in r.items()]

family_router = APIRouter(
    prefix=version + '/families',
    tags=['family']
)


@family_router.get(path="",
                   response_model=schemas.PagedFamilies,
                   response_class=JSONResponse,
                   summary=default_route_summary)
def families(habitat: str = None,
             sample: str = None,
             microbial_source: str = None,
             page_size: int = 5,
             page: int = 0):
    return crud.get_families(
        page=page, page_size=page_size,
        habitat=habitat, microbial_source=microbial_source, sample=sample)


@family_router.get(path="/{accession}",
                   response_model=schemas.Family,
                   response_class=JSONResponse,
                   summary=default_route_summary)
def family(accession: str):
    return crud.get_family(accession)



@family_router.get(path="/{accession}/features",
                   response_model=Dict[str, schemas.AMPFeatures],
                   response_class=JSONResponse,
                   summary=default_route_summary)
def fam_features(accession: str = 'SPHERE-III.001_396'):
    return crud.get_fam_features(accession)


@family_router.get(path="/{accession}/distributions",
                   response_model=schemas.Distributions,
                   response_class=JSONResponse,
                   summary=default_route_summary)
def fam_distributions(accession: str = 'SPHERE-III.001_396'):
    return crud.get_distributions(accession=accession)


@family_router.get("/{accession}/downloads",
                   response_model=schemas.FamilyDownloads,
                   response_class=JSONResponse,
                   summary=default_route_summary)
def fam_downloads(accession: str):
    return crud.get_fam_downloads(accession=accession)


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
def get_statistics():
    return crud.get_statistics()


@default_router.get(path='/all_available_options',
                    # response_model=schemas.Filters,
                    response_class=JSONResponse,
                    summary=default_route_summary
                    )
def get_all_options():
    return crud.get_all_options()


@default_router.get(path='/in_db/{entity_type}/{accession}',
                    summary=default_route_summary)
def in_db(entity_type: str = 'family',
          accession: str = 'SPHERE-III.000_428'):
    return crud.entity_in_db(entity_type=entity_type, accession=accession)



_DOWNLOAD_FILE_PATH = {
        'AMPSphere_latest.sqlite':      'ampsphere_main_db/AMPSphere_latest.sqlite',
        'AMPSphere_latest.mmseqsdb':    'mmseqs_db/AMPSphere_latest.mmseqsdb',
        'AMPSphere_latest.hmm':         'hmmprofile_db/AMPSphere_latest.hmm',

        'AMP.tsv':                      'tables/AMP.tsv',
        'GMSCMetadata.tsv':             'tables/GMSCMetadata.tsv',
}

@default_router.get(path="/downloads",
                    response_model=List[str],
                    response_class=JSONResponse,
                    summary=default_route_summary)
def get_downloads():
    return list(sorted(_DOWNLOAD_FILE_PATH.keys()))


@default_router.get(path="/downloads/{file}",
                    response_class=FileResponse,
                    summary=default_route_summary)
def download_file(file: str):
    fpath = f'data/{_DOWNLOAD_FILE_PATH[file]}'
    return FileResponse(fpath)


@default_router.get(path="/search/mmseqs",
                    response_model=List[schemas.mmSeqsSearchResult],
                    response_class=JSONResponse,
                    summary=default_route_summary)
def mmseqs_search(query: str = 'KKVKSIFKKALAMMGENEVKAWGIGIK'):
    return crud.mmseqs_search(query)


@default_router.get(path="/search/hmmer",
                    response_model=List[schemas.HMMERSearchResult],
                    response_class=JSONResponse,
                    summary=default_route_summary)
def hmmscan_search(query: str = 'KKVKSIFKKALAMMGENEVKAWGIGIK'):
    return crud.hmmscan_search(query)

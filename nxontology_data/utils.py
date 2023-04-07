import gzip
import json
import logging
import sys
from pathlib import Path
from typing import Any

import bioregistry.resolve
import pandas as pd
from bioregistry.resource_manager import _safe_curie_to_str
from networkx.readwrite.json_graph import node_link_data
from nxontology import NXOntology
from rdflib.plugins.sparql.processor import SPARQLResult

logger = logging.getLogger(__name__)


def get_output_dir() -> Path:
    """Local output directory in this repository."""
    return Path(__file__).parent.parent.joinpath("output")


def get_source_output_dir(source: str) -> Path:
    output_dir = get_output_dir().joinpath(source)
    output_dir.mkdir(parents=True, exist_ok=True)
    return output_dir


def write_ontology(
    nxo: NXOntology[Any], output_dir: Path, compression_threshold_mb: float = 10.0
) -> Path:
    data = node_link_data(nxo.graph)
    json_bytes = json.dumps(data, indent=2, ensure_ascii=False).encode()
    json_size_mb = sys.getsizeof(json_bytes) / 1_000_000
    path = output_dir.joinpath(f"{nxo.name}.json")
    if json_size_mb > compression_threshold_mb:
        json_bytes = gzip.compress(json_bytes, mtime=0)
        path = path.with_name(f"{path.name}.gz")
        logger.info(
            f"{path.name}: gzip reduced size from {json_size_mb:.1f} to {sys.getsizeof(json_bytes) / 1_000_000:.1f} MB"
        )
    path.write_bytes(json_bytes)
    logger.info(f"Wrote ontology to {path}")
    # ensure JSON is valid and check_is_dag
    nxo.read_node_link_json(path.as_posix())
    return path


def write_dataframe(df: pd.DataFrame, path: Path) -> None:
    df.to_json(
        path,
        orient="records",
        compression={"method": "gzip", "mtime": 0},
        indent=2,
        date_format="iso",
    )


def sparql_results_to_df(results: SPARQLResult) -> pd.DataFrame:
    """
    Export results from an rdflib SPARQL query into a `pandas.DataFrame`,
    using Python types. See https://github.com/RDFLib/rdflib/issues/1179
    and https://github.com/RDFLib/sparqlwrapper/issues/205.
    """
    return pd.DataFrame(
        data=([None if x is None else x.toPython() for x in row] for row in results),
        columns=[str(x) for x in results.vars],
    )


def normalize_parsed_curie(xref_prefix: str, xref_accession: str) -> str | None:
    """
    Normalize a parsed CURIE according to Bioregistry.
    Return a string using preferred prefix capitalization.
    https://github.com/biopragmatics/bioregistry/issues/790
    """
    prefix, accession = bioregistry.resolve.normalize_parsed_curie(
        xref_prefix, xref_accession, use_preferred=True
    )
    return _safe_curie_to_str(prefix, accession)

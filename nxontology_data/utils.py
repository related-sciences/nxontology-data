import gzip
import json
import logging
import sys
from pathlib import Path
from typing import Any

import pandas as pd
from networkx.readwrite.json_graph import node_link_data
from nxontology import NXOntology
from rdflib.plugins.sparql.processor import SPARQLResult

logger = logging.getLogger(__name__)


def get_output_dir() -> Path:
    """Local output directory in this repository."""
    return Path(__file__).parent.parent.joinpath("output")


def write_ontology(nxo: NXOntology[Any], output_dir: Path) -> Path:
    data = node_link_data(nxo.graph)
    json_bytes = json.dumps(data, indent=2, ensure_ascii=False).encode()
    json_size_mb = sys.getsizeof(json_bytes) / 1_000_000
    path = output_dir.joinpath(f"{nxo.graph.graph['name']}.json")
    if json_size_mb > 10.0:
        json_bytes = gzip.compress(json_bytes, mtime=0)
        path = path.with_name(f"{path.name}.gz")
        logger.info(
            f"{path.name}: gzip reduced size from {json_size_mb:.1f} to {sys.getsizeof(json_bytes) / 1_000_000:.1f} MB"
        )
    path.write_bytes(json_bytes)
    logger.info(f"Wrote ontology to {path}")
    return path


def write_dataframe(df: pd.DataFrame, path: Path) -> None:
    df.to_json(
        path,
        orient="records",
        compression={"method": "gzip", "mtime": 0},
        indent=2,
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

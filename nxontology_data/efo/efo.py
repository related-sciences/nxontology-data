import functools
import json
import logging
import shutil
from pathlib import Path
from typing import Any

import fsspec
import networkx as nx
import pandas as pd
import rdflib
from nxontology import NXOntology

from nxontology_data.utils import (
    get_source_output_dir,
    normalize_parsed_curie,
    sparql_results_to_df,
    write_dataframe,
    write_ontology,
)

logger = logging.getLogger(__name__)


class EfoProcessor:
    name: str
    version: str
    EFO_REPO = "https://github.com/EBISPOT/efo"

    def __init__(
        self, name: str = "efo_otar_profile", version: str = "v3.52.0"
    ) -> None:
        self.name = name
        # TODO: support bioversions lookup
        self.version = version

    @property
    def owl_url(self) -> str:
        """
        Get download URL from an EFO GitHub release for a specific version and filename.
        Example return URL:
        https://github.com/EBISPOT/efo/releases/download/v3.52.0/efo_otar_profile.owl
        """
        return f"{self.EFO_REPO}/releases/download/{self.version}/{self.name}.owl"

    @property
    def owl_path(self) -> Path:
        return get_source_output_dir("efo").joinpath(f"{self.name}.owl.xz")

    def download_owl(self) -> Path:
        """
        Download an EFO release file from GitHub with compression.
        """
        logger.info(f"Downloading {self.owl_url} to {self.owl_path}")
        with fsspec.open(self.owl_url, mode="rb") as src, fsspec.open(
            self.owl_path, mode="wb", compression="infer"
        ) as dst:
            shutil.copyfileobj(src, dst)
        return self.owl_path

    @functools.cache  # noqa: B019
    def load_rdf(self) -> rdflib.graph.Graph:
        """
        Read raw EFO ontology as RDF graph.
        """
        logger.info(f"Loading {self.owl_path} into rdflib")
        rdf = rdflib.Graph()
        with fsspec.open(self.owl_path, "rt", compression="infer") as read_file:
            rdf.parse(source=read_file, format="xml")
        logger.info("Loading complete.")
        return rdf

    @staticmethod
    def _get_query(name: str) -> str:
        """
        Read SPARQL query from text file.
        """
        return Path(__file__).parent.joinpath(f"queries/{name}.rq").read_text()

    def run_query(self, name: str, cache: bool = False) -> pd.DataFrame:
        """
        Run SPARQL query on an rdflib.Graph instance of th4e MeSH RDF,
        returning the results as a pandas.DataFrame.
        Enable cache to cache results by the query text (not name/path).
        """
        rdf = self.load_rdf()
        if cache and not hasattr(rdf, "cached_query"):
            rdf.cached_query = functools.cache(rdf.query)
        query = self._get_query(name)
        results = (rdf.cached_query if cache else rdf.query)(query)
        return sparql_results_to_df(results)

    def get_terms_df(self) -> pd.DataFrame:
        return self.run_query("terms", cache=True)

    def get_subclass_df(self) -> pd.DataFrame:
        return self.run_query("subclasses", cache=True)

    def get_obsolete_df(self) -> pd.DataFrame:
        return self.run_query("terms_obsolete", cache=True)

    def get_alt_id_df(self) -> pd.DataFrame:
        return self.run_query("alt_id", cache=True)

    def get_xrefs_df(self) -> pd.DataFrame:
        xref_df = self.run_query("xrefs", cache=True)
        xref_df["xref_bioregistry"] = xref_df.apply(
            lambda row: normalize_parsed_curie(row.xref_prefix, row.xref_accession),
            axis="columns",
        )
        return xref_df

    def get_replaced_terms(self) -> dict[str, list[str]]:
        """Get a mapping from current EFO terms to their alternative IDs or replaced obsolete terms."""
        logger.info("Generating replaced terms")
        old_to_new = dict(
            self.get_obsolete_df()[["efo_id", "replaced_by_efo_id"]]
            .dropna()
            .query("efo_id != replaced_by_efo_id")
            .to_dict("split")["data"]
        )
        logger.info(
            f"Loaded {len(old_to_new):,} old term to new term mappings from obsolete replacements."
        )
        current_terms = set(self.get_terms_df()["efo_id"])
        for row in self.get_alt_id_df().itertuples():
            if row.alt_id in current_terms:
                continue
            if row.alt_id in row.efo_id:
                continue
            if row.efo_id not in current_terms:
                continue
            old_to_new[row.alt_id] = row.efo_id
        logger.info(
            f"Loaded alternative IDs: old_to_new now contains {len(old_to_new):,} items."
        )

        def update_term(old_term: str) -> str:
            new_term = old_to_new.get(old_term)
            if new_term is None:
                return old_term
            return update_term(new_term)

        current_to_old: dict[str, set[str]] = {}
        for old_term, new_term in old_to_new.items():
            if old_term in current_terms:
                # not expected to happen
                continue
            newest_term = update_term(new_term)
            if newest_term in current_terms:
                current_to_old.setdefault(new_term, set()).add(old_term)
        logger.info(
            f"{len(current_to_old):,} current terms have 1 or more replaced/alternative terms."
        )
        return {k: sorted(v) for k, v in current_to_old.items()}

    def get_nodes(self) -> list[dict[str, Any]]:
        logger.info("Generating nodes")
        node_df = self.get_terms_df()
        node_df["replaces"] = node_df.efo_id.map(self.get_replaced_terms())
        node_df["xrefs"] = node_df.efo_id.map(
            self.get_xrefs_df()
            .query("efo_id != xref_bioregistry")
            .groupby("efo_id")
            .apply(lambda df: sorted(set(df.xref_bioregistry.dropna())))
        )
        # Use .to_json and not .to_dict to convert NaN to None
        return json.loads(node_df.to_json(orient="records"))  # type: ignore [no-any-return]

    def create_nxo(self) -> NXOntology[str]:
        nxo: NXOntology[str] = NXOntology()
        nxo.graph.graph["name"] = self.name
        nxo.graph.graph["version"] = self.version
        nxo.graph.graph["source_url"] = self.owl_url
        nxo.set_graph_attributes(
            node_name_attribute="efo_label",
            node_identifier_attribute="{node}",
            node_url_attribute="efo_uri",
        )
        for data in self.get_nodes():
            nxo.add_node(data["efo_id"], **data)
        for edge in self.get_subclass_df().to_dict(orient="records"):
            source = edge["efo_id"]
            target = edge["child_efo_id"]
            try:
                nxo.add_edge(source, target)
            except nx.NodeNotFound as e:
                logger.warning(f"Skipping edge {source}, {target} due to: {e}")
        logging.info(
            f"Created {nxo.__class__.__name__} for {self.name} {self.version} with "
            f"{nxo.n_nodes:,} nodes and {nxo.graph.number_of_edges():,} edges."
        )
        return nxo

    def write_outputs(self) -> None:
        output_dir = get_source_output_dir("efo")
        nxo = self.create_nxo()
        write_ontology(nxo, output_dir)
        write_dataframe(self.get_xrefs_df(), output_dir.joinpath("efo_xrefs.json.gz"))
        write_dataframe(
            self.get_obsolete_df(), output_dir.joinpath("efo_obsolete.json")
        )


def process_efo() -> None:
    processor = EfoProcessor()
    processor.download_owl()
    processor.write_outputs()

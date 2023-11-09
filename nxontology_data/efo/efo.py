import functools
import json
import logging
import re
import shutil
from pathlib import Path
from typing import Any

import bioversions
import curies
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
        self, name: str = "efo_otar_profile", version: str | None = "current"
    ) -> None:
        """
        name: variant of efo. Valid options include 'efo', 'efo_otar_profile', and 'efo_otar_slim'.
              Note that efo_otar_slim gets computed by nxontology-data from efo_otar_profile.
        version: EFO version to use like 'v3.52.0'. If None, use the latest version from bioversions.
              If version='current', download current tag from GitHub <https://github.com/EBISPOT/efo/releases/tag/current>.
              See owl_version property for version extracted from the OWL download.
        """
        self.name = name
        if version is None:
            # WARNING: Bioregistry version is out of date
            version = bioversions.get_version("efo")
            if not version.startswith("v"):
                version = f"v{version}"
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

    @functools.cached_property
    def owl_version(self) -> str | None:
        """
        Get the version from the versionIRI from the EFO OWL file.
        https://github.com/EBISPOT/efo/issues/1972
        """
        with fsspec.open(self.owl_path, "rt", compression="infer") as read_file:
            # <owl:versionIRI rdf:resource="http://www.ebi.ac.uk/efo/releases/v3.52.0/efo.owl"/>
            for line in read_file:
                line = line.strip()
                if line.startswith("<owl:versionIRI rdf:resource"):
                    break
            else:
                logger.warning(f"Could not find a versionIRI line in {self.owl_path}")
                # efo_otar_profile.owl and efo_otar_slim.owl missing versionIRI
                # https://github.com/EBISPOT/efo/issues/1972
                return None
        logger.info(f"Found versionIRI line: {line}")
        match = re.match(
            r"<owl:versionIRI rdf:resource=\"http://www.ebi.ac.uk/efo/releases/(?P<version>v\d+\.\d+\.\d+)/\w+.owl\"/>",
            line,
        )
        assert match is not None
        version = match.group("version")
        logger.info(f"Extracted version from versionIRI: {version!r}")
        return version

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

    def get_xref_sources_df(self) -> pd.DataFrame:
        return self.run_query("xref_sources", cache=True)

    def get_mapping_properties_df(self) -> pd.DataFrame:
        converter = curies.get_bioregistry_converter()

        converter.add_prefix(
            "icd10cm-missing-prefix", "http://purl.bioontology.org/ontology/ICD10CM/"
        )

        df = (
            self.run_query("mapping_properties", cache=True)
            .assign(
                xref_id=lambda df: df["xref_id"].apply(
                    lambda xref: converter.compress(xref)
                )
            )
            .dropna()
            .assign(
                xref_id=lambda df: df["xref_id"]
                .str.replace("icd10cm-missing-prefix:", "icd10cm:")
                .str.replace("obo:Orphanet_", "Orphanet:")
                .str.split(":", expand=True)
                .apply(
                    lambda row: normalize_parsed_curie(
                        xref_prefix=row[0],
                        xref_accession=row[1],
                        collapse_orphanet=True,
                    ),
                    axis="columns",
                )
            )
        )

        return df

    def get_synonyms(self) -> dict[str, dict[str, str]]:
        synonym_scopes = {
            "hasExactSynonym": "exact",
            "hasNarrowSynonym": "narrow",
            "hasBroadSynonym": "broad",
            "hasRelatedSynonym": "related",
        }
        df = self.run_query("synonyms", cache=True)
        df["scope"] = df.predicate_id.map(synonym_scopes)
        df = (
            df.rename(columns={"synonym": "name"})[["efo_id", "name", "scope"]]
            .dropna()
            .drop_duplicates()
            .sort_values(["efo_id", "name", "scope"])
        )
        return {
            k: v[["name", "scope"]].to_dict(orient="records")
            for k, v in df.groupby("efo_id")
        }

    def get_subsets(self) -> dict[str, list[str]]:
        df = self.run_query("subsets", cache=True)
        return {
            k: sorted(set(v["subset_id"].dropna())) for k, v in df.groupby("efo_id")
        }

    def get_xrefs_df(self) -> pd.DataFrame:
        xref_df = self.run_query("xrefs", cache=True)
        xref_df["xref_bioregistry"] = xref_df.apply(
            lambda row: normalize_parsed_curie(
                row.xref_prefix, row.xref_accession, collapse_orphanet=True
            ),
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

    def get_xref_details(self) -> dict[str, dict[str, str | list[str] | None]]:
        xrefs = self.get_xrefs_df()[["efo_id", "xref_bioregistry"]].rename(
            columns={"xref_bioregistry": "xref_id"}
        )

        xref_sources = (
            self.get_xref_sources_df()
            .assign(
                xref_id=lambda df: df["xref"]
                .str.split(":", expand=True)
                .apply(
                    lambda row: normalize_parsed_curie(
                        xref_prefix=row[0],
                        xref_accession=row[1],
                        collapse_orphanet=True,
                    ),
                    axis="columns",
                )
            )
            .groupby(["efo_id", "xref_id"])["axiom_source"]
            .apply(list)
            .reset_index()
            .rename(columns={"axiom_source": "sources"})
        )

        def get_relation(x: list[str]) -> str | None:
            if "skos:exactMatch" in x or "mondo:exactMatch" in x:
                return "skos:exactMatch"
            if "skos:closeMatch" in x or "mondo:closeMatch" in x:
                return "skos:closeMatch"
            return None

        mapping_properties = (
            self.get_mapping_properties_df()
            .groupby(["efo_id", "xref_id"])["mapping_property_id"]
            .apply(list)
            .reset_index()
            .rename(columns={"mapping_property_id": "mapping_properties"})
            .assign(
                relation=lambda x: x["mapping_properties"].apply(get_relation),
            )
        )

        xref_details = (
            xrefs.merge(
                mapping_properties,
                how="outer",
                on=["efo_id", "xref_id"],
            )
            .merge(
                xref_sources,
                how="outer",
                on=["efo_id", "xref_id"],
            )
            .query("efo_id != xref_id")
        )

        return {
            k: v[["xref_id", "sources", "relation"]].to_dict(orient="records")
            for k, v in xref_details.groupby("efo_id")
        }

    def get_nodes(self) -> list[dict[str, Any]]:
        logger.info("Generating nodes")
        node_df = self.get_terms_df()
        node_df["synonyms"] = node_df.efo_id.map(self.get_synonyms())
        node_df["replaces"] = node_df.efo_id.map(self.get_replaced_terms())
        node_df["xrefs"] = node_df.efo_id.map(
            self.get_xrefs_df()
            .query("efo_id != xref_bioregistry")
            .groupby("efo_id")
            .apply(lambda df: sorted(set(df.xref_bioregistry.dropna())))
        )
        node_df["subsets"] = node_df.efo_id.map(self.get_subsets())
        node_df["xref_details"] = node_df.efo_id.map(self.get_xref_details())
        # Use .to_json and not .to_dict to convert NaN to None
        return json.loads(node_df.to_json(orient="records"))  # type: ignore [no-any-return]

    def create_nxo(self) -> NXOntology[str]:
        nxo: NXOntology[str] = NXOntology()
        nxo.graph.graph["name"] = self.name
        nxo.graph.graph["version"] = self.owl_version
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
            f"Created {nxo.__class__.__name__} for {self.name} {self.owl_version} with "
            f"{nxo.n_nodes:,} nodes and {nxo.graph.number_of_edges():,} edges."
        )
        return nxo

    def write_outputs(self) -> None:
        output_dir = get_source_output_dir("efo")
        nxo = self.create_nxo()
        write_ontology(nxo, output_dir)
        write_dataframe(
            self.get_xrefs_df(), output_dir.joinpath(f"{self.name}_xrefs.json.gz")
        )
        write_dataframe(
            self.get_obsolete_df(), output_dir.joinpath(f"{self.name}_obsolete.json.gz")
        )
        write_dataframe(
            self.get_mapping_properties_df(),
            output_dir.joinpath(f"{self.name}_mapping_properties.json.gz"),
        )
        write_dataframe(
            self.get_xref_sources_df(),
            output_dir.joinpath(f"{self.name}_xref_sources.json.gz"),
        )
        if nxo.name == "efo_otar_profile":
            nxo_slim = self.create_slim_nxo(nxo)
            write_ontology(nxo_slim, output_dir, compression_threshold_mb=30.0)

    @staticmethod
    def create_slim_nxo(nxo: NXOntology[str]) -> NXOntology[str]:
        """
        EFO OTAR Slim is created by pruning EFO OTAR Profile to only include therapeutic area terms and their descendants.
        https://github.com/EBISPOT/efo/issues/926
        """
        logger.info("Creating EFO OTAR slim")
        assert nxo.name == "efo_otar_profile"
        otar_slim_nodes = set()
        for node, data in nxo.graph.nodes(data=True):
            if data.get("therapeutic_area"):
                otar_slim_nodes |= nxo.node_info(node).descendants
        nxo_slim: NXOntology[str] = NXOntology(
            nxo.graph.subgraph(otar_slim_nodes).copy()
        )
        nxo_slim.graph.graph["name"] = "efo_otar_slim"
        nxo_slim.graph.graph[
            "note"
        ] = "EFO OTAR Slim was created from EFO OTAR Profile by nxontology-data."
        return nxo_slim


def process_efo(
    name: str = "efo_otar_profile", version: str | None = "current"
) -> None:
    processor = EfoProcessor(name=name, version=version)
    processor.download_owl()
    processor.write_outputs()


def process_efo_all(version: str | None = "current") -> None:
    for name in "efo", "efo_otar_profile":
        process_efo(name=name, version=version)

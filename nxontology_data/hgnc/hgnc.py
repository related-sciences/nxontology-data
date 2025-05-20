import json
import logging
import zipfile
from pathlib import Path
from typing import Any

import pandas as pd
import requests
from nxontology import NXOntology

from nxontology_data.utils import get_source_output_dir, write_ontology

logger = logging.getLogger(__name__)


def get_hgnc_output_dir() -> Path:
    return get_source_output_dir(source="hgnc")


class HgncGeneGroupDownloader:
    BASE_URL = "https://storage.googleapis.com/public-download-files/hgnc/csv/csv/genefamily_db_tables/"
    FILENAMES = [
        "README.txt",
        "external_resource.csv",
        "family.csv",
        "family_alias.csv",
        "family_has_external_resource.csv",
        "gene_has_family.csv",
        "hierarchy.csv",
        "hierarchy_closure.csv",
    ]
    OUTPUT_FILENAME = "hgnc_gene_group_db_tables.zip"
    SYMBOL_URL = "https://www.genenames.org/cgi-bin/download/custom?col=gd_hgnc_id&col=gd_app_sym&status=Approved&hgnc_dbtag=on&order_by=gd_hgnc_id&format=text&submit=submit"
    """Custom download of gene symbols from HGNC."""

    @classmethod
    def download_zip(cls) -> Path:
        """
        Download all files in genefamily_db_tables to a zip archive.
        Also downloads a custom query for gene symbols.
        """
        zip_path = get_hgnc_output_dir().joinpath(cls.OUTPUT_FILENAME)
        # download all files to a zip archive
        with zipfile.ZipFile(
            zip_path, mode="w", compression=zipfile.ZIP_DEFLATED
        ) as zip_file:
            # genefamily_db_tables downloads
            for filename in cls.FILENAMES:
                url = cls.BASE_URL + filename
                logger.info(f"Downloading {url}")
                response = requests.get(url)
                response.raise_for_status()
                zip_file.writestr(filename, response.content)
            # custom gene symbol download
            logger.info(f"Downloading {cls.SYMBOL_URL}")
            response = requests.get(cls.SYMBOL_URL)
            response.raise_for_status()
            zip_file.writestr("gene_symbols.csv", response.text.replace("\t", ","))
        return zip_path


class HgncGeneGroupNxoLoader:
    @classmethod
    def load_tables(cls) -> dict[str, pd.DataFrame]:
        zip_path = HgncGeneGroupDownloader.download_zip()
        with zipfile.ZipFile(zip_path, mode="r") as zip_file:
            return {
                name.removesuffix(".csv"): pd.read_csv(zip_file.open(name))
                for name in zip_file.namelist()
                if name.endswith(".csv")
            }

    @classmethod
    def export_hgnc_outputs(cls) -> None:
        tables = HgncGeneGroupNxoLoader.load_tables()
        nxo = cls._create_nxo_from_tables(tables)
        # set a higher compression threshold, because the git diff will help monitor for changes.
        write_ontology(
            nxo=nxo, output_dir=get_hgnc_output_dir(), compression_threshold_mb=25.0
        )

    @classmethod
    def _create_nxo_from_tables(
        cls, tables: dict[str, pd.DataFrame]
    ) -> NXOntology[int]:
        symbol_map = (
            tables["gene_symbols"]
            .copy()
            .set_index("HGNC ID")
            .to_dict()["Approved symbol"]
        )
        logger.info("Creating HGNC Gene Group NXOntology")
        nxo: NXOntology[int] = NXOntology()
        nxo.graph.graph["name"] = "hgnc_gene_group"
        nxo.graph.graph["description"] = "HGNC Gene Group / Family Ontology"
        nxo.graph.graph["data_license_source"] = (
            "https://www.genenames.org/about/license/"
        )
        nxo.graph.graph["data_license_url"] = (
            "https://creativecommons.org/publicdomain/zero/1.0/"
        )
        nxo.graph.graph["data_license_spdx_id"] = "CC0-1.0"
        nxo.set_graph_attributes(
            node_identifier_attribute="id", node_name_attribute="name"
        )
        for node_data in cls._get_nodes(tables):
            nxo.add_node(node_data["id"], **node_data)
        for edge in tables["hierarchy"].itertuples(index=False):
            nxo.add_edge(edge.parent_fam_id, edge.child_fam_id)
        for node in nxo.graph.nodes:
            node_info = nxo.node_info(node)
            genes_closure = set()
            for descendant in node_info.descendants:
                genes_direct = nxo.node_info(descendant).data["genes_direct"]
                if genes_direct:
                    genes_closure |= set(genes_direct)
                node_info.data["genes_closure"] = sorted(genes_closure)
            node_info.data["genes_direct_count"] = len(node_info.data["genes_direct"])
            node_info.data["genes_closure_count"] = len(node_info.data["genes_closure"])
        for node in nxo.graph.nodes:
            # expand genes_direct/genes_closure from strs to dicts to include symbols
            node_info = nxo.node_info(node)
            node_info.data["genes_direct"] = cls._add_symbols(
                hgnc_ids=node_info.data["genes_direct"], symbol_map=symbol_map
            )
            node_info.data["genes_closure"] = cls._add_symbols(
                hgnc_ids=node_info.data["genes_closure"], symbol_map=symbol_map
            )
        return nxo

    @staticmethod
    def _add_symbols(
        hgnc_ids: list[str], symbol_map: dict[str, str]
    ) -> list[dict[str, str | None]]:
        genes = [
            {
                "hgnc_id": hgnc_id,
                "symbol": symbol_map.get(hgnc_id),
            }
            for hgnc_id in hgnc_ids
        ]
        genes.sort(key=lambda x: x["symbol"] or "")
        return genes

    @classmethod
    def _get_nodes(cls, tables: dict[str, pd.DataFrame]) -> list[dict[str, Any]]:
        df = (
            tables["family"]
            .sort_values("id")
            .rename(columns={"abbreviation": "root_symbol"})
        )
        desc_source = df["desc_source"].str.split("|", n=1, regex=False, expand=True)
        df["desc_source"] = desc_source[0]
        df["desc_source_url"] = desc_source[1]
        df["pubmed_ids"] = df["pubmed_ids"].str.split(",")
        df["name_aliases"] = df["id"].map(cls._get_aliases(tables))
        df["external_resources"] = df["id"].map(cls._get_external_resources(tables))
        df["genes_direct"] = (
            df["id"]
            .map(cls._get_gene_assignments(tables))
            .apply(lambda x: x if isinstance(x, list) else [])
        )
        # reorder columns
        df = df[
            [
                "id",
                "name",
                "name_aliases",
                "root_symbol",
                "typical_gene",
                "desc_label",
                "desc_comment",
                "desc_source",
                "desc_source_url",
                "desc_go",
                "pubmed_ids",
                "external_note",
                "external_resources",
                "genes_direct",
            ]
        ]
        # Use .to_json and not .to_dict to convert NaN to None
        return json.loads(df.to_json(orient="records"))  # type: ignore [no-any-return]

    @staticmethod
    def _get_external_resources(
        tables: dict[str, pd.DataFrame]
    ) -> dict[int, list[dict[str, Any]]]:
        return (  # type: ignore [no-any-return]
            tables["family_has_external_resource"]
            .merge(tables["external_resource"].add_prefix("ext_"))
            .sort_values(["family_id", "ext_id"])
            .groupby("family_id")
            .apply(lambda df: df.drop(columns=["family_id"]).to_dict(orient="records"))
            .to_dict()
        )

    @staticmethod
    def _get_aliases(tables: dict[str, pd.DataFrame]) -> dict[int, list[str]]:
        return (  # type: ignore [no-any-return]
            tables["family_alias"]
            .sort_values(["family_id", "id"])
            .groupby("family_id")
            .apply(lambda df: df["alias"].tolist())
            .to_dict()
        )

    @staticmethod
    def _get_gene_assignments(tables: dict[str, pd.DataFrame]) -> dict[int, list[int]]:
        return (  # type: ignore [no-any-return]
            tables["gene_has_family"]
            .sort_values(["family_id", "hgnc_id"])
            .groupby("family_id")
            .apply(lambda df: [f"HGNC:{x}" for x in df["hgnc_id"]])
        )

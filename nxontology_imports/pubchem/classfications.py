import gzip
import json
import logging
import re
import sys
from pathlib import Path
from typing import Any, Counter

import requests
from networkx.readwrite.json_graph import node_link_data
from nxontology import NXOntology

logger = logging.getLogger(__name__)


root_dir = Path(__file__).parent.parent.parent
output_dir = root_dir.joinpath("output", "pubchem")


class PubchemClassificationApi:
    """
    For a JSON index of all PubChem classifications:
    <https://pubchem.ncbi.nlm.nih.gov/classification/cgi/classifications.fcgi?hid=index&format=json>

    For a web interface to PubChem classification hierarchies:
    <https://pubchem.ncbi.nlm.nih.gov/classification/#hid=87>

    Possibly unrelated PubChem API docs:
    <https://pubchemdocs.ncbi.nlm.nih.gov/pug-rest$classification_nodes>
    """

    REST_API = (
        "https://pubchem.ncbi.nlm.nih.gov/classification/cgi/classifications.fcgi"
    )

    @staticmethod
    def _check_hierarchy_type(hierarchy: Any) -> dict[str, Any]:
        assert isinstance(hierarchy, dict)
        for key in hierarchy:
            assert isinstance(key, str)
        return hierarchy

    @classmethod
    def get_hierarchy_catalog(cls) -> list[dict[str, Any]]:
        """Return API index of all PubChem classification heirarchies."""
        params = {
            "format": "json",
            "hid": "index",
        }
        response = requests.get(cls.REST_API, params=params)
        logger.info(f"Queried for the pubchem hierarchy catalog at {response.url}")
        hierarchies = response.json()["Hierarchies"]["Hierarchy"]
        assert isinstance(hierarchies, list)
        for hierarchy in hierarchies:
            cls._check_hierarchy_type(hierarchy)
        return hierarchies

    @classmethod
    def get_hierarchy(cls, hierarchy_id: int) -> dict[str, Any]:
        """
        For chembl protein class, hierarchy_id=87.
        """
        params = {
            "format": "json",
            "hid": hierarchy_id,
            # "depth": 30,  # max depth
            "start": "root",
        }
        response = requests.get(cls.REST_API, params=params)
        if not response.ok:
            logger.info(
                f"Response failed with code {response.status_code} for {response.url}"
            )
        response.raise_for_status()
        logger.info(f"Queried for pubchem hierarchy {hierarchy_id} at {response.url}")
        hierarchy = response.json()["Hierarchies"]["Hierarchy"][0]
        cls._check_hierarchy_type(hierarchy)
        return hierarchy  # type: ignore [no-any-return]

    @classmethod
    def get_metadata(cls, hierarchy: dict[str, Any]) -> dict[str, Any]:
        """
        Extract a dictionary of hierarchy metadata from the Information field.
        """
        info = hierarchy["Information"]
        return {
            "name": cls._get_simple_name(hierarchy),
            "pubchem_hierarch_id": int(hierarchy["HID"]),
            "pubchem_source_id": hierarchy["SourceID"],
            "pubchem_source_name": hierarchy["SourceName"],
            "pubchem_description": info.get("Description"),
            "pubchem_comments": info.get("Comments"),
            "source_url": info.get("URL"),
        }

    @staticmethod
    def convert_node_id(node: str) -> int:
        if node == "root":
            raise ValueError("Not intended to be used on the root node.")
        return int(node.removeprefix("node_"))

    @classmethod
    def create_nxo(cls, hierarchy_id: int) -> NXOntology[int]:
        hierarchy = cls.get_hierarchy(hierarchy_id=hierarchy_id)
        nxo: NXOntology[int] = NXOntology()
        nxo.graph.graph.update(cls.get_metadata(hierarchy))
        # sort nodes for cleaner output
        nodes = hierarchy["Node"]
        nodes.sort(key=lambda node: cls.convert_node_id(node["NodeID"]))
        # add nodes
        for node in nodes:
            node_id = cls.convert_node_id(node["NodeID"])
            info = node["Information"]
            if description := info.get("Description"):
                if isinstance(description, list):
                    description = description[0]
                assert isinstance(description, str)
            nxo.add_node(
                node_id,
                name=info.get("Name"),
                description=description,
                pubchem_hnid=info["HNID"],
                url=info.get("URL"),
            )
        # add edges
        for node in nodes:
            node_id = cls.convert_node_id(node["NodeID"])
            for parent in node["ParentID"]:
                if parent == "root":
                    # root appears to be a placeholder node added by PubChem
                    continue
                parent_id = cls.convert_node_id(parent)
                nxo.add_edge(parent_id, node_id)
        return nxo

    @classmethod
    def _get_simple_name(cls, hierarchy: dict[str, Any]) -> str:
        sep = "_"
        name = "{HID:03d} {SourceName} {SourceID}".format(**hierarchy)
        name = name.lower()
        # replace non-alphanumeric characters with the separator
        name = "".join(c if c.isalnum() else sep for c in name)
        # collapse sequences of the separator
        name = re.sub(f"{re.escape(sep)}+", re.escape(sep), name)
        # remove duplicate words
        return sep.join(Counter(name.split(sep)))

    @classmethod
    def write_hierarchy_catalog(cls) -> list[dict[str, Any]]:
        hierarchies = cls.get_hierarchy_catalog()
        hierarchies.sort(key=lambda h: h["HID"])  # type: ignore [no-any-return]
        for hierarchy in hierarchies:
            simple_name = cls._get_simple_name(hierarchy)
            hierarchy["nxo_name"] = simple_name
        json_str = json.dumps(hierarchies, indent=2, ensure_ascii=False)
        output_dir.joinpath("catalog.json").write_text(json_str)
        return hierarchies


skip_hierarchy_ids = [
    2,  # 002_chebi_obo fails with status 500
]


def write_ontology(nxo: NXOntology[Any]) -> Path:
    data = node_link_data(nxo.graph)
    json_bytes = json.dumps(data, indent=2, ensure_ascii=False).encode()
    json_size_mb = sys.getsizeof(json_bytes) / 1_000_000
    path = output_dir.joinpath(f"{nxo.graph.graph['name']}.json")
    if json_size_mb > 10.0:
        json_bytes = gzip.compress(json_bytes, mtime=0)
        path = path.with_name(f"{path.name}.gz")
        logging.info(
            f"{path.name}: gzip reduced size from {json_size_mb:.1f} to {sys.getsizeof(json_bytes) / 1_000_000:.1f} MB"
        )
    path.write_bytes(json_bytes)
    return path


def export_all_heirarchies() -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    hierarchies = PubchemClassificationApi.write_hierarchy_catalog()
    for hierarchy in hierarchies:
        nxo_name = hierarchy["nxo_name"]
        if hierarchy["HID"] in skip_hierarchy_ids:
            logging.info(f"Skipping {nxo_name} since it's in the skip list.")
            continue
        if list(output_dir.glob(f"{nxo_name}.*")):
            logging.info(f"Skipping {nxo_name} since output file exists.")
            continue
        logging.info(f"Beginning create_nxo for {nxo_name}.")
        try:
            nxo = PubchemClassificationApi.create_nxo(hierarchy_id=hierarchy["HID"])
        except requests.HTTPError:
            logging.info(f"Skipping {nxo_name} because request failed.")
            continue
        write_ontology(nxo)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    export_all_heirarchies()

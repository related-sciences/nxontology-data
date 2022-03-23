from pathlib import Path
from typing import Any

import requests
from nxontology import NXOntology


class PubchemClassificationApi:
    """
    https://pubchemdocs.ncbi.nlm.nih.gov/pug-rest$classification_nodes
    """

    REST_API = (
        "https://pubchem.ncbi.nlm.nih.gov/classification/cgi/classifications.fcgi"
    )

    @classmethod
    def get_results(cls, hierarchy_id: int = 87) -> dict[str, Any]:
        """
        For chembl protein class, hierarchy_id=87.
        """
        params = {
            "format": "json",
            "hid": hierarchy_id,
            "depth": 30,  # max depth
            "start": "root",
        }
        response = requests.get(cls.REST_API, params=params)
        hierarchy = response.json()["Hierarchies"]["Hierarchy"][0]
        assert isinstance(hierarchy, dict)
        for key in hierarchy:
            assert isinstance(key, str)
        return hierarchy

    @classmethod
    def get_metadata(cls, hierarchy: dict[str, Any]) -> dict[str, Any]:
        """
        Extract a dictionary of hierarchy metadata from the Information field.
        """
        info = hierarchy["Information"]
        ["SourceName", "SourceID", "RootID", "HID", "Information", "Node"]
        return {
            "name": "{SourceName} {SourceID}".format(**hierarchy),
            "pubchem_hierarch_id": int(hierarchy["HID"]),
            "description": info["Description"],
        }


class ProteinClassImporter:
    PUBCHEM_HIERARCHY_ID = 87

    @staticmethod
    def convert_node_id(node: str) -> int:
        if node == "root":
            raise ValueError("Not intended to be used on the root node.")
        return int(node.removeprefix("node_"))

    @classmethod
    def get_pubchem_hierarchy(cls) -> dict[str, Any]:
        return PubchemClassificationApi.get_results(
            hierarchy_id=cls.PUBCHEM_HIERARCHY_ID
        )

    @classmethod
    def create_nxo(cls) -> NXOntology[int]:
        hierarchy = cls.get_pubchem_hierarchy()
        nxo: NXOntology[int] = NXOntology()
        nxo.graph.graph.update(PubchemClassificationApi.get_metadata(hierarchy))
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
                name=info["Name"],
                pubchem_hnid=info["HNID"],
                description=description,
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


if __name__ == "__main__":
    nxo = ProteinClassImporter.create_nxo()
    root_dir = Path(__file__).parent.parent.parent
    nxo.write_node_link_json(
        root_dir.joinpath("output", "chembl-protein-class-node-link.json").as_posix()
    )

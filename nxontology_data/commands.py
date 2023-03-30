import logging

import fire
from nxontology import NXOntology

from nxontology_data.hgnc.hgnc import HgncGeneGroupNxoLoader
from nxontology_data.mesh.mesh import MeshLoader
from nxontology_data.pubchem.classifications import export_all_heirarchies
from nxontology_data.utils import get_source_output_dir, write_ontology


def write_test_output() -> None:
    """Export an empty test ontology, useful for testing CI deployment."""
    output_dir = get_source_output_dir("test")
    nxo: NXOntology[str] = NXOntology()
    nxo.graph.graph["name"] = "test"
    write_ontology(nxo, output_dir)


def cli() -> None:
    """
    Run like `poetry run nxontology_data`
    """
    logging.basicConfig()
    logging.getLogger().setLevel(logging.INFO)
    commands = {
        "mesh": MeshLoader.export_mesh_outputs,
        "pubchem": export_all_heirarchies,
        "hgnc": HgncGeneGroupNxoLoader.export_hgnc_outputs,
        "test": write_test_output,
    }
    fire.Fire(commands)

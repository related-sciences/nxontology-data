import logging

import fire
from nxontology import NXOntology

from nxontology_data.mesh.mesh import MeshLoader
from nxontology_data.pubchem.classifications import export_all_heirarchies
from nxontology_data.utils import get_output_dir, write_ontology


def write_test_output() -> None:
    """Export an empty test ontology, useful for testing CI deployment."""
    output_dir = get_output_dir().joinpath("test")
    output_dir.mkdir(parents=True, exist_ok=True)
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
        "test": write_test_output,
    }
    fire.Fire(commands)

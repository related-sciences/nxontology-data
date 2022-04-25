import logging

import fire

from nxontology_data.mesh.mesh import MeshLoader
from nxontology_data.pubchem.classifications import export_all_heirarchies


def cli() -> None:
    """
    Run like `poetry run nxontology_data`
    """
    logging.basicConfig()
    logging.getLogger().setLevel(logging.INFO)
    commands = {
        "mesh": MeshLoader.export_mesh_outputs,
        "pubchem": export_all_heirarchies,
    }
    fire.Fire(commands)

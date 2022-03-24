import logging

import fire

from nxontology_imports.pubchem.classifications import export_all_heirarchies


def cli() -> None:
    """
    Run like `poetry run nxontology_import`
    """
    logging.basicConfig()
    logging.getLogger().setLevel(logging.INFO)
    commands = {
        "pubchem": export_all_heirarchies,
    }
    fire.Fire(commands)

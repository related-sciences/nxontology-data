from nxontology_imports.pubchem.classifications import PubchemClassificationApi


def test_get_hierarchy_metadata() -> None:
    hierarchy = {
        "SourceName": "ChEMBL",
        "SourceID": "Target Tree",
        "HID": 87,
        "Information": {
            "Name": "Target Tree",
            "Description": [
                "The ChEMBL Protein Target Tree is a structured classification of the protein target entities contained with the ChEMBL resource.",
                "Author: ChEMBL curation team",
                "ChEMBL Release version 29",
                "For any queries contact chembl-help@ebi.ac.uk",
                "Created on 01/31/2022 09:16:50",
            ],
            "URL": "https://www.ebi.ac.uk/chembl/target/browser",
        },
    }
    metadata = PubchemClassificationApi.get_metadata(hierarchy)
    assert metadata["name"] == "087_chembl_target_tree"
    assert metadata["pubchem_source_id"] == "Target Tree"
    assert metadata["pubchem_source_name"] == "ChEMBL"
    assert metadata["pubchem_hierarchy_id"] == 87
    assert metadata["pubchem_description"][0].startswith(
        "The ChEMBL Protein Target Tree"
    )
    assert metadata["pubchem_comments"] is None
    assert metadata["source_url"] == "https://www.ebi.ac.uk/chembl/target/browser"

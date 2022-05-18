import pathlib

import fsspec
import networkx as nx
import pytest
import rdflib
from nxontology import NXOntology

from nxontology_data.mesh.mesh import MeshLoader

test_data_dir = pathlib.Path(__file__).parent.joinpath("rdf-2020-subset")


def get_testing_mesh_rdf() -> rdflib.Graph:
    return MeshLoader._read_mesh_rdf(test_data_dir.as_posix(), "mesh2020-subset.nt")


@pytest.fixture
def rdf() -> rdflib.Graph:
    return get_testing_mesh_rdf()


def test_get_mesh_rdf(rdf: rdflib.Graph) -> None:
    assert "meshv" in dict(rdf.namespaces())


def test_create_nxo(rdf: rdflib.Graph) -> None:
    nxo, id_df = MeshLoader.create_nxo(rdf, year_yyyy="2020")
    assert nxo.graph.nodes["C000598941"]["mesh_label"] == "Keratoactinomycosis"
    assert nxo.graph.nodes["C000598941"]["mesh_class"] == "SCR_Disease"
    assert (
        nxo.graph.nodes["C000598941"]["mesh_uri"]
        == "http://id.nlm.nih.gov/mesh/2020/C000598941"
    )
    assert nxo.graph.nodes["C000598941"]["mesh_id"] == "C000598941"
    path = pathlib.Path(__file__).parent.joinpath("rdf-2020-subset/nxo-node-link.json")
    # uncomment to regenerate nxo-node-link.json
    # nxo.write_node_link_json(path.as_posix())
    expected: NXOntology[str] = NXOntology.read_node_link_json(path.as_posix())
    assert nx.is_isomorphic(nxo.graph, expected.graph)


def test_create_vocab_digraph(rdf: rdflib.Graph) -> None:
    vocab_nxo = MeshLoader.create_vocab_digraph(rdf)
    assert isinstance(vocab_nxo, nx.DiGraph)
    assert vocab_nxo.number_of_nodes() == 18


def create_testing_nt(full_nxo: NXOntology[str], full_nt_path: str) -> None:
    """
    Regenerate testing mesh2020-subset.nt from full mesh release.
    """
    source_mesh_id = "C000598941"
    nodes = full_nxo.node_info(source_mesh_id).ancestors
    with fsspec.open(full_nt_path, "rt", compression="infer") as rf, open(
        f"{test_data_dir}/mesh2020-subset.nt", "wt"
    ) as wf:
        for line in rf:
            for node in nodes:
                if f"/{node}>" in line:
                    wf.write(line)

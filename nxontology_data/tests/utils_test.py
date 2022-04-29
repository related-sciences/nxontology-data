import pytest
import rdflib

from nxontology_data.utils import get_output_dir, sparql_results_to_df


def test_get_output_dir() -> None:
    output_dir = get_output_dir()
    root = output_dir.parent
    # file that should only exist in root directory
    assert root.joinpath("pyproject.toml").exists()


@pytest.fixture
def rdflib_foaf_graph() -> rdflib.Graph:
    """
    FOAF (Friend of a Friend) testing graph from rdflib.
    """
    graph = rdflib.Graph()
    return graph.parse(
        source="https://github.com/RDFLib/rdflib/raw/56dc4207ce6e7b11ed7b45fb4fd4020ba548e718/examples/foaf.n3",
        format="n3",
    )


_foaf_sparql = """\
SELECT
  ?subject
  ?subject_is_tim
  (COUNT(*) AS ?n_triples)
  (MIN(?predicate) AS ?sample_predicate)
  (SAMPLE(?missing) AS ?missing) 
WHERE {
  ?subject ?predicate ?object.
  BIND(?subject = <http://www.w3.org/People/Berners-Lee/card#i> AS ?subject_is_tim)
  OPTIONAL {?subject <this_predicate_does_not_exist> ?missing .}
}
GROUP BY ?subject ?subject_is_tim
ORDER BY DESC(?n_triples) ?subject
LIMIT 10
"""


def test_sparql_results_to_df(rdflib_foaf_graph: rdflib.Graph) -> None:
    results = rdflib_foaf_graph.query(_foaf_sparql)
    df = sparql_results_to_df(results)
    assert len(df) == 10
    # test column values (no ? prefix), type (as strings), and order
    assert list(df.columns) == [
        "subject",
        "subject_is_tim",
        "n_triples",
        "sample_predicate",
        "missing",
    ]
    first_row = next(df.itertuples())
    # test value of subject, ensuring type conversion to str
    assert first_row.subject == "http://www.w3.org/People/Berners-Lee/card#i"
    # test value of subject_is_tim, ensuring type conversion to bool
    assert first_row.subject_is_tim is True
    # test value of n_triples, ensuring type conversion to int
    assert first_row.n_triples == 61
    # test value of sample_predicate, ensuring type conversion to str
    assert (
        first_row.sample_predicate == "http://www.w3.org/1999/02/22-rdf-syntax-ns#type"
    )
    # test value of missing, ensuring it's None
    assert first_row.missing is None

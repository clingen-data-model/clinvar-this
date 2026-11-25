"""Submission conversion must preserve its source statements."""

from pathlib import Path

from ga4gh.cat_vrs.models import CategoricalVariant, DefiningAlleleConstraint
from ga4gh.core.models import iriReference
from ga4gh.va_spec.base import EvidenceLine
from ga4gh.vrs.models import Allele

from clinvar_api.models import CollectionMethod
from clinvar_this.io.gkm_json import read_gkm_json_file
from clinvar_this.io.gkm_json.base import BatchMetadata


def test_conversion_preserves_defining_allele_expressions() -> None:
    """Member expressions must not accumulate on the defining allele."""
    path = (
        Path(__file__).parents[2]
        / "data/io_gkm_json/native/oncogenicity_civic_assertions.json"
    )
    transformer, statements = read_gkm_json_file(str(path))
    statement = statements[0]
    variant = statement.proposition.subject
    assert isinstance(variant, CategoricalVariant)
    assert variant.constraints
    constraint = variant.constraints[0].root
    assert isinstance(constraint, DefiningAlleleConstraint)
    assert isinstance(constraint.allele, Allele)
    assert constraint.allele.expressions

    # Keep a cDNA member so conversion extends the defining protein expressions.
    variant_data = variant.model_dump(mode="json")
    variant_data["members"] = [variant_data["members"][0]]
    statement.proposition.subject = CategoricalVariant.model_validate(variant_data)
    original = statement.model_dump(mode="json")
    metadata = BatchMetadata(collection_method=CollectionMethod.CURATION)

    first = transformer.records_to_submission_container([statement], metadata)
    assert statement.model_dump(mode="json") == original
    second = transformer.records_to_submission_container([statement], metadata)

    assert first.oncogenicity_submission
    assert (
        first.oncogenicity_submission[0].variant_set.variant[0].hgvs
        == "NM_020975.4:c.2753T>C"
    )
    assert second.model_dump() == first.model_dump()
    assert statement.model_dump(mode="json") == original


def test_conversion_preserves_statement_citations() -> None:
    """Evidence citations must not accumulate in statement-level references."""
    path = (
        Path(__file__).parents[2]
        / "data/io_gkm_json/native/oncogenicity_civic_assertions.json"
    )
    transformer, statements = read_gkm_json_file(str(path))
    statement = statements[0]
    statement.reportedIn = [iriReference("https://example.org/statement-source")]
    assert statement.hasEvidenceLines
    evidence_line = statement.hasEvidenceLines[0]
    assert isinstance(evidence_line, EvidenceLine)
    evidence_line.reportedIn = [iriReference("https://example.org/evidence-source")]
    original = statement.model_dump(mode="json")
    metadata = BatchMetadata(collection_method=CollectionMethod.CURATION)

    first = transformer.records_to_submission_container([statement], metadata)
    assert statement.model_dump(mode="json") == original
    second = transformer.records_to_submission_container([statement], metadata)

    citations = first.oncogenicity_submission[0].oncogenicity_classification.citation
    urls = {citation.url for citation in citations}
    assert "https://example.org/statement-source" in urls
    assert "https://example.org/evidence-source" in urls
    assert second.model_dump() == first.model_dump()
    assert statement.model_dump(mode="json") == original

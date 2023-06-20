from typing import Optional
import requests

from pydantic import BaseModel
import pandas as pd

from multiprocessing import Pool


class GenomicVariant(BaseModel):

    assembly: str
    chr: str
    position: int

    ref: str
    alt: str

class Prediction(BaseModel):

    variant: GenomicVariant
    score_value: str
    score_label: Optional[str] = None
    score_name: str

    score_source: str


HEADERS = {
    "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/114.0.0.0 Safari/537.36"
}


def get_cadd_results(variant: GenomicVariant):
    CADD_URL = "https://cadd.gs.washington.edu/snv/{assembly}-v1.6_inclAnno/{chr}:{position}_{ref}_{alt}"
    url = CADD_URL.format(assembly=variant.assembly, chr=variant.chr, position=variant.position, ref=variant.ref, alt=variant.alt)

    response = requests.get(url, headers=HEADERS)
    response.raise_for_status()
    tables = pd.read_html(response.text)

    cadd_data_table = tables[1].T

    cadd_data_table.columns = cadd_data_table.iloc[0]
    cadd_data_table = cadd_data_table[1:]

    predictions = [
        Prediction(
            variant=variant,
            score_value=cadd_data_table["PHRED"].max(),
            score_name="CADD_PHRED",
            score_source="CADD",
        )
    ]
    return predictions


def get_gnomad_results(variant: GenomicVariant):
    GNOMAD_URL = "https://gnomad.broadinstitute.org/api/"

    json_data = {
        'operationName': 'GnomadVariant',
        'query':"""
query GnomadVariant($variantId: String!, $datasetId: DatasetId!, $referenceGenome: ReferenceGenomeId!, $includeLocalAncestry: Boolean!, $includeLiftoverAsSource: Boolean!, $includeLiftoverAsTarget: Boolean!) {
  variant(variantId: $variantId, dataset: $datasetId) {
    variant_id
    reference_genome
    chrom
    pos
    ref
    alt
    caid
    colocated_variants
    coverage {
      exome {
        mean
      }
      genome {
        mean
      }
    }
    multi_nucleotide_variants {
      combined_variant_id
      changes_amino_acids
      n_individuals
      other_constituent_snvs
    }
    exome {
      ac
      an
      ac_hemi
      ac_hom
      faf95 {
        popmax
        popmax_population
      }
      filters
      populations {
        id
        ac
        an
        ac_hemi
        ac_hom
      }
      local_ancestry_populations @include(if: $includeLocalAncestry) {
        id
        ac
        an
      }
      age_distribution {
        het {
          bin_edges
          bin_freq
          n_smaller
          n_larger
        }
        hom {
          bin_edges
          bin_freq
          n_smaller
          n_larger
        }
      }
      quality_metrics {
        allele_balance {
          alt {
            bin_edges
            bin_freq
            n_smaller
            n_larger
          }
        }
        genotype_depth {
          all {
            bin_edges
            bin_freq
            n_smaller
            n_larger
          }
          alt {
            bin_edges
            bin_freq
            n_smaller
            n_larger
          }
        }
        genotype_quality {
          all {
            bin_edges
            bin_freq
            n_smaller
            n_larger
          }
          alt {
            bin_edges
            bin_freq
            n_smaller
            n_larger
          }
        }
        site_quality_metrics {
          metric
          value
        }
      }
    }
    genome {
      ac
      an
      ac_hemi
      ac_hom
      faf95 {
        popmax
        popmax_population
      }
      filters
      populations {
        id
        ac
        an
        ac_hemi
        ac_hom
      }
      local_ancestry_populations @include(if: $includeLocalAncestry) {
        id
        ac
        an
      }
      age_distribution {
        het {
          bin_edges
          bin_freq
          n_smaller
          n_larger
        }
        hom {
          bin_edges
          bin_freq
          n_smaller
          n_larger
        }
      }
      quality_metrics {
        allele_balance {
          alt {
            bin_edges
            bin_freq
            n_smaller
            n_larger
          }
        }
        genotype_depth {
          all {
            bin_edges
            bin_freq
            n_smaller
            n_larger
          }
          alt {
            bin_edges
            bin_freq
            n_smaller
            n_larger
          }
        }
        genotype_quality {
          all {
            bin_edges
            bin_freq
            n_smaller
            n_larger
          }
          alt {
            bin_edges
            bin_freq
            n_smaller
            n_larger
          }
        }
        site_quality_metrics {
          metric
          value
        }
      }
    }
    flags
    lof_curations {
      gene_id
      gene_symbol
      verdict
      flags
      project
    }
    rsids
    transcript_consequences {
      domains
      gene_id
      gene_version
      gene_symbol
      hgvs
      hgvsc
      hgvsp
      is_canonical
      is_mane_select
      is_mane_select_version
      lof
      lof_flags
      lof_filter
      major_consequence
      polyphen_prediction
      sift_prediction
      transcript_id
      transcript_version
    }
    in_silico_predictors {
      id
      value
      flags
    }
    non_coding_constraint {
      start
      stop
      possible
      observed
      expected
      oe
      z
    }
  }

  clinvar_variant(variant_id: $variantId, reference_genome: $referenceGenome) {
    clinical_significance
    clinvar_variation_id
    gold_stars
    last_evaluated
    review_status
    submissions {
      clinical_significance
      conditions {
        name
        medgen_id
      }
      last_evaluated
      review_status
      submitter_name
    }
  }

  liftover(source_variant_id: $variantId, reference_genome: $referenceGenome) @include(if: $includeLiftoverAsSource) {
    liftover {
      variant_id
      reference_genome
    }
    datasets
  }

  liftover_sources: liftover(liftover_variant_id: $variantId, reference_genome: $referenceGenome) @include(if: $includeLiftoverAsTarget) {
    source {
      variant_id
      reference_genome
    }
    datasets
  }

  meta {
    clinvar_release_date
  }
}
""",
        'variables': {
            'datasetId': 'gnomad_r2_1',
            'includeLocalAncestry': False,
            'includeLiftoverAsSource': True,
            'includeLiftoverAsTarget': False,
            'referenceGenome': variant.assembly,
            'variantId': f'{variant.chr}-{variant.position}-{variant.ref}-{variant.alt}',
        },
    }

    response = requests.post(GNOMAD_URL, json=json_data, headers=HEADERS)
    response.raise_for_status()
    resp_data = response.json()
    predictions = []
    if variant_data := resp_data["data"]["variant"]:
        exome_data = variant_data["exome"]
        genome_data = variant_data["genome"]

        total_an = exome_data["an"] + genome_data["an"]
        total_ac = exome_data["ac"] + genome_data["ac"]
        total_hom = exome_data["ac_hom"] + exome_data["ac_hemi"] + genome_data["ac_hom"] + genome_data["ac_hemi"]

        subpopulations = [
            { "id": p["id"] , **{
                key: sum(pp[key] for pp in exome_data["populations"] + genome_data["populations"] if pp["id"] == p["id"] and pp[key])
                for key in ("ac", "an", "ac_hemi", "ac_hom")
            }} for p in exome_data["populations"] if "_" not in p["id"]
        ]

        popmax = max(subpopulations, key=lambda p: p["ac"] / p["an"])

        predictions.append(Prediction(
            variant=variant,
            score_value=f"{(total_ac / total_an) * 100:.5f}%",
            score_label="",
            score_name="MAF",
            score_source="gnomAD 2.1",
        ))
        predictions.append(Prediction(
            variant=variant,
            score_value=f"{total_hom}",
            score_label="",
            score_name="HOM",
            score_source="gnomAD 2.1",
        ))
        predictions.append(Prediction(
            variant=variant,
            score_value=f"{(popmax['ac'] / popmax['an']) * 100:.5f}%",
            score_label=f"{popmax['id'].upper()}",
            score_name="POPMAX AF",
            score_source="gnomAD 2.1",
        ))

    return predictions


def get_mt21_results(variant: GenomicVariant):
    MT21_URL = "https://www.genecascade.org/MTc2021/ChrPos102.cgi"
    response = requests.get(MT21_URL, params={"chromosome": variant.chr, "position": variant.position, "ref": variant.ref, "alt": variant.alt}, headers=HEADERS)
    response.raise_for_status()
    tables = pd.read_html(response.text)
    predictions = [
        Prediction(
            variant=variant,
            score_value=pred["Tree vote"],
            score_label=pred["Prediction"],
            score_name=f"{pred['Transcript']}_{pred['Model']}",
            score_source="MutationTaster2021",
        )
        for pred in
        tables[0].to_dict("records")
    ]
    return predictions


def get_franklin_predictions(variant: GenomicVariant):
    FRANKLIN_URL = "https://franklin.genoox.com/api/fetch_predictions"
    json_data = {
        "chr": f"{variant.chr}",
        "pos": f"{variant.position}",
        "ref": variant.ref,
        "alt": variant.alt,
        "version":"",
        "analysis_id":"",
        "reference_version":variant.assembly
    }
    response = requests.post(FRANKLIN_URL, json=json_data, headers=HEADERS)
    response.raise_for_status()
    resp_data = response.json()

    if "prediction" in resp_data:
        predictions = [
            Prediction(
                score_name=algo_name,
                score_value=algo_results["score"],
                score_label=algo_results.get("prediction"),
                score_source="Franklin", variant=variant
            ) for algo_name, algo_results in resp_data["prediction"].items()
        ]
        return predictions
    return []


VARIANT_SOURCES = [
    get_cadd_results,
    get_mt21_results,
    get_franklin_predictions,
    get_gnomad_results,
]


v = GenomicVariant(assembly="GRCh37", chr="19", position=11233940, ref="G", alt="A")

def run(fun, var):
    return fun(var)

with Pool(10) as p:
    result = p.starmap(run, [(f, v) for f in VARIANT_SOURCES])

for r in result:
    print(r)

from typing import List, Optional, Dict, Tuple
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


class Scraper:

    @classmethod
    @property
    def name(cls):
        return cls.__name__


class CADD(Scraper):
    CADD_URL = "https://cadd.gs.washington.edu/snv/{assembly}-v1.6_inclAnno/{chr}:{position}_{ref}_{alt}"

    @classmethod
    def url(cls, variant: GenomicVariant):
        var_url = cls.CADD_URL.format(
            assembly=variant.assembly,
            chr=variant.chr,
            position=variant.position,
            ref=variant.ref,
            alt=variant.alt
        )
        return var_url

    @classmethod
    def query(cls, variant: GenomicVariant):
        response = requests.get(cls.url(variant), headers=HEADERS)
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
                score_source=cls.name,
            )
        ]
        return predictions


class GNOMAD(Scraper):

    @classmethod
    def url(cls, variant: GenomicVariant) -> str:
        GNOMAD_VARIANT_URL = "https://gnomad.broadinstitute.org/variant/{chr}-{position}-{ref}-{alt}?dataset=gnomad_r2_1"
        url = GNOMAD_VARIANT_URL.format(chr=variant.chr, position=variant.position, ref=variant.ref, alt=variant.alt)
        return url

    @classmethod
    def query(cls, variant: GenomicVariant):
        GNOMAD_URL = "https://gnomad.broadinstitute.org/api/"

        json_data = {
            'operationName': 'GnomadVariant',
            'query':"""
    query GnomadVariant($variantId: String!, $datasetId: DatasetId!, $includeLocalAncestry: Boolean!) {
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
        }
      }
    }
    """,
            'variables': {
                'datasetId': 'gnomad_r2_1',
                'includeLocalAncestry': False,
                'variantId': f'{variant.chr}-{variant.position}-{variant.ref}-{variant.alt}',
            },
        }

        response = requests.post(GNOMAD_URL, json=json_data, headers=HEADERS)
        response.raise_for_status()
        resp_data = response.json()
        predictions = []
        if variant_data := resp_data["data"]["variant"]:
            print(variant_data)
            exome_data = variant_data["exome"] or {"an": 0, "ac": 0, "ac_hom": 0, "ac_hemi": 0, "populations": []}
            genome_data = variant_data["genome"] or {"an": 0, "ac": 0, "ac_hom": 0, "ac_hemi": 0, "populations": []}

            total_an = exome_data["an"] + genome_data["an"]  # type: ignore
            total_ac = exome_data["ac"] + genome_data["ac"]  # type: ignore
            total_hom = exome_data["ac_hom"] + exome_data["ac_hemi"] + genome_data["ac_hom"] + genome_data["ac_hemi"]  # type: ignore

            subpopulations = [
                { "id": p["id"] , **{
                    key: sum(pp[key] for pp in exome_data["populations"] + genome_data["populations"] if pp["id"] == p["id"] and pp[key])  # type: ignore
                    for key in ("ac", "an", "ac_hemi", "ac_hom")
                }} for p in exome_data["populations"] if "_" not in p["id"]  # type: ignore
            ]

            popmax = max(subpopulations, key=lambda p: p["ac"] / p["an"])

            predictions.append(Prediction(
                variant=variant,
                score_value=f"{(total_ac / total_an) * 100:.5f}%",
                score_label="",
                score_name="MAF",
                score_source=cls.name,
            ))
            predictions.append(Prediction(
                variant=variant,
                score_value=f"{total_hom}",
                score_label="",
                score_name="HOM",
                score_source=cls.name,
            ))
            predictions.append(Prediction(
                variant=variant,
                score_value=f"{(popmax['ac'] / popmax['an']) * 100:.5f}%",
                score_label=f"{popmax['id'].upper()}",
                score_name="POPMAX AF",
                score_source=cls.name,
            ))

        return predictions


class MutationTaster2021(Scraper):

    @classmethod
    def url(cls, variant: GenomicVariant) -> str:
        MT21_VARIANT_URL = "https://www.genecascade.org/MTc2021/ChrPos102.cgi?chromosome={chr}&position={position}&ref={ref}&alt={alt}"
        url = MT21_VARIANT_URL.format(chr=variant.chr, position=variant.position, ref=variant.ref, alt=variant.alt)
        return url

    @classmethod
    def query(cls, variant: GenomicVariant):
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
                score_source=cls.name,
            )
            for pred in
            tables[0].to_dict("records")
        ]
        return predictions


class Franklin(Scraper):

    @classmethod
    def url(cls, variant: GenomicVariant) -> str:
        FRANKLIN_URL = "https://franklin.genoox.com/clinical-db/variant/snp/chr{chr}-{position}-{ref}-{alt}"
        url = FRANKLIN_URL.format(chr=variant.chr, position=variant.position, ref=variant.ref, alt=variant.alt)
        return url

    @classmethod
    def query(cls, variant: GenomicVariant):
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
                    score_source=cls.name, variant=variant
                ) for algo_name, algo_results in resp_data["prediction"].items()
            ]
            return predictions
        return []


VARIANT_SOURCES = [
    CADD, MutationTaster2021, Franklin, GNOMAD
]

NUM_PROCESSES = len(VARIANT_SOURCES)


def __run(cls, var) -> Tuple[str, List[Prediction]]:
    return (cls.name, cls.query(var))


def queryUrls(variant: GenomicVariant) -> Dict[str, str]:
    return {
        cls.name: cls.url(variant) for cls in VARIANT_SOURCES
    }

def queryVariant(variant: GenomicVariant) -> Dict[str, List[Prediction]]:
    with Pool(NUM_PROCESSES) as pool:
        predictions = dict(pool.starmap(__run, [(f, variant) for f in VARIANT_SOURCES]))
        return predictions

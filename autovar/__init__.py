import re
from typing import Any, Callable, List, Optional, Dict, Tuple
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

    def __str__(self):
        return f"chr{self.chr}:g.{self.position}{self.ref}>{self.alt}"


class Prediction(BaseModel):

    variant: GenomicVariant
    score_value: str
    score_label: Optional[str] = None
    score_name: str

    score_source: str


def get_prediction(predictions: List[Prediction], name=None, label=None, source=None) -> List[Prediction]:
    """Returns predictions for which all of the submitted criteria match. Otherwise return nothing.
    """
    matched = []
    for prediction in predictions:
        matched_name = name is None or prediction.score_name == name
        matched_label = label is None or prediction.score_label == label
        matched_source = source is None or prediction.score_source == source
        if matched_name and matched_label and matched_source:
            matched.append(prediction)
    return matched


HEADERS = {
    "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/114.0.0.0 Safari/537.36"
}


def get_spliceai_results(variant: GenomicVariant):
    SPLICE_AI_URL = "https://spliceailookup-api.broadinstitute.org/spliceai/?hg=37&distance=500&mask=1&variant=19-11233940-G-A&raw=19-11233940%20-%20G-A"

    if variant.assembly in ("GRCh37", "hg19"):
        assemblyint = 37
    elif variant.assembly in ("GRCh38", "hg38"):
        assemblyint = 38
    else:
        raise RuntimeError(f"Unsupported assembly string {variant.assembly}")

    params = {
        "hg": assemblyint,
        "distance": 500,
        "mask": 1,
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


class Metadome(Scraper):

    @classmethod
    def score_to_label(cls, score) -> str:
        """Taken from: https://github.com/cmbi/metadome/blob/84b281bf621a492a8cb75c04293c8631dbd3fede/metadome/presentation/web/static/js/dashboard/visualization.js#L1277

        """
        if score:
            if score <= 0.175:
                return "highly intolerant"
            elif score <= 0.525:
                return "intolerant"
            elif score <= 0.7:
                return "slightly intolerant"
            elif score <= 0.875:
                return "neutral"
            elif score <= 1.025:
                return "slightly tolerant"
            elif score <= 1.375:
                return "tolerant"
            elif score > 1.375:
                return "highly tolerant"
        return ""

    @classmethod
    def url(cls, variant: GenomicVariant) -> str:
        return ""

    @classmethod
    def query(cls, variant: GenomicVariant):
        effects = MyVariant.get_effect(variant)
        if not effects:
            return []

        def match_nm(nm_a, nm_b):
            nm_a_core, *_ = nm_a.split(".")
            nm_b_core, *_ = nm_b.split(".")
            return nm_a_core == nm_b_core

        if type(effects) is dict:
            effects= [effects]

        preds = []
        for effect in effects:
            gene_name = effect.get("genename")
            if effect and gene_name:
                protein_data = effect.get("prot", {})
                protein_pos = int(protein_data.get("protpos") or -1)
                transcripts = Metadome.get_transcripts(gene_name)
                for transcript in transcripts:
                    if transcript.get("has_protein_data"):
                        effect_enst = effect.get("feature_id")
                        transcript_enst = transcript.get("gencode_id")
                        transcript_nms = transcript.get("refseq_nm_numbers")
                        if match_nm(effect_enst, transcript_enst):
                            metadome_data = Metadome.get_transcript_predictions(transcript_enst)
                            metadome_effect = metadome_data["positional_annotation"][protein_pos]
                            metadome_effect["domains"]
                            score = metadome_effect["sw_dn_ds"]
                            label = cls.score_to_label(score)
                            pred = Prediction(variant=variant, score_name=f"{transcript_enst}({transcript_nms})", score_label=label, score_value=f"{score:.3f}", score_source=cls.name)
                            preds.append(pred)
        return preds

    @classmethod
    def get_transcripts(cls, gene_name: str) -> list:
        url = f"https://stuart.radboudumc.nl/metadome/api/get_transcripts/{gene_name}"
        resp = requests.get(url)
        resp.raise_for_status()
        data = resp.json()
        print(data)

        d = data.get("transcript_ids")
        if not d:
            return data.get("trancript_ids", [])
        return d

    @classmethod
    def get_transcript_predictions(cls, enst_id: str) -> dict:
        url = f"https://stuart.radboudumc.nl/metadome/api/result/{enst_id}/"
        resp = requests.get(url)
        resp.raise_for_status()
        data = resp.json()
        return data


def default_converter(results: dict):
    return ", ".join(str(r) for r in results)


def convert_list(results):
    if results:
        if type(results[0]) in (float, int):
            return results[0]
        elif len(results[0]) == 0:
            return ""
        maximum = max(results[0])
        return maximum
    return ""

def fmt_hgvs(results):
    if results:
        annotations = results[0]
        all_lines = []
        if type(annotations) is dict:
            annotations = [annotations]

        for annot in annotations:
            line = f"{annot['genename']}({annot['feature_id']}):{annot['hgvs_c']}"
            if "hgvs_p" in annot:
                line += f" {annot['hgvs_p']}"
            line += f" ({annot['effect']})"
            all_lines.append(line)
        return "|".join(all_lines)
    return ""


class Mapping(BaseModel):
    name: str
    keys: List[tuple]
    converter: Callable = default_converter


def fetch_mapped(entry, mapping):
    values = []
    for keys in mapping.keys:
        value = entry
        for key in keys:
            value = value.get(key, "")
            if not value:
                break
        values.append(value)
    fmt_value = mapping.converter(values)
    return fmt_value



class MyVariant(Scraper):
    BASE_URL = "https://myvariant.info/v1/variant"
    BASE_PARAMS = {
        "email": "max.zhao@charite.de",
        "size": 3,
    }

    MAPPINGS = [
        Mapping(name="CADD_PHRED", keys=[("cadd", "phred")]),
        Mapping(name="REVEL", keys=[("dbnsfp", "revel", "score")], converter=convert_list),
        Mapping(name="Polyphen2", keys=[("dbnsfp", "polyphen2", "hdiv", "rankscore")]),
        Mapping(name="Effect", keys=[("snpeff", "ann"),], converter=fmt_hgvs)
    ]

    @classmethod
    def url(cls, variant: GenomicVariant) -> str:
        return f"{cls.BASE_URL}/{str(variant)}"

    @classmethod
    def get_effect(cls, variant: GenomicVariant) -> dict:
        params = {**cls.BASE_PARAMS, **{
        }}
        variant_id = str(variant)
        resp = requests.post(cls.BASE_URL, params=params, json={"ids": [variant_id]})
        resp.raise_for_status()
        annotation_data = resp.json()

        if len(annotation_data) < 1 or annotation_data[0].get("notfound", False):
            return {}
        annotation_entry = annotation_data[0]

        effects, = fetch_mapped(annotation_entry, Mapping(name="Effect", keys=[("cadd", "gene")], converter=lambda d: d))
        return effects


    @classmethod
    def query(cls, variant: GenomicVariant):
        params = {**cls.BASE_PARAMS, **{
        }}
        variant_id = str(variant)
        resp = requests.post(cls.BASE_URL, params=params, json={"ids": [variant_id]})
        resp.raise_for_status()
        annotation_data = resp.json()

        if len(annotation_data) < 1 or annotation_data[0].get("notfound", False):
            return []

        annotation_entry = annotation_data[0]
        predictions = []
        for prediction_mapping in cls.MAPPINGS:
            result = fetch_mapped(annotation_entry, prediction_mapping)
            pred = Prediction(variant=variant, score_value=result, score_name=prediction_mapping.name, score_source=cls.name)
            if result:
                predictions.append(pred)

        return predictions



VARIANT_SOURCES = [
    CADD, MutationTaster2021, Franklin, GNOMAD, MyVariant, Metadome
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

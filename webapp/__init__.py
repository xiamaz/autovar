import re
from typing import Optional
from autovar import queryVariant, queryUrls, GenomicVariant
from flask import Flask, render_template, request

app = Flask(__name__)


def parse(variant_string: str) -> Optional[GenomicVariant]:
    variant_string = variant_string.replace(",", "")
    if m := re.search(r"(\d+)[ :\-_)]+(\d+)[ \-_:]*([ATCGatcg])[ >\-_:]([ATCGatcg])", variant_string):
        chrom = m.group(1)
        pos = int(m.group(2))
        ref = m.group(3)
        alt = m.group(4)
        return GenomicVariant(assembly="GRCh37", chr=chrom, position=pos, ref=ref, alt=alt)
    return None


def first_entry(name="", source=""):
    def wraps(predictions):
        for pred in predictions:
            if (not name or pred.score_name == name) and (not source or pred.score_source == source):
                return pred.score_value
        return "NOT FOUND"
    return wraps


def highest_value_label(name="", source=""):
    def wraps(predictions):
        candidates = []
        cand_labels = []
        for pred in predictions:
            if (not name or pred.score_name == name) and (not source or pred.score_source == source):
                score_value = pred.score_value
                if "|" in score_value:
                    score_path, score_ben = score_value.split("|")
                    print(score_path, score_ben, score_value)
                    score_value = float(score_path) / (float(score_ben) + float(score_path))
                    candidates.append(score_value)
                    cand_labels.append(pred.score_label)
                elif score_value:
                    score_value = float(score_value)
                    candidates.append(score_value)
                    cand_labels.append(pred.score_label)
        if candidates:
            return cand_labels[candidates.index(max(candidates))]
        return "NOT FOUND"
    return wraps


PRETTY_KEYS = {
    "MutationTaster": highest_value_label(source="MutationTaster2021"),
    "Polyphen": first_entry("Polyphen2"),
    "SIFT": first_entry("sift"),
    "CADD": first_entry("CADD_PHRED"),
    "REVEL": first_entry("REVEL"),
    "SpliceAI": first_entry("splice_ai_genome"),
    "GNOMAD MAF": first_entry("MAF"),
    "GNOMAD HOM": first_entry("HOM"),
}

def format_entry(name, mapper, predictions):
    all_preds = [p for pp in predictions.values() for p in pp]
    entry_line = {
        "name": name,
        "value": mapper(all_preds),
    }
    return entry_line

@app.route("/", methods=["POST", "GET"])
def main():
    predictions = {}
    sel_preds = []
    effects = []
    urls = {}
    raw_variant = ""
    error = ""
    if request.method == "POST":
        raw_variant = request.form["variant"]

        if variant := parse(raw_variant):
            predictions = queryVariant(variant)
            urls = queryUrls(variant)
        else:
            error = f"Failed to parse {raw_variant} into a genomic variant"

    if predictions:
        for name, mapper in PRETTY_KEYS.items():
            fmt_key = format_entry(name, mapper, predictions)
            sel_preds.append(fmt_key)

        for preds in predictions["MyVariant"]:
            if preds.score_name == "Effect":
                effects.append(preds.score_value)

    return render_template("main.html", predictions=predictions, urls=urls, query=raw_variant, error=error, selected_predictions=sel_preds, effects=effects)

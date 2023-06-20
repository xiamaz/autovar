import re
from typing import Optional
from autovar import queryVariant, queryUrls, GenomicVariant
from flask import Flask, render_template, request

app = Flask(__name__)


def parse(variant_string: str) -> Optional[GenomicVariant]:
    if m := re.search(r"(\d+)[ :\-_)]+(\d+)[ \-_:]*([ATCGatcg])[ >\-_:]([ATCGatcg])", variant_string):
        chrom = m.group(1)
        pos = int(m.group(2))
        ref = m.group(3)
        alt = m.group(4)
        return GenomicVariant(assembly="GRCh37", chr=chrom, position=pos, ref=ref, alt=alt)
    return None

@app.route("/", methods=["POST", "GET"])
def main():
    predictions = {}
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

    return render_template("main.html", predictions=predictions, urls=urls, query=raw_variant, error=error)

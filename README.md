# Autovar

Automatically query genomic variants.

## Setup

Just have docker and docker compose installed, run as a service using the following command:

```
docker compose up -d
```

Afterwards the service will be exposed on Port 5001. This can be configured by changing the Port in the `Dockerfile` and the `docker-compose.yml`.


## Usage

Just query a genomic position in hg19 or grch37 syntax including ref and alt allele. Currently only single nucleotide exchanges (both coding and non-coding) are supported.

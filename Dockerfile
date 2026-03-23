FROM python:3.11-slim AS builder

RUN apt-get update && \
    apt-get install -y --no-install-recommends build-essential && \
    rm -rf /var/lib/apt/lists/*

ARG VERSION
RUN pip install --no-cache-dir afquery==${VERSION}


FROM python:3.11-slim

COPY --from=builder /usr/local/lib/python3.11/site-packages /usr/local/lib/python3.11/site-packages
COPY --from=builder /usr/local/bin/afquery /usr/local/bin/afquery

ENTRYPOINT ["afquery"]

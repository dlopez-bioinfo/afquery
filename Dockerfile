FROM python:3.11-slim

ARG VERSION
RUN pip install --no-cache-dir afquery==${VERSION}

ENTRYPOINT ["afquery"]

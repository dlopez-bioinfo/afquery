# Contributing

## Documentation deployment

The published site at <https://dlopez-bioinfo.github.io/afquery/> is built with [MkDocs](https://www.mkdocs.org/) (Material theme) and versioned with [mike](https://github.com/jimporter/mike). Source content lives under `docs/` and is configured by `mkdocs.yml`.

### CI workflows

| Workflow | Trigger | What it does |
|---|---|---|
| `.github/workflows/docs.yml` | push to `master` (paths: `docs/**`, `mkdocs.yml`, `src/**`) and manual dispatch | `mike deploy --push --update-aliases dev` — publishes the working master under the `dev` alias |
| `.github/workflows/release.yml` (job `docs`) | push of a `v*` tag (non-`rc`) | `mike deploy --push --update-aliases <version> latest` followed by `mike set-default --push latest` — publishes the tagged version, points `latest` at it, and makes `latest` the site root |

### Cutting a release

Tag the commit with a [PEP 440](https://peps.python.org/pep-0440/)-compatible version prefixed with `v` (e.g. `v0.3.0`) and push the tag. The `release.yml` workflow handles PyPI, the GitHub release, the Docker image, and the versioned docs deploy. Pre-release tags (anything containing `rc`) skip Docker and docs publishing.

### Local preview

```bash
pip install -e ".[docs]"

# Plain build — no versioning, fastest iteration on content
mkdocs serve

# Versioned layout — only after running mike deploy locally at least once
mike deploy 0.0.0-test dev      # writes to local gh-pages branch (no push)
mike serve                      # serves the gh-pages branch with version selector
```

To discard local mike state: `git branch -D gh-pages`.

### Migrating from a flat `mkdocs gh-deploy` (one-time)

When `gh-pages` still contains content from the previous flat deploy, the first `mike deploy` will leave the old root-level files in place. To wipe `gh-pages` clean before deploying, manually trigger the *Deploy Documentation* workflow with `bootstrap: true`. This runs `mike delete --all --push --allow-empty` before deploying, leaving only the versioned layout (`versions.json` + version subdirectories + redirector at root).

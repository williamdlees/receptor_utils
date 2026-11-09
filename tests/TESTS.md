## Running tests

```
pip install -e ".[test]"
pytest                      # runs the suite (network tests are skipped by default)
pytest --cov                # with coverage
pytest -m network           # optional tests that contact live services
```


The suite runs on every push and pull request via GitHub Actions (`.github/workflows/tests.yml`).

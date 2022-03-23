# Import ontologies to NXOntology

## Development

```shell
# Install the environment
poetry install --no-root

# Update the lock file
poetry update

# Run tests
pytest

# Set up the git pre-commit hooks.
# `git commit` will now trigger automatic checks including linting.
pre-commit install

# Run all pre-commit checks (CI will also run this).
pre-commit run --all
```

## License

This repository is released under an Apache License 2.0 License (see [LICENSE.md](LICENSE.md)).
Please see the source licenses of ontologies imported into NXOntology objects.

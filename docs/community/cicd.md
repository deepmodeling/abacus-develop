# Lists of continuous integration (CI) pipelines

The directory `./github/workflows` contains the continuous integration (CI) pipelines for the project. The pipelines are written in YAML format and are executed by GitHub Actions. Check the [Actions page](https://github.com/deepmodeling/abacus-develop/actions) of the repo for the status of the pipelines.

> The tests without mentioning are not executed.

## On Pull Request (PR)

The following CI pipelines are triggered on pull request (PR) creation or update (the user pushes a new commit to the incoming branch):

- Integration test and unit tests (`test.yml`): This pipeline builds ABACUS with all available features, runs integration tests, and runs unit tests.

- Building tests with CMake (`build_test_cmake.yml`): This pipeline builds ABACUS with each feature separately, ensuring: (i) there are no conflicts between features, (ii) the compilation is successful whether the feature is enabled, and (iii) it works well on multiple platforms, i.e. with GNU+OpenBLAS toolchain and Intel+MKL toolchain.

- [Rerender the docs site](https://readthedocs.org/projects/abacus-rtd/builds/): This pipeline rerenders the documentation site on Read the Docs. It is automatically triggered when the documentation is updated.

- Testing GPU features (`cuda.yml`): This pipeline builds ABACUS with GPU support and runs several tests on the GPU. **Currently disabled for the lack of GPU resource.**

After the PR merges into the main branch, the following pipelines are triggered:

- Building Docker images(`devcontainer.yml`): This pipeline builds the Docker images with latest codes and executables. The images are tagged as `abacus-gnu:latest`, `abacus-intel:latest`, and `abacus-cuda:latest`, and then pushed to the GitHub Container Registry (ghcr.io/deepmodeling) and AliCloud Container mirror(registry.dp.tech/deepmodeling). For example: `docker pull ghcr.io/deepmodeling/abacus-intel:latest`.

## On Routine

- Dynamic analysis (`dynamic_analysis.yml`): This pipeline runs integration tests with [AddressSanitizer](https://github.com/google/sanitizers/wiki/AddressSanitizer) to detect memory errors. The pipeline is scheduled to run **every Sunday**. The results are published on [GitHub Pages](https://deepmodeling.github.io/abacus-develop/).

## On Release

- Coverage test (`coverage.yml`): This pipeline builds ABACUS with all available features, runs integration tests, and runs unit tests. It also measures the code coverage of the tests. The results are published at [codecov.io](https://app.codecov.io/gh/deepmodeling/abacus-develop).
- Building tagged Docker images (`image.yml`): The built image is tagged in the pattern of `abacus:3.5.0`, and pushed to the GitHub Container Registry (ghcr.io/deepmodeling) and AliCloud Container mirror(registry.dp.tech/deepmodeling). Use `abacus:latest` to fetch the latest image. For example: `docker pull ghcr.io/deepmodeling/abacus:latest`.

# Lists of continuous integration (CI) actions

The directory `.github/workflows` contains the continuous integration (CI) actions for the project. The actions are written in YAML format and are executed by GitHub Actions. Check the [Actions page](https://github.com/deepmodeling/abacus-develop/actions) of the repo for the status of the actions.

## On Pull Request (PR)

The following CI actions are triggered on pull request (PR) creation or update (the user pushes a new commit to the incoming branch):

- Integration test and unit tests (`test.yml`): This action builds ABACUS with all available features, runs integration tests, and runs unit tests. It also performs [static analysis](../CONTRIBUTING.md#code-formatting-style) with <pre-commit.ci>.
- Building tests with CMake (`build_test_cmake.yml`) and Makefile (`build_test_makefile.yml`): This action builds ABACUS with each feature separately, ensuring: (i) there are no conflicts between features, (ii) the compilation is successful whether the feature is enabled, and (iii) it works well on multiple platforms, i.e. with GNU+OpenBLAS toolchain and Intel+MKL toolchain.
- [Rerender the docs site](https://readthedocs.org/projects/abacus-rtd/builds/): This action rerenders the documentation site on Read the Docs. It is automatically triggered when the documentation is updated.
- Testing GPU features (`cuda.yml`): This action builds ABACUS with GPU support and runs several tests on the GPU.

> Some tests are executed on self-hosted runners, which are maintained by the deepmodeling community. To save resources, the actions triggered by new commits may cancel the previous actions in the same PR.

## On PR Merge

After the PR merges into the main branch, the following actions are triggered:

- Building Docker images (`devcontainer.yml`): This action builds the Docker images with latest codes and executables. The images are tagged as `abacus-gnu:latest`, `abacus-intel:latest`, and `abacus-cuda:latest`, and then pushed to the GitHub Container Registry (`ghcr.io/deepmodeling`) and AliCloud Container Registry (`registry.dp.tech/deepmodeling`, recommended for one having issue connecting to GitHub). For example: `docker pull ghcr.io/deepmodeling/abacus-intel:latest`.
- Generate doxygen site (`doxygen.yml`): This action generates the Doxygen site for the project. The site is published on [GitHub Pages](https://deepmodeling.github.io/abacus-develop/).
- Mirror the repo to Gitee (`mirror.yml`): This action mirrors the repo to [Gitee](https://gitee.com/deepmodeling/abacus-develop).

## On Routine

- Dynamic analysis (`dynamic_analysis.yml`): This action runs integration tests with [AddressSanitizer](https://github.com/google/sanitizers/wiki/AddressSanitizer) to detect memory errors. The action is scheduled to run **every Sunday**. The results are published to [the dashboard branch](https://github.com/deepmodeling/abacus-develop/blob/dashboard/README.md).

## On Release

- Coverage test (`coverage.yml`): This action builds ABACUS with all available features, runs integration tests, and runs unit tests. It also measures the code coverage of the tests. The results are published at [codecov.io](https://app.codecov.io/gh/deepmodeling/abacus-develop).
- Building tagged Docker images (`devcontainer.yml`): Same as that action above; in addition the built image is tagged in the pattern of `abacus-intel:3.6.0`. For example: `docker pull ghcr.io/deepmodeling/abacus-intel:3.6.0`.

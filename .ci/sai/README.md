# SAI GPU CI manual

## Purpose and scope

This workflow rebuilds ABACUS at one exact commit on the SAI GPU cluster and
runs the matrix in `gpu-matrix.ini` through Slurm. The current matrix has 49
cases. It reports build, resource, and case states and saves the raw logs.

This is validation, not a benchmark or a general remote-shell service. The
client requires a committed checkout: it bundles the checked-out `HEAD` and
checks that it equals `--source-sha`. Uncommitted files are not sent.

**Warning:** approved candidate code is compiled and executed as the configured
SAI SSH user. Use an account and project root with the intended permissions.

## GitHub setup

Create these GitHub Environments for the `SAI GPU validation` workflow:

- `sai-ssh-scheduled`: no required reviewers. The daily schedule uses this
  environment.
- `sai-ssh-manual`: require the reviewers who should approve a remote run.
  Manual dispatches and the pull-request comment trigger use this environment.

In each environment, configure the SSH secret and connection variables:

- Secret: `SAI_SSH_PRIVATE_KEY` (the private key used by the runner).
- Variables: `SAI_PROJECT_ROOT`, `SAI_SSH_HOST`, `SAI_SSH_PORT`, and
  `SAI_SSH_USER`.

`SAI_PROJECT_ROOT` should be an absolute path below the SAI user's home, for
example `/home/<group>/<user>/abacus_sai_gpu_ci`. A manual dispatch may provide
the optional `project_root` input to override this variable for that run.

## Triggers

- **Manual:** open Actions, choose `SAI GPU validation`, and select
  **Run workflow**. `source_sha` is required and must be a 40-character commit
  SHA. The workflow must be dispatched from the repository default branch.
  `project_root` is optional. The `sai-ssh-manual` approval is requested before
  the build job starts.
- **Daily:** the schedule is `30 20 * * *` (20:30 UTC). It tests the default
  branch of `deepmodeling/abacus-develop` and uses `sai-ssh-scheduled`.
- **Pull request:** comment exactly `/abacus-ci sai-gpu` on an open pull
  request. The commenter needs Triage, Write, Maintain, or Admin permission.
  The pull request head commit is tested after `sai-ssh-manual` approval.

Other issue comments are ignored.

## Run locally

Run from the repository checkout whose `HEAD` is the commit to test. The
following uses the committed `HEAD` SHA and the same required options as CI:

```bash
SOURCE_SHA="$(git rev-parse HEAD)"
RUN_ID="$(date +%s)"
python3 .ci/sai/sai.py local \
  --ssh-config "$HOME/.ssh/sai-config" \
  --target sai-ci \
  --project-root /home/<group>/<user>/abacus_sai_gpu_ci \
  --source-repository "$PWD" \
  --source-sha "$SOURCE_SHA" \
  --namespace manual \
  --run-id "$RUN_ID" \
  --run-attempt 1 \
  --artifacts ./sai-artifacts
```

The SSH config must define the `sai-ci` target (host, port, user, key, and
known-host checking). `--source-repository` must be the checkout used to
compute `SOURCE_SHA`. The command waits for Slurm completion and exits with a
non-zero status when a case or infrastructure result is not fully successful.

For all available subcommands and options, run:

```bash
python3 .ci/sai/sai.py --help
python3 .ci/sai/sai.py local --help
```

## Matrix configuration

`gpu-matrix.ini` is validated before jobs are submitted.

- `[cluster]`: Slurm `partition`, absolute `mapping_root` for the MPI mapping
  script, `disable_nccl_ib` (`true` or `false`), and `poll_seconds` (1-300).
- `[build]`: build-job `qos`, `nodes`, `tasks_per_node`, `gpus_per_node`, and
  `time_seconds`.
- `[resource.NAME]`: the same allocation fields plus `parallelism`, the
  maximum number of array tasks running at once. Each resource must have a
  case. There is one rank per GPU and no resource may exceed 16 GPUs.
- `[case.NNN]`: contiguous, zero-padded sections with `suite`, `name`,
  `resource`, and `runner` (`autotest` or `cusolvermp`).

Resource component labels are generated, not configured separately. A
single-node resource is shown as `N GPU` or `N GPUs`; a multi-node resource is
shown as `N nodes / M GPUs`. Thus `gpu1`, `gpu2`, and `gpu4` display `1 GPU`,
`2 GPUs`, and `4 GPUs`; `gpu8x2` displays `2 nodes / 16 GPUs`.
`case.049` is `15_rtTDDFT_GPU/19_NO_Si48_CUSOLVERMP_TDDFT_GPU`; it uses
`gpu8x2` and the `cusolvermp` runner.

## Results and retention

On SAI, a run is created below:

```
<SAI_PROJECT_ROOT>/runs/<namespace>/<run-id>-<attempt>/
```

Its `results/` directory contains `result.json`, `summary.md`, build and case
logs, Slurm output, module/tool records, and status files. The coordinator and
working data are alongside it while the run is active. After results are
collected, CI archives `results/` and `jobs/` as:

```
<SAI_PROJECT_ROOT>/archives/<namespace>/<run-id>-<attempt>.tar.gz
```

The client removes archived files older than 72 hours when preparing a later
run, and removes the active run after archiving. On the GitHub runner,
`ARTIFACT_ROOT` is `${runner.temp}/sai-artifacts`; it contains the collected
`results/`, `run.json`, and `client.log`. CI uploads that directory as
`sai-gpu-<run-id>-<attempt>` and retains it for 30 days. A pull-request comment
links to the Actions run and the uploaded raw files.

## Troubleshooting

**SSH fails.** Check `SAI_SSH_HOST`, `SAI_SSH_PORT`, and `SAI_SSH_USER`, that
the key matches the account, and that the target is reachable. CI uses the
committed `.ci/sai/known_hosts` with strict host-key checking. Test the same
target with the SSH config before retrying; do not disable host-key checking.

**A module cannot be loaded.** `modules.sh` sources Lmod, purges modules, and
loads `cmake/3.31.6` and the configured ABACUS dependency module. Ask the SAI
administrator to provide or update those modules. Modules provide the
compiler, CUDA, MPI, and library dependencies; do not add library paths to CI
(`LD_LIBRARY_PATH`, `CPATH`, or `CMAKE_PREFIX_PATH`) or hard-code site paths.

**CMake or linking fails.** Inspect `results/configure.log`, `build.log`,
`install.log`, `CMakeCache.txt`, `tools.txt`, and `ldd.txt`. The build uses Unix
Makefiles, CUDA architecture 70, CUDA MPI, cuSOLVERMP, cuBLASMP, and
NCCL parallel-device options. A missing runtime library causes the `ldd` check
to fail; fix the module environment rather than adding a CI path.

**A job stays pending or times out.** Check the selected partition and QoS,
GPU availability, node and task limits, and the `time_seconds` value for that
resource. Slurm output is in `results/`; an allocation or queue delay is an
infrastructure issue, not a case failure.

**MPI/PMIx initialization fails.** Autotest retries once after a recognized
`PMIX_ERR_FILE_OPEN_FAILURE` or `PMIX_ERR_OUT_OF_RESOURCE` initialization
failure. If it persists, inspect the attempt logs and the loaded MPI module,
Slurm allocation, and mapping file. The `cusolvermp` runner does not use this
retry.

# Remote GPU validation

## Purpose and scope

This workflow rebuilds ABACUS at one exact commit on a remote GPU cluster and
runs the matrix in `config.ini` through Slurm. The configuration contains the
remote connection, Slurm resources, and 49 test cases. The workflow reports
build, resource, and case states and saves the raw logs.

This is validation, not a benchmark or a general remote-shell service. The
client requires a committed checkout: it bundles the checked-out `HEAD` and
checks that it equals `--source-sha`. Uncommitted files are not sent.

**Warning:** approved candidate code is compiled and executed as the configured
remote SSH user. Use an account and project root with the intended permissions.

## Reference deployment

The maintained deployment currently runs on the SAI cluster. Its queue,
module, MPI mapping, and SSH host-key settings are kept in the configuration
files in this directory and can be adapted for another Slurm site.

## GitHub setup

Create these GitHub Environments for the `GPU validation` workflow:

- `gpu-ci-scheduled`: no required reviewers. The daily schedule uses this
  environment.
- `gpu-ci-manual`: require the reviewers who should approve a remote run.
  Manual dispatches and the pull-request comment trigger use this environment.

In each environment, configure one secret:

- Secret: `REMOTE_SSH_PRIVATE_KEY` (the private key used by the runner).

The two environments have different approval rules, so the same key must be
added to each one. Host, port, user, and the default project root come from the
trusted `[remote]` section in `config.ini`. A manual dispatch may provide the
optional `project_root` input to override the configured root for that run.

## Triggers

- **Manual:** open Actions, choose `GPU validation`, and select
  **Run workflow**. `source_sha` is required and must be a 40-character commit
  SHA. The workflow must be dispatched from the repository default branch.
  `project_root` is optional. The `gpu-ci-manual` approval is requested before
  the build job starts.
- **Daily:** the schedule is `30 20 * * *` (20:30 UTC). It tests the default
  branch of `deepmodeling/abacus-develop` and uses `gpu-ci-scheduled`.
- **Pull request:** comment exactly `/abacus-ci gpu` on an open pull
  request. The commenter needs Triage, Write, Maintain, or Admin permission.
  The pull request head commit is tested after `gpu-ci-manual` approval.

Other issue comments are ignored.

## Run locally

Run from the repository checkout whose committed `HEAD` is the revision to
test. If the SSH config defines the default `gpu-ci` target, no options are
required:

```bash
python3 .ci/slurm/runner.py run
```

By default, the command uses `~/.ssh/config`, the `gpu-ci` host alias, the
current Git checkout and its `HEAD`, and `~/abacus_gpu_ci` below the remote
user's home. It writes downloaded results to `./gpu-ci-artifacts`. Override any
value shown by `run --help`; for example, use `--target my-cluster` for a
different local host alias. The command waits for Slurm completion and exits
with a non-zero status when a case or infrastructure result is not fully
successful.

Source is sent as one compressed Git bundle split into eight parallel rsync
transfers. The remote side checks the merged SHA-256 and the Git bundle before
updating its cache. Later runs send only commits after the cached SHA; a run at
the same SHA sends no source data.

For the public commands and all `run` options, use:

```bash
python3 .ci/slurm/runner.py --help
python3 .ci/slurm/runner.py run --help
```

## Configuration

`config.ini` is validated before jobs are submitted.

- `[remote]`: SSH `host`, `port`, `user`, and `project_root`. The project root
  must be an absolute path or use `~/` relative to the remote user's home.
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

On the remote cluster, a run is created below:

```
<project-root>/runs/<namespace>/<run-id>-<attempt>/
```

Its `results/` directory contains `result.json`, `summary.md`, build and case
logs, Slurm output, module/tool records, and status files. The coordinator and
working data are alongside it while the run is active. After results are
collected, CI archives `results/` and `jobs/` as:

```
<project-root>/archives/<namespace>/<run-id>-<attempt>.tar.gz
```

The client removes archived files older than 72 hours when preparing a later
run, and removes the active run after archiving. On the GitHub runner,
`ARTIFACT_ROOT` is `${runner.temp}/gpu-ci-artifacts`; it contains the collected
`results/`, `run.json`, and `client.log`. CI uploads that directory as
`gpu-validation-<run-id>-<attempt>` and retains it for 30 days. A pull-request comment
links to the Actions run and the uploaded raw files.

## Troubleshooting

**SSH fails.** Check the `[remote]` values in `config.ini`, that the key matches
the configured account, and that the target is reachable. CI uses the committed
`.ci/slurm/known_hosts` with strict host-key checking. Test the same target with
the SSH config before retrying; do not disable host-key checking.

**A module cannot be loaded.** `modules.sh` sources Lmod, purges modules, and
loads `cmake/3.31.6` and the configured ABACUS dependency module. Ask the site
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

**MPI/PMIx initialization fails.** Both runners retry once after a recognized
MPI startup failure. If it persists, inspect the attempt logs and the loaded
MPI module, Slurm allocation, and mapping file.

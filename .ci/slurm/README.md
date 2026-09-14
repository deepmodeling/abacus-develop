# Remote GPU validation

This workflow rebuilds one committed ABACUS revision on a remote GPU cluster, runs the test matrix in `config.ini` through Slurm, and reports each build and test group separately. It is a functional test, not a benchmark.

The maintained setup runs at the [Open Source Supercomputing Center of SAI](https://www.open-sai.com/). The same client can be configured for another Slurm cluster.

**Trust boundary:** the selected commit is compiled and executed as the remote SSH user. Approve only code that may run with that account's permissions.

## Set up GitHub

A repository administrator performs these steps once. Forks start disabled because GitHub does not copy variables or secrets from the parent repository.

1. Open **Settings > Secrets and variables > Actions > Variables**, choose **New repository variable**, and set:

   - Name: `GPU_VALIDATION_ENABLED`
   - Value: `true`

   The value is the lowercase text `true`. If this variable is absent or has another value, the workflow skips all remote work.

2. Open **Settings > Environments** and create:

   - `gpu-ci-scheduled`, with no required reviewers, for daily tests.
   - `gpu-ci-manual`, with the maintainers who may approve PR and manual tests listed as required reviewers.

3. Open each environment, choose **Environment secrets > Add environment secret**, and add:

   - Name: `REMOTE_SSH_PRIVATE_KEY`
   - Value: the complete private key, including its `BEGIN` and `END` lines.

   Install the matching public key in `authorized_keys` for the remote account named in `config.ini`. Add the private key to both environments because they have different approval rules. Do not create a repository-level SSH secret; the workflow reads this environment secret only after entering the selected environment.

Host, port, user, and the normal remote project directory are read from the trusted `[remote]` section of `config.ini`.

## Run validation

All GitHub methods below require `GPU_VALIDATION_ENABLED=true`. The workflow must already be present on the repository's default branch.

### Test a pull request

On an open pull request, add this exact comment:

```text
/abacus-ci gpu
```

The comment cannot contain other text. The author of the comment needs Triage, Write, Maintain, or Admin permission. The bot immediately posts a link to the queued Actions run. After a reviewer approves the `gpu-ci-manual` environment, the workflow tests the PR head commit and updates that same bot comment with the result and raw-file link.

### Run the daily test

No manual action is needed. The workflow is scheduled every day at 20:30 UTC (`30 20 * * *`). It tests the current default branch of `deepmodeling/abacus-develop` and uses `gpu-ci-scheduled`, so it does not wait for approval.

### Start a run from Actions

1. Open **Actions > GPU validation > Run workflow**.
2. Select the repository default branch under **Use workflow from**.
3. Enter the full, lowercase 40-character commit SHA in `source_sha`.
4. Leave `project_root` empty to use `config.ini`, or enter another permitted remote directory.
5. Start the run and approve the `gpu-ci-manual` environment when prompted.

The commit must exist in the repository where the workflow is running. For an external contributor's pull request, use the PR comment command instead.

### Run from a local checkout

Create an SSH host entry. The default alias is `gpu-ci`; use the host, port, user, and key for your account:

```sshconfig
Host gpu-ci
    HostName <host>
    Port <port>
    User <user>
    IdentityFile ~/.ssh/<private-key>
```

Then run this command from a committed ABACUS checkout:

```bash
python3 .ci/slurm/runner.py run
```

By default, the command uses the checkout's committed `HEAD`, `~/.ssh/config`, the `gpu-ci` alias, and the remote directory from `config.ini`. Uncommitted candidate-source changes are not sent. The local command does use the current `.ci/slurm` control files, including local changes to its scripts and templates. The command waits for Slurm, prints live build and test progress, downloads the results, and exits nonzero if validation fails.

Local results go to `/tmp/abacus_gpu_ci_<uid>/<namespace>/<run_id>_<attempt>/`. Use `--artifacts` for a permanent local directory or `--target my-cluster` for another SSH alias. All available options and defaults are shown by:

```bash
python3 .ci/slurm/runner.py --help
python3 .ci/slurm/runner.py run --help
```

## Configuration

`config.ini` is validated before jobs are submitted.

- `[site]`: the resource acknowledgement, site name, and public URL shown at the end of result reports. Change these values for another cluster.
- `[remote]`: SSH `host`, `port`, `user`, `project_root`, and comma-separated `allowed_project_roots`. The project root may use `~/` or an absolute path, but its remote resolved path must be below one of the allowed roots. Prefixes are also resolved remotely, so aliases such as `/home` pointing to `/org` are accepted.
- `[cluster]`: Slurm `partition`, absolute `mapping_root` for the MPI mapping script, `disable_nccl_ib` (`true` or `false`), and `poll_seconds` (1-300).
- `[build]`: build-job `qos`, `nodes`, `tasks_per_node`, `gpus_per_node`, and `time_seconds`.
- `[resource.NAME]`: the same allocation fields plus `parallelism`, the maximum number of array tasks running at once. Each resource must have a case. There is one rank per GPU and no resource may exceed 16 GPUs.
- `[case.NNN]`: contiguous, zero-padded sections with `suite`, `name`, `resource`, and `runner` (`autotest` or `cusolvermp`).

Resource labels are generated, not configured separately. A single-node resource is shown as `N GPU` or `N GPUs`; a multi-node resource is shown as `N nodes / M GPUs`. Thus `gpu1`, `gpu2`, and `gpu4` display `1 GPU`, `2 GPUs`, and `4 GPUs`; `gpu4x2` displays `2 nodes / 8 GPUs`. Test results are reported by their `tests/` folder even though Slurm arrays remain grouped by resource. `case.049` is `15_rtTDDFT_GPU/19_NO_Si48_CUSOLVERMP_TDDFT_GPU`; it uses `gpu4x2` and the `cusolvermp` runner.

## Results and retention

On the remote cluster, a run is created below:

```
<project-root>/runs/<namespace>/<run-id>-<attempt>/
```

Its `results/` directory contains `result.json`, `summary.md`, build and case logs, Slurm output, module/tool records, and status files. The coordinator and working data are alongside it while the run is active. After results are collected, the client archives `results/` and `jobs/` as:

```
<project-root>/archives/<namespace>/<run-id>-<attempt>.tar.gz
```

The client removes archived files older than 72 hours when preparing a later run, and removes the active run after archiving. On the GitHub runner, `ARTIFACT_ROOT` is `${runner.temp}/gpu-ci-artifacts`; it contains the collected `results/`, `run.json`, and `client.log`. CI uploads that directory as `gpu-validation-<run-id>-<attempt>` and retains it for 30 days. A pull-request comment links to the Actions run and the uploaded raw files. If the client stops before completion, the remote run is left in place so that its detached coordinator and Slurm jobs are not interrupted.

Source is sent as a compressed Git bundle. The remote Git cache keeps the three most recent PR or manual revisions, the latest two daily dates, the first daily revision of every UTC month, and one weekly revision for the current month. Weekly revisions from earlier months are removed. Concurrent runs reserve the cache revisions they use, so another run cannot remove their transfer base.

## Troubleshooting

**SSH fails.** Check the `[remote]` values in `config.ini`, that the key matches the configured account, and that the target is reachable. CI uses the committed `.ci/slurm/known_hosts` with strict host-key checking. Test the same target with the SSH config before retrying; do not disable host-key checking.

**A module cannot be loaded.** `modules.sh` sources Lmod, purges modules, and loads `cmake/3.31.6` and the configured ABACUS dependency module. Ask the site administrator to provide or update those modules. Modules provide the compiler, CUDA, MPI, and library dependencies; do not add library paths to CI (`LD_LIBRARY_PATH`, `CPATH`, or `CMAKE_PREFIX_PATH`) or hard-code site paths.

**CMake or linking fails.** Inspect `results/configure.log`, `build.log`, `install.log`, `CMakeCache.txt`, `tools.txt`, and `ldd.txt`. The build uses Unix Makefiles, CUDA architecture 70, CUDA MPI, cuSOLVERMP, cuBLASMP, and NCCL parallel-device options. A missing runtime library causes the `ldd` check to fail; fix the module environment rather than adding a CI path.

**A job stays pending or times out.** Check the selected partition and QoS, GPU availability, node and task limits, and the `time_seconds` value for that resource. Slurm output is in `results/`; an allocation or queue delay is an infrastructure issue, not a case failure.

**MPI/PMIx initialization fails.** Both runners retry once after a recognized MPI startup failure. If it persists, inspect the attempt logs and the loaded MPI module, Slurm allocation, and mapping file.

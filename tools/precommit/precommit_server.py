#!/usr/bin/env python3

# Based on CP2K's tools/precommit/precommit_server.py, adapted for ABACUS.

import logging
import os
from os import path
import subprocess
from subprocess import PIPE, STDOUT
import tempfile
from time import time

from flask import Flask, request

app = Flask(__name__)
app.config["MAX_CONTENT_LENGTH"] = 1024 * 1024  # 1MiB
app.logger.setLevel(logging.INFO)
app.logger.info("ABACUS Precommit Server is up and running :-)")


# ======================================================================================
@app.route("/")
def hello():
    return "abacus precommit server revision: " + os.environ.get("REVISION", "unknown")


# ======================================================================================
@app.route("/black", methods=["POST"])
def black():
    return run_tool(["black"])


# ======================================================================================
@app.route("/shfmt", methods=["POST"])
def shfmt():
    return run_tool(["shfmt", "-i=2", "-ci", "-sr", "-w"])


# ======================================================================================
@app.route("/shellcheck", methods=["POST"])
def shellcheck():
    return run_tool(["shellcheck"], timeout=30)


# ======================================================================================
@app.route("/mdformat", methods=["POST"])
@app.route("/markdownlint", methods=["POST"])
def mdformat():
    return run_tool(["mdformat", "--wrap=100"], timeout=30)


# ======================================================================================
@app.route("/clangformat", methods=["POST"])
def clangformat():
    return run_tool(["clang_format_wrapper.sh"])


# ======================================================================================
@app.route("/cmakeformat", methods=["POST"])
def cmakeformat():
    return run_tool(["cmake-format", "-i"], timeout=30)


# ======================================================================================
def run_tool(cmd, timeout=3):
    assert len(request.files) == 1
    orig_fn = list(request.files.keys())[0]
    data_before = request.files[orig_fn].read()
    data_kb = len(data_before) / 1024.0
    fn = path.basename(orig_fn)

    with tempfile.TemporaryDirectory() as workdir:
        abs_fn = path.join(workdir, fn)
        with open(abs_fn, "wb") as f:
            f.write(data_before)

        t1 = time()
        try:
            p = subprocess.run(
                cmd + [fn], cwd=workdir, timeout=timeout, stdout=PIPE, stderr=STDOUT
            )
        except subprocess.TimeoutExpired:
            app.logger.info(f"Timeout of {cmd[0]} on {data_kb:.1f}KB after {timeout}s.")
            return f"Timeout while running {cmd[0]} - please try again.", 408
        t2 = time()
        app.logger.info(f"Ran {cmd[0]} on {data_kb:.1f}KB in {t2-t1:.1f}s.")

        if p.returncode != 0:
            return p.stdout, 422
        with open(abs_fn, "rb") as f:
            data_after = f.read()
        if data_after == data_before:
            return "Not Modified", 304
        return data_after, 200


# ======================================================================================
if __name__ == "__main__":
    app.run()

# EOF

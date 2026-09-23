# -*- coding: utf-8 -*-
# pytest configuration for the standalone pipeline tests.
#
# paths.py resolves its roots from DHSVM_* environment variables and
# creates the output folder at import time, so the environment has to
# point at a scratch folder before any pipeline module is imported. A
# test that runs a stage sets DHSVM_OUT to its own folder and reloads
# paths; this file only makes the first import harmless and puts
# pipeline/ on sys.path so the flat modules import by name.

import os
import sys
import tempfile
from pathlib import Path

PIPELINE = Path(__file__).resolve().parent.parent / "pipeline"
if str(PIPELINE) not in sys.path:
    sys.path.insert(0, str(PIPELINE))

_scratch = tempfile.mkdtemp(prefix="dhsvm_tests_")
os.environ.setdefault("DHSVM_OUT", _scratch)
os.environ.setdefault("DHSVM_INPUTS", _scratch)
os.environ.setdefault("DHSVM_REF", _scratch)
os.environ.setdefault("DHSVM_EPSG", "32617")

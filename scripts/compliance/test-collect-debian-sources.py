#!/usr/bin/env python3

"""Test Debian Sources control parsing without network access."""

from __future__ import annotations

import importlib.util
from pathlib import Path

script = Path(__file__).with_name("collect-debian-sources.py")
spec = importlib.util.spec_from_file_location("debian_sources", script)
assert spec and spec.loader
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)

text = """Package: example
Version: 1.0-1
Checksums-Sha256:
 aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa 10 example.dsc
 bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb 20 example.tar.xz
Directory: pool/main/e/example

"""
stanzas = module.parse_stanzas(text)
assert len(stanzas) == 1
files = module.stanza_files(stanzas[0])
assert files == [
    {
        "sha256": "a" * 64,
        "size": "10",
        "name": "example.dsc",
    },
    {
        "sha256": "b" * 64,
        "size": "20",
        "name": "example.tar.xz",
    },
]
print("Debian source collector parser test passed.")

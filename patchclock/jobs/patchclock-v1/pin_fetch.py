#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Munim, Inc. Written by AI agents working for MV2 of Munim, Inc.
"""Fetch the frozen patchclock 0.1.0 files at their pinned commits and verify
each file's git blob sha1. Any mismatch or download failure exits non-zero, so
the job never runs code other than what PREREG-PATCH-v0.1 froze."""
import hashlib
import sys
import time
import urllib.request

BASE = "https://raw.githubusercontent.com/munimtechnologies/psp-five-minute-test/"
PINS = (
    ("1ede6866ca98fc1d750ec4ded50d4aa055210b80", "patchclock/patchclock.py",
     "patchclock.py", "baa3cc16783c7e62ce7efe9c8ff7812ccc39bffe"),
    ("1ede6866ca98fc1d750ec4ded50d4aa055210b80", "patchclock/test_patchclock.py",
     "test_patchclock.py", "a93722d3ff6e9e8fa56517b14c0dc539fdd7ddf0"),
    ("8913a2d50735288b3651531ce9ea86d161b1d270", "patchclock/PREREG-PATCH-v0.1.md",
     "PREREG-PATCH-v0.1.md", "ce060d8375928972f51400c89bb24bfb4d397f0e"),
)


def blob_sha1(data: bytes) -> str:
    return hashlib.sha1(b"blob %d\0" % len(data) + data).hexdigest()


def get(url: str) -> bytes:
    last = None
    for attempt in range(4):
        try:
            with urllib.request.urlopen(url, timeout=60) as r:
                return r.read()
        except Exception as exc:  # noqa: BLE001
            last = exc
            time.sleep(5 * (attempt + 1))
    raise RuntimeError(f"download failed: {url}: {last}")


def main() -> int:
    rc = 0
    for commit, path, local, want in PINS:
        data = get(BASE + commit + "/" + path)
        got = blob_sha1(data)
        ok = got == want
        print(f"PIN {local} commit={commit[:10]} blob={got} want={want} {'OK' if ok else 'MISMATCH'}")
        if not ok:
            rc = 1
            continue
        with open(local, "wb") as fh:
            fh.write(data)
    print(f"PIN_RESULT rc={rc}")
    return rc


if __name__ == "__main__":
    sys.exit(main())

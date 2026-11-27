"""Run every verification script and report a combined result."""

import subprocess, sys, os, time

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = ["check_A_gaussian.py", "check_B_smooth.py", "check_C_rough.py", "check_D_combinatorics.py"]

t0 = time.time()
failed = []
total_pass = total_fail = 0
for s in SCRIPTS:
    print("=" * 70)
    print(s)
    print("=" * 70)
    r = subprocess.run([sys.executable, "-I", os.path.join(HERE, s)], capture_output=True, text=True)
    print(r.stdout, end="")
    if r.stderr:
        print(r.stderr, end="")
    last = [l for l in r.stdout.splitlines() if "passed" in l and "failed" in l]
    if last:
        parts = last[-1].split(":")[1].split(",")
        total_pass += int(parts[0].split()[0])
        total_fail += int(parts[1].split()[0])
    if r.returncode != 0:
        failed.append(s)
print("=" * 70)
print(f"TOTAL: {total_pass} checks passed, {total_fail} failed, {time.time() - t0:.1f} s")
if failed:
    print("FAILED scripts:", ", ".join(failed))
    sys.exit(1)
print("All verification scripts passed.")

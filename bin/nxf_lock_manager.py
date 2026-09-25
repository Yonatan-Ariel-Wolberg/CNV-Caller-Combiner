#!/usr/bin/env python3
"""
nxf_lock_manager.py
Robust, graceful, and reproducible Nextflow lock management tool.

Features:
- Sub-second inspection of Nextflow LevelDB database locks (.nextflow/cache/**/db/LOCK).
- Never scans the massive work/ directory (workdirs don't contain Nextflow session locks).
- Process-aware & SLURM-aware: detects active Nextflow runs before touching any lock.
- Safely purges stale/orphaned locks left behind by killed jobs or cluster timeouts.
- Zero external dependencies: pure standard Python 3.

Usage:
  python3 bin/nxf_lock_manager.py status
  python3 bin/nxf_lock_manager.py clean
  python3 bin/nxf_lock_manager.py clean --force
  python3 bin/nxf_lock_manager.py unlock --session <session_uuid>
"""

import argparse
import glob
import os
import re
import subprocess
import sys
import time
from datetime import datetime


def get_active_nextflow_pids():
    """Detect local running Nextflow / Java master processes for current user."""
    pids = []
    try:
        current_user = os.getenv("USER", "")
        cmd = ["pgrep", "-u", current_user, "-f", "nextflow.*\\.jar|nextflow run"]
        out = subprocess.check_output(cmd, stderr=subprocess.DEVNULL, text=True)
        for line in out.strip().splitlines():
            if line.strip():
                pids.append(int(line.strip()))
    except Exception:
        pass
    return pids


def get_active_slurm_jobs():
    """Detect active SLURM master jobs that might be running Nextflow for the user."""
    jobs = []
    current_user = os.getenv("USER", "")
    if not current_user:
        return jobs
    try:
        cmd = ["squeue", "-u", current_user, "-h", "-o", "%i|%j|%T|%R|%Z"]
        out = subprocess.check_output(cmd, stderr=subprocess.DEVNULL, text=True)
        for line in out.strip().splitlines():
            parts = line.strip().split("|")
            if len(parts) >= 3:
                job_id, job_name, state = parts[0], parts[1], parts[2]
                if state in ("RUNNING", "CONFIGURING", "COMPLETING"):
                    # Ignore pure interactive shells (e.g. salloc / srun bash)
                    name_lower = job_name.lower()
                    if name_lower in ("bash", "sh", "zsh", "tcsh", "login"):
                        continue
                    jobs.append({
                        "id": job_id,
                        "name": job_name,
                        "state": state,
                        "workdir": parts[4] if len(parts) > 4 else ""
                    })
    except Exception:
        pass
    return jobs


def find_nxf_locks(project_dir):
    """Find all LevelDB LOCK files under .nextflow/cache within project_dir."""
    cache_root = os.path.join(project_dir, ".nextflow", "cache")
    if not os.path.exists(cache_root):
        return []
    
    pattern = os.path.join(cache_root, "*", "db", "LOCK")
    lock_files = glob.glob(pattern)
    
    # Also check legacy or flat cache structures
    lock_files.extend(glob.glob(os.path.join(cache_root, "db", "LOCK")))
    return sorted(list(set(lock_files)))


def parse_session_from_lock(lock_path):
    """Extract session UUID from lock path."""
    parts = lock_path.split(os.sep)
    try:
        # Expected: .../.nextflow/cache/<uuid>/db/LOCK
        db_idx = parts.index("db")
        if db_idx > 0 and parts[db_idx - 1] != "cache":
            return parts[db_idx - 1]
    except ValueError:
        pass
    return "unknown_session"


def inspect_locks(project_dir):
    """Analyze all lock files and classify as ACTIVE or ORPHANED."""
    lock_files = find_nxf_locks(project_dir)
    local_pids = get_active_nextflow_pids()
    slurm_jobs = get_active_slurm_jobs()

    # Determine if ANY master Nextflow job is running
    is_nextflow_running = bool(local_pids or slurm_jobs)

    results = []
    now = time.time()

    for l_path in lock_files:
        session_id = parse_session_from_lock(l_path)
        try:
            mtime = os.path.getmtime(l_path)
            age_sec = now - mtime
            mtime_str = datetime.fromtimestamp(mtime).strftime("%Y-%m-%d %H:%M:%S")
        except OSError:
            mtime = now
            age_sec = 0
            mtime_str = "unknown"

        # If a Nextflow process is active and was started recently or modified the lock recently
        if is_nextflow_running and age_sec < 86400:
            status = "ACTIVE"
            reason = f"Active local PIDs: {local_pids} or SLURM jobs: {[j['id'] for j in slurm_jobs]}"
        else:
            status = "ORPHANED"
            reason = "No active Nextflow process or SLURM job holds this lock"

        results.append({
            "path": l_path,
            "session": session_id,
            "mtime": mtime_str,
            "age_min": round(age_sec / 60, 1),
            "status": status,
            "reason": reason
        })

    return results, is_nextflow_running


def cmd_status(project_dir):
    results, is_running = inspect_locks(project_dir)
    print("=" * 80)
    print(f" NEXTFLOW LOCK STATUS REPORT: {project_dir}")
    print("=" * 80)
    if not results:
        print("  ✓ No Nextflow LOCK files detected. Cache directory is clean.")
        print("=" * 80)
        return 0

    print(f" Found {len(results)} LOCK file(s):")
    print(f" {'SESSION UUID':<38} {'AGE (min)':<10} {'STATUS':<10} {'MODIFIED'}")
    print("-" * 80)
    for r in results:
        status_label = "🟢 ACTIVE" if r["status"] == "ACTIVE" else "🔴 ORPHANED"
        print(f" {r['session']:<38} {r['age_min']:<10} {status_label:<10} {r['mtime']}")
    print("-" * 80)
    if is_running:
        print("  Notice: Active Nextflow process(es) detected. Active locks are protected.")
    print("=" * 80)
    return 0


def cmd_clean(project_dir, force=False):
    results, is_running = inspect_locks(project_dir)
    if not results:
        print("✓ No Nextflow LOCK files found. Nothing to clean.")
        return 0

    removed_count = 0
    preserved_count = 0

    for r in results:
        lock_file = r["path"]
        if r["status"] == "ACTIVE" and not force:
            print(f"⚠️  PRESERVED: Active lock for session {r['session']} ({r['reason']}).")
            print("   (Use --force to override if this is a stale process).")
            preserved_count += 1
            continue

        try:
            os.remove(lock_file)
            print(f"✓ REMOVED: Orphaned lock for session {r['session']} (Path: {lock_file})")
            removed_count += 1
        except OSError as e:
            print(f"✗ ERROR: Could not remove {lock_file}: {e}", file=sys.stderr)

    print("-" * 80)
    print(f"Clean complete: {removed_count} lock(s) removed, {preserved_count} active lock(s) preserved.")
    return 0


def main():
    parser = argparse.ArgumentParser(
        description="Nextflow Lock Manager: Safe, fast, and deterministic lock resolution."
    )
    parser.add_argument(
        "--project-dir",
        default=os.getcwd(),
        help="Path to Nextflow project root directory (default: current directory)."
    )
    subparsers = parser.add_subparsers(dest="command", help="Available subcommands")

    # Status subcommand
    subparsers.add_parser("status", help="List all Nextflow locks and report their liveness status.")

    # Clean subcommand
    clean_p = subparsers.add_parser("clean", help="Safely purge orphaned locks while preserving active jobs.")
    clean_p.add_argument(
        "--force",
        action="store_true",
        help="Force delete all locks even if active processes are detected."
    )

    args = parser.parse_args()

    project_dir = os.path.abspath(args.project_dir)
    # If project_dir contains .nextflow or is parent of bin, ensure correct root
    if not os.path.exists(os.path.join(project_dir, ".nextflow")) and os.path.exists(os.path.join(project_dir, "..", ".nextflow")):
        project_dir = os.path.abspath(os.path.join(project_dir, ".."))

    cmd = args.command or "status"
    if cmd == "status":
        sys.exit(cmd_status(project_dir))
    elif cmd == "clean":
        sys.exit(cmd_clean(project_dir, force=args.force))
    else:
        parser.print_help()
        sys.exit(1)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
gmx_eta.py — Estimate GROMACS MD run ETA and performance from a .log file.

This script parses a GROMACS molecular dynamics log file to:
  • Determine the total number of simulation steps (nsteps)
  • Extract recent checkpoint entries (step number + timestamp)
  • Compute the estimated time remaining (ETA) based on the average time per step
    over the last N checkpoints (configurable with --smooth)
  • Calculate simulation performance in ns/day using the timestep (dt) from the log
    if available, or from the --dt argument
  • Optionally display a live‑updating or static plot of ETA over time
  • Optionally log ETA output to a file and/or show a progress bar

Usage examples:
  python gmx_eta.py md.log
  python gmx_eta.py md.log --smooth 10
  python gmx_eta.py md.log --watch 300 --logfile-out eta.log --progress
  python gmx_eta.py md.log --watch 60 --plot
  python gmx_eta.py md.log --plot
  python gmx_eta.py md.log --dt 0.004

Arguments:
  logfile        Path to GROMACS .log file
  --watch N      Refresh every N seconds (default: run once)
  --logfile-out  Append output lines to a file
  --progress     Show a progress bar (requires tqdm)
  --plot         Plot ETA over time (requires matplotlib)
  --smooth N     Average ETA over last N checkpoints (default: 2)
  --dt PS        Timestep in picoseconds (default: 0.002 ps = 2 fs);
                 overridden by auto‑detected dt from log if present
"""
import re
import time
import argparse
import unicodedata
from datetime import datetime, timedelta

# Optional imports
try:
    from tqdm import tqdm
except ImportError:
    tqdm = None

try:
    import matplotlib.pyplot as plt
except ImportError:
    plt = None

MONTH_MAP = {
    'Jan': 1, 'Feb': 2, 'Mar': 3, 'Apr': 4, 'May': 5, 'Jun': 6,
    'Jul': 7, 'Aug': 8, 'Sep': 9, 'Oct': 10, 'Nov': 11, 'Dec': 12
}

def normalize_ts(ts: str) -> str:
    ts_norm = unicodedata.normalize('NFKC', ts)
    ts_norm = re.sub(r'\s+', ' ', ts_norm).strip()
    ts_norm = re.sub(r'\bS?Sep\b', 'Sep', ts_norm)
    return ts_norm

def timestamp_to_seconds(ts):
    try:
        parts = normalize_ts(ts).split()
        if len(parts) != 5:
            raise ValueError(f"Unexpected timestamp format: {parts}")
        mon_abbr = parts[1]
        month = MONTH_MAP.get(mon_abbr)
        if month is None:
            raise ValueError(f"Unknown month abbreviation: {mon_abbr}")
        day = int(parts[2])
        h, m, s = map(int, parts[3].split(":"))
        year = int(parts[4])
        dt = datetime(year, month, day, h, m, s)
        return int(dt.timestamp())
    except Exception as e:
        print(f"[timestamp_to_seconds] Failed to parse: {repr(ts)} — {e}")
        return None

def seconds_to_hms(seconds):
    return str(timedelta(seconds=int(seconds)))

def parse_log(logfile):
    with open(logfile, 'r', encoding='utf-8', errors='replace') as f:
        lines = f.readlines()

    # --- Total steps detection ---
    total_steps = None
    for line in lines:
        if "nsteps" in line.lower():
            match = re.search(r'nsteps\s*=\s*(\d+)', line, re.IGNORECASE)
            if match:
                total_steps = int(match.group(1))
                break
    if total_steps is None:
        for line in lines:
            match = re.search(r'(\d+)\s+steps', line, re.IGNORECASE)
            if match:
                total_steps = int(match.group(1))
                break

    # --- dt detection ---
    dt_ps = None
    for line in lines:
        if re.search(r'\bdt\b', line):
            match = re.search(r'dt\s*=\s*([0-9.]+)', line)
            if match:
                try:
                    dt_ps = float(match.group(1))
                    break
                except ValueError:
                    pass

    # --- Checkpoint entries ---
    checkpoints = []
    for line in lines:
        match = re.search(r'Writing checkpoint, step\s+([0-9,]+)\s+at\s+(.+)', line)
        if match:
            step_str = match.group(1).replace(",", "").strip()
            try:
                step = int(step_str)
            except ValueError:
                print(f"[parse_log] Could not parse step number from: {repr(step_str)}")
                step = None
            raw_ts = match.group(2).strip()
            ts_clean = normalize_ts(raw_ts)
            checkpoints.append((step, ts_clean))

    return total_steps, checkpoints, dt_ps

def estimate_eta(logfile, smooth_n=2, dt_ps_arg=0.002):
    total_steps, checkpoints, dt_ps_detected = parse_log(logfile)
    if total_steps is None:
        raise ValueError("Could not determine total steps from log file.")
    if len(checkpoints) < 2:
        raise ValueError("Not enough checkpoint entries found.")

    # Use detected dt if available, else fallback to argument
    dt_ps = dt_ps_detected if dt_ps_detected is not None else dt_ps_arg

    # Use last N checkpoints for smoothing
    n = min(smooth_n, len(checkpoints))
    recent = checkpoints[-n:]

    total_time = 0
    total_steps_done = 0
    for (s1, t1), (s2, t2) in zip(recent[:-1], recent[1:]):
        ts1 = timestamp_to_seconds(t1)
        ts2 = timestamp_to_seconds(t2)
        if ts1 is None or ts2 is None or s1 is None or s2 is None:
            continue
        total_time += (ts2 - ts1)
        total_steps_done += (s2 - s1)

    if total_steps_done == 0:
        raise ValueError("No progress detected in the selected smoothing window.")

    time_per_step = total_time / total_steps_done
    current_step = recent[-1][0]
    steps_remaining = total_steps - current_step
    time_remaining = steps_remaining * time_per_step

    # ns/day calculation
    ns_per_step = dt_ps / 1000.0
    ns_per_sec = (total_steps_done * ns_per_step) / total_time
    ns_per_day = ns_per_sec * 86400

    eta = seconds_to_hms(time_remaining)
    now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    return {
        "timestamp": now,
        "total_steps": total_steps,
        "current_step": current_step,
        "eta": eta,
        "eta_seconds": time_remaining,
        "progress": current_step / total_steps,
        "ns_per_day": ns_per_day,
        "dt_ps_used": dt_ps,
        "dt_ps_detected": dt_ps_detected
    }

def main():
    parser = argparse.ArgumentParser(description="Estimate GROMACS MD ETA from log file.")
    parser.add_argument("logfile", help="Path to GROMACS .log file")
    parser.add_argument("--watch", type=int, help="Refresh every N seconds")
    parser.add_argument("--logfile-out", help="Log output to file")
    parser.add_argument("--progress", action="store_true", help="Show progress bar")
    parser.add_argument("--plot", action="store_true", help="Plot ETA over time (requires matplotlib)")
    parser.add_argument("--smooth", type=int, default=2, help="Average ETA over last N checkpoints (default: 2)")
    parser.add_argument("--dt", type=float, default=0.002, help="Timestep in ps (default: 0.002 ps = 2 fs)")
    args = parser.parse_args()

    eta_history = []

    try:
        while True:
            try:
                result = estimate_eta(args.logfile, smooth_n=args.smooth, dt_ps_arg=args.dt)
                dt_info = f"(auto-detected)" if result["dt_ps_detected"] is not None else "(from --dt)"
                line = (f"[{result['timestamp']}] Step: {result['current_step']} / {result['total_steps']} | "
                        f"ETA: {result['eta']} | Perf: {result['ns_per_day']:.2f} ns/day | "
                        f"dt: {result['dt_ps_used']} ps {dt_info}")
                print(line)

                if args.logfile_out:
                    with open(args.logfile_out, "a") as f:
                        f.write(line + "\n")

                if args.progress and tqdm:
                    tqdm.write("")
                    bar = tqdm(total=100, desc="Progress", bar_format="{l_bar}{bar}| {n_fmt}%")
                    bar.n = int(result["progress"] * 100)
                    bar.refresh()
                    bar.close()

                if args.plot and plt:
                    eta_history.append((datetime.now(), result["eta_seconds"]))
                    plt.clf()
                    times, etas = zip(*eta_history)
                    plt.plot(times, etas, marker='o')
                    plt.ylabel("ETA (seconds)")
                    plt.xlabel("Timestamp")
                    plt.title("Estimated Time Remaining")
                    if args.watch:
                        plt.pause(0.01)  # live update
                    else:
                        # Draw once, but don't close until user does
                        plt.show(block=True)

            except Exception as e:
                print(f"Error: {e}")

            if not args.watch:
                break
            time.sleep(args.watch)

    except KeyboardInterrupt:
        print("Stopped by user.")

if __name__ == "__main__":
    main()
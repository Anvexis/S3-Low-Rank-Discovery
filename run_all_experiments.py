# run_all_experiments.py
"""
Orchestrator: runs all empirical tests sequentially with progress tracking.
Outputs are logged to results/logs/ and CSVs are saved to results/.
"""
import subprocess, sys, os
from tqdm import tqdm
from pathlib import Path
from datetime import datetime

EXPERIMENTS = [
    ("S3 Verified Streaming Rank", "experiments/streaming_rank_verified.py"),
    ("S5 Rank Test", "experiments/test_S5_rank.py"),
    ("S6 Rank Test", "experiments/test_S6_rank.py"),
    ("Sparse Large P (p≈2^80)", "experiments/sparse_large_p.py"),
    ("Fast Streaming Solver", "experiments/fast_streaming_solver.py"),
]

LOG_DIR = Path("results/logs")
LOG_DIR.mkdir(parents=True, exist_ok=True)

def run_experiment(name, script_path):
    log_file = LOG_DIR / f"{name.replace(' ', '_').lower()}.log"
    print(f"\n🚀 Running: {name}")
    print(f"📝 Log: {log_file}")
    
    with open(log_file, "w") as f:
        f.write(f"# {name} | {datetime.now().isoformat()}\n")
        f.write(f"# Command: python3 {script_path}\n\n")
        try:
            result = subprocess.run(
                [sys.executable, script_path],
                stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                text=True, check=True, timeout=3600
            )
            f.write(result.stdout)
            print(f"✅ {name} completed successfully.")
            return True
        except subprocess.TimeoutExpired:
            f.write("❌ TIMEOUT (>1 hour)\n")
            print(f"⏳ {name} timed out.")
            return False
        except subprocess.CalledProcessError as e:
            f.write(f"❌ ERROR (exit {e.returncode})\n{e.stdout}")
            print(f"❌ {name} failed. Check log.")
            return False

def main():
    print("="*60)
    print("🔬 Semaev Low-Rank Structure: Full Experiment Suite")
    print("="*60)
    
    success_count = 0
    for name, script in tqdm(EXPERIMENTS, desc="Experiments", unit="test"):
        if not Path(script).exists():
            print(f"⚠️ Skipping {name}: {script} not found.")
            continue
        if run_experiment(name, script):
            success_count += 1
            
    print("\n" + "="*60)
    print(f"🏁 Finished: {success_count}/{len(EXPERIMENTS)} experiments succeeded.")
    print(f"📂 Logs: {LOG_DIR}")
    print(f"📊 Data: results/*.csv")
    print("="*60)

if __name__ == "__main__":
    main()
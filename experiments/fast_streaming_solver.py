# fast_streaming_solver.py
import random, time, csv, os, psutil
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed

try:
    import gmpy2
    HAS_GMPY2 = True
    mpz = gmpy2.mpz
    mod_inv = lambda a, p: int(gmpy2.invert(a, mpz(p)))
    mod_mul = lambda a, b, p: int((mpz(a) * mpz(b)) % mpz(p))
    mod_sub = lambda a, b, p: int((mpz(a) - mpz(b)) % mpz(p))
    print("✅ gmpy2 active")
except ImportError:
    HAS_GMPY2 = False
    mpz = lambda x: x
    mod_inv = lambda a, p: pow(a, -1, p)
    mod_mul = lambda a, b, p: (a * b) % p
    mod_sub = lambda a, b, p: (a - b) % p
    print("⚠️ gmpy2 not found")

# =============================================================================
# CORE MATH
# =============================================================================
def S3(x1, x2, x3, a, b, p):
    e1 = (x1 + x2 + x3) % p
    e2 = (x1*x2 + x2*x3 + x3*x1) % p
    e3 = (x1*x2*x3) % p
    return (pow(e2 - a, 2, p) - 4 * e1 * (e3 + b)) % p

def find_prime(bits):
    target = 1 << bits
    start = target + (1 - target % 6) % 6
    for p in range(start, start + 2000000, 6):
        if pow(2, p-1, p) == 1: return p
    return None

def generate_valid_B(p, n, a, b, max_att=3000000):
    B, seen = [], set()
    for _ in range(max_att):
        if len(B) >= n: break
        x = random.randint(0, p-1)
        if x in seen: continue
        y2 = (x**3 + a*x + b) % p
        if y2 == 0 or pow(y2, (p-1)//2, p) == 1:
            B.append(x); seen.add(x)
    return B

# =============================================================================
# PARALLEL WORKER
# =============================================================================
def generate_chunk_rows(args):
    """Worker: generates and pre-reduces rows for a chunk of indices"""
    idx_start, idx_end, B, xR, a, b, p, basis = args
    rows = []
    for i in range(idx_start, idx_end):
        row = [S3(B[i], B[j], xR, a, b, p) for j in range(len(B))]
        # Reduce against current basis
        for piv_col, b_row in basis:
            if row[piv_col] != 0:
                f = mod_mul(row[piv_col], mod_inv(b_row[piv_col], p), p)
                row = [mod_sub(row[k], mod_mul(f, b_row[k], p), p) for k in range(len(B))]
        rows.append((i, row))
    return rows

# =============================================================================
# STREAMING SOLVER WITH EARLY EXIT
# =============================================================================
def streaming_rank_parallel(n, p, a, b, xR, B, chunk_size=100, stable_window=80, max_workers=None):
    basis = []  # List of (pivot_col, row)
    rows_since_pivot = 0
    total_processed = 0
    
    workers = max_workers or max(1, os.cpu_count() - 1)
    
    with ProcessPoolExecutor(max_workers=workers) as executor:
        for start in range(0, n, chunk_size):
            end = min(start + chunk_size, n)
            future = executor.submit(generate_chunk_rows, (start, end, B, xR, a, b, p, basis))
            
            for i, row in future.result():
                total_processed += 1
                pivot = next((k for k in range(n) if row[k] != 0), None)
                
                if pivot is not None:
                    basis.append((pivot, row))
                    rows_since_pivot = 0
                    print(f"  📈 Pivot found at col {pivot} | Rank = {len(basis)} | Rows processed: {total_processed}")
                else:
                    rows_since_pivot += 1
                    
                if rows_since_pivot >= stable_window:
                    print(f"  ⏹️ Early exit: rank stabilized at {len(basis)} after {stable_window} zero-reduced rows.")
                    return len(basis)
                    
            if total_processed >= n: break
            
    return len(basis)

# =============================================================================
# MAIN
# =============================================================================
def main():
    random.seed(1337)
    a, b = 0, 7
    Path("results").mkdir(exist_ok=True)
    
    p = find_prime(64)
    if not p: print("❌ Prime ~2^64 not found"); return
    print(f"⚡ Field: p ≈ 2^64 ({p})\n")
    
    results = []
    for n in [2000, 3000, 5000]:
        B = generate_valid_B(p, n, a, b)
        xR = random.randint(0, p-1)
        
        mem_before = psutil.Process(os.getpid()).memory_info().rss / 1024**2
        t0 = time.time()
        r = streaming_rank_parallel(n, p, a, b, xR, B, chunk_size=200, stable_window=100)
        dt = time.time() - t0
        mem_after = psutil.Process(os.getpid()).memory_info().rss / 1024**2
        
        print(f"  ✅ |B|={n:4d} → Rank = {r} | Time: {dt:.2f}s | ΔMem: {mem_after-mem_before:.1f}MB\n")
        results.append({"n": n, "rank": r, "time_s": f"{dt:.2f}", "mem_delta_MB": f"{mem_after-mem_before:.1f}"})
        
    csv_path = Path("results/fast_streaming_results.csv")
    with open(csv_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=results[0].keys())
        w.writeheader(); w.writerows(results)
    print(f"💾 Saved → {csv_path}")

if __name__ == "__main__":
    main()
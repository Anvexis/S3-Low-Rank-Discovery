# sparse_large_p.py
import random, time, csv, sys, os, psutil
from pathlib import Path

try:
    import gmpy2
    HAS_GMPY2 = True
    mpz = gmpy2.mpz
    mod_inv = lambda a, p: int(gmpy2.invert(a, mpz(p)))
    mod_mul = lambda a, b, p: int((mpz(a) * mpz(b)) % mpz(p))
    mod_sub = lambda a, b, p: int((mpz(a) - mpz(b)) % mpz(p))
except ImportError:
    HAS_GMPY2 = False
    mpz = lambda x: x
    mod_inv = lambda a, p: pow(a, -1, p)
    mod_mul = lambda a, b, p: (a * b) % p
    mod_sub = lambda a, b, p: (a - b) % p

print(f"✅ {'gmpy2' if HAS_GMPY2 else 'Python'} arithmetic loaded")

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

def generate_valid_B(p, n, a, b, max_att=2000000):
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
# STREAMING LOW-RANK SOLVER (Matrix-Free)
# =============================================================================
def streaming_low_rank(n, p, a, b, xR, B, max_rank=10):
    basis = []  # List of (row, pivot_col)
    for i in range(n):
        # Generate row i on-the-fly
        row = [S3(B[i], B[j], xR, a, b, p) for j in range(n)]
        
        # Reduce against existing basis
        for b_row, piv_col in basis:
            if row[piv_col] != 0:
                factor = mod_mul(row[piv_col], mod_inv(b_row[piv_col], p), p)
                for k in range(n):
                    row[k] = mod_sub(row[k], mod_mul(factor, b_row[k], p), p)
                    
        # Find new pivot
        pivot = next((k for k in range(n) if row[k] != 0), None)
        if pivot is not None:
            basis.append((row, pivot))
            if len(basis) >= max_rank:
                return max_rank
    return len(basis)

# =============================================================================
# MAIN EXPERIMENT
# =============================================================================
def main():
    random.seed(1337)
    a, b = 0, 7
    Path("results").mkdir(exist_ok=True)
    
    p = find_prime(80)
    if not p: print("❌ Prime ~2^80 not found"); return
    print(f"⚡ Field: p ≈ 2^80 ({p})\n")
    
    results = []
    for n in [200, 500, 1000]:
        B = generate_valid_B(p, n, a, b)
        xR = random.randint(0, p-1)
        
        mem_before = psutil.Process(os.getpid()).memory_info().rss / 1024**2
        t0 = time.time()
        r = streaming_low_rank(n, p, a, b, xR, B, max_rank=8)
        dt = time.time() - t0
        mem_after = psutil.Process(os.getpid()).memory_info().rss / 1024**2
        
        print(f"  |B|={n:4d} → Rank = {r} | Time: {dt:.2f}s | ΔMem: {mem_after-mem_before:.1f}MB")
        results.append({"n": n, "rank": r, "time_s": f"{dt:.2f}", "mem_delta_MB": f"{mem_after-mem_before:.1f}"})
        
    csv_path = Path("results/sparse_large_p.csv")
    with open(csv_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=results[0].keys())
        w.writeheader(); w.writerows(results)
    print(f"\n💾 Saved → {csv_path}")

if __name__ == "__main__":
    main()
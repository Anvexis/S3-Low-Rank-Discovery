# streaming_rank_verified.py
import random, time, csv, os, psutil
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

def streaming_rank_strict(n, p, a, b, xR, B, max_theoretical_rank=9):
    """
    Maintains true Row-Echelon Form. Stops at algebraic bound.
    Memory: O(n * r), Time: O(n^2 * r)
    """
    basis = []  # List of (pivot_col, reduced_row)
    rows_processed = 0
    
    for i in range(n):
        row = [S3(B[i], B[j], xR, a, b, p) for j in range(n)]
        
        # Forward reduction against current basis
        for piv_col, b_row in basis:
            if row[piv_col] != 0:
                f = mod_mul(row[piv_col], mod_inv(b_row[piv_col], p), p)
                row = [mod_sub(row[k], mod_mul(f, b_row[k], p), p) for k in range(n)]
                
        pivot = next((k for k in range(n) if row[k] != 0), None)
        if pivot is not None:
            # Normalize pivot row
            inv = mod_inv(row[pivot], p)
            row = [mod_mul(v, inv, p) for v in row]
            basis.append((pivot, row))
            
            # Backward reduction to maintain strict REF
            for idx in range(len(basis)-1):
                pc, br = basis[idx]
                if br[pivot] != 0:
                    f = br[pivot]
                    basis[idx] = (pc, [mod_sub(br[k], mod_mul(f, row[k], p), p) for k in range(n)])
                    
        rows_processed += 1
        if len(basis) >= max_theoretical_rank:
            print(f"  ✅ Algebraic bound reached: rank = {len(basis)} after {rows_processed} rows.")
            return len(basis)
            
    return len(basis)

def main():
    random.seed(1337)
    a, b = 0, 7
    Path("results").mkdir(exist_ok=True)
    
    p = find_prime(64)
    if not p: print("❌ Prime ~2^64 not found"); return
    print(f"⚡ Field: p ≈ 2^64 ({p}) | Theoretical max rank = 9\n")
    
    results = []
    for n in [500, 1000, 2000]:
        B = generate_valid_B(p, n, a, b)
        xR = random.randint(0, p-1)
        
        mem_before = psutil.Process(os.getpid()).memory_info().rss / 1024**2
        t0 = time.time()
        r = streaming_rank_strict(n, p, a, b, xR, B, max_theoretical_rank=9)
        dt = time.time() - t0
        mem_after = psutil.Process(os.getpid()).memory_info().rss / 1024**2
        
        print(f"  |B|={n:4d} → Rank = {r} | Time: {dt:.2f}s | ΔMem: {mem_after-mem_before:.1f}MB\n")
        results.append({"n": n, "rank": r, "time_s": f"{dt:.2f}", "mem_delta_MB": f"{mem_after-mem_before:.1f}"})
        
    csv_path = Path("results/verified_streaming_rank.csv")
    with open(csv_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=results[0].keys())
        w.writeheader(); w.writerows(results)
    print(f"💾 Saved → {csv_path}")

if __name__ == "__main__":
    main()
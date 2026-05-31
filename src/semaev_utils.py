# src/semaev_utils.py
"""
Shared utilities for Semaev summation polynomial experiments.
Provides modular arithmetic, prime generation, factor base sampling,
and exact rank computation over finite fields.
"""
import random, os, psutil
from pathlib import Path

# =============================================================================
# 1. MODULAR ARITHMETIC (gmpy2 accelerated with native fallback)
# =============================================================================
try:
    import gmpy2
    HAS_GMPY2 = True
    _mpz = gmpy2.mpz
    mod_inv = lambda a, p: int(gmpy2.invert(_mpz(a), _mpz(p)))
    mod_mul = lambda a, b, p: int((_mpz(a) * _mpz(b)) % _mpz(p))
    mod_sub = lambda a, b, p: int((_mpz(a) - _mpz(b)) % _mpz(p))
except ImportError:
    HAS_GMPY2 = False
    mod_inv = lambda a, p: pow(a, -1, p)
    mod_mul = lambda a, b, p: (a * b) % p
    mod_sub = lambda a, b, p: (a - b) % p

# =============================================================================
# 2. CORE POLYNOMIALS
# =============================================================================
def S3(x1, x2, x3, a, b, p):
    """Semaev's 3rd summation polynomial: S3=0 ⇔ P1+P2+P3=O"""
    e1 = (x1 + x2 + x3) % p
    e2 = (x1*x2 + x2*x3 + x3*x1) % p
    e3 = (x1*x2*x3) % p
    return (pow(e2 - a, 2, p) - 4 * e1 * (e3 + b)) % p

# =============================================================================
# 3. FIELD & FACTOR BASE UTILS
# =============================================================================
def find_prime(bits, mod_base=6, mod_target=1, search_range=2000000):
    """Find prime p ≡ mod_target (mod mod_base) near 2^bits"""
    target = 1 << bits
    offset = (mod_target - target % mod_base) % mod_base
    start = target + offset
    for p in range(start, start + search_range, mod_base):
        if pow(2, p-1, p) == 1:  # Fermat primality check
            return p
    return None

def generate_valid_B(p, n, a, b, max_att=3000000):
    """Sample n unique x-coordinates where x^3+ax+b is a quadratic residue"""
    B, seen = [], set()
    for _ in range(max_att):
        if len(B) >= n: break
        x = random.randint(0, p-1)
        if x in seen: continue
        y2 = (x**3 + a*x + b) % p
        if y2 == 0 or pow(y2, (p-1)//2, p) == 1:
            B.append(x); seen.add(x)
    return B

def get_memory_mb():
    """Current process RSS memory in MB"""
    return psutil.Process(os.getpid()).memory_info().rss / 1024**2

# =============================================================================
# 4. EXACT RANK SOLVER (Gaussian Elimination over F_p)
# =============================================================================
def gf_rank(matrix, p):
    """Compute exact matrix rank over F_p. Returns integer rank."""
    M = [row[:] for row in matrix]
    rows, cols = len(M), len(M[0])
    rank = r = 0
    for c in range(cols):
        if r >= rows: break
        pivot = next((i for i in range(r, rows) if M[i][c] % p != 0), -1)
        if pivot == -1: continue
        M[r], M[pivot] = M[pivot], M[r]
        inv = mod_inv(M[r][c], p)
        M[r] = [mod_mul(v, inv, p) for v in M[r]]
        for i in range(rows):
            if i != r and M[i][c]:
                f = M[i][c]
                M[i] = [mod_sub(M[i][j], mod_mul(f, M[r][j], p), p) for j in range(cols)]
        rank += 1; r += 1
    return rank

def save_results(results, filename, folder="results"):
    """Save list of dicts to CSV"""
    Path(folder).mkdir(exist_ok=True)
    import csv
    csv_path = Path(folder) / filename
    with open(csv_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=results[0].keys())
        writer.writeheader(); writer.writerows(results)
    return csv_path
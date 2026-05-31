# test_S6_rank.py
import random, time, csv
from pathlib import Path

try:
    import gmpy2
    HAS_GMPY2 = True
    print("✅ gmpy2 active (accelerating rank solver)")
except ImportError:
    HAS_GMPY2 = False
    print("⚠️ gmpy2 not found (using native Python)")

# =============================================================================
# 1. MODULAR LINEAR ALGEBRA
# =============================================================================
def det_mod(M, p):
    """Determinant of square matrix over F_p via Gaussian elimination"""
    n = len(M)
    M = [row[:] for row in M]
    det = 1
    for i in range(n):
        pivot = next((k for k in range(i, n) if M[k][i] % p != 0), None)
        if pivot is None: return 0
        if pivot != i:
            M[i], M[pivot] = M[pivot], M[i]
            det = (-det) % p
        det = (det * M[i][i]) % p
        inv = pow(M[i][i], -1, p)
        for k in range(i+1, n):
            if M[k][i] != 0:
                f = (M[k][i] * inv) % p
                for j in range(i, n):
                    M[k][j] = (M[k][j] - f * M[i][j]) % p
    return det % p

def poly_resultant(c1, c2, p):
    """Resultant of two polynomials c1, c2 (coeffs [c0, c1, ..., cd])"""
    d1, d2 = len(c1)-1, len(c2)-1
    n = d1 + d2
    if n == 0: return 1
    syl = [[0]*n for _ in range(n)]
    for i in range(d2):
        for j in range(d1+1): syl[i][i+j] = c1[j] % p
    for i in range(d1):
        for j in range(d2+1): syl[d2+i][i+j] = c2[j] % p
    return det_mod(syl, p)

# =============================================================================
# 2. SEMAEV POLYNOMIALS (Recursive via Resultants)
# =============================================================================
def S3_coeffs_X(x1, x2, a, b, p):
    """S3(x1,x2,X) = A*X^2 + B*X + C → returns [C, B, A]"""
    s = (x1 + x2) % p
    prod = (x1 * x2) % p
    A = (s*s - 4*prod) % p
    B = (2*s*(prod - a) - 4*(s*prod + b)) % p
    C = ((prod - a)**2 - 4*s*b) % p
    return [C, B, A]

def interpolate_poly(eval_func, xs, p):
    """Lagrange interpolation over F_p for points (x, eval_func(x))"""
    n = len(xs)
    coeffs = [0]*n
    for i in range(n):
        xi, yi = xs[i], eval_func(xs[i])
        L = 1
        for j in range(n):
            if i == j: continue
            xj = xs[j]
            L = (L * (-xj) * pow(xi - xj, -1, p)) % p
        for k in range(n):
            coeffs[k] = (coeffs[k] + yi * L) % p
    return coeffs

def S4_eval(x3, x4, x5, X, a, b, p):
    return poly_resultant(S3_coeffs_X(x3, x4, a, b, p), S3_coeffs_X(x5, X, a, b, p), p)

def S4_coeffs_X(x3, x4, x5, a, b, p):
    return interpolate_poly(lambda X: S4_eval(x3, x4, x5, X, a, b, p), range(5), p)

def S5_eval(x3, x4, x5, x6, X, a, b, p):
    return poly_resultant(S3_coeffs_X(x3, x4, a, b, p), S4_coeffs_X(x5, x6, X, a, b, p), p)

def S5_coeffs_X(x3, x4, x5, x6, a, b, p):
    return interpolate_poly(lambda X: S5_eval(x3, x4, x5, x6, X, a, b, p), range(9), p)

def S6_eval(x1, x2, x3, x4, x5, x6, a, b, p):
    """S6 = Res_X(S3(x1,x2,X), S5(x3,x4,x5,x6,X))"""
    return poly_resultant(S3_coeffs_X(x1, x2, a, b, p), S5_coeffs_X(x3, x4, x5, x6, a, b, p), p)

def S3_eval(x1, x2, x3, a, b, p):
    e1 = (x1 + x2 + x3) % p
    e2 = (x1*x2 + x2*x3 + x3*x1) % p
    e3 = (x1*x2*x3) % p
    return (pow(e2 - a, 2, p) - 4 * e1 * (e3 + b)) % p

# =============================================================================
# 3. RANK SOLVER (gmpy2 accelerated)
# =============================================================================
def gf_rank(matrix, p):
    M = [row[:] for row in matrix]
    rows, cols = len(M), len(M[0])
    rank = r = 0
    if HAS_GMPY2:
        p_mpz = gmpy2.mpz(p)
        M = [[gmpy2.mpz(v) for v in row] for row in M]
        inv_fn = lambda x: int(gmpy2.invert(x, p_mpz))
        mul_fn = lambda a, b: int((a * b) % p_mpz)
        sub_fn = lambda a, b: int((a - b) % p_mpz)
    else:
        inv_fn = lambda x: pow(x, -1, p)
        mul_fn = lambda a, b: (a * b) % p
        sub_fn = lambda a, b: (a - b) % p

    for c in range(cols):
        if r >= rows: break
        pivot = next((i for i in range(r, rows) if M[i][c] % p != 0), -1)
        if pivot == -1: continue
        M[r], M[pivot] = M[pivot], M[r]
        inv = inv_fn(M[r][c])
        M[r] = [mul_fn(v, inv) for v in M[r]]
        for i in range(rows):
            if i != r and M[i][c]:
                f = M[i][c]
                M[i] = [sub_fn(M[i][j], mul_fn(f, M[r][j])) for j in range(cols)]
        rank += 1; r += 1
    return rank

# =============================================================================
# 4. HELPERS & MAIN
# =============================================================================
def generate_valid_B(p, n, a, b, max_att=1000000):
    B, seen = [], set()
    for _ in range(max_att):
        if len(B) >= n: break
        x = random.randint(0, p-1)
        if x in seen: continue
        y2 = (x**3 + a*x + b) % p
        if y2 == 0 or pow(y2, (p-1)//2, p) == 1:
            B.append(x); seen.add(x)
    return B

def find_prime_mod1(bits):
    target = 1 << bits
    start = target + (1 - target % 6) % 6
    for p in range(start, start + 500000, 6):
        if pow(2, p-1, p) == 1: return p
    return None

def main():
    random.seed(42)
    a, b = 0, 7
    Path("results").mkdir(exist_ok=True)
    
    p = find_prime_mod1(26)
    if not p: print("❌ Prime not found"); return
    print(f"🔬 Field: p ≈ 2^26 ({p})\n")
    
    results = []
    for n in [6, 8, 10]:
        B = generate_valid_B(p, n, a, b)
        xR = [random.randint(0, p-1) for _ in range(4)]
        
        t0 = time.time()
        M = [[S6_eval(B[i], B[j], xR[0], xR[1], xR[2], xR[3], a, b, p) for j in range(n)] for i in range(n)]
        r = gf_rank(M, p)
        dt = time.time() - t0
        results.append({"poly": "S6", "n": n, "rank": r, "time_s": f"{dt:.3f}"})
        print(f"  |B|={n:2d} → Rank(S6) = {r} ({dt:.3f}s)")
        
    B_ref = generate_valid_B(p, 10, a, b)
    xR3 = random.randint(0, p-1)
    M3 = [[S3_eval(B_ref[i], B_ref[j], xR3, a, b, p) for j in range(10)] for i in range(10)]
    r3 = gf_rank(M3, p)
    results.append({"poly": "S3", "n": 10, "rank": r3, "time_s": "-"})
    print(f"  |B|=10 → Rank(S3) = {r3} (reference)")
    
    csv_path = Path("results/S6_rank_test.csv")
    with open(csv_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=results[0].keys())
        w.writeheader(); w.writerows(results)
    print(f"\n💾 Saved → {csv_path}")

if __name__ == "__main__":
    main()
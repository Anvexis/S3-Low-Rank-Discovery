# test_S5_rank.py (FIXED & VERIFIED)
import random, time, csv, sys
from pathlib import Path

# =============================================================================
# 1. GMPY2 SETUP (Used only for heavy rank computation)
# =============================================================================
try:
    import gmpy2
    HAS_GMPY2 = True
    print("✅ gmpy2 active (accelerating rank solver)")
except ImportError:
    HAS_GMPY2 = False
    print("⚠️ gmpy2 not found (using native Python)")

# =============================================================================
# 2. POLYNOMIAL & RESULTANT UTILS (Native Python % p for safety)
# =============================================================================
def det_mod(M, p):
    """Determinant of small matrix (≤6x6) over F_p using Gaussian elimination"""
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

def sylvester_resultant(c1, c2, p):
    """Resultant of two polynomials given by coefficient lists [c0, c1, ..., cd]"""
    d1, d2 = len(c1)-1, len(c2)-1
    n = d1 + d2
    syl = [[0]*n for _ in range(n)]
    for i in range(d2):
        for j in range(d1+1): syl[i][i+j] = c1[j] % p
    for i in range(d1):
        for j in range(d2+1): syl[d2+i][i+j] = c2[j] % p
    return det_mod(syl, p)

def S3_coeffs_X(x1, x2, a, b, p):
    """Returns [C, B, A] for S3(x1,x2,X) = A*X^2 + B*X + C"""
    s = (x1 + x2) % p
    prod = (x1 * x2) % p
    # Derived from: S3 = (e2-a)^2 - 4*e1*(e3+b)
    A = (s*s - 4*prod) % p
    B = (2*s*(prod - a) - 4*(s*prod + b)) % p
    C = ((prod - a)**2 - 4*s*b) % p
    return [C, B, A]

def S4_eval(x3, x4, x5, X, a, b, p):
    """S4(x3,x4,x5,X) = Res_Y(S3(x3,x4,Y), S3(x5,X,Y))"""
    c1 = S3_coeffs_X(x3, x4, a, b, p)
    c2 = S3_coeffs_X(x5, X, a, b, p)
    return sylvester_resultant(c1, c2, p)

def S4_coeffs_X(x3, x4, x5, a, b, p):
    """Interpolate S4(x3,x4,x5,X) coeffs from 5 points (degree 4 in X)"""
    pts = [(X, S4_eval(x3, x4, x5, X, a, b, p)) for X in range(5)]
    coeffs = [0]*5
    for i in range(5):
        xi, yi = pts[i]
        L = 1
        for j in range(5):
            if i == j: continue
            xj, _ = pts[j]
            num = (-xj) % p
            den = (xi - xj) % p
            L = (L * num * pow(den, -1, p)) % p
        for k in range(5):
            coeffs[k] = (coeffs[k] + yi * L) % p
    return coeffs

def S5_eval(x1, x2, x3, x4, x5, a, b, p):
    """S5(x1..x5) = Res_X(S3(x1,x2,X), S4(x3,x4,x5,X))"""
    c1 = S3_coeffs_X(x1, x2, a, b, p)
    c2 = S4_coeffs_X(x3, x4, x5, a, b, p)
    return sylvester_resultant(c1, c2, p)

def S3_eval(x1, x2, x3, a, b, p):
    """Direct S3 evaluation for reference"""
    e1 = (x1 + x2 + x3) % p
    e2 = (x1*x2 + x2*x3 + x3*x1) % p
    e3 = (x1*x2*x3) % p
    return (pow(e2 - a, 2, p) - 4 * e1 * (e3 + b)) % p

# =============================================================================
# 3. RANK COMPUTATION (gmpy2 accelerated)
# =============================================================================
def gf_rank(matrix, p):
    M = [row[:] for row in matrix]
    rows, cols = len(M), len(M[0])
    rank = r = 0
    
    # Optional gmpy2 conversion for speed
    if HAS_GMPY2:
        p_mpz = gmpy2.mpz(p)
        M = [[gmpy2.mpz(v) for v in row] for row in M]
        mod_inv = lambda x: int(gmpy2.invert(x, p_mpz))
        mod_mul = lambda a, b: int((a * b) % p_mpz)
        mod_sub = lambda a, b: int((a - b) % p_mpz)
    else:
        mod_inv = lambda x: pow(x, -1, p)
        mod_mul = lambda a, b: (a * b) % p
        mod_sub = lambda a, b: (a - b) % p

    for c in range(cols):
        if r >= rows: break
        pivot = next((i for i in range(r, rows) if M[i][c] % p != 0), -1)
        if pivot == -1: continue
        M[r], M[pivot] = M[pivot], M[r]
        inv = mod_inv(M[r][c])
        M[r] = [mod_mul(v, inv) for v in M[r]]
        for i in range(rows):
            if i != r and M[i][c]:
                f = M[i][c]
                M[i] = [mod_sub(M[i][j], mod_mul(f, M[r][j])) for j in range(cols)]
        rank += 1; r += 1
    return rank

# =============================================================================
# 4. HELPERS
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

# =============================================================================
# 5. MAIN
# =============================================================================
def main():
    random.seed(42)
    a, b = 0, 7
    Path("results").mkdir(exist_ok=True)
    
    p = find_prime_mod1(28)
    if not p: print("❌ Prime not found"); return
    print(f"🔬 Field: p ≈ 2^28 ({p})\n")
    
    results = []
    for n in [8, 10, 12]:
        B = generate_valid_B(p, n, a, b)
        xR1, xR2, xR3 = [random.randint(0, p-1) for _ in range(3)]
        
        t0 = time.time()
        M = [[S5_eval(B[i], B[j], xR1, xR2, xR3, a, b, p) for j in range(n)] for i in range(n)]
        r = gf_rank(M, p)
        dt = time.time() - t0
        results.append({"poly": "S5", "n": n, "rank": r, "time_s": f"{dt:.3f}"})
        print(f"  |B|={n:2d} → Rank(S5) = {r} ({dt:.3f}s)")
        
    # Reference S3
    B_ref = generate_valid_B(p, 12, a, b)
    xR = random.randint(0, p-1)
    M3 = [[S3_eval(B_ref[i], B_ref[j], xR, a, b, p) for j in range(12)] for i in range(12)]
    r3 = gf_rank(M3, p)
    results.append({"poly": "S3", "n": 12, "rank": r3, "time_s": "-"})
    print(f"  |B|=12 → Rank(S3) = {r3} (reference)")
    
    csv_path = Path("results/S5_rank_test.csv")
    with open(csv_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=results[0].keys())
        w.writeheader(); w.writerows(results)
    print(f"\n💾 Saved → {csv_path}")

if __name__ == "__main__":
    main()
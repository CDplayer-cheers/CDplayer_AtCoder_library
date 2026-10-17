ll mod(ll a, ll m)
{
    ll ans = a % m;
    if (ans < 0)
        return ans + m;
    else
        return ans;
}
// ax+by=gcd(a,b)をみたす(x,y). O(log(min(a,b)))
pair<ll, ll> ext_euc(ll a, ll b)
{
    ll x1 = 1, y1 = 0, x2 = 0, y2 = 1, d1 = a, d2 = b;
    while (d2 != 0)
    {
        ll t = d1 / d2;
        x1 -= x2 * t;
        swap(x1, x2);
        y1 -= y2 * t;
        swap(y1, y2);
        d1 -= d2 * t;
        swap(d1, d2);
    }
    if (d1 < 0)
    {
        x1 *= -1;
        y1 *= -1;
    }
    return make_pair(x1, y1);
}
ll modinv(ll a, ll m)
{
    return mod(ext_euc(a, m).first, m);
}
template <ll MOD>
struct Modint
{
    ll n;
    Modint(ll input = 0)
    {
        n = mod(input, MOD);
    }

    Modint &operator=(ll input)
    {
        n = mod(input, MOD);
        return *this;
    }

    Modint &operator+=(const Modint &other)
    {
        n = (n + other.n) % MOD;
        return *this;
    }
    Modint &operator-=(const Modint &other)
    {
        n = mod(n - other.n, MOD);
        return *this;
    }
    Modint &operator*=(const Modint &other)
    {
        n = (n * other.n) % MOD;
        return *this;
    }
    Modint &operator/=(const Modint &other)
    {
        n = (n * modinv(other.n, MOD)) % MOD;
        return *this;
    }

    Modint operator-() { return Modint(-n); }

    friend Modint operator+(Modint lhs, const Modint &rhs) { return lhs += rhs; }
    friend Modint operator-(Modint lhs, const Modint &rhs) { return lhs -= rhs; }
    friend Modint operator*(Modint lhs, const Modint &rhs) { return lhs *= rhs; }
    friend Modint operator/(Modint lhs, const Modint &rhs) { return lhs /= rhs; }

    friend bool operator==(const Modint &lhs, const Modint &rhs) { return (lhs.n == rhs.n); }
    friend bool operator!=(const Modint &lhs, const Modint &rhs) { return (lhs.n != rhs.n); }

    friend istream &operator>>(istream &is, Modint &m)
    {
        ll val;
        is >> val;
        m = val;
        return is;
    }
    friend ostream &operator<<(ostream &os, const Modint &m)
    {
        return os << m.n;
    }

    Modint pow(ll p)
    {
        if (p < 0)
            return Modint<MOD>(1) / this->pow(-p);

        Modint ans = 1, base = *this;
        while (p > 0)
        {
            if (p & 1)
                ans *= base;
            base *= base;
            p >>= 1;
        }
        return ans;
    }
};
template <ll MOD>
struct Combination
{
    vector<Modint<MOD>> fact, inv_fact;

    Combination(int max_n) : fact(max_n + 1), inv_fact(max_n + 1)
    {
        fact.at(0) = 1;
        repb(i, 1, max_n + 1)
        {
            fact.at(i) = fact.at(i - 1) * i;
        }
        inv_fact.at(max_n) = Modint<MOD>(1) / fact.at(max_n);
        reps(i, max_n - 1, -1, -1)
        {
            inv_fact.at(i) = inv_fact.at(i + 1) * (i + 1);
        }
    }

    Modint<MOD> c(int n, int r)
    {
        if (r < 0 || r > n)
            return 0;
        return fact.at(n) * inv_fact.at(r) * inv_fact.at(n - r);
    }

    Modint<MOD> p(int n, int r)
    {
        if (r < 0 || r > n)
            return 0;
        return fact.at(n) * inv_fact.at(n - r);
    }

    Modint<MOD> h(int n, int r)
    {
        if (r < 0 || r > n)
            return 0;
        if (n == 0 && r == 0)
            return 1;
        return c(n + r - 1, r);
    }
};
using mint = Modint<998244353LL>;
using mcomb = Combination<998244353LL>;
// nが大きくrが小さいとき.O(r)
mint mcomb_rsmall(ll n, ll r)
{
    mint over = 1;
    mint under = 1;
    repb(i, 1, r + 1)
    {
        over *= n - i + 1;
        under *= i;
    }
    return over / under;
}

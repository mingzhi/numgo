package random

import (
	"math"
	"math/rand"
)

// Rand embeds math/rand.Rand.
// So it has all functions available in math/rand.Rand.
type Rand struct {
	*rand.Rand
}

// New returns a new Rand that uses random values from src
// to generate other random values.
func New(src rand.Source) *Rand {
	s := &Rand{}
	s.Rand = rand.New(src)
	return s
}

// standardExpFloat64 returns a random number
// from a standard exponential distribution with scale = 1.
func (r *Rand) standardExpFloat64() float64 {
	return -math.Log(1.0 - r.Float64())
}

// ExpFloat64 returns a random number of a exponential distribution.
func (r *Rand) ExpFloat64(scale float64) float64 {
	return scale * r.standardExpFloat64()
}

// UniformFloat64 samples a number from a uniform distribution.
func (r *Rand) UniformFloat64(loc, scale float64) float64 {
	return loc + scale*r.Float64()
}

// PoissonInt64 samples a number from a Poisson distribution.
func (r *Rand) PoissonInt64(lam float64) int64 {
	if lam >= 10 {
		return r.poissonPtrs(lam)
	} else if lam == 0 {
		return 0
	} else {
		return r.poissonMult(lam)
	}
}

func (r *Rand) poissonMult(lam float64) int64 {
	var X int64
	U := 0.0
	enlam := math.Exp(-lam)
	prod := 1.0
	for {
		U = r.Float64()
		prod *= U
		if prod > enlam {
			X++
		} else {
			return X
		}
	}
}

func (r *Rand) poissonPtrs(lam float64) int64 {
	var k int64
	var U, V, slam, loglam, a, b, invalpha, vr, us float64
	slam = math.Sqrt(lam)
	loglam = math.Log(lam)
	b = 0.931 + 2.53*slam
	a = -0.059 + 0.02483*b
	invalpha = 1.1239 + 1.1328/(b-3.4)
	vr = 0.9277 - 3.6224/(b-2)
	for {
		U = r.Float64() - 0.5
		V = r.Float64()
		us = 0.5 - math.Abs(U)
		k = int64((2*a/us+b)*U + lam + 0.43)
		if (us >= 0.07) && (V <= vr) {
			return k
		}
		if (k < 0) || ((us < 0.013) && (V > us)) {
			continue
		}
		if (math.Log(V) + math.Log(invalpha) - math.Log(a/(us*us)+b)) <= (-lam + float64(k)*loglam - loggam(float64(k+1))) {
			return k
		}
	}
}

func loggam(x float64) float64 {
	var (
		x0, x2, xp, gl, gl0 float64
		k, n                int64
	)

	a := []float64{8.333333333333333e-02, -2.777777777777778e-03,
		7.936507936507937e-04, -5.952380952380952e-04,
		8.417508417508418e-04, -1.917526917526918e-03,
		6.410256410256410e-03, -2.955065359477124e-02,
		1.796443723688307e-01, -1.39243221690590e+00}

	x0 = x
	n = 0
	if (x == 1.0) || (x == 2.0) {
		return 0.0
	} else if x <= 7.0 {
		n = int64(7 - x)
		x0 = x + float64(n)
	}

	x2 = 1.0 / (x0 * x0)
	xp = 2 * math.Pi
	gl0 = a[9]
	for k = 8; k >= 0; k-- {
		gl0 *= x2
		gl0 += a[k]
	}
	gl = gl0/x0 + 0.5*math.Log(xp) + (x0-0.5)*math.Log(x0) - x0
	if x <= 7.0 {
		for k = 1; k <= n; k++ {
			gl -= math.Log(x0 - 1.0)
			x0 -= 1.0
		}
	}

	return gl
}

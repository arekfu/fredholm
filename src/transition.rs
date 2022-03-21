use crate::options::Params;
use crate::state::State;
use rand::Rng;
use rgsl::lambert_w::{lambert_W0, lambert_Wm1};

pub trait Transition<T>: Send + Sync
where
    T: Rng,
{
    fn sample(&self, from: &State, cutoff: f64, rng: &mut T) -> Option<State>;
}

pub struct PNTransition {
    alpha: f64,
    lambda: f64,
    sigma: f64,
}

impl PNTransition {
    pub fn new(params: &Params) -> PNTransition {
        PNTransition {
            alpha: params.alpha,
            lambda: params.lambda,
            sigma: params.sigma,
        }
    }

    fn displace_positive(&self, xi: f64) -> f64 {
        (-self.alpha + self.sigma - self.alpha * lambert_Wm1((-1.0_f64 + xi) / (1.0_f64.exp())))
            / (self.alpha * self.sigma)
    }

    fn displace_negative(&self, xi: f64, exp: f64) -> f64 {
        (-self.alpha + self.sigma
            - self.alpha
                * lambert_W0(
                    (-(self.alpha * xi) + (-1.0_f64 + xi) * (self.alpha - self.sigma) / exp)
                        / (1.0_f64.exp() * self.alpha),
                ))
            / (self.alpha * self.sigma)
    }
}

impl<T: Rng> Transition<T> for PNTransition {
    /// Sample a transition using the positive-negative decomposition
    ///
    /// The kernel that we are sampling here is
    ///
    ///     f(d) = λ Σ Exp[-Σ d] (α d - 1)
    ///
    /// λ is just treated as an absorption probability, so we discard it for the
    /// rest of the discussion. The rest of the kernel is split into a negative
    /// (0 < d < 1/a) and a positive part (d > 1/a), which we sample separately.
    ///
    /// Mathematica says that the integral of the kernel from 0 to d is
    ///
    ///     F(d) = (α - Σ + Exp[-Σ d] (Σ - α (Σ d + 1))) / Σ
    ///
    /// We will need an implementation of the Lambert W function (provided by
    /// [rgsl]) to sample this by the CDF inversion method.
    ///
    /// The norm of the negative part is
    ///
    ///     N_- = ( Σ - α (1 - Exp[-Σ/α]) ) / Σ
    ///         = 1 - (α/Σ) (1 - Exp[-Σ/α])
    ///         = 1 - (α/Σ) + (α/Σ) Exp[-Σ/α]
    ///
    /// and the norm of the positive part is
    ///
    ///     N_+ = (α/Σ) Exp[-Σ/α]
    ///
    /// so
    ///
    ///     N_- + N_+ = 1 - (α/Σ) + 2 (α/Σ) Exp[-Σ/α]
    ///
    /// The sampling formula for the negative part (0<d<1/α) is
    ///
    ///     d = (Σ - α (1 + W[(- α ξ + (ξ - 1) (α - Σ) Exp[Σ/α] ) / (e α)]))
    ///         / (α Σ)
    ///
    ///     [CForm] (-alpha + sigma - alpha*ProductLog((-(alpha*xi) +
    ///              exp(sigma/alpha)*(-1 + xi)*(alpha -
    ///              sigma))/(e*alpha)))/(alpha*sigma)
    ///
    /// and the sampling formula for the positive part (d>1/α) is
    ///
    ///     d = (Σ - α (1 + W_{-1}[ (ξ-1)/e ])) / (α Σ)
    ///
    ///     [CForm] (-alpha + sigma - alpha*W_{-1}((-1 + xi)/E))/(alpha*sigma)
    ///
    /// Here W is the principal branch of the Lambert W function (ProductLog[x]
    /// in Mathematica) and W1 is the secondary real-valued branch of the
    /// Lambert function (ProductLog[-1,x]). The use of the secondary branch is
    /// important because otherwise the function will not have the right limits
    /// as ξ approaches 0 and 1.
    ///
    fn sample(&self, from: &State, cutoff: f64, rng: &mut T) -> Option<State>
    where
        T: Rng,
    {
        let exp = (-self.sigma / self.alpha).exp();
        let norm_plus = self.alpha * exp / self.sigma;
        let norm_minus = 1.0_f64 - self.alpha / self.sigma + norm_plus;
        debug_assert!(norm_plus >= 0.0);
        debug_assert!(norm_minus >= 0.0);
        let norm_tot = norm_plus + norm_minus;
        let mut weight_multiplier = self.lambda * norm_tot;
        let prob_plus = norm_plus / norm_tot;

        let xi: f64 = rng.gen();
        let xi2: f64 = rng.gen();
        let displacement: f64 = if xi < prob_plus {
            // positive branch was selected
            self.displace_positive(xi2)
        } else {
            // negative branch was selected
            weight_multiplier = -weight_multiplier;
            self.displace_negative(xi2, exp)
        };

        let new_position = from.position + displacement;
        if new_position > cutoff {
            return None;
        }

        let new_weight = from.weight * weight_multiplier;
        Some(State {
            position: new_position,
            weight: new_weight,
        })
    }
}

pub struct ELTransition {
    alpha: f64,
    lambda: f64,
    sigma: f64,
}

impl ELTransition {
    pub fn new(params: &Params) -> ELTransition {
        ELTransition {
            alpha: params.alpha,
            lambda: params.lambda,
            sigma: params.sigma,
        }
    }
}

impl<T: Rng> Transition<T> for ELTransition {
    /// Sample a transition using the exp-linear decomposition
    ///
    /// The kernel that we are sampling here is
    ///
    ///     f(d) = λ Σ Exp[-Σ d] (α d - 1)
    ///
    /// We decompose the kernel into a linear-exponential part and an
    /// exponential part:
    ///
    ///     f(d) = λ (1 + α/Σ) ((1-p) Σ² d Exp[-Σ d] - p Σ Exp[-Σ d])
    ///
    /// The probability p reads
    ///
    ///     p = α / (α + Σ)
    ///
    /// The sampling algorithm is the following:
    ///
    /// - The factor λ (1 + α/Σ) is treated as a weight multiplier.
    /// - With probability p: we sample an exponential displacement from the
    ///   distribution Σ Exp[-Σ d] and we flip the sign of the weight;
    /// - With probability 1-p: we sample *two* exponential displacements from
    ///  the distribution Σ Exp[-Σ d]; the actual particle displacement is the
    ///  sum of the sampled displacements. Here we use the fact that the sum of
    ///  two exponential variates has a linear-exponential distribution.
    ///
    fn sample(&self, from: &State, cutoff: f64, rng: &mut T) -> Option<State>
    where
        T: Rng,
    {
        let norm_tot = 1.0_f64 + self.alpha / self.sigma;
        let mut weight_multiplier = self.lambda * norm_tot;
        let prob_exp = 1.0_f64 / norm_tot;

        let xi: f64 = rng.gen();
        let xi2: f64 = rng.gen();
        let exp2 = -xi2.ln() / self.sigma;
        let displacement: f64 = if xi < prob_exp {
            // exponential branch was selected
            weight_multiplier = -weight_multiplier;
            exp2
        } else {
            // linear-exponential branch was selected
            let xi3: f64 = rng.gen();
            let exp3 = -xi3.ln() / self.sigma;
            exp2 + exp3
        };

        let new_position = from.position + displacement;
        if new_position > cutoff {
            return None;
        }

        let new_weight = from.weight * weight_multiplier;
        Some(State {
            position: new_position,
            weight: new_weight,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::PNTransition;

    const TOLERANCE: f64 = 1e-10;

    #[test]
    fn positive() {
        let transition = PNTransition {
            alpha: 0.375,
            lambda: 0.3435,
            sigma: 1.234,
        };
        let d0 = transition.displace_positive(0.0_f64);
        let d1 = transition.displace_positive(0.99999_f64);
        assert!(
            (d0 - 1.0_f64 / transition.alpha).abs() < TOLERANCE,
            "d0 = {} < {} = 1/α",
            d0,
            1.0_f64 / transition.alpha
        );
        assert!(
            d1 >= 1.0_f64 / transition.alpha,
            "d1 = {} < {} = 1/α",
            d1,
            1.0_f64 / transition.alpha
        );
    }

    #[test]
    fn negative() {
        let transition = PNTransition {
            alpha: 0.375,
            lambda: 0.3435,
            sigma: 1.234,
        };
        let exp = (-transition.sigma / transition.alpha).exp();
        let d0 = transition.displace_negative(0.0_f64, exp);
        let d1 = transition.displace_negative(1.0_f64, exp);
        assert!(d0.abs() < TOLERANCE, "d0 = {}", d0);
        assert!(
            (d1 - 1.0_f64 / transition.alpha).abs() < TOLERANCE,
            "d1 = {} > {} = 1/α",
            d1,
            1.0_f64 / transition.alpha
        );
    }
}

#include "mc.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace pp {

namespace {

double box_volume(
    const vec& lower,
    const vec& upper
)
{
    if (lower.size() != upper.size()) {
        throw std::invalid_argument(
            "Monte Carlo: lower and upper must have the same dimension"
        );
    }

    if (lower.size() == 0) {
        throw std::invalid_argument(
            "Monte Carlo: integration region must have at least one dimension"
        );
    }

    double volume = 1.0;

    for (int i = 0; i < lower.size(); ++i) {
        if (upper[i] <= lower[i]) {
            throw std::invalid_argument(
                "Monte Carlo: upper[i] must be larger than lower[i]"
            );
        }

        volume *= upper[i] - lower[i];
    }

    return volume;
}


double sample_variance(
    double sum,
    double sum_squared,
    std::size_t count
)
{
    if (count < 2) {
        return 0.0;
    }

    const double n =
        static_cast<double>(count);

    const double mean =
        sum / n;

    const double mean_squared =
        sum_squared / n;

    return std::max(
        0.0,
        mean_squared - mean * mean
    );
}


mc_result stratified_recursive(
    const mc_function& f,
    const vec& lower,
    const vec& upper,
    std::size_t N,
    const random_function& random_double,
    std::size_t nmin
)
{

    if (N <= nmin || N - nmin < 2) {
        return plainmc(
            f,
            lower,
            upper,
            N,
            random_double
        );
    }

    const int dim = lower.size();

    vec midpoint(dim);
    vec x(dim);

    for (int i = 0; i < dim; ++i) {
        midpoint[i] =
            0.5 * (lower[i] + upper[i]);
    }


    std::vector<double> sum_left(dim, 0.0);
    std::vector<double> sum_right(dim, 0.0);

    std::vector<double> sum2_left(dim, 0.0);
    std::vector<double> sum2_right(dim, 0.0);

    std::vector<std::size_t> count_left(dim, 0);
    std::vector<std::size_t> count_right(dim, 0);

    for (std::size_t k = 0; k < nmin; ++k) {

        for (int i = 0; i < dim; ++i) {
            const double u = random_double();

            if (!std::isfinite(u) || u < 0.0 || u >= 1.0) {
                throw std::runtime_error(
                    "stratifiedmc: random generator must return numbers in [0,1)"
                );
            }

            x[i] =
                lower[i]
                + u * (upper[i] - lower[i]);
        }

        const double fx = f(x);

        if (!std::isfinite(fx)) {
            throw std::runtime_error(
                "stratifiedmc: integrand returned a non-finite value"
            );
        }

        for (int i = 0; i < dim; ++i) {

            if (x[i] < midpoint[i]) {
                ++count_left[i];

                sum_left[i] += fx;
                sum2_left[i] += fx * fx;
            }
            else {
                ++count_right[i];

                sum_right[i] += fx;
                sum2_right[i] += fx * fx;
            }
        }
    }

    std::vector<double> variance_left(dim, 0.0);
    std::vector<double> variance_right(dim, 0.0);

    int split_dimension = 0;
    double largest_subvariance = -1.0;


    for (int i = 0; i < dim; ++i) {

        variance_left[i] =
            sample_variance(
                sum_left[i],
                sum2_left[i],
                count_left[i]
            );

        variance_right[i] =
            sample_variance(
                sum_right[i],
                sum2_right[i],
                count_right[i]
            );


        if (
            count_left[i] < 2
            || count_right[i] < 2
        ) {
            continue;
        }

        const double combined_subvariance =
            variance_left[i]
            + variance_right[i];

        if (combined_subvariance > largest_subvariance) {
            largest_subvariance =
                combined_subvariance;

            split_dimension = i;
        }
    }


    vec left_upper = upper;
    vec right_lower = lower;

    left_upper[split_dimension] =
        midpoint[split_dimension];

    right_lower[split_dimension] =
        midpoint[split_dimension];

    const std::size_t remaining_points =
        N - nmin;

    const double left_weight =
        variance_left[split_dimension];

    const double right_weight =
        variance_right[split_dimension];

    std::size_t N_left;

    if (left_weight + right_weight > 0.0) {
        const double left_fraction =
            left_weight
            / (left_weight + right_weight);

        N_left = static_cast<std::size_t>(
            std::llround(
                left_fraction
                * static_cast<double>(remaining_points)
            )
        );
    }
    else {

        N_left =
            remaining_points / 2;
    }


    N_left = std::max<std::size_t>(
        1,
        std::min(
            remaining_points - 1,
            N_left
        )
    );

    const std::size_t N_right =
        remaining_points - N_left;

    const mc_result left_result =
        stratified_recursive(
            f,
            lower,
            left_upper,
            N_left,
            random_double,
            nmin
        );

    const mc_result right_result =
        stratified_recursive(
            f,
            right_lower,
            upper,
            N_right,
            random_double,
            nmin
        );


    return {
        left_result.integral
            + right_result.integral,

        std::hypot(
            left_result.error,
            right_result.error
        )
    };
}

} // namespace


double lcg::next()
{
    if (m_ <= 1) {
        throw std::runtime_error(
            "lcg: modulus must be larger than one"
        );
    }

    seed_ =
        (a_ * seed_ + c_) % m_;


    return
        (static_cast<double>(seed_) + 1.0)
        /
        (static_cast<double>(m_) + 1.0);
}


double corput(
    std::size_t n,
    unsigned base
)
{
    if (base < 2) {
        throw std::invalid_argument(
            "corput: base must be at least two"
        );
    }

    double q = 0.0;
    double base_factor =
        1.0 / static_cast<double>(base);

    while (n > 0) {
        const std::size_t digit =
            n % base;

        q +=
            static_cast<double>(digit)
            * base_factor;

        n /= base;
        base_factor /=
            static_cast<double>(base);
    }

    return q;
}


std::vector<unsigned> prime_numbers(
    std::size_t count
)
{
    std::vector<unsigned> primes;
    primes.reserve(count);

    unsigned candidate = 2;

    while (primes.size() < count) {

        bool candidate_is_prime = true;

        for (const unsigned prime : primes) {

  
            if (prime > candidate / prime) {
                break;
            }

            if (candidate % prime == 0) {
                candidate_is_prime = false;
                break;
            }
        }

        if (candidate_is_prime) {
            primes.push_back(candidate);
        }

        ++candidate;
    }

    return primes;
}


Halton::Halton(
    std::size_t dimension,
    std::size_t prime_offset
)
{
    if (dimension == 0) {
        throw std::invalid_argument(
            "Halton: dimension must be positive"
        );
    }

    const std::vector<unsigned> primes =
        prime_numbers(
            dimension + prime_offset
        );

    bases_.assign(
        primes.begin()
            + static_cast<std::ptrdiff_t>(prime_offset),

        primes.begin()
            + static_cast<std::ptrdiff_t>(
                prime_offset + dimension
            )
    );
}


vec Halton::x(std::size_t n) const
{
    const int dim =
        static_cast<int>(bases_.size());

    vec point(dim);

    for (int i = 0; i < dim; ++i) {
        point[i] =
            corput(
                n,
                bases_[static_cast<std::size_t>(i)]
            );
    }

    return point;
}


mc_result plainmc(
    const mc_function& f,
    const vec& lower,
    const vec& upper,
    std::size_t N,
    const random_function& random_double
)
{
    if (N == 0) {
        throw std::invalid_argument(
            "plainmc: N must be positive"
        );
    }

    const double volume =
        box_volume(lower, upper);

    const int dim =
        lower.size();

    vec x(dim);

    double sum = 0.0;
    double sum_squared = 0.0;

    for (std::size_t k = 0; k < N; ++k) {

        for (int i = 0; i < dim; ++i) {
            const double u =
                random_double();

            if (!std::isfinite(u) || u < 0.0 || u >= 1.0) {
                throw std::runtime_error(
                    "plainmc: random generator must return numbers in [0,1)"
                );
            }

            x[i] =
                lower[i]
                + u * (upper[i] - lower[i]);
        }

        const double fx = f(x);

        if (!std::isfinite(fx)) {
            throw std::runtime_error(
                "plainmc: integrand returned a non-finite value"
            );
        }

        sum += fx;
        sum_squared += fx * fx;
    }

    const double number_of_points =
        static_cast<double>(N);

    const double mean =
        sum / number_of_points;

    const double variance =
        sample_variance(
            sum,
            sum_squared,
            N
        );

    const double integral =
        volume * mean;

    const double error =
        volume
        * std::sqrt(
            variance / number_of_points
        );

    return {
        integral,
        error
    };
}


mc_result quasimc(
    const mc_function& f,
    const vec& lower,
    const vec& upper,
    std::size_t points_per_sequence,
    const Halton& sequence1,
    const Halton& sequence2
)
{
    if (points_per_sequence == 0) {
        throw std::invalid_argument(
            "quasimc: points_per_sequence must be positive"
        );
    }

    const double volume =
        box_volume(lower, upper);

    const int dim =
        lower.size();

    if (
        sequence1.dimension()
            != static_cast<std::size_t>(dim)
        ||
        sequence2.dimension()
            != static_cast<std::size_t>(dim)
    ) {
        throw std::invalid_argument(
            "quasimc: Halton dimensions must match the integration dimension"
        );
    }

    vec point1(dim);
    vec point2(dim);

    double sum1 = 0.0;
    double sum2 = 0.0;

    for (
        std::size_t k = 0;
        k < points_per_sequence;
        ++k
    ) {
    
        const vec unit_point1 =
            sequence1.x(k + 1);

        const vec unit_point2 =
            sequence2.x(k + 1);

        for (int i = 0; i < dim; ++i) {
            point1[i] =
                lower[i]
                + unit_point1[i]
                * (upper[i] - lower[i]);

            point2[i] =
                lower[i]
                + unit_point2[i]
                * (upper[i] - lower[i]);
        }

        const double f1 =
            f(point1);

        const double f2 =
            f(point2);

        if (
            !std::isfinite(f1)
            || !std::isfinite(f2)
        ) {
            throw std::runtime_error(
                "quasimc: integrand returned a non-finite value"
            );
        }

        sum1 += f1;
        sum2 += f2;
    }

    const double number_of_points =
        static_cast<double>(
            points_per_sequence
        );

    const double integral1 =
        volume * sum1 / number_of_points;

    const double integral2 =
        volume * sum2 / number_of_points;


    const double integral =
        0.5 * (integral1 + integral2);

    const double error =
        std::abs(integral1 - integral2);

    return {
        integral,
        error
    };
}


mc_result stratifiedmc(
    const mc_function& f,
    const vec& lower,
    const vec& upper,
    std::size_t N,
    const random_function& random_double,
    std::size_t nmin
)
{
    if (N == 0) {
        throw std::invalid_argument(
            "stratifiedmc: N must be positive"
        );
    }

    if (nmin < 4) {
        throw std::invalid_argument(
            "stratifiedmc: nmin must be at least four"
        );
    }

    box_volume(lower, upper);

    return stratified_recursive(
        f,
        lower,
        upper,
        N,
        random_double,
        nmin
    );
}

} // namespace pp
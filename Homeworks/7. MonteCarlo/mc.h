#pragma once

#include "vector.h"

#include <cstddef>
#include <cstdint>
#include <functional>
#include <vector>

namespace pp {


struct mc_result {
    double integral;
    double error;
};



using mc_function =
    std::function<double(const vec&)>;


using random_function =
    std::function<double()>;



class lcg {
private:
    std::uint64_t seed_;
    std::uint64_t a_;
    std::uint64_t c_;
    std::uint64_t m_;

public:
    explicit lcg(
        std::uint64_t seed = 0,
        std::uint64_t a = 1664525ULL,
        std::uint64_t c = 1013904223ULL,
        std::uint64_t m = (1ULL << 32)
    )
        : seed_(seed),
          a_(a),
          c_(c),
          m_(m)
    {
    }

    double next();
};



double corput(
    std::size_t n,
    unsigned base
);



std::vector<unsigned> prime_numbers(
    std::size_t count
);



class Halton {
private:
    std::vector<unsigned> bases_;

public:
    explicit Halton(
        std::size_t dimension,
        std::size_t prime_offset = 0
    );

    vec x(std::size_t n) const;

    std::size_t dimension() const
    {
        return bases_.size();
    }
};



mc_result plainmc(
    const mc_function& f,
    const vec& lower,
    const vec& upper,
    std::size_t N,
    const random_function& random_double
);



mc_result quasimc(
    const mc_function& f,
    const vec& lower,
    const vec& upper,
    std::size_t points_per_sequence,
    const Halton& sequence1,
    const Halton& sequence2
);



mc_result stratifiedmc(
    const mc_function& f,
    const vec& lower,
    const vec& upper,
    std::size_t N,
    const random_function& random_double,
    std::size_t nmin = 64
);


} // namespace pp
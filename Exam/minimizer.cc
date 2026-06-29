#include "minimizer.h"

#include <algorithm>
#include <stdexcept>
#include <utility>
#include <vector>

// Debugging mainly helped with by ChatGPT, entire minimizer here inspired by homework 9.

namespace pp {

    namespace {

        void sort_simple(
            std::vector<vec>& points,
            std::vector<double>& values
        ) {
            const int count =
                static_cast<int>(points.size());

            for (int i = 0; i < count; ++i) {
                for (int j = i + 1; j < count; ++j) {
                    if (values[j] < values[i]) {
                        std::swap(values[i], values[j]);
                        std::swap(points[i], points[j]);
                    }
                }
            }
        }


        double simple_size(
            const std::vector<vec>& points
        ) {
            double largest_dist = 0.0;

            for (
                int i = 1;
                i < static_cast<int>(points.size());
                ++i
            ) {
                const double dist =
                    norm(points[i] - points[0]);

                largest_dist =
                    std::max(largest_dist, dist);
            }

            return largest_dist;
        }

    } // anonymous namespace


    simple_result down_simple(
        const std::function<double(const vec&)>& f,
        const vec& start,
        double step,
        double tol,
        int max_iterations
    ) {
        const int n = start.size();

        if (n < 1) {
            throw std::invalid_argument(
                "Simplex minimizer needs at least one variable"
            );
        }

        if (step <= 0.0) {
            throw std::invalid_argument(
                "Simplex step must be positive"
            );
        }

        if (tol <= 0.0) {
            throw std::invalid_argument(
                "Simplex tolerance must be positive"
            );
        }


        std::vector<vec> points;
        points.reserve(n + 1);

        points.push_back(start);


        for (int i = 0; i < n; ++i) {
            vec point = start;
            point[i] += step;
            points.push_back(point);
        }


        std::vector<double> values(n + 1);

        for (int i = 0; i < n + 1; ++i) {
            values[i] = f(points[i]);
        }

        const double reflection_factor = 1.0;
        const double expansion_factor = 2.0;
        const double contraction_factor = 0.5;
        const double reduction_factor = 0.5;

        int iteration = 0;
        bool converged = false;

        for (
            iteration = 0;
            iteration < max_iterations;
            ++iteration
        ) {

            sort_simple(points, values);

            if (simple_size(points) < tol) {
                converged = true;
                break;
            }


            vec centroid(n);

            for (int i = 0; i < n; ++i) {
                centroid += points[i];
            }

            centroid /= static_cast<double>(n);

  
            const vec reflected =
                centroid
                + reflection_factor
                * (centroid - points[n]);

            const double reflected_value =
                f(reflected);

            if (reflected_value < values[0]) {
                const vec expanded =
                    centroid
                    + expansion_factor
                    * (reflected - centroid);

                const double expanded_value =
                    f(expanded);

                if (expanded_value < reflected_value) {
                    points[n] = expanded;
                    values[n] = expanded_value;
                }
                else {
                    points[n] = reflected;
                    values[n] = reflected_value;
                }

                continue;
            }

            if (reflected_value < values[n - 1]) {
                points[n] = reflected;
                values[n] = reflected_value;
                continue;
            }


            vec contracted(n);
            double comparison_value = values[n];

            if (reflected_value < values[n]) {

                contracted =
                    centroid
                    + contraction_factor
                    * (reflected - centroid);

                comparison_value =
                    reflected_value;
            }
            else {

                contracted =
                    centroid
                    + contraction_factor
                    * (points[n] - centroid);
            }

            const double contracted_value =
                f(contracted);

            if (contracted_value < comparison_value) {
                points[n] = contracted;
                values[n] = contracted_value;
                continue;
            }

    
            for (int i = 1; i < n + 1; ++i) {
                points[i] =
                    points[0]
                    + reduction_factor
                    * (points[i] - points[0]);

                values[i] =
                    f(points[i]);
            }
        }


        sort_simple(points, values);

        return {
            points[0],
            values[0],
            iteration,
            converged
        };
    }

} // namespace pp
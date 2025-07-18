/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2024 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This file is part of Psi4.
 *
 * Psi4 is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * Psi4 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along
 * with Psi4; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

 /**
 *  Borrowed RangeStep implementation from @kb1rd
 *  lcom implementation.
 */

#include "rangestep.h"

#include <tuple>
#include <vector>
#include <algorithm>
#include "psi4/libpsi4util/exception.h"

namespace psi_lcom {

RangeStep::RangeStep(): global(0.0), coefs() {}
RangeStep::RangeStep(double global, const std::vector<RangeCutCoef>& coefs): global(global), coefs(coefs) {
    this->prune();
}
RangeStep::RangeStep(const RangeStep& other): global(other.global), coefs(other.coefs) {}

RangeStep RangeStep::fromRangeValues(double local, const std::vector<RangeCutCoef>& valuesAt) {
    RangeStep this_;
    this_.global = local;
    this_.coefs = std::vector<RangeCutCoef>(valuesAt.size());

    double coef = local;
    for(size_t i = 0; i < valuesAt.size(); i++) {
        coef = std::get<1>(valuesAt[i]) - coef;
        this_.coefs[i] = std::make_tuple(std::get<0>(valuesAt[i]), coef);
    }
    this_.prune();
    return this_;
}
RangeStep RangeStep::fromAlphaOmegaBeta(double alpha, double omega, double beta) {
    RangeStep this_;
    this_.global = alpha;
    this_.coefs = std::vector<RangeCutCoef>(
        omega != 0.0 && beta != 0.0 ? 1 : 0,
        std::make_tuple(omega, beta)
    );
    return this_;
}
RangeStep RangeStep::sum(std::vector<RangeStep> rangesteps) {
    RangeStep this_;
    size_t n_total;
    for (const auto& rc : rangesteps) {
        this_.global += rc.global;
        n_total += rc.coefs.size();
    }
    this_.coefs.reserve(n_total); // Pre-reserve the total size required; This is faster
    for (const auto& rc : rangesteps) {
        this_.coefs.insert(this_.coefs.end(), rc.coefs.begin(), rc.coefs.end());
    }
    this_.prune(); // Now, sort & remove duplicates
    return this_;
}
RangeStep RangeStep::scaled(RangeStep other, double scalef) {
    if (scalef == 0) { return RangeStep(); }

    RangeStep this_(other);
    this_.global *= scalef;
    for(size_t i = 0; i < this_.coefs.size(); i++) {
        this_.coefs[i] = std::make_tuple(
            std::get<0>(this_.coefs[i]),
            std::get<1>(this_.coefs[i]) * scalef
        );
    }
    return this_;
}

void RangeStep::prune() {
    // Must be ordered for this algorithm to work
    this->ensureOrdered();

    double omega_last = 0.0;
    for(size_t i = 0; i < this->coefs.size(); i++) {
        if (std::get<0>(this->coefs[i]) == omega_last) { // Duplicate omegas
            if (omega_last == 0.0) {
                this->global += std::get<1>(this->coefs[i]);
            } else {
                double newcoef = std::get<1>(this->coefs[i-1]);
                newcoef += std::get<1>(this->coefs[i]);
                this->coefs[i-1] = std::make_tuple(std::get<0>(this->coefs[i-1]), newcoef);
            }
            this->coefs.erase(this->coefs.begin() + i);
            i--;
        } else if (std::get<1>(this->coefs[i]) == 0) { // Zero betas
            this->coefs.erase(this->coefs.begin() + i);
            i--;
        }
        omega_last = std::get<0>(this->coefs[i+1]);
    }

    this->coefs.shrink_to_fit(); // Free up unused space
}
void RangeStep::ensureOrdered() {
    std::sort(this->coefs.begin(), this->coefs.end());
}

std::tuple<double, double, double> RangeStep::toAlphaOmegaBeta() const {
    double alpha, beta, omega;

    alpha = this->global;
    if (this->n_w() > 1) {
        throw psi::PSIEXCEPTION(
            "Unable to convert hybrid range step into a simple alpha, omega, beta. "
            "You're probably running code that doesn't yet support multi-RSH functionals with such a functional"
        );
    } else if (this->n_w() > 0) {
        omega = std::get<0>(this->coefficient(0));
        beta = std::get<1>(this->coefficient(0));
    } else {
        omega = beta = 0;
    }
    return std::make_tuple(alpha, omega, beta);
}
RangeCutCoef RangeStep::value(size_t n) const {
    double coef = this->global;
    for (size_t i = 0; i < n; i++) {
        coef += std::get<1>(this->coefficient(i));
    }
    return std::make_tuple(std::get<0>(this->coefficient(n)), coef);
}

} // namespace psi
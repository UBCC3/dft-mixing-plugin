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

#ifndef RANGESTEP_H
#define RANGESTEP_H

/**
 *  Borrowed RangeStep implementation from @kb1rd
 *  lcom implementation.
 */

#include <tuple>
#include <vector>
#include <cstddef>

namespace psi_lcom {

/** Tuple of (omega, beta) */
typedef std::tuple<double, double> RangeCutCoef;

/**
 * RangeStep: A utility class for handling range separation.
 *
 * Originally designed to replace the Hartree-Fock exchange "alpha", "omega", and "beta" coefficients, this class is
 * designed to provide coefficients for a series of range seperated calculations (i,e, calculations that switch on at a
 * particular range, "omega").
 *
 * ^
 * | *------------------------------------*
 * | |         GLOBAL TERM                | <- Stored in `global`
 * | *-------------------------------------
 * |         *----------------------------*
 * |         |          b_1 TERM          | <- Stored in `coefs[0]`
 * |         *----------------------------*
 * |               *----------------------*
 * |               |       b_2 TERM       | <- Stored in `coefs[1]`
 * |               *----------------------*
 * |                         *------------*
 * |                         |  b_n TERM  | <- Stored in `coefs[n-1]`
 * |                         *------------*
 * |
 * *---------------------------------------> r
 *   0      w_1   w_2 ....  w_n
 *
 * There are two types of numbers you could get at a particular cutoff: The coefficient and the value. The coefficient
 * is the `b_n` in the diagram above for a single calculation (a single block above). The value is the sum of
 * coefficients at that point. So, it would be the sum of blocks at a particular r value in the diagram.
 *
 * One final note: In reality, the range separation cutoffs are implemented with error functions, so the above graphic
 * is just for illustration since it's a bit easier to think about blocks than discrete functions ;)
 *
 */
class RangeStep {
protected:
    double global;
    std::vector<RangeCutCoef> coefs;
public:
    RangeStep();
    RangeStep(double global, const std::vector<RangeCutCoef>& coefs);
    RangeStep(const RangeStep& other);

    /** Create a RangeStep from the desired values at certain cutoffs. */
    static RangeStep fromRangeValues(double local, const std::vector<RangeCutCoef>& valuesAt);
    /** Create a RangeStep from the legacy alpha, omega, and beta coefficients. */
    static RangeStep fromAlphaOmegaBeta(double alpha, double omega, double beta);
    /** Create a RangeStep by summing other RangeSteps. */
    static RangeStep sum(std::vector<RangeStep> rangesteps);
    /** Create a RangeStep by scaling another RangeStep. */
    static RangeStep scaled(RangeStep other, double scalef);

    /**
     * For the RangeStep to work, the coefficients must be in order. This method ensures this is the case.
     * You shouldn't need to call this manually.
     */
    void ensureOrdered();
    /**
     * Simplify this RangeStep if possible by merging and eliminating terms.
     * You shouldn't need to call this manually.
     */
    void prune();

    inline double globalCoefficient() const { return this->global; }
    inline RangeCutCoef coefficient(size_t n) const { return this->coefs.at(n); }
    inline const std::vector<RangeCutCoef>& coefficients() const { return this->coefs; }
    /** Get the number of cutoffs. */
    inline size_t n_w() const { return coefs.size(); }
    /** Convert this RangeStep to alpha, omega, beta coefs if possible. Throws if not possible. */
    std::tuple<double, double, double> toAlphaOmegaBeta() const;

    /** True if this RangeStep has a global term. */
    inline bool has_global() const { return this->global != 0.0; }
    /** True if this RangeStep has ranged terms. */
    inline bool is_lrc() const { return coefs.size() > 0; }
    /** True if this RangeStep has any nonzero terms. */
    inline bool is_nonzero() const { return this->has_global() || this->is_lrc(); }

    /** Get the value (not coef) for this RangeStep at the local level. */
    inline double localValue() const { return this->global; } // Only the global term influences the coef at the local level
    /** Get the value (not coef) for this RangeStep at a particular cutoff indexed by n. */
    RangeCutCoef value(size_t n) const;
};

} // namespace psi

#endif
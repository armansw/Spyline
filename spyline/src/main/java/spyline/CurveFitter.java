package spyline;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;

import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.CholeskyDecomposition;
import org.apache.commons.math3.linear.MatrixUtils;


public class CurveFitter {
    CurveFitter(Spline spline) {
        this.mDelta = 0.0;
        this.p = 0.0;
        this.mError = 0.0;
        this.mPenalty = this.penalty(spline);
    }

    public double penalty(Spline spline) {
        ArrayList<Double> knots = spline.getKnots();
        int g = spline.getInternalKnotsNum();
        this.mPenalty = 0.0;
        for (int i = 0; i < g + 1; ++i) {
            this.mPenalty += 1.0 / (knots.get(i + 1) - knots.get(i));
        }
        this.mError = this.mDelta + this.p * this.mPenalty;
        return this.mPenalty;
    }

    public double delta(Spline spline, Points points, double smoothingWeight) {
        this.mDelta = 0;
        for (int i = 0; i < points.len(); ++i) {
            double e = points.getItem(i).getW() * (points.getItem(i).getY() - spline.getValue(points.getItem(i).getX()));
            this.mDelta += e * e;
        }
        double nu = 0.0;
        if (smoothingWeight > 0) {
            int k = spline.getDegree();
            int g = spline.getInternalKnotsNum();
            ArrayList<Double> coefficients = spline.getCoefficients();
            for (int q = k + 1; q < g + k + 1; ++q) {
                double e = 0.0;
                for (int i = q - k - 1; i < q + 1; ++i) {
                    e += coefficients.get(i) * spline.getLeadDerivativeDifference(i, q);
                    nu += e * e;
                }
            }
            nu *= smoothingWeight;
        }
        this.mDelta += nu;
        this.mError = this.mDelta + this.p * this.mPenalty;
        return this.mDelta;
    }

    // Checked twice!

    public static double penaltyDerivative(Spline spline, int knotId) {
        ArrayList<Double> knots = spline.getKnots();
        double a = knots.get(knotId - 1);
        double b = knots.get(knotId);
        double c = knots.get(knotId + 1);
        double bma = b - a;
        double cmb = c - b;
        return 1.0 / (cmb * cmb) - 1.0 / (bma * bma);
    }

    // Checked!
    public double error(Spline spline, Points points, double sw) {
        return this.delta(spline, points, sw) + this.p * this.penalty(spline);
    }

    // Checked twice!
    public double errorGradient(Spline spline, Points points, double smoothingWeight, int knotId) {
        double gradError = 0.0;
        for (int i = 0; i < points.len(); ++i) {
            double wSq = points.getItem(i).getW() * points.getItem(i).getW();
            double diff = points.getItem(i).getY() - spline.getValue(points.getItem(i).getX());
            gradError -= spline.getValueDerivativeKnot(points.getItem(i).getX(), knotId + spline.getDegree()) * wSq * diff;
        }

        if (this.p > 0) {
            gradError += 0.5 * this.p * CurveFitter.penaltyDerivative(spline, knotId);
        }

        if (smoothingWeight > 0) {
            double smError = 0.0;
            int k = spline.getDegree();
            int g = spline.getInternalKnotsNum();
            ArrayList<Double> coefficients = spline.getCoefficients();
            for (int q = k + 1; q < g + k + 1; ++q) {
                double sum1 = 0.0;
                double sum2 = 0.0;
                for (int i = q - k - 1; i < q + 1; ++i) {
                    double ci = coefficients.get(i);
                    double leadDerDiff = spline.getLeadDerivativeDifference(i, q);
                    sum1 += ci * leadDerDiff;
                    sum2 += ci * spline.getLeadDerDiffDerKnot(leadDerDiff, i, q, knotId + spline.getDegree());
                }
                smError += sum1 * sum2;
            }
            gradError += smoothingWeight * smError;
        }

        return 2 * gradError;
    }

    // Checked twice!
    public double theta(Spline spline, Points points, double sw, double alpha, ArrayList<Double> direction,
            ArrayList<Double> fixedKnots) {
        int g = spline.getInternalKnotsNum();
        ArrayList<Double> knots = new ArrayList<>(Collections.nCopies(g + 2, 0.0));
        knots.set(0, spline.getLeftBound());
        knots.set(g + 1, spline.getRightBound());
        for (int i = 0; i < g; ++i) {
            knots.set(i + 1, fixedKnots.get(i + 1) + alpha * direction.get(i));
        }
        spline.setKnots(knots);

        if (CurveFitter.approximate(spline, points, sw)) {
            return this.error(spline, points, sw);
        }
        return -1.0;
    }

    // Checked!
    public static double norm(ArrayList<Double> v) {
        double sum = 0;
        for (Double x : v) {
            sum += x.doubleValue() * x.doubleValue();
        }
        return sum;
    }

    // Checked!
    public static boolean approximate(Spline spline, Points points, double smoothingWeight) {
        int k = spline.getDegree();
        int g = spline.getInternalKnotsNum();
        int n = points.len();
        ArrayList<Double> coefficients = new ArrayList<>(Collections.nCopies(g + k + 1, 0.0));
        ArrayList<ArrayList<Double>> A = new ArrayList<ArrayList<Double>>(g + k + 1);
        for (int i = 0; i < g + k + 1; ++i) {
            A.add(new ArrayList<Double>(Collections.nCopies(g + k + 1, 0.0)));
        }

        int l = 0;
        for (int r = 0; r < n; ++r) {
            l = spline.getLeftNodeIndex(points.getItem(r).getX(), l);
            if (l < 0) {
                return false;
            }
            ArrayList<Double> bSplines = spline.bSplines(points.getItem(r).getX(), k);
            for (int i = 0; i < k + 1; ++i) {
                double wSq = points.getItem(r).getW() * points.getItem(r).getW();
                for (int j = 0; j < i + 1; ++j) {
                    A.get(i + l - k).set(j + l - k,
                            A.get(i + l - k).get(j + l - k) + wSq * bSplines.get(i) * bSplines.get(j));
                }
                coefficients.set(i + l - k, coefficients.get(i + l - k) + wSq * points.getItem(r).getY() * bSplines.get(i));
            }
        }

        if (smoothingWeight > 0) {
            for (int q = 0; q < g; ++q) {
                for (int i = q; i < q + k + 2; ++i) {
                    double ai = spline.getLeadDerivativeDifference(i, q + k + 1);
                    for (int j = q; j < i + 1; ++j) {
                        A.get(i).set(j, A.get(i).get(j)
                                + smoothingWeight * ai * spline.getLeadDerivativeDifference(j, q + k + 1));
                    }
                }
            }
        }

        for (int i = 0; i < g + k + 1; ++i) {
            for (int j = 0; j < i; ++j) {
                A.get(j).set(i, A.get(i).get(j));
            }
        }

        double[][] ei = new double[A.size()][A.size()];
        for (int i = 0; i < A.size(); ++i) {
            for (int j = 0; j < A.size(); ++j) {
                ei[i][j] = A.get(i).get(j);
            }
        }
        RealMatrix AR = MatrixUtils.createRealMatrix(ei);
        // RealMatrix AR = new Array2DRowRealMatrix();
        // AR.createMatrix(A.size(), A.size());
        // AR.setSubMatrix(ei, 0, 0);
        
        CholeskyDecomposition decomposition = new CholeskyDecomposition(AR);
        RealMatrix Lu = decomposition.getL();
        double[][] Lb = Lu.getData();
        ArrayList<ArrayList<Double>> L = new ArrayList<ArrayList<Double>>(A.size());
        for (int i = 0; i < A.size(); ++i) {
            L.add(new ArrayList<Double>(Collections.nCopies(A.size(), 0.0)));
        }
        for (int i = 0; i < A.size(); ++i) {
            for (int j = 0; j < A.size(); ++j) {
                L.get(i).set(j, Lb[i][j]);
            }
        }

        for (int i = 0; i < g + k + 1; ++i) {
            for (int j = 0; j < i; ++j) {
                coefficients.set(i, coefficients.get(i) - L.get(i).get(j) * coefficients.get(j));
            }
            coefficients.set(i, coefficients.get(i) / L.get(i).get(i));
        }

        for (int i = g + k; i >= 0; --i) {
            for (int j = g + k; j >= i + 1; --j) {
                coefficients.set(i, coefficients.get(i) - L.get(j).get(i) * coefficients.get(j));
            }
            coefficients.set(i, coefficients.get(i) / L.get(i).get(i));
        }

        spline.setCoefficients(coefficients);

        return true;
    }

    // Checked twice!
    public boolean specDimensionalMinimization(Spline spline, Points points, double sw,
            ArrayList<Double> direction, ArrayList<Double> errorDerivative, ArrayList<Double> fixedKnots) {
        int g = spline.getInternalKnotsNum();
        ArrayList<Double> knots = spline.getKnots();
        double alphaMax = Double.POSITIVE_INFINITY;
        double a = spline.getLeftBound();
        double b = spline.getRightBound();

        if (direction.get(0) < 0) {
            alphaMax = (a - knots.get(1)) / direction.get(0);
        }
        for (int i = 0; i < g - 1; ++i) {
            if (direction.get(i) > direction.get(i + 1)) {
                alphaMax = Math.min(alphaMax,
                        (knots.get(i + 2) - knots.get(i + 1)) / (direction.get(i) - direction.get(i + 1)));
            }
        }
        if (direction.get(g - 1) > 0) {
            alphaMax = Math.min(alphaMax, (b - knots.get(g)) / direction.get(g - 1));
        }

        double theta0 = this.mError;
        double theta0Der = 0.0;
        Iterator<Double> iter1 = direction.iterator();
        Iterator<Double> iter2 = errorDerivative.iterator();
        while (iter1.hasNext() && iter2.hasNext()) {
            theta0Der += iter1.next() * iter2.next();
        }
        double alpha0 = 0.0;
        double alpha2 = alphaMax / (1 - theta0 / alphaMax / theta0Der);
        double alpha1 = 0.5 * alpha2;
        double q0 = this.mDelta;
        double r0 = this.mPenalty;
        double theta1 = this.theta(spline, points, sw, alpha1, direction, fixedKnots);
        if (theta1 < 0) {
            return false;
        }
        double q1 = this.mDelta;
        double r1 = this.mPenalty;

        int iteration = 0;
        int maxNumOfIterations = 10;
        while (theta1 >= theta0 && iteration < maxNumOfIterations) {
            double alphaTilde = -0.5 * theta0Der * alpha1 * alpha1 / (theta1 - theta0 - theta0Der * alpha1);
            alpha1 = Math.max(0.1 * alpha1, alphaTilde);
            theta1 = this.theta(spline, points, sw, alpha1, direction, fixedKnots);
            if (theta1 < 0) {
                return false;
            }
            q1 = this.mDelta;
            r1 = this.mPenalty;
            iteration += 1;
        }

        if (iteration > 0) {
            if (theta1 > theta0) {
                this.theta(spline, points, sw, alpha0, direction, fixedKnots);
            }
            return true;
        }

        double theta2 = this.theta(spline, points, sw, alpha2, direction, fixedKnots);
        if (theta2 < 0) {
            return false;
        }
        double q2 = this.mDelta;
        double r2 = this.mPenalty;

        while (theta2 < theta1) {
            alpha0 = alpha1;
            q0 = q1;
            r0 = r1;
            alpha1 = alpha2;
            theta1 = theta2;
            q1 = q2;
            r1 = r2;
            alpha2 = Math.min(2 * alpha1, 0.5 * (alphaMax + alpha1));
            theta2 = this.theta(spline, points, sw, alpha2, direction, fixedKnots);
            if (theta2 < 0) {
                return false;
            }
            q2 = this.mDelta;
            r2 = this.mPenalty;
        }

        double a0 = q0;
        double diff1 = alpha1 - alpha0;
        double diff2 = alpha2 - alpha0;
        double a2 = (q1 - q0) / diff1;
        a2 -= (q2 - q0) / diff2;
        a2 /= alpha1 - alpha2;
        double a1 = (q1 - a0) / diff1;
        a1 -= a2 * diff1;

        double fraction = diff1 / diff2;
        double numerator = r1 - r0 - fraction * (r2 - r0);
        double temp = Math.log((alphaMax - alpha1) / (alphaMax - alpha0));
        double denominator = temp - fraction * Math.log((alphaMax - alpha2) / (alphaMax - alpha0));
        double b2 = numerator / denominator;
        double b1 = (r1 - r0 - b2 * temp) / diff1;

        a = -2 * a2;
        b = -a * (alphaMax + alpha0) - this.p * b1 - a1;
        double c = (this.p * b1 + a1 + a * alpha0) * alphaMax - this.p * b2;

        double root1 = -0.5 * (b + Math.sqrt(b * b - 4 * a * c)) / a;
        double root2 = -b / a - root1;

        double alphaRes = 0.0;
        if (0 < root1 && root1 < alphaMax) {
            alphaRes = root1;
        } else if (0 < root2 && root2 < alphaMax) {
            alphaRes = root2;
        }

        double thetaRes = this.theta(spline, points, sw, alphaRes, direction, fixedKnots);

        return thetaRes >= 0;
    }

    // Checked!
    public static boolean initiateGrid(Spline spline, Points points) {
        int k = spline.getDegree();
        int g = spline.getInternalKnotsNum();
        ArrayList<Double> knots = new ArrayList<>(Collections.nCopies(g + 2, 0.0));
        knots.set(0, spline.getLeftBound());
        knots.set(g + 1, spline.getRightBound());
        int n = points.len();

        int uniqueSize = 0;
        int index = 0;

        while (index < n && points.getItem(index).getX() < knots.get(g + 1)) {
            if (index != 0 && points.getItem(index).getX() != points.getItem(index - 1).getX()) {
                uniqueSize += 1;
            }
            index += 1;
        }

        if (uniqueSize <= 0) {
            return false;
        }

        if (uniqueSize < g + k + 1) {
            return false;
        }

        double pointsPerKnot = (double)uniqueSize / (g + 1);
        int knotIndex = 1;
        int i = 1;
        int counter = 0;

        while (knotIndex < g + 1) {
            while (counter < knotIndex * pointsPerKnot || points.getItem(i).getX() == points.getItem(i - 1).getX()) {
                if (points.getItem(i).getX() != points.getItem(i - 1).getX()) {
                    counter += 1;
                }
                i += 1;
            }
            knots.set(knotIndex, 0.5 * (points.getItem(i).getX() + points.getItem(i - 1).getX()));
            knotIndex += 1;
        }

        spline.setKnots(knots);
        return true;
    }

    // Checked!
    public boolean approximateWithOptimalGrid(Spline spline, Points points, double smoothingWeight,
            double eps1, double eps2) {
        if (!CurveFitter.initiateGrid(spline, points)) {
            return false;
        }

        int g = spline.getInternalKnotsNum();

        if (!CurveFitter.approximate(spline, points, smoothingWeight)) {
            return false;
        }

        ArrayList<Double> direction = new ArrayList<>(Collections.nCopies(g, 0.0));
        ArrayList<Double> errorDerivative = new ArrayList<>(Collections.nCopies(g, 0.0));

        this.delta(spline, points, smoothingWeight);
        this.p = eps1 * this.mDelta * (spline.getRightBound() - spline.getLeftBound()) / (g + 1) / (g + 1);
        this.mError = this.mDelta + this.p * this.penalty(spline);

        for (int i = 0; i < g; ++i) {
            errorDerivative.set(i, this.errorGradient(spline, points, smoothingWeight, i + 1));
            direction.set(i, -errorDerivative.get(i));
        }

        double oldNorm = CurveFitter.norm(direction);
        double criteria1 = eps1 + eps2;
        double criteria2 = criteria1;
        int maxNumOfIter = 1000;

        int iteration = 0;
        double eps2Seq = eps2 * eps2;

        while ((criteria1 >= eps1 || criteria2 >= eps2Seq) && iteration < maxNumOfIter) {
            ArrayList<Double> copy = new ArrayList<>(spline.getKnots());
            ArrayList<Double> fixedKnots = copy;
            double oldError = this.mError;
            if (!this.specDimensionalMinimization(spline, points, smoothingWeight, direction, errorDerivative,
                    fixedKnots)) {
                spline.setKnots(fixedKnots);
                return CurveFitter.approximate(spline, points, smoothingWeight);
            }

            for (int i = 0; i < g; ++i) {
                errorDerivative.set(i, this.errorGradient(spline, points, smoothingWeight, i + 1));
            }

            double newNorm = CurveFitter.norm(errorDerivative);

            if (iteration % g == 0) {
                for (int i = 0; i < g; ++i) {
                    direction.set(i, -errorDerivative.get(i));
                }
            } else {
                double temp = newNorm / oldNorm;
                for (int i = 0; i < g; ++i) {
                    direction.set(i, direction.get(i) * temp);
                    direction.set(i, direction.get(i) - errorDerivative.get(i));
                }
            }

            double numerator = 0.0;
            double denominator = 0.0;

            ArrayList<Double> knots = spline.getKnots();
            for (int i = 1; i < knots.size() - 1; ++i) {
                double temp = knots.get(i) - fixedKnots.get(i);
                numerator += temp * temp;
                denominator += fixedKnots.get(i) * fixedKnots.get(i);
            }

            criteria1 = Math.abs(oldError - this.mError) / oldError;
            criteria2 = numerator / denominator;

            oldNorm = newNorm;
        }

        return true;
    }

    // Checked!

    private double mDelta;
    private double p;
    private double mError;
    private double mPenalty;
}
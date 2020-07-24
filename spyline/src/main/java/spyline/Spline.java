package spyline;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

public class Spline {
    Spline(ArrayList<Double> coefficients, ArrayList<Double> horizontalKnots, int splineDegree) {
        this.k = splineDegree;
        this.kFact = factorial(this.k);
        this.coefficients = new ArrayList<>(coefficients);
        Collections.sort(horizontalKnots);
        this.g = horizontalKnots.size() - 2;
        this.a = horizontalKnots.get(0);
        this.b = horizontalKnots.get(this.g + 1);
        this.knots = new ArrayList<Double>();
        for (int i = 0; i < this.k + 1; ++i) {
            this.knots.add(this.a);
        }
        for (int i = 0; i < this.g; ++i) {
            this.knots.add(0.0);
        }
        for (int i = 0; i < this.k + 1; ++i) {
            this.knots.add(this.b);
        }
        for (int i = 0; i < this.g; ++i) {
            this.knots.set(i + this.k + 1, horizontalKnots.get(i + 1));
        }
    }

    private int factorial(int n) {
        if (n <= 2) {
            return n;
        }
        return n * factorial(n - 1);
    }

    public int getLeftNodeIndex(double point, int minId) {
        if (point < this.a || point > this.b) {
            return -1;
        }
        int l = minId;
        while ((l < this.g + this.k) && (this.knots.get(l) > point || this.knots.get(l + 1) <= point)) {
            l += 1;
        }
        return l;
    }

    public int getInternalKnotsNum() {
        return this.g;
    }

    public int getDegree() {
        return this.k;
    }

    public double getLeftBound() {
        return this.a;
    }

    public double getRightBound() {
        return this.b;
    }

    public int getDegreeFactorial() {
        return this.kFact;
    }

    public ArrayList<Double> getCoefficients() {
        return this.coefficients;
    }

    public void setCoefficients(List<Double> coefficients) {
        this.coefficients = new ArrayList<>(coefficients);
    }

    public void setLeftEdge(double leftEdge) {
        this.a = leftEdge;
    }

    public void setRightEdge(double rightEdge) {
        this.b = rightEdge;
    }

    public ArrayList<Double> getKnots() {
        ArrayList<Double> slicedKnots = new ArrayList<>();
        for (int i = this.k; i < this.k + this.g + 2; ++i) {
            slicedKnots.add(0.0);
            slicedKnots.set(i - this.k, this.knots.get(i));
        }
        return slicedKnots;
    }

    public void setKnots(ArrayList<Double> horizontalKnots) {
        Collections.sort(horizontalKnots);
        this.g = horizontalKnots.size() - 2;
        this.a = horizontalKnots.get(0);
        this.b = horizontalKnots.get(this.g + 1);
        ArrayList<Double> localKnots = new ArrayList<>();
        for (int i = 0; i < this.k + 1; ++i) {
            localKnots.add(this.a);
        }
        for (int i = 0; i < this.g; ++i) {
            localKnots.add(0.0);
        }
        for (int i = 0; i < this.k + 1; ++i) {
            localKnots.add(this.b);
        }
        this.knots = localKnots;
        for (int i = 0; i < this.g; ++i) {
            this.knots.set(i + this.k + 1, horizontalKnots.get(i + 1));
        }
    }

    public double bSpline(double point, int deg, int knotId) {
        if (point < this.knots.get(knotId) || point > this.knots.get(knotId + deg + 1)) {
            return 0.0;
        }

        if (deg == 0) {
            return (point != this.knots.get(knotId + deg + 1)) ? 1.0 : 0.0;
        }

        if (this.knots.get(knotId + deg) < this.knots.get(knotId + deg + 1)) {
            int j = 0;
            while (j < deg && this.knots.get(knotId + j).equals(knots.get(knotId + j + 1))) {
                j += 1;
            }
            if (j == deg) {
                return Math.pow((this.knots.get(knotId + deg + 1) - point)
                        / (this.knots.get(knotId + deg + 1) - this.knots.get(knotId)), deg);
            }
        }

        int l = this.getLeftNodeIndex(point, knotId);
        ArrayList<Double> buff = new ArrayList<>(Collections.nCopies(deg + 1, 0.0));
        buff.set(knotId - 1 + deg, 1.0);

        for (int j = 1; j < deg + 1; ++j) {
            for (int i = l; i >= l - deg + j; --i) {
                double alpha = (point - this.knots.get(i)) / (this.knots.get(i + 1 + deg - j) - this.knots.get(i));
                buff.set(i - l + deg, alpha * buff.get(i - l + deg) + (1 - alpha) * buff.get(i - 1 - l + deg));
            }
        }

        return buff.get(deg);
    }

    // Until this Checked!
    ArrayList<Double> bSplines(double point, int deg) {
        ArrayList<Double> buff = new ArrayList<>(Collections.nCopies(this.k + 1, 0.0));
        if (deg > this.k) {
            return buff;
        }

        int l = this.getLeftNodeIndex(point, 0);
        buff.set(deg, 1.0);

        for (int r = 1; r < deg + 1; ++r) {
            int v = l - r + 1;
            double w2 = (this.knots.get(v + r) - point) / (this.knots.get(v + r) - this.knots.get(v));
            buff.set(deg - r, w2 * buff.get(deg - r + 1));
            for (int i = deg - r + 1; i < deg; ++i) {
                double w1 = w2;
                v += 1;
                w2 = (this.knots.get(v + r) - point) / (this.knots.get(v + r) - this.knots.get(v));
                buff.set(i, (1 - w1) * buff.get(i) + w2 * buff.get(i + 1));
            }
            buff.set(deg, (1 - w2) * buff.get(deg));
        }
        return buff;
    }

    // Checked!
    public void insertNode(int coordinate) {
        int j = this.getLeftNodeIndex(coordinate, 0);

        if (j < 0 || this.knots.get(j) == coordinate) {
            return;
        }

        this.knots.add((double) coordinate);
        Collections.sort(this.knots);
        this.coefficients.add(this.coefficients.get(this.g + this.k));

        for (int i = this.g + this.k; i >= j + 1; --i) {
            this.coefficients.set(i, this.coefficients.get(i - 1));
        }

        for (int i = j; i >= j - this.k + 1; --i) {
            double ri = (coordinate - this.knots.get(i)) / (this.knots.get(i + this.k + 1) - this.knots.get(i));
            this.coefficients.set(i, ri * this.coefficients.get(i) + (1 - ri) * this.coefficients.get(i - 1));
        }

        this.g += 1;
    }

    // Checked!
    public double bSplineDerivative(double point, int l, int i, int derDegree) {
        if (derDegree == 0) {
            return this.bSpline(point, l, i);
        }

        if (l == 0) {
            return 0.0;
        }

        double spline = 0.0;
        double c1 = this.knots.get(i + l) - this.knots.get(i);
        double c2 = this.knots.get(i + l + 1) - this.knots.get(i + 1);
        if (c1 != 0) {
            spline += this.bSplineDerivative(point, l - 1, i, derDegree - 1) / c1;
        }
        if (c2 != 0) {
            spline -= this.bSplineDerivative(point, l - 1, i + 1, derDegree - 1) / c2;
        }
        return l * spline;
    }

    // Checked!
    public double getValue(double point) {
        if (point < this.a || point > this.b) {
            return 0.0;
        }

        int l = this.getLeftNodeIndex(point, 0);
        if (l < 0) {
            return 0;
        }

        ArrayList<Double> buff = new ArrayList<>(Collections.nCopies(this.k + 1, 0.0));

        for (int i = 0; i < this.k + 1; ++i) {
            buff.set(i, this.coefficients.get(i + l - this.k));
        }

        for (int j = 1; j < this.k + 1; ++j) {
            for (int i = l; i >= l - this.k + j; --i) {
                double alpha = (point - this.knots.get(i)) / (this.knots.get(i + 1 + this.k - j) - this.knots.get(i));
                buff.set(i - l + this.k, alpha * buff.get(i - l + this.k) + (1 - alpha) * buff.get(i - 1 - l + this.k));
            }
        }

        return buff.get(this.k);
    }
    // Checked!

    public double getValueDerivative(double point, int derDegree) {
        if (derDegree == 0) {
            return this.getValue(point);
        }

        if (derDegree > this.k) {
            return 0.0;
        }

        int l = this.getLeftNodeIndex(point, 0);

        if (l < 0) {
            if (derDegree > 1) {
                return 0.0;
            }
            return -1.0;
        }

        double alpha = 1.0;
        double spline = 0.0;

        for (int i = 0; i < derDegree; ++i) {
            alpha *= (this.k - i);
        }

        ArrayList<Double> buff = new ArrayList<>(Collections.nCopies(this.k + 1, 0.0));
        for (int i = 0; i < this.k + 1; ++i) {
            buff.set(i, this.coefficients.get(i + l - this.k));
        }

        for (int j = 1; j < derDegree; ++j) {
            for (int i = j; i >= l - this.k + j; --i) {
                buff.set(i - l + this.k, (buff.get(i - l + this.k) - buff.get(i - l - 1 + this.k))
                        / (this.knots.get(i + 1 + this.k - j) - this.knots.get(i)));
            }
        }

        for (int i = derDegree; i < this.k + 1; ++i) {
            spline += buff.get(i) * this.bSpline(point, this.k - derDegree, l + i - this.k);
        }

        return alpha * spline;
    }

    // Checked!
    public double getLeadDerivativeDifference(int i, int q) {
        if (i < q - this.k - 1 || i > q) {
            return 0.0;
        }
        double numerator = (2 * (this.k % 2) - 1) * this.kFact * (this.knots.get(i + this.k + 1) - this.knots.get(i));
        double denominator = 1.0;
        for (int j = i; j < i + this.k + 2; ++j) {
            if (j != q) {
                denominator *= this.knots.get(q) - this.knots.get(j);
            }
        }
        return numerator / denominator;
    }

    // Checked!
    public double getLeadDerDiffDerKnot(double lddk, int i, int q, int l) {
        if (l < i || l > i + this.k + 1) {
            return 0.0;
        }

        if (l != i && l != q && l != i + this.k + 1) {
            return lddk / (this.knots.get(q) - this.knots.get(l));
        }

        double c = (double) (2 * (this.k % 2) - 1) * this.kFact;
        double product = c / lddk;
        double totalSum = 0.0;

        if (q != i && q != i + this.k + 1) {
            if (l == i) {
                double temp = (this.knots.get(i + this.k + 1) - this.knots.get(i))
                        * (this.knots.get(i + this.k + 1) - this.knots.get(q));
                return lddk / (this.knots.get(q) - this.knots.get(i) / temp);
            }
            if (l == i + this.k + 1) {
                double temp = (this.knots.get(i + this.k + 1) - this.knots.get(i))
                        * (this.knots.get(q) - this.knots.get(i));
                return lddk / (this.knots.get(q) - this.knots.get(i + this.k + 1)) / temp;
            }
            if (q == l) {
                product *= this.knots.get(i + this.k + 1) - this.knots.get(i);
                for (int j = i; j < i + this.k + 2; ++j) {
                    if (j != q) {
                        totalSum += product / (this.knots.get(q) - this.knots.get(j));
                    }
                }
                product *= product;
                return -c * (this.knots.get(i + this.k + 1) - this.knots.get(i)) * totalSum / product;
            }
        } else {
            for (int j = i + 1; j < i + this.k + 1; ++j) {
                totalSum += product / (this.knots.get(q) - this.knots.get(j));
            }
            return -c * totalSum / (product * product);
        }
        return 0.0;
    }

    // Checked!
    public double getLeadDerDiffDerKnot(int i, int q, int l) {
        if (l < i || l > i + this.k + 1) {
            return 0.0;
        }
        return this.getLeadDerDiffDerKnot(this.getLeadDerivativeDifference(i, q), i, q, l);
    }

    // Checked!
    public double getValueDerivativeKnot(double point, int knotId) {
        if (point < this.a) {
            if (knotId == this.k + 1) {
                double numerator = -this.k * (this.coefficients.get(1) - this.coefficients.get(0)) * (point - this.a);
                double denominator = (this.knots.get(knotId) - this.a) * (this.knots.get(knotId) - this.a);
                return numerator / denominator;
            }
            return 0.0;
        }

        if (point > this.b) {
            if (knotId == this.g + this.k) {
                double numerator = -this.k
                        * (this.coefficients.get(this.g + this.k) - this.coefficients.get(this.g + this.k - 1))
                        * (point - this.b);
                double denominator = (this.knots.get(knotId) - this.b) * (this.knots.get(knotId) - this.b);
                return numerator / denominator;
            }
            return 0.0;
        }

        if (point <= this.knots.get(knotId - this.k) || point >= this.knots.get(knotId + this.k)) {
            return 0.0;
        }

        int l = this.getLeftNodeIndex(point, 0);
        if (l < 0) {
            return 0.0;
        }

        if (l >= knotId) {
            l += 1;
        }

        ArrayList<Double> buff = new ArrayList<>(Collections.nCopies(this.k + 1, 0.0));

        for (int i = 0; i < this.k + 1; ++i) {
            if (i < knotId - l || i > knotId - l + this.k) {
                buff.set(i, 0.0);
            } else {
                buff.set(i, this.coefficients.get(i + l - this.k - 1) - this.coefficients.get(i + l - this.k));
                if (i + l + 1 <= knotId) {
                    buff.set(i, buff.get(i) / (this.knots.get(i + l + 1) - this.knots.get(i + l - this.k)));
                } else if (i <= knotId) {
                    buff.set(i, buff.get(i) / (this.knots.get(i + l) - this.knots.get(i + l - this.k)));
                } else {
                    buff.set(i, buff.get(i) / (this.knots.get(i + l) - this.knots.get(i + l - this.k - 1)));
                }
            }
        }
        double alpha;
        for (int j = 1; j < this.k + 1; ++j) {
            for (int i = l; i >= l - this.k + j; --i) {
                if (i + 1 + this.k - j <= knotId) {
                    alpha = (point - this.knots.get(i)) / (this.knots.get(i + 1 + this.k - j) - this.knots.get(i));
                } else if (i <= knotId) {
                    alpha = (point - this.knots.get(i)) / (this.knots.get(i + this.k - j) - this.knots.get(i));
                } else {
                    alpha = (point - this.knots.get(i - 1)) / (this.knots.get(i + this.k - j) - this.knots.get(i - 1));
                }
                buff.set(i - l + this.k, alpha * buff.get(i - l + this.k) + (1 - alpha) * buff.get(i - 1 - l + this.k));
            }
        }
        return buff.get(this.k);
    }

    private int k = 3;
    private int kFact;
    private ArrayList<Double> coefficients;
    private int g;
    private double a;
    private double b;
    private ArrayList<Double> knots;
}

// Everything checked!
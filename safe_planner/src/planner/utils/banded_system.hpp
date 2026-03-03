#pragma once
#include <algorithm>
namespace safe_planner::planner::trajectory::utils{
class BandedSystem {
public:
    BandedSystem() = default;

    ~BandedSystem() {
        destroy();
    }

    // 拷贝构造
    BandedSystem(const BandedSystem& other)
        : N(other.N), lowerBw(other.lowerBw), upperBw(other.upperBw)
    {
        if (other.ptrData) {
            int size = N * (lowerBw + upperBw + 1);
            ptrData = new double[size];
            std::copy(other.ptrData, other.ptrData + size, ptrData);
        }
    }

    // 拷贝赋值
    BandedSystem& operator=(const BandedSystem& other) {
        if (this != &other) {
            destroy();
            N = other.N;
            lowerBw = other.lowerBw;
            upperBw = other.upperBw;

            if (other.ptrData) {
                int size = N * (lowerBw + upperBw + 1);
                ptrData = new double[size];
                std::copy(other.ptrData, other.ptrData + size, ptrData);
            }
        }
        return *this;
    }

    // 移动构造
    BandedSystem(BandedSystem&& other) noexcept
        : N(other.N), lowerBw(other.lowerBw), upperBw(other.upperBw), ptrData(other.ptrData)
    {
        other.ptrData = nullptr;
        other.N = other.lowerBw = other.upperBw = 0;
    }

    // 移动赋值
    BandedSystem& operator=(BandedSystem&& other) noexcept {
        if (this != &other) {
            destroy();
            N = other.N;
            lowerBw = other.lowerBw;
            upperBw = other.upperBw;
            ptrData = other.ptrData;

            other.ptrData = nullptr;
            other.N = other.lowerBw = other.upperBw = 0;
        }
        return *this;
    }

    void create(int n, int p, int q) {
        destroy();
        N = n;
        lowerBw = p;
        upperBw = q;
        int size = N * (lowerBw + upperBw + 1);
        ptrData = new double[size];
        std::fill(ptrData, ptrData + size, 0.0);
    }

    void destroy() {
        delete[] ptrData;
        ptrData = nullptr;
    }

    private:
        int N;
        int lowerBw;
        int upperBw;
        // Compulsory nullptr initialization here
        double *ptrData = nullptr;

    public:
        // Reset the matrix to zero
        inline void reset(void)
        {
            std::fill_n(ptrData, N * (lowerBw + upperBw + 1), 0.0);
            return;
        }

        // The band matrix is stored as suggested in "Matrix Computation"
        inline const double &operator()(const int &i, const int &j) const
        {
            return ptrData[(i - j + upperBw) * N + j];
        }

        inline double &operator()(const int &i, const int &j)
        {
            return ptrData[(i - j + upperBw) * N + j];
        }

        // This function conducts banded LU factorization in place
        // Note that NO PIVOT is applied on the matrix "A" for efficiency!!!
        inline void factorize_lu()
        {
            int iM, jM;
            double cVl;
            for (int k = 0; k <= N - 2; k++)
            {
                iM = std::min(k + lowerBw, N - 1);
                cVl = operator()(k, k);
                for (int i = k + 1; i <= iM; i++)
                {
                    if (operator()(i, k) != 0.0)
                    {
                        operator()(i, k) /= cVl;
                    }
                }
                jM = std::min(k + upperBw, N - 1);
                for (int j = k + 1; j <= jM; j++)
                {
                    cVl = operator()(k, j);
                    if (cVl != 0.0)
                    {
                        for (int i = k + 1; i <= iM; i++)
                        {
                            if (operator()(i, k) != 0.0)
                            {
                                operator()(i, j) -= operator()(i, k) * cVl;
                            }
                        }
                    }
                }
            }
            return;
        }

        // This function solves Ax=b, then stores x in b
        // The input b is required to be N*m, i.e.,
        // m vectors to be solved.
        template <typename EIGENMAT>
        inline void solve(EIGENMAT &b) const
        {
            int iM;
            for (int j = 0; j <= N - 1; j++)
            {
                iM = std::min(j + lowerBw, N - 1);
                for (int i = j + 1; i <= iM; i++)
                {
                    if (operator()(i, j) != 0.0)
                    {
                        b.row(i) -= operator()(i, j) * b.row(j);
                    }
                }
            }
            for (int j = N - 1; j >= 0; j--)
            {
                b.row(j) /= operator()(j, j);
                iM = std::max(0, j - upperBw);
                for (int i = iM; i <= j - 1; i++)
                {
                    if (operator()(i, j) != 0.0)
                    {
                        b.row(i) -= operator()(i, j) * b.row(j);
                    }
                }
            }
            return;
        }

        // This function solves ATx=b, then stores x in b
        // The input b is required to be N*m, i.e.,
        // m vectors to be solved.
        template <typename EIGENMAT>
        inline void solve_adj(EIGENMAT &b) const
        {
            int iM;
            for (int j = 0; j <= N - 1; j++)
            {
                b.row(j) /= operator()(j, j);
                iM = std::min(j + upperBw, N - 1);
                for (int i = j + 1; i <= iM; i++)
                {
                    if (operator()(j, i) != 0.0)
                    {
                        b.row(i) -= operator()(j, i) * b.row(j);
                    }
                }
            }
            for (int j = N - 1; j >= 0; j--)
            {
                iM = std::max(0, j - lowerBw);
                for (int i = iM; i <= j - 1; i++)
                {
                    if (operator()(j, i) != 0.0)
                    {
                        b.row(i) -= operator()(j, i) * b.row(j);
                    }
                }
            }
            return;
        }
    };
}
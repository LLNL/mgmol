#ifndef MGMOL_LOCALVECTOR
#define MGMOL_LOCALVECTOR

#include <vector>

template<class T>
class LocalVector
{
private:
    std::vector<T> data_;

public:
    LocalVector(const int n)
        : data_(n)
    {}

    LocalVector(const std::vector<T>& v) :
        data_(v)
    {}

    T* data(){ return data_.data(); }
    const T* data()const{ return data_.data(); }

    void normalize()
    {
        double norm = Tnrm2(data_.size(), data_.data());
        Tscal(data_.size(), 1. / norm, data_.data());
    }

    T dot(const LocalVector<T>& v)const
    {
        return Tdot(data_.size(), data_.data(), v.data_.data());
    }

    void swap(LocalVector<T>& v)
    {
        data_.swap(v.data_);
    }

    T scaledDiff2(const LocalVector<T>& v, const double theta)
    {
        double diff2 = 0.;
        int n        = v.data_.size();

        for(int i=0;i<n;i++)
        {
            double tmp = data_[i] - theta * v.data_[i];
            diff2 += tmp * tmp;
         }
         return diff2;
    }
};

#endif

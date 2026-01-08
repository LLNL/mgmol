#ifndef MGMOL_LOCALVECTOR
#define MGMOL_LOCALVECTOR

#include "memory_space.h"
#include "mputils.h"

#include <vector>

template <typename DataType, typename MemorySpaceType>
class LocalVector
{
private:
    std::vector<DataType> data_;

public:
    LocalVector(const int n) : data_(n) { }

    LocalVector(const std::vector<DataType>& v) : data_(v) { }

    DataType* data() { return data_.data(); }
    const DataType* data() const { return data_.data(); }

    void normalize()
    {
        double norm = Tnrm2(data_.size(), data_.data());
        Tscal(data_.size(), 1. / norm, data_.data());
    }

    DataType dot(const LocalVector<DataType, MemorySpaceType>& v) const
    {
        return Tdot(data_.size(), data_.data(), v.data_.data());
    }

    void swap(LocalVector<DataType, MemorySpaceType>& v)
    {
        data_.swap(v.data_);
    }
    void swap(std::vector<DataType>& v) { data_.swap(v); }

    DataType scaledDiff2(
        const LocalVector<DataType, MemorySpaceType>& v, const double theta)
    {
        double diff2 = 0.;
        int n        = v.data_.size();

        for (int i = 0; i < n; i++)
        {
            double tmp = data_[i] - theta * v.data_[i];
            diff2 += tmp * tmp;
        }
        return diff2;
    }
};

#endif

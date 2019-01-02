#include "LDAonGrid.h"
#include "PBEonGrid.h"
#include "PBEonGridSpin.h"
#include "Potentials.h"

#include <iostream>

template <class T>
class XCfunctionalFactory
{
public:
    static XConGrid* create(const int xctype, const int nspin,
                            Rho<T>& rho,
                            Potentials& pot)
    {
        if (xctype == 0)
        {
            return new LDAonGrid<T>(rho, pot);
        }
        else if (xctype == 2)
        {
            if (nspin > 1)
                return new PBEonGridSpin<T>(rho, pot);
            else
                return new PBEonGrid<T>(rho, pot);
        }
        else
        {
            std::cerr << "Invalid XC option" << std::endl;
            return 0;
        }
    }
};

